
#include "sas.h"
#include "preprocess.h"
#include "parser.h"
#include "plans.h"
#include "forward.h"
#include "nodeset.h"
#include "sas_heuristic.h"
#include "cost_table.h"
#include "bfs.h"

#define VALIDATE_ABSTRACT_PLAN
#define VALIDATE_FINAL_PLAN
#define USE_RESET_LIST

void print_sas_short(const HSPS::SASInstance& ins, std::ostream& to)
{
  for (HSPS::index_type k = 0; k < ins.n_variables(); k++) {
    ins.write_variable(to, ins.variables[k]);
    to << std::endl;
  }
  to << "init: ";
  ins.write_partial_state(to, ins.init_state);
  to << std::endl << "goal: ";
  ins.write_partial_state(to, ins.goal_state);
  to << std::endl;
}

BEGIN_HSPS_NAMESPACE

void compute_dtg
(const SASInstance& ins, index_type var, weighted_graph& g)
{
  assert(var < ins.n_variables());
  g.init(ins.variables[var].n_values());
  for (index_type k = 0; k < ins.n_actions(); k++) {
    index_type v_post = ins.actions[k].post.value_of(var);
    if (v_post != no_such_index) {
      assert(v_post < ins.variables[var].n_values());
      index_type v_pre = ins.actions[k].pre.value_of(var);
      assert(v_pre != no_such_index);
      assert(v_pre < ins.variables[var].n_values());
      if (g.adjacent(v_pre, v_post)) {
	if (ins.actions[k].cost < g.weight(v_pre, v_post))
	  g.set_weight(v_pre, v_post, ins.actions[k].cost);
      }
      else {
	g.add_edge(v_pre, v_post, ins.actions[k].cost);
      }
    }
  }
}

void make_value_map_from_dvs
(const SASInstance& ins,
 const index_set_vec& dvs,
 value_map& map)
{
  assert(dvs.size() == ins.n_variables());
  map.set_length(ins.n_variables());
  for (index_type k = 0; k < ins.n_variables(); k++) {
    map[k].assign_identity(ins.variables[k].n_values());
    assert(ins.init_state.value_of(k) != no_such_index);
    if (!dvs[k].empty()) {
      index_type other = no_such_index;
      if (!dvs[k].contains(ins.init_state.value_of(k)))
	other = ins.init_state.value_of(k);
      for (index_type i = 0; i < ins.variables[k].n_values(); i++) {
	if (!dvs[k].contains(i)) {
	  if (other == no_such_index)
	    other = i;
	  map[k][i] = other;
	}
      }
    }
    else {
      index_type other = ins.init_state.value_of(k);
      for (index_type i = 0; i < ins.variables[k].n_values(); i++)
	map[k][i] = other;
    }
  }
}

// create an map from the internal state representation of the abstraction
// defined by m1 to the internal state representation of the abstraction
// defined by m0; m1 must be a refinement of m0.
void make_internal_map
(const value_map& m0,
 const value_map& m1,
 value_map& map)
{
  assert(m0.size() == m1.size());
  map.set_length(m0.size());
  for (index_type k = 0; k < m0.size(); k++) {
    mapping c0;
    mapping::compact(m0[k], c0);
    mapping c1;
    mapping::compact(m1[k], c1);
    map[k].assign_identity(c1.range());
    for (index_type i = 0; i < map[k].size(); i++) {
      index_set a;
      mapping::inverse_map_image(c1, i, a);
      index_set b;
      mapping::map_image(c0, a, b);
      assert(b.size() == 1);
      map[k][i] = b[0];
    }
    //std::cerr << "m0 = " << m0[k] << ", m1 = " << m1[k]
    //      << ", result = " << map[k] << std::endl;
  }
}

void compute_distances
(const weighted_graph& g, bool reverse, index_type src, cost_vec& dist)
{
  dist.assign_value(POS_INF, g.size());
  dist[src] = 0;
  bool done = false;
  while (!done) {
    done = true;
    for (index_type i = 0; i < g.size(); i++)
      if (FINITE(dist[i])) {
	if (reverse) {
	  for (index_type j = 0; j < g.predecessors(i).size(); j++)
	    if ((dist[i] + g.weight(g.predecessors(i)[j], i))
		< dist[g.predecessors(i)[j]]) {
	      dist[g.predecessors(i)[j]] =
		dist[i] + g.weight(g.predecessors(i)[j], i);
	      done = false;
	    }
	}
	else {
	  for (index_type j = 0; j < g.successors(i).size(); j++)
	    if ((dist[i] + g.weight(i, g.successors(i)[j]))
		< dist[g.successors(i)[j]]) {
	      dist[g.successors(i)[j]] =
		dist[i] + g.weight(i, g.successors(i)[j]);
	      done = false;
	    }
	}
      }
  }
}

void connected_without
(const weighted_graph& g, index_type n0, const bool_vec& fs, bool_vec& rs)
{
  rs.assign_value(false, g.size());
  rs[n0] = true;
  bool done = false;
  while (!done) {
    done = true;
    for (index_type i = 0; i < g.size(); i++)
      if (rs[i]) {
	for (index_type j = 0; j < g.bidirectional(i).size(); j++)
	  if (!rs[g.bidirectional(i)[j]] && !fs[g.bidirectional(i)[j]]) {
	    rs[g.bidirectional(i)[j]] = true;
	    done = false;
	  }
      }
  }
}

void make_partition
(const SASInstance& ins,
 index_type var,
 const weighted_graph& dtg,
 index_type val0,
 index_type val1,
 equivalence& eq)
{
  assert(var < ins.n_variables());
  index_type n_val = ins.variables[var].n_values();
  assert(dtg.size() == n_val);
  assert(val0 < n_val);
  assert(val1 < n_val);
  cost_vec fwd_dist;
  compute_distances(dtg, false, val0, fwd_dist);
  cost_vec bwd_dist;
  compute_distances(dtg, true, val1, bwd_dist);
  //std::cerr << "dtg = " << dtg << std::endl;
  //std::cerr << "val0 = " << val0 << ", val1 = " << val1 << std::endl;
  //std::cerr << "fwd = " << fwd_dist << std::endl;
  //std::cerr << "bwd = " << bwd_dist << std::endl;
  assert(fwd_dist[val1] == bwd_dist[val0]);
  NTYPE opt = fwd_dist[val1];
  bool_vec on_opt_path(false, n_val);
  for (index_type i = 0; i < n_val; i++)
    if ((fwd_dist[i] + bwd_dist[i]) == opt)
      on_opt_path[i] = true;
  assert(on_opt_path[val0] && on_opt_path[val1]);
  eq.reset(n_val);
  for (index_type i = 0; i < n_val; i++)
    if (on_opt_path[i]) {
      for (index_type j = i + 1; j < n_val; j++)
	if (on_opt_path[j])
	  if (fwd_dist[i] == fwd_dist[j])
	    if ((i != val1) && (j != val1))
	      eq.merge(i, j);
    }
    else {
      bool_vec rsi;
      connected_without(dtg, i, on_opt_path, rsi);
      for (index_type j = i + 1; j < n_val; j++)
	if (!on_opt_path[j])
	  if (rsi[j] ||
	      (((fwd_dist[i] + bwd_dist[j]) > opt) &&
	       ((fwd_dist[j] + bwd_dist[i]) > opt)))
	    eq.merge(i, j);
    }
  //std::cerr << "eq = " << eq << std::endl;
}

void refine_value_map(const equivalence& eq, mapping& m)
{
  assert(eq.size() == m.size());
  equivalence ceq(m.size());
  for (index_type i = 0; i < m.size(); i++)
    for (index_type j = i + 1; j < m.size(); j++)
      if ((m[i] == m[j]) && eq(i, j))
	ceq.merge(i, j);
  //std::cerr << "m = " << m << std::endl;
  //std::cerr << "eq = " << eq << std::endl;
  //std::cerr << "c = " << ceq << std::endl;
  ceq.make_map(m);
  //std::cerr << "m' = " << m << std::endl;
}

NTYPE validate_sas_plan
(const SASInstance& ins, const index_vec& plan)
{
  NTYPE plancost = 0;
  partial_state s(ins.init_state);
  for (index_type k = 0; k < plan.size(); k++) {
    assert(plan[k] < ins.n_actions());
    if (!s.implies(ins.actions[plan[k]].pre)) {
      std::cerr << "precondition ";
      ins.write_partial_state(std::cerr, ins.actions[plan[k]].pre);
      std::cerr << " of action " << ins.actions[plan[k]].name
		<< " not met in ";
      ins.write_partial_state(std::cerr, s);
      std::cerr << std::endl;
      return POS_INF;
    }
    if (!s.implies(ins.actions[plan[k]].prv)) {
      std::cerr << "prevail-condition ";
      ins.write_partial_state(std::cerr, ins.actions[plan[k]].prv);
      std::cerr << " of action " << ins.actions[plan[k]].name
		<< " not met in ";
      ins.write_partial_state(std::cerr, s);
      std::cerr << std::endl;
      return POS_INF;
    }
    s.assign(ins.actions[plan[k]].post);
    plancost += ins.actions[plan[k]].cost;
  }
  if (!s.implies(ins.goal_state)) {
    std::cerr << "goal ";
    ins.write_partial_state(std::cerr, ins.goal_state);
    std::cerr << " not met in ";
    ins.write_partial_state(std::cerr, s);
    std::cerr << std::endl;
    return POS_INF;
  }
  return plancost;
}

class AbstractionLevel {
public:
  typedef enum {fwd_blind, bwd_blind, fwd_ha_light, bwd_ha_light,
		fwd_ha, bwd_ha,	fwd_lmcut, bwd_lmcut, switchback}
    solver_type;

  static count_type fwd_successors;
  static count_type fwd_expansions;
  static count_type bwd_successors;
  static count_type bwd_expansions;

private:
  index_type level;
  solver_type solver;
  const SASInstance& baseins;
  const lvector<weighted_graph>& dtgs;
  const index_set_vec& dvs;
  const value_map& vm;
  Statistics& stats;
  int verbose_level;

  SASInstance absins;
  NTYPE mincost;

  AbstractionLevel* above;
  value_map   lvm;
  HashNodeSet nodes;
  NodeQueue   open;
#ifdef USE_RESET_LIST
  node_vec    reset;
#endif

  class ALState : public ProgressionState {
    AbstractionLevel& level; // points to *this* level
    partial_state state;
    index_type    act;

    // successor constructor
    ALState(ALState& p, index_type a);

  public:
    ALState(AbstractionLevel& l, const partial_state& s0)
      : level(l), state(s0), act(no_such_index) { };
    ALState(const ALState& s)
      : level(s.level), state(s.state), act(s.act) { };

    virtual Transition* transition();
    virtual NTYPE delta_cost();
    virtual NTYPE est_cost();
    virtual bool  is_final();
    virtual bool  is_max();
    virtual NTYPE expand(Search& s, NTYPE bound);
    virtual int compare(const State& s);
    virtual index_type hash();
    virtual State* copy();
    virtual void insert(Plan& p);
    virtual void write(::std::ostream& s);
  };

  class ALNewState : public Search {
    AbstractionLevel& level; // points to *this* level
    Node* current_node;
    count_type n_new;
  public:
    ALNewState(AbstractionLevel& l, Node* n)
      : level(l), current_node(n), n_new(0) { };
    ~ALNewState()
    {
      if (current_node && !level.stats.break_signal_raised()) {
	if (level.reverse()) {
	  AbstractionLevel::bwd_successors += n_new;
	  AbstractionLevel::bwd_expansions += 1;
	}
	else {
	  AbstractionLevel::fwd_successors += n_new;
	  AbstractionLevel::fwd_expansions += 1;
	}
      }
    };
    virtual NTYPE new_state(State& s, NTYPE bound);
    virtual bool  solved() const { return false; };
    virtual bool  optimal() const { return false; };
    virtual bool  done() const { return false; };
  };

  friend class ALState;
  friend class ALNewState;

  bool reverse() const {
    switch (solver) {
    case fwd_blind:
    case fwd_ha_light:
    case fwd_ha:
    case fwd_lmcut:
      return false;
    case bwd_blind:
    case bwd_ha_light:
    case bwd_ha:
    case bwd_lmcut:
      return true;
    case switchback:
      return ((level % 2) == 0);
    default:
      assert(0);
    }
  };

  NTYPE myopic(const partial_state& s);
  void  clear_search_info(node_vec& ns);
  void  clear_search_info();
  void cache_pg(NTYPE sc, node_vec& ns);
  void  cache_pg(NTYPE sc);
  Node* lookup(const partial_state& s);

  Node* switch_back(const partial_state& g);
  NTYPE HA(const partial_state& g);
  NTYPE base_level_search(ActionSequence& plan);
  void  init_open();
  void  dealloc();

  NTYPE solve_blind(ActionSequenceSet& plans);

  bool check(const ActionSequence& ap,
	     const partial_state state,
	     index_type step,
	     ActionSequence& plan,
	     index_set_vec& cs);
  bool check(const ActionSequence& ap,
	     ActionSequence& plan,
	     index_set_vec& dvs);

  bool check2(const ActionSequence& ap,
	      const partial_state state,
	      index_type step,
	      ActionSequence& plan,
	      lvector<pair_set>& cs);
  bool check2(const ActionSequence& ap,
	      ActionSequence& plan,
	      index_set_vec& dvs);

  AbstractionLevel(index_type l,
		   const SASInstance& b,
		   const lvector<weighted_graph>& gs,
		   const index_set_vec& d,
		   const value_map& m,
		   const value_map& lm,
		   AbstractionLevel* a,
		   Statistics& s,
		   int v = 0);

public:
  AbstractionLevel(index_type l,
		   solver_type so,
		   const SASInstance& b,
		   const lvector<weighted_graph>& gs,
		   const index_set_vec& d,
		   const value_map& m,
		   Statistics& st,
		   int v = 0);

  NTYPE cegar(ActionSequence& plan, NTYPE ub, index_type db);
};

count_type AbstractionLevel::fwd_successors = 0;
count_type AbstractionLevel::fwd_expansions = 0;
count_type AbstractionLevel::bwd_successors = 0;
count_type AbstractionLevel::bwd_expansions = 0;

AbstractionLevel::AbstractionLevel
(index_type l,
 solver_type so,
 const SASInstance& b,
 const lvector<weighted_graph>& gs,
 const index_set_vec& d,
 const value_map& m,
 Statistics& s,
 int v)
  : level(l), solver(so), baseins(b), dtgs(gs), dvs(d), vm(m), stats(s),
    verbose_level(v), absins(b, m), mincost(POS_INF), above(0)
{
  absins.cross_reference();
  if (verbose_level > 0) {
    std::cerr << "abstraction level " << level << ": signature = "
	      << absins.signature << " (" << absins.n_actions()
	      << " actions)" << std::endl;
  }
  for (index_type k = 0; k < absins.n_actions(); k++)
    mincost = MIN(absins.actions[k].cost, mincost);
}

AbstractionLevel::AbstractionLevel
(index_type l,
 const SASInstance& b,
 const lvector<weighted_graph>& gs,
 const index_set_vec& d,
 const value_map& m,
 const value_map& lm,
 AbstractionLevel* a,
 Statistics& s,
 int v)
  : level(l), solver(a->solver), baseins(b), dtgs(gs), dvs(d), vm(m), stats(s),
    verbose_level(v), absins(b, m), mincost(a->mincost), above(a), lvm(lm)
{
  absins.cross_reference();
  if (verbose_level > 0) {
    std::cerr << "level " << level << ": signature = "
	      << absins.signature << " (" << absins.n_actions()
	      << " actions)" << std::endl;
  }
}

NTYPE AbstractionLevel::myopic(const partial_state& s)
{
  if (reverse()) {
    if (s.implies(absins.init_state)) return 0;
  }
  else {
    if (s.implies(absins.goal_state)) return 0;
  }
  NTYPE c_min = POS_INF;
  for (index_type k = 0; k < absins.n_actions(); k++) {
    bool app = s.implies(absins.actions[k].prv);
    if (reverse())
      app = (app && s.implies(absins.actions[k].post));
    else
      app = (app && s.implies(absins.actions[k].pre));
    if (app) {
      c_min = MIN(c_min, absins.actions[k].cost);
    }
  }
  return c_min;
}

AbstractionLevel::ALState::ALState
(ALState& s, index_type a)
  : level(s.level), state(s.state), act(a)
{
  if (level.reverse())
    state.assign(level.absins.actions[a].pre);
  else
    state.assign(level.absins.actions[a].post);
  set_predecessor(&s);
}

Transition* AbstractionLevel::ALState::transition()
{
  if (act != no_such_index)
    return new SeqProgTrans(act, 1);
  else
    return 0;
}

NTYPE AbstractionLevel::ALState::delta_cost()
{
  if (act != no_such_index)
    return level.absins.actions[act].cost;
  else
    return 0;
}

NTYPE AbstractionLevel::ALState::est_cost()
{
  if (level.solver == fwd_lmcut) {
    SASLMCutBase lmc(level.absins);
    NTYPE v = lmc.compute(state, level.absins.goal_state,
			  SASCostACF(level.absins));
    return v;
  }
  if (level.solver == bwd_lmcut) {
    SASLMCutBase lmc(level.absins);
    NTYPE v = lmc.compute(level.absins.init_state, state,
			  SASCostACF(level.absins));
    return v;
  }
  if (level.solver == switchback) {
    if (level.above) {
      partial_state s_above(state, level.lvm);
      Node* n = level.above->switch_back(s_above);
      if (n)
	return n->acc;
      else if (level.stats.break_signal_raised())
	return 0;
      else
	return POS_INF;
    }
  }
  if ((level.solver == fwd_ha) || (level.solver == bwd_ha)) {
    if (level.above) {
      partial_state s_above(state, level.lvm);
      return level.above->HA(s_above);
    }
  }
  if ((level.solver == fwd_ha_light) || (level.solver == bwd_ha_light)) {
    if (level.above) {
      partial_state s_above(state, level.lvm);
      Node* n_above = level.above->lookup(s_above);
      if (n_above) {
	if (FINITE(n_above->opt)) return n_above->opt;
	return n_above->est;
      }
    }
  }
  // blind search, or at the top level
  if (is_final())
    return 0;
  else
    return level.mincost;
}

bool  AbstractionLevel::ALState::is_final()
{
  if (level.reverse())
    return state.implies(level.absins.init_state);
  else
    return state.implies(level.absins.goal_state);
}

bool  AbstractionLevel::ALState::is_max()
{
  return false;
}

// note: this does not implement the correct semantics of the
// expand function; it may call search.new_state with successor
// states that exceed the bound.
NTYPE AbstractionLevel::ALState::expand(Search& s, NTYPE bound)
{
  NTYPE c_min = POS_INF;
  for (index_type k = 0; k < level.absins.n_actions(); k++) {
    bool app = state.implies(level.absins.actions[k].prv);
    if (level.reverse())
      app = (app && state.implies(level.absins.actions[k].post));
    else
      app = (app && state.implies(level.absins.actions[k].pre));
    if (app) {
      ALState* s_new = new ALState(*this, k);
      NTYPE c_new = s.new_state(*s_new, bound - s_new->delta_cost());
      c_min = MIN(c_min, c_new);
      delete s_new;
      if (s.done()) {
	return c_new;
      }
    }
  }
  return c_min;
}

int AbstractionLevel::ALState::compare(const State& s)
{
  ALState& sbs = (ALState&)s;
  if (state < sbs.state)
    return -1;
  else if (state > sbs.state)
    return 1;
  else
    return 0;
}

index_type AbstractionLevel::ALState::hash()
{
  return level.absins.state_hash_function.index(state);
}

State* AbstractionLevel::ALState::copy()
{
  return new ALState(*this);
}

void AbstractionLevel::ALState::insert(Plan& p)
{
  if (act != no_such_index) {
    p.insert(act);
    p.advance(1);
  }
}

void AbstractionLevel::ALState::write(::std::ostream& s)
{
  level.absins.write_partial_state(s, state);
}

NTYPE AbstractionLevel::ALNewState::new_state
(State& s, NTYPE bound)
{
  n_new += 1;
  NTYPE s_delta = s.delta_cost();
  Node* n = level.nodes.insert_node(s);
  if (n->state) {
    assert(current_node);
    if ((current_node->acc + s_delta) < n->acc) {
      // if (n->exp > 0) {
      // 	std::cerr << "error: reached ";
      // 	n->write(std::cerr);
      // 	std::cerr << " from ";
      // 	current_node->write(std::cerr);
      // 	std::cerr << " with delta cost = " << s_delta << std::endl;
      // }
      if (level.solver == switchback) {
	assert(n->exp == 0);
	assert(n->pos != no_such_index);
      }
      n->bp_pre = current_node;
      if (n->bp_trans) delete n->bp_trans;
      n->bp_trans = s.transition();
      n->bp_trans->set_predecessor(current_node->bp_trans);
      n->bp_delta = s_delta;
      n->acc = current_node->acc + s_delta;
      n->val = n->acc + n->est;
      if (n->pos != no_such_index)
	level.open.shift_up(n->pos);
      else
	level.open.enqueue(n);
    }
  }
  else {
    n->state = s.copy();
    n->state->set_predecessor(0);
    if (current_node) {
      n->bp_pre = current_node;
      n->bp_trans = s.transition();
      n->bp_trans->set_predecessor(current_node->bp_trans);
      n->bp_delta = s_delta;
      n->acc = current_node->acc + s_delta;
    }
    else {
      n->bp_pre = 0;
      n->bp_trans = 0;
      n->bp_delta = s_delta;
      n->acc = s_delta;
    }
    n->est = s.est_cost();
    n->val = n->acc + n->est;
    n->exp = 0;
    level.open.enqueue(n);
  }
#ifdef USE_RESET_LIST
  if ((level.solver == fwd_ha) || (level.solver == bwd_ha))
    level.reset.append(n);
#endif
  return POS_INF;
}

Node* AbstractionLevel::switch_back(const partial_state& g)
{
  ALState sg(*this, g);
  Node* ng = nodes.find_node(sg);
  if (ng) {
    if (ng->closed)
      return ng;
  }
  while (!open.empty()) {
    if (stats.break_signal_raised()) return 0;
    Node* current_node = open.dequeue();
    // node should only ever be expanded once
    assert(current_node->exp == 0);
    if (FINITE(current_node->val)) {
      stats.expand_node(*(current_node->state));
      ALNewState ns(*this, current_node);
      current_node->state->expand(ns, POS_INF);
      current_node->exp += 1;
    }
    current_node->closed = true;
    if (current_node->state->compare(sg) == 0)
      return current_node;
  }
  return 0;
}

Node* AbstractionLevel::lookup(const partial_state& s)
{
  ALState ss(*this, s);
  Node* n = nodes.find_node(ss);
  return n;
}

void AbstractionLevel::clear_search_info(node_vec& ns)
{
  for (index_type k = 0; k < ns.length(); k++) {
    ns[k]->closed = false;
    ns[k]->acc = POS_INF;
    ns[k]->val = POS_INF;
    ns[k]->pos = no_such_index;
    ns[k]->exp = 0;
    ns[k]->bp_pre = 0;
    if (ns[k]->bp_trans) delete ns[k]->bp_trans;
    ns[k]->bp_trans = 0;
    ns[k]->bp_delta = 0;
  }
}

void AbstractionLevel::clear_search_info()
{
  node_vec v(0, 0);
  nodes.collect_nodes(v);
  clear_search_info(v);
}

void AbstractionLevel::cache_pg(NTYPE sc, node_vec& ns)
{
  for (index_type k = 0; k < ns.length(); k++)
    if (ns[k]->closed) {
      if (INFINITE(ns[k]->opt))
	if ((sc - ns[k]->acc) > ns[k]->est)
	  ns[k]->est = (sc - ns[k]->acc);
    }
}

void AbstractionLevel::cache_pg(NTYPE sc)
{
  node_vec ns(0, 0);
  nodes.collect_nodes(ns);
  cache_pg(sc, ns);
}

NTYPE AbstractionLevel::HA(const partial_state& s)
{
  ALState s0(*this, s);
  Node* n0 = nodes.insert_node(s0);
  if (n0->state) {
    if (FINITE(n0->opt))
      return n0->opt;
    if (INFINITE(n0->est))
      return n0->est;
  }
  if (!(n0->state)) {
    n0->state = s0.copy();
    n0->state->set_predecessor(0);
    n0->est = s0.est_cost();
  }
#ifdef USE_RESET_LIST
  clear_search_info(reset);
  reset.clear();
  reset.append(n0);
#else
  clear_search_info();
#endif
  n0->bp_pre = 0;
  n0->bp_trans = 0;
  n0->bp_delta = 0;
  n0->acc = 0;
  n0->val = n0->acc + n0->est;
  open.clear();
  open.enqueue(n0);
  while (!open.empty()) {
    if (stats.break_signal_raised()) return 0;
    Node* current_node = open.dequeue();
    // we've found an already solved node
    if (current_node->solved()) {
      current_node->cache_optimal_path(current_node->opt);
      assert(n0->opt == (current_node->acc + current_node->opt));
#ifdef USE_RESET_LIST
      cache_pg(n0->opt, reset);
#else
      cache_pg(n0->opt);
#endif
      return n0->opt;
    }
    // we've found a new goal node
    if (current_node->state->is_final()) {
      current_node->cache_optimal_path(0);
      assert(n0->opt == current_node->acc);
#ifdef USE_RESET_LIST
      cache_pg(n0->opt, reset);
#else
      cache_pg(n0->opt);
#endif
      return n0->opt;
    }
    if (FINITE(current_node->val)) {
      stats.expand_node(*(current_node->state));
      ALNewState ns(*this, current_node);
      current_node->state->expand(ns, POS_INF);
      current_node->exp += 1;
    }
    current_node->closed = true;
  }
  n0->est = POS_INF;
#ifdef USE_RESET_LIST
  cache_pg(n0->opt, reset);
#else
  cache_pg(POS_INF);
#endif
  return POS_INF;
}

NTYPE AbstractionLevel::base_level_search(ActionSequence& plan)
{
  NTYPE f_max = 0;
  while (!open.empty()) {
    if (stats.break_signal_raised()) {
      std::cerr << "level " << level <<": highest f-value = "
		<< f_max << std::endl;
      return f_max;
    }
    Node* current_node = open.dequeue();
    if (current_node->val > f_max) {
      f_max = current_node->val;
      if (verbose_level > 1) {
	std::cerr << "f = " << f_max << " (" << stats << ")" << std::endl;
      }
    }
    // if the heuristic is consistent, a node should only ever be
    // expanded once, but we can't guarantee that when using HA*
    //assert(current_node->exp == 0);
    if (FINITE(current_node->val)) {
      stats.expand_node(*(current_node->state));
      ALNewState ns(*this, current_node);
      current_node->state->expand(ns, POS_INF);
      current_node->exp += 1;
    }
    current_node->closed = true;
    if (current_node->state->is_final()) {
      plan.clear();
      //std::cerr << "solution found:" << std::endl;
      //current_node->write_back_path(std::cerr);
      current_node->bp_trans->insert_path(plan);
      if (reverse()) {
	for (index_type i = 0; i < plan.size()/2; i++) {
	  assert(i < plan.size() - (i + 1));
	  index_type a = plan[plan.size() - (i + 1)];
	  plan[plan.size() - (i + 1)] = plan[i];
	  plan[i] = a;
	}
      }
      current_node->cache_optimal_path(0);
#ifdef USE_RESET_LIST
      if ((solver == fwd_ha) || (solver == bwd_ha))
	cache_pg(current_node->acc, reset);
      else if ((solver == fwd_ha_light) || (solver == bwd_ha_light))
	cache_pg(current_node->acc);
#else
      if ((solver == fwd_ha) || (solver == bwd_ha) ||
	  (solver == fwd_ha_light) || (solver == bwd_ha_light))
	cache_pg(current_node->acc);
#endif
      return current_node->acc;
    }
  }
  return POS_INF;
}

void AbstractionLevel::init_open()
{
  ALNewState ns(*this, 0);
  if (reverse()) {
    PartialStateEnumerator e(absins.signature, absins.goal_state);
    bool more = e.first();
    while (more) {
      ALState s0(*this, e.current_state());
      ns.new_state(s0, POS_INF);
      more = e.next();
    }
  }
  else {
    ALState s0(*this, absins.init_state);
    ns.new_state(s0, POS_INF);
  }
  if (reverse() && (verbose_level > 0)) {
    std::cerr << "level " << level << ": open initialised with "
	      << open.size() << " nodes" << std::endl;
  }
}

void AbstractionLevel::dealloc()
{
  open.clear();
  nodes.clear();
}

NTYPE AbstractionLevel::solve_blind
(ActionSequenceSet& plans)
{
  SASCostACF cost(absins);
  SASHeuristic h;
  SASSeqProgState s0(absins, cost, h, absins.init_state);
  Result search_res(&plans);
  search_res.set_stop_condition(Result::stop_at_first);
  Statistics my_stats(&stats);
  BFS search(my_stats, search_res, 1000007);
  NTYPE v = search.start(s0);
  return v;
}

// a slightly refined version of the basic first-fail check
bool AbstractionLevel::check
(const ActionSequence& ap,
 const partial_state state,
 index_type step,
 ActionSequence& plan,
 index_set_vec& new_dvs)
{
  if (step < ap.size()) {
    assert(ap[step] < absins.n_actions());
    index_set cands;
    mapping::inverse_map_image(absins.action_reduce_map, ap[step], cands);
    //std::cerr << "step " << step << ": abs = " << ap[step]
    //      << ": cands = " << cands << std::endl;
    index_set failed;
    for (index_type k = 0; k < cands.size(); k++) {
      assert(cands[k] < baseins.n_actions());
      assert(baseins.actions[cands[k]].cost ==
	     absins.actions[ap[step]].cost);
      bool ok = true;
      if (state.implies(baseins.actions[cands[k]].pre) &&
	  state.implies(baseins.actions[cands[k]].prv)) {
	partial_state s_next(state);
	s_next.assign(baseins.actions[cands[k]].post);
	bool ok = check(ap, s_next, step + 1, plan, new_dvs);
	if (ok) {
	  ((index_vec&)plan).insert(cands[k], 0);
	  return true;
	}
      }
      else {
	failed.insert(cands[k]);
      }
    }
    lvector<pair_set> fcs;
    for (index_type k = 0; k < failed.size(); k++) {
      pair_set fc;
      for (index_type i = 0; i < baseins.actions[failed[k]].pre.size(); i++)
	if (!state.contains(baseins.actions[failed[k]].pre[i])) {
	  index_type var = baseins.actions[failed[k]].pre[i].first;
	  index_type val = baseins.actions[failed[k]].pre[i].second;
	  assert(!dvs[var].contains(val));
	  fc.insert(baseins.actions[failed[k]].pre[i]);
	}
      for (index_type i = 0; i < baseins.actions[failed[k]].prv.size(); i++)
	if (!state.contains(baseins.actions[failed[k]].prv[i])) {
	  index_type var = baseins.actions[failed[k]].prv[i].first;
	  index_type val = baseins.actions[failed[k]].prv[i].second;
	  assert(!dvs[var].contains(val));
	  fc.insert(baseins.actions[failed[k]].prv[i]);
	}
      assert(!fc.empty());
      bool found = false;
      for (index_type i = 0; (i < fc.size()) && !found; i++)
	if (new_dvs[fc[i].first].contains(fc[i].second))
	  found = true;
      if (!found) fcs.append(fc);
      bool_vec rem_to_hit(true, fcs.size());
      for (index_type k = 0; k < fcs.size(); k++)
	if (rem_to_hit[k]) {
	  assert(!fcs[k].empty());
	  index_type i_best = no_such_index;
	  index_type n_best = 0;
	  for (index_type i = 0; i < fcs[k].size(); i++) {
	    index_type n = 1;
	    for (index_type j = k + 1; j < fcs.size(); j++)
	      if (rem_to_hit[j] && fcs[j].contains(fcs[k][i])) n += 1;
	    if (n > n_best) {
	      i_best = i;
	      n_best = n;
	    }
	  }
	  assert((i_best != no_such_index) && (i_best < fcs[k].size()));
	  index_type var = fcs[k][i_best].first;
	  index_type val = fcs[k][i_best].second;
	  assert(!dvs[var].contains(val));
	  new_dvs[var].insert(val);
	  rem_to_hit[k] = false;
	  for (index_type j = k + 1; j < fcs.size(); j++)
	    if (rem_to_hit[j] && fcs[j].contains(fcs[k][i_best]))
	      rem_to_hit[j] = false;
	}
    }
    return false;
  }
  else if (state.implies(baseins.goal_state)) {
    plan.clear();
    return true;
  }
  for (index_type i = 0; i < baseins.goal_state.size(); i++)
    if (!state.contains(baseins.goal_state[i])) {
      index_type var = baseins.goal_state[i].first;
      index_type val = baseins.goal_state[i].second;
      assert(!dvs[var].contains(val));
      new_dvs[var].insert(val);
      return false;
    }
  assert(0);
  return false;
}

bool AbstractionLevel::check
(const ActionSequence& ap,
 ActionSequence& plan,
 index_set_vec& new_dvs)
{
  if (verbose_level > 0) {
    std::cerr << "level " << level << ": checking abstract plan of "
	      << ap.size() << " actions..." << std::endl;
  }
  new_dvs.assign_value(EMPTYSET, baseins.n_variables());
  return check(ap, baseins.init_state, 0, plan, new_dvs);
}

// bool AbstractionLevel::check2
// (const ActionSequence& ap,
//  const partial_state state,
//  index_type step,
//  ActionSequence& plan,
//  lvector<pair_set>& cs)
// {
//   if (step < ap.size()) {
//     assert(ap[step] < absins.n_actions());
//     index_set cands;
//     mapping::inverse_map_image(absins.action_reduce_map, ap[step], cands);
//     //std::cerr << "step " << step << ": abs = " << ap[step]
//     //      << ": cands = " << cands << std::endl;
//     assert(cs.size() > 0);
//     pair_set prev_cs(cs[cs.size() - 1]);
//     for (index_type k = 0; k < cands.size(); k++) {
//       assert(cands[k] < baseins.n_actions());
//       assert(baseins.actions[cands[k]].cost ==
// 	     absins.actions[ap[step]].cost);
//       bool ok_now = true;
//       for (index_type i = 0; i < baseins.actions[cands[k]].pre.size(); i++)
// 	if (!state.contains(baseins.actions[cands[k]].pre[i])) {
// 	  cs[cs.size() - 1].insert(baseins.actions[cands[k]].pre[i]);
// 	  ok_now = false;
// 	}
//       for (index_type i = 0; i < baseins.actions[cands[k]].prv.size(); i++)
// 	if (!state.contains(baseins.actions[cands[k]].prv[i])) {
// 	  cs[cs.size() - 1].insert(baseins.actions[cands[k]].prv[i]);
// 	  ok_now = false;
// 	}
//       partial_state s_next(state);
//       s_next.assign(baseins.actions[cands[k]].post);
//       bool ok_rest = check2(ap, s_next, step + 1, plan, cs);
//       if (ok_now && ok_rest) {
// 	((index_vec&)plan).insert(cands[k], 0);
// 	return true;
//       }
//       if ((k + 1) < cands.size()) {
// 	cs.append(prev_cs);
//       }
//     }
//     return false;
//   }
//   else if (state.implies(baseins.goal_state)) {
//     plan.clear();
//     return true;
//   }
//   return false;
// }

// bool AbstractionLevel::check2
// (const ActionSequence& ap,
//  ActionSequence& plan,
//  index_set_vec& new_dvs)
// {
//   if (verbose_level > 0) {
//     std::cerr << "level " << level << ": checking abstract plan of "
// 	      << ap.size() << " actions..." << std::endl;
//   }
//   pair_set empty_cs;
//   lvector<pair_set> cs(empty_cs, 1);
//   bool ok = check2(ap, baseins.init_state, 0, plan, cs);
//   new_dvs.assign_value(EMPTYSET, baseins.n_variables());
//   if (ok) return true;
//   if (verbose_level > 0) {
//     std::cerr << "level " << level << ": examined " << cs.size()
// 	      << " executions" << std::endl;
//   }
//   bool_vec rem_to_hit(true, cs.size());
//   for (index_type k = 0; k < cs.size(); k++)
//     if (rem_to_hit[k]) {
//       assert(!cs[k].empty());
//       index_type i_best = no_such_index;
//       index_type n_best = 0;
//       for (index_type i = 0; i < cs[k].size(); i++) {
// 	index_type n = 1;
// 	for (index_type j = k + 1; j < cs.size(); j++)
// 	  if (rem_to_hit[j] && cs[j].contains(cs[k][i])) n += 1;
// 	if (n > n_best) {
// 	  i_best = i;
// 	  n_best = n;
// 	}
//       }
//       assert((i_best != no_such_index) && (i_best < cs[k].size()));
//       index_type var = cs[k][i_best].first;
//       index_type val = cs[k][i_best].second;
//       assert(!dvs[var].contains(val));
//       new_dvs[var].insert(val);
//       rem_to_hit[k] = false;
//       for (index_type j = k + 1; j < cs.size(); j++)
// 	if (rem_to_hit[j] && cs[j].contains(cs[k][i_best]))
// 	  rem_to_hit[j] = false;
//     }
//   return false;
// }

NTYPE AbstractionLevel::cegar(ActionSequence& plan, NTYPE ub, index_type db)
{
  stats.start();
  ActionSequence abs_plan;
  if (verbose_level > 0)
    std::cerr << "level " << level << ": searching (" << solver
	      << ", " << reverse() << ")..." << std::endl;
  init_open();
  NTYPE v = base_level_search(abs_plan);
  if (stats.break_signal_raised()) {
    stats.stop();
    return v;
  }
  std::cerr << "level " << level << ": cost = " << PRINT_NTYPE(v)
	    << " (" << stats.total_nodes() << " nodes, "
	    << stats.total_time() << " sec., "
	    << stats.peak_memory() << "k)" << std::endl;
  if (verbose_level > 1) {
    std::cerr << "abstract plan:" << std::endl;
    PrintSASActions printer(absins, std::cerr);
    abs_plan.output(printer);
    std::cerr << std::endl;
  }
  if (INFINITE(v) || (v >= ub)) {
    stats.stop();
    return v;
  }
#ifdef VALIDATE_ABSTRACT_PLAN
  NTYPE vv = validate_sas_plan(absins, abs_plan);
  assert(vv == v);
#endif
  index_set_vec new_dvs;
  bool ok = check(abs_plan, plan, new_dvs);
  if (ok) {
    if (verbose_level > 0) {
      std::cerr << "level " << level << ": abstract plan works!" << std::endl;
    }
#ifdef VALIDATE_FINAL_PLAN
    NTYPE vv = validate_sas_plan(baseins, plan);
    if (vv != v) {
      std::cerr << "vv = " << vv << ", v = " << v << std::endl;
    }
    assert(vv == v);
#endif
    stats.stop();
    return v;
  }
  else if (level < db) {
    index_type n_new = 0;
    for (index_type k = 0; k < baseins.n_variables(); k++) {
      assert(!new_dvs[k].have_common_element(dvs[k]));
      n_new += new_dvs[k].size();
    }
    if (verbose_level > 0) {
      std::cerr << "level " << level << ": abstract plan failed, "
		<< n_new << " new distinguished values" << std::endl;
      if (verbose_level > 1) {
	for (index_type k = 0; k < baseins.n_variables(); k++)
	  for (index_type i = 0; i < new_dvs[k].size(); i++)
	    std::cerr << baseins.variables[k].name << ":"
		      << baseins.variables[k].domain[new_dvs[k][i]]
		      << std::endl;
      }
    }
    assert(n_new > 0);
    value_map refined_vm(vm);
    for (index_type k = 0; k < baseins.n_variables(); k++) {
      if (new_dvs[k].size() > 0) {
	index_type v0 = baseins.init_state.value_of(k);
	for (index_type i = 0; i < new_dvs[k].size(); i++) {
	  equivalence eq;
	  make_partition(baseins, k, dtgs[k], v0, new_dvs[k][i], eq);
	  refine_value_map(eq, refined_vm[k]);
	}
      }
      new_dvs[k].insert(dvs[k]);
    }
    if (verbose_level > 2) {
      std::cerr << "next level dvs: " << new_dvs << std::endl;
      std::cerr << "next level map: " << refined_vm << std::endl;
    }
    value_map internal_vm;
    make_internal_map(vm, refined_vm, internal_vm);
    AbstractionLevel* next =
      new AbstractionLevel(level + 1, baseins, dtgs, new_dvs, refined_vm,
			   internal_vm, this, stats, verbose_level);
    if ((solver == fwd_blind) || (solver == bwd_blind))
      dealloc();
    if ((solver == fwd_ha_light) || (solver == bwd_ha_light))
      if (above)
	above->dealloc();
    NTYPE v_next = next->cegar(plan, ub, db);
    v = MAX(v, v_next);
    delete next;
  }
  else {
    std::cerr << "level " << level << ": depth limit reached" << std::endl;
  }
  stats.stop();
  return v;
}

END_HSPS_NAMESPACE

HSPS::index_pair most_expensive_goal
(HSPS::Instance& ins, HSPS::SASInstance& sins, HSPS::Statistics& stats)
{
  HSPS::CostTable h1(ins, stats);
  h1.compute_H1(HSPS::CostACF(ins));
  HSPS::ToSASHAdapter sh1(sins, h1);
  HSPS::index_type g_max = HSPS::no_such_index;
  NTYPE c_max = 0;
  for (HSPS::index_type i = 0; i < sins.goal_state.size(); i++) {
    HSPS::partial_state sgi;
    sgi.insert(sins.goal_state[i]);
    NTYPE c = sh1.eval(sgi);
    if ((c > c_max) || (g_max == HSPS::no_such_index)) {
      g_max = i;
      c_max = c;
    }
  }
  assert((g_max != HSPS::no_such_index) && (g_max < sins.goal_state.size()));
  return sins.goal_state[g_max];
}

int main(int argc, char *argv[]) {
  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  int  verbose_level = 1;
  int  opt_prep = 1;
  bool opt_rel = true;
  int  opt_find = 0;
  int  opt_verify = 0;
  bool opt_sas_sel = false;
  bool opt_unit_cost = false;
  bool opt_print_plan = true;
  int  opt_solver = HSPS::AbstractionLevel::switchback;
  bool opt_init_all = false;
  bool opt_flat = false;
  NTYPE ub = POS_INF;
  HSPS::index_type db = HSPS::index_type_max;

  long          time_limit = 0;
  unsigned long memory_limit = 0;

  HSPS::Statistics  stats;
  stats.enable_interrupt();
  stats.start();
  HSPS::Statistics  read_stats(&stats);
  HSPS::Statistics  prep_stats(&stats);
  HSPS::Statistics  cegar_stats(&stats);

  HSPS::Parser* reader = new HSPS::Parser(symbols);

  for (int k = 1; k < argc; k++) {
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      verbose_level = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-find") == 0) && (k < argc - 1)) {
      opt_find = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-verify") == 0) && (k < argc - 1)) {
      opt_verify = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-prep") == 0) && (k < argc - 1)) {
      opt_prep = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-no-rel") == 0) {
      opt_rel = false;
    }
    else if (strcmp(argv[k],"-sas-sel") == 0) {
      opt_sas_sel = true;
    }
    else if (strcmp(argv[k],"-u") == 0) {
      opt_unit_cost = true;
    }
    else if (strcmp(argv[k],"-p") == 0) {
      opt_print_plan = false;
    }
    else if ((strcmp(argv[k],"-b") == 0) && (k < argc - 1)) {
      ub = A_TO_N(argv[++k]);
    }
    else if ((strcmp(argv[k],"-d") == 0) && (k < argc - 1)) {
      db = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-s") == 0) && (k < argc - 1)) {
      opt_solver = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-a") == 0) {
      opt_init_all = true;
    }
    else if (strcmp(argv[k],"-F") == 0) {
      opt_flat = true;
    }
    else if ((strcmp(argv[k],"-t") == 0) && (k < argc - 1)) {
      time_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-y") == 0) && (k < argc - 1)) {
      memory_limit = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-dba-semantics") == 0) {
      HSPS::PDDL_Base::del_before_add_semantics = true;
    }
    else if (*argv[k] != '-') {
      read_stats.start();
      reader->read(argv[k], false);
      read_stats.stop();
    }
  }

  reader->post_process();

  if (time_limit > 0) stats.enable_time_out(time_limit);
  if (memory_limit > 0) stats.enable_memory_limit(memory_limit);

  prep_stats.start();
  if (opt_find == 1) {
    std::cerr << "searching for invariants..." << std::endl;
    reader->find_cc();
  }
  std::cerr << reader->dom_sc_invariants.size() << " invariants" << std::endl;

  std::cerr << "instantiating..." << std::endl;
  HSPS::PDDL_Base::name_instance_by_problem_file = true;
  HSPS::Instance instance;
  reader->instantiate(instance);
  if (stats.break_signal_raised()) return 0;

  HSPS::Preprocessor prep(instance, prep_stats);
  if (opt_prep > 0) {
    std::cerr << "preprocessing..." << std::endl;
    prep.preprocess(opt_prep > 1);
    if (stats.break_signal_raised()) return 0;
    if (opt_rel) {
      prep.compute_irrelevant_atoms();
      prep.remove_irrelevant_atoms();
      if (stats.break_signal_raised()) return 0;
      if (!instance.cross_referenced()) {
	std::cerr << "re-cross referencing..." << std::endl;
	instance.cross_reference();
      }
    }
  }
  else {
    instance.cross_reference();
  }
  if (stats.break_signal_raised()) return 0;

  if (opt_find == 2) {
    std::cerr << "searching for invariants..." << std::endl;
    HSPS::graph* g_inc = prep.inconsistency_graph();
    prep.find_inconsistent_set_invariants(*g_inc);
    instance.add_missing_negation_invariants();
  }
  else if (opt_find == 3) {
    std::cerr << "searching for invariants..." << std::endl;
    prep.bfs_find_invariants();
  }
  if (stats.break_signal_raised()) return 0;
  if (opt_verify == 2) {
    prep.verify_invariants(*(prep.inconsistency()));
    prep.remove_unverified_invariants();
  }
  else if (opt_verify == 1) {
    instance.verify_invariants();
    prep.remove_unverified_invariants();
  }
  if (stats.break_signal_raised()) return 0;

  std::cerr << instance.n_invariants() << " invariants" << std::endl;
  if (verbose_level > 1) {
    for (HSPS::index_type k = 0; k < instance.n_invariants(); k++)
      instance.print_invariant(std::cerr, instance.invariants[k]);
  }
  std::cerr << "constructing SAS instance..." << std::endl;
  HSPS::SASInstance sasp(instance, opt_sas_sel, false, true);
  if (stats.break_signal_raised()) return 0;
  std::cerr << "cross-referencing..." << std::endl;
  sasp.cross_reference();
  if (stats.break_signal_raised()) return 0;
  prep_stats.stop();
  std::cerr << "SAS instance " << instance.name << " built in "
	    << prep_stats.total_time() << " seconds" << std::endl;
  std::cout << ";; "<< sasp.name << ": "
	    << sasp.n_variables() << " variables, "
	    << sasp.n_actions() << " actions, "
	    << std::endl;

  if (verbose_level > 0)
    std::cerr << "base level signature: " << sasp.signature
	      << " (" << sasp.n_actions() << " actions)"
	      << std::endl;
  if (verbose_level > 2) {
    print_sas_short(sasp, std::cerr);
  }
  if (db == 0) exit(0);

  HSPS::index_type n_dvs = 0;
  HSPS::index_set_vec dvs(HSPS::EMPTYSET, sasp.n_variables());

#ifdef USE_DVS
  if (opt_init_all) {
    for (HSPS::index_type k = 0; k < sasp.n_variables(); k++) {
      assert(sasp.init_state.defines(k));
      if (sasp.goal_state.defines(k)) {
    	dvs[k].insert(sasp.init_state.value_of(k));
    	dvs[k].insert(sasp.goal_state.value_of(k));
    	n_dvs += dvs[k].size();
      }
    }
  }
  else {
    HSPS::index_pair g_max = most_expensive_goal(instance, sasp, prep_stats);
    HSPS::index_type var = g_max.first;
    HSPS::index_type val_g = g_max.second;
    HSPS::index_type val_i = sasp.init_state.value_of(var);
    dvs[var].insert(val_g);
    dvs[var].insert(val_i);
  }

  if (verbose_level > 0) {
    std::cerr << n_dvs << " initial distinguished values" << std::endl;
    if (verbose_level > 1) {
      for (HSPS::index_type k = 0; k < sasp.n_variables(); k++)
	for (HSPS::index_type i = 0; i < dvs[k].size(); i++)
	  std::cerr << sasp.variables[k].name << ":"
		    << sasp.variables[k].domain[dvs[k][i]]
		    << std::endl;
    }
    if (verbose_level > 2)
      std::cerr << "initial dvs: " << dvs << std::endl;
  }

  HSPS::value_map initial_vm(HSPS::mapping(), sasp.n_variables());
  make_value_map_from_dvs(sasp, dvs, initial_vm);
  if (stats.break_signal_raised()) return 0;
  if (verbose_level > 2)
    std::cerr << "initial map: " << initial_vm << std::endl;

#else

  HSPS::lvector<HSPS::weighted_graph>
    dtgs(HSPS::weighted_graph(), sasp.n_variables());
  for (HSPS::index_type k = 0; k < sasp.n_variables(); k++)
    HSPS::compute_dtg(sasp, k, dtgs[k]);

  HSPS::value_map initial_vm(HSPS::mapping(), sasp.n_variables());

  if (opt_flat) {
    for (HSPS::index_type k = 0; k < sasp.n_variables(); k++)
      initial_vm[k].assign_identity(sasp.variables[k].n_values());
  }
  else if (opt_init_all) {
    for (HSPS::index_type k = 0; k < sasp.n_variables(); k++) {
      if (sasp.goal_state.defines(k)) {
	HSPS::index_type v0 = sasp.init_state.value_of(k);
	assert(v0 != HSPS::no_such_index);
	HSPS::index_type vg = sasp.goal_state.value_of(k);
	assert(vg != HSPS::no_such_index);
	HSPS::equivalence deq;
	HSPS::make_partition(sasp, k, dtgs[k], v0, vg, deq);
	//std::cerr << sasp.variables[k].name
	//	  << ": " << v0 << "." << sasp.variables[k].domain[v0]
	//	  << " -> " << vg << "." << sasp.variables[k].domain[vg]
	//	  << ", deq = " << deq << std::endl;
	deq.make_map(initial_vm[k]);
	dvs[k].insert(v0);
	dvs[k].insert(vg);
      }
      else {
	HSPS::index_type v0 = sasp.init_state.value_of(k);
	assert(v0 != HSPS::no_such_index);
	initial_vm[k].assign_value(v0, sasp.variables[k].n_values());
      }
    }
  }
  else {
    HSPS::index_pair g_max = most_expensive_goal(instance, sasp, prep_stats);
    HSPS::index_type var = g_max.first;
    HSPS::index_type val_g = g_max.second;
    HSPS::index_type val_i = sasp.init_state.value_of(var);
    HSPS::equivalence deq;
    HSPS::make_partition(sasp, var, dtgs[var], val_i, val_g, deq);
    deq.make_map(initial_vm[var]);
    dvs[var].insert(val_g);
    dvs[var].insert(val_i);
    for (HSPS::index_type k = 0; k < sasp.n_variables(); k++)
      if (k != var) {
	HSPS::index_type v0 = sasp.init_state.value_of(k);
	assert(v0 != HSPS::no_such_index);
	initial_vm[k].assign_value(v0, sasp.variables[k].n_values());
      }
  }

  if (stats.break_signal_raised()) return 0;
  if (verbose_level > 2)
    std::cerr << "initial map: " << initial_vm << std::endl;

#endif

  HSPS::AbstractionLevel* l0 =
    new HSPS::AbstractionLevel
    (0, (HSPS::AbstractionLevel::solver_type)opt_solver,
     sasp, dtgs, dvs, initial_vm, cegar_stats, verbose_level);
  if (stats.break_signal_raised()) return 0;
  HSPS::ActionSequence plan;
  NTYPE v_max = l0->cegar(plan, ub, db);

  std::cout << ";; highest lower bound = " << PRINT_NTYPE(v_max) << std::endl;
  std::cout << ";; nodes expanded = " << cegar_stats.total_nodes() << std::endl;
  std::cout << ";; branching factor: forward = "
	    << (HSPS::AbstractionLevel::fwd_expansions > 0 ?
		(HSPS::AbstractionLevel::fwd_successors/
		 (double)HSPS::AbstractionLevel::fwd_expansions) : 0)
	    << ", backward = "
	    << (HSPS::AbstractionLevel::bwd_expansions > 0 ?
		(HSPS::AbstractionLevel::bwd_successors/
		 (double)HSPS::AbstractionLevel::bwd_expansions) : 0)
	    << std::endl;
  std::cout << ";; time (reading input) = " << read_stats.total_time()
	    << std::endl;
  std::cout << ";; time (preprocessing) = " << prep_stats.total_time()
	    << std::endl;
  std::cout << ";; time (CEGAR) = " << cegar_stats.total_time()
	    << std::endl;
  std::cout << ";; time (total) = " << stats.total_time()
	    << std::endl;
  std::cout << ";; peak memory = " << stats.peak_memory() << "k" << std::endl;
  if (FINITE(v_max) && !stats.break_signal_raised() && opt_print_plan) {
    HSPS::PrintSASActions printer(sasp, std::cout);
    plan.output(printer);
    std::cout << std::endl;
  }

  return 0;
}
