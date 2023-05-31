
#include "soft.h"

#include "seq_reg.h"
#include "ida.h"

BEGIN_HSPS_NAMESPACE

SoftInstance::SoftInstance()
  : soft(SoftGoal(), 0), null_value(0)
{
  // done
}

SoftInstance::SoftInstance(const Name* n)
  : Instance(n), soft(SoftGoal(), 0), null_value(0)
{
  // done
}

SoftInstance::SoftGoal& SoftInstance::new_soft_goal()
{
  soft.inc_length();
  return soft[soft.length() - 1];
}

bool SoftInstance::empty_plan_valid()
{
  for (index_type k = 0; k < hard.length(); k++)
    if (!atoms[hard[k]].init) return false;
  return true;
}

NTYPE SoftInstance::empty_plan_value()
{
  if (trace_level > 1) {
    std::cerr << "computing value of empty plan..." << std::endl;
  }
  if (!empty_plan_valid()) {
    if (trace_level > 1) {
      std::cerr << "empty plan not valid, value = -INF" << std::endl;
    }
    return NEG_INF;
  }
  NTYPE v = null_value;
  for (index_type k = 0; k < soft.length(); k++)
    if (soft[k].is_sat_init(*this)) {
      if (trace_level > 1) {
	std::cerr << "soft goal " << soft[k].name << ": ";
	write_atom_set(std::cerr, soft[k].atoms);
	std::cerr << " is satisfied in init state, value +"
		  << soft[k].weight << std::endl;
      }
      v += soft[k].weight;
    }
  if (trace_level > 1) {
    std::cerr << "empty plan value = " << v << std::endl;
  }
  return v;
}

void SoftInstance::remap_hard_goals(const index_vec& atom_map)
{
  remap_set(hard, atom_map);
}

void SoftInstance::remap_soft_goals(const index_vec& atom_map)
{
  for (index_type k = 0; k < n_soft(); k++) {
    remap_set(soft[k].atoms, atom_map);
  }
}

long SoftInstance::integrify_weights()
{
#ifdef NTYPE_RATIONAL
  long common_m = 1;
  for (index_type k = 0; k < n_soft(); k++) {
    common_m = lcm(common_m, soft[k].weight.divisor());
  }

  std::cerr << "common lcm = " << common_m << std::endl;

  for (index_type k = 0; k < n_soft(); k++) {
    NTYPE wi(soft[k].weight.numerator() * (common_m / soft[k].weight.divisor()), 1);
    std::cerr << "replacing " << soft[k].weight << " by " << wi << std::endl;
    soft[k].weight = wi;
  }
  return common_m;
#else
  std::cerr << "error: integrify not implemented for NTYPE_FLOAT" << std::endl;
  exit(255);
#endif
}

void SoftInstance::compile_KG(Instance& ins)
{
  std::cerr << "compiling soft goals with Keyder & Geffner method..."
	    << std::endl;
  ins.copy(*this);
  Atom& a_mod = ins.new_atom(new StringName("normal-mode"));
  a_mod.init = true;
  a_mod.goal = false;
  for (index_type k = 0; k < ins.n_actions(); k++)
    ins.actions[k].pre.insert(a_mod.index);
  NTYPE residual_value = 0;
  index_set g(hard);
  index_type prev = no_such_index;
  for (index_type k = 0; k < n_soft(); k++) {
    const Name* n =
      (soft[k].name ? soft[k].name : new EnumName("preference", k));
    Atom& a = ins.new_atom(n);
    a.init = false;
    a.goal = true;
    g.insert(a.index);
    Action& a_sat = ins.new_action(new ModName(n, "sat"));
    a_sat.pre.insert(soft[k].atoms);
    if (prev != no_such_index)
      a_sat.pre.insert(prev);
    a_sat.add.insert(a.index);
    a_sat.del.insert(a_mod.index);
    a_sat.cost = 0;
    Action& a_unsat = ins.new_action(new ModName(n, "unsat"));
    if (prev != no_such_index)
      a_unsat.pre.insert(prev);
    a_unsat.add.insert(a.index);
    a_unsat.del.insert(a_mod.index);
    a_unsat.cost = soft[k].weight;
    residual_value += (-1 * soft[k].weight);
    prev = a.index;
  }
  residual_value -= null_value;
  if (residual_value != 0) {
    std::cerr << "warning: compiled problem metric differs by "
	      << residual_value
	      << std::endl;
  }
  ins.set_goal(g);
}

bool SoftInstance::compile_direct_cost(Instance& ins)
{
  std::cerr << "compiling soft goals with direct action costs..." << std::endl;
  bool_vec is_dec(false, n_soft());

  // check applicability and find missing negations
  index_set missing_neg;
  for (index_type k = 0; k < n_soft(); k++) {
    if (soft[k].atoms.size() != 1) {
      std::cerr << "warning: cannot do direct cost compilation for soft goal ";
      write_atom_set(std::cerr, soft[k].atoms);
      std::cerr << std::endl;
      return false;
    }
    index_type p = soft[k].atoms[0];
    if (atoms[p].add_by.empty()) {
      is_dec[k] = true;
    }
    else if (atoms[p].del_by.empty()) {
      is_dec[k] = false;
    }
    else {
      std::cerr << "warning: cannot do direct cost compilation for soft goal ";
      write_atom_set(std::cerr, soft[k].atoms);
      std::cerr << " (it is non-monotonic)" << std::endl;
      return false;
    }
    if (atoms[p].neg == no_such_index)
      missing_neg.insert(p);
  }

  ins.copy(*this);
  std::cerr << missing_neg.size() << " missing negtions" << std::endl;
  ins.complete_atom_negations(missing_neg);
  if (!ins.cross_referenced())
    ins.cross_reference();
  index_set g(hard);

  for (index_type k = 0; k < n_soft(); k++) {
    index_type p = soft[k].atoms[0];
    index_type not_p = ins.atoms[p].neg;
    assert(not_p != no_such_index);
    std::cerr << "compiling soft goal #" << k << ": " << ins.atoms[p].name
	      << std::endl;
    // goal can only be destroyed; first action to do it gets the penalty
    if (is_dec[k]) {
      index_type n_act = ins.n_actions();
      for (index_type i = 0; i < n_act; i++)
	if (ins.actions[i].del.contains(p)) {
	  Action& a2 = ins.copy_action(i);
	  Action& a1 = ins.actions[i];
	  a1.pre.insert(not_p);
	  a2.pre.insert(p);
	  a2.cost += soft[k].weight;
	}
//       for (index_type a = 0; a < atoms[p].del_by.size(); a++) {
// 	Action& a2 = ins.copy_action(atoms[p].del_by[a]);
// 	Action& a1 = ins.actions[atoms[p].del_by[a]];
// 	a1.pre.insert(not_p);
// 	a2.pre.insert(p);
// 	a2.cost += soft[k].weight;
//       }
    }
    // goal cannot be destroyed; add a "fake" option and make it a hard goal
    else {
      const Name* n = atoms[p].name;
      Action& a_fake = ins.new_action(new ModName(n, "fake"));
      a_fake.pre.insert(not_p);
      a_fake.add.insert(p);
      a_fake.del.insert(not_p);
      a_fake.cost = soft[k].weight;
      g.insert(p);
    }
    Stopwatch::seconds();
    std::cerr << "now " << ins.n_actions() << " actions, "
	      << Stopwatch::peak_memory() << "k" << std::endl;
  }
  ins.set_goal(g);
  return true;
}

void SoftInstance::create_decision_problem
(const bool_vec& sel, NTYPE nb_min, Instance& ins)
{
  index_set s(hard);
  NTYPE v = null_value;
  for (index_type k = 0; k < n_soft(); k++) if (sel[k]) {
    s.insert(soft[k].atoms);
    v += soft[k].weight;
  }
  ins.copy(*this);
  ins.set_goal(s);
  ins.set_cost_bound(v - nb_min);
}

void SoftInstance::create_decision_problem
(const bool_vec& sel, Instance& ins)
{
  index_set s(hard);
  for (index_type k = 0; k < n_soft(); k++) if (sel[k]) {
    s.insert(soft[k].atoms);
  }
  ins.copy(*this);
  ins.set_goal(s);
}

NTYPE SoftInstance::compute_epsilon()
{
#ifdef NTYPE_RATIONAL
  bool first = true;
  NTYPE e;
  for (index_type k = 0; k < n_actions(); k++) if (actions[k].cost != 0) {
    if (first) {
      e = actions[k].cost;
      first = false;
    }
    else {
      e = rational::rgcd(e, actions[k].cost);
    }
  }
  for (index_type k = 0; k < n_soft(); k++) if (soft[k].weight != 0) {
    if (first) {
      e = soft[k].weight;
      first = false;
    }
    else {
      e = rational::rgcd(e, soft[k].weight);
    }
  }
  return e;
#else
  std::cerr << "error: epsilon calculation not implemented for NTYPE_FLOAT"
	    << std::endl;
  exit(255);
#endif
}

NTYPE SoftInstance::eval_goal_state(const index_set& s)
{
  NTYPE v = null_value;
  for (index_type k = 0; k < n_soft(); k++)
    if (s.contains(soft[k].atoms)) v += soft[k].weight;
  return v;
}

NTYPE SoftInstance::eval_goal_state(const index_set& s, index_set& g)
{
  NTYPE v = null_value;
  g.clear();
  for (index_type k = 0; k < n_soft(); k++)
    if (s.contains(soft[k].atoms)) {
      v += soft[k].weight;
      g.insert(k);
    }
  return v;
}

NTYPE SoftInstance::eval_plan(Schedule& s)
{
  if (trace_level > 1) {
    std::cerr << "evaluating plan";
    if (s.plan_name()) std::cerr << " " << s.plan_name();
    std::cerr << "..." << std::endl;
  }
  index_set f;
  bool ok = s.simulate(f, 0);
  if (!ok) {
    if (trace_level > 1) {
      std::cerr << " - plan is not executable" << std::endl;
    }
    return NEG_INF;
  }
  if (trace_level > 1) {
    std::cerr << " - achieved atoms: ";
    write_atom_set(std::cerr, f);
    std::cerr << std::endl;
  }
  if (!f.contains(hard)) {
    if (trace_level > 1) {
      std::cerr << " - hard goals not satisfied" << std::endl;
    }
    return NEG_INF;
  }
  NTYPE v = eval_goal_state(f);
  if (trace_level > 1) {
    std::cerr << " - value is " << v << std::endl;
  }
  return v;
}

void SoftInstance::eval_plan_set(ScheduleSet& s, cost_vec& v)
{
  v.assign_value(0, s.length());
  for (index_type k = 0; k < s.length(); k++) {
    assert(s[k]);
    v[k] = eval_plan(*s[k]);
  }
}

void SoftInstance::write_problem_goal(std::ostream& s) const
{
  if (write_PDDL3) {
    if ((n_hard() + n_soft()) > 0) {
      s << " (:goal";
      if ((n_hard() + n_soft()) > 1) s << " (and";
      for (index_type k = 0; k < n_hard(); k++) {
	s << " (";
	atoms[hard[k]].name->write(s, Name::NC_INSTANCE);
	s << ")";
      }
      for (index_type k = 0; k < n_soft(); k++) {
	s << " (preference";
	if (soft[k].name) {
	  s << " ";
	  soft[k].name->write(s, Name::NC_INSTANCE);
	}
	if (soft[k].atoms.length() > 1) s << " (and";
	for (index_type i = 0; i < soft[k].atoms.length(); i++) {
	  s << " (";
	  atoms[soft[k].atoms[i]].name->write(s, Name::NC_INSTANCE);
	  s << ")";
	}
	if (soft[k].atoms.length() > 1) s << ")";
	s << ")";
      }
      if ((n_hard() + n_soft()) > 1) s << ")";
      s << ")" << std::endl;
    }
  }
  else {
    Instance::write_problem_goal(s);
  }
}

void SoftInstance::write_goal_value_expression(std::ostream& s) const
{
  if (n_soft() > 0) {
    if (null_value != 0) {
      s << "(+ " << PRINT_NTYPE(null_value) << " ";
    }
    for (index_type k = 0; k < n_soft() - 1; k++) {
      s << "(+ (* " << PRINT_NTYPE(soft[k].weight)
	<< " (- 1 (is-violated " << soft[k].name << "))) ";
    }
    s << "(* " << PRINT_NTYPE(soft[n_soft() - 1].weight)
      << " (- 1 (is-violated " << soft[n_soft() - 1].name << ")))";
    for (index_type k = 0; k < n_soft() - 1; k++) {
      s << ")";
    }
    if (null_value != 0) {
      s << ")";
    }
  }
  else if (null_value != 0) {
    s << PRINT_NTYPE(null_value);
  }
  else {
    s << "0";
  }
}

void SoftInstance::write_problem_metric(std::ostream& s) const
{
  if (!write_metric) return;
  if ((max_cost > 0) && (n_soft() > 0)) {
    if (!write_PDDL2 || !write_PDDL3) {
      std::cerr << "warning: can't write correct :metric without PDDL2/PDDL3"
		<< std::endl;
      return;
    }
    s << " (:metric maximize (- ";
    write_goal_value_expression(s);
    s << " (total-cost)))" << std::endl;
  }
  else if ((max_cost == 0) && (n_soft() > 0)) {
    if (!write_PDDL3) {
      std::cerr << "warning: can't write correct :metric without PDDL3"
		<< std::endl;
      return;
    }
    s << " (:metric maximize ";
    write_goal_value_expression(s);
    s << ")" << std::endl;
  }
  else {
    if (!write_PDDL2) {
      std::cerr << "warning: can't write correct :metric without PDDL2"
		<< std::endl;
      return;
    }
    s << " (:metric minimize (total-cost))" << std::endl;
  }
}

void SoftInstance::write_soft_goal_set
(std::ostream& s, const index_set& set) const
{
  s << "{";
  for (index_type k = 0; k < set.length(); k++) {
    assert(set[k] < n_soft());
    if (k > 0) s << ",";
    if (soft[set[k]].name) s << soft[set[k]].name;
    else s << "#" << set[k];
  }
  s << "}";
}

void SoftInstance::write_soft_goal_set
(std::ostream& s, const bool_vec& set) const
{
  s << "{";
  bool first = true;
  for (index_type k = 0; k < n_soft(); k++) if (set[k]) {
    if (!first) s << ",";
    first = false;
    if (soft[k].name) s << soft[k].name;
    else s << "#" << k;
  }
  s << "}";
}

void SoftInstance::print(std::ostream& s) const
{
  Instance::print(s);
  s << "hard goals: " << hard << " = ";
  write_atom_set(s, hard);
  s << std::endl << "soft goals:" << std::endl;
  for (index_type k = 0; k < n_soft(); k++) {
    s << " " << soft[k].name << ": " << soft[k].atoms << " = ";
    write_atom_set(s, soft[k].atoms);
    s << " (weight = " << PRINT_NTYPE(soft[k].weight) << ")" << std::endl;
  }
  s << "metric null value: " << PRINT_NTYPE(null_value) << std::endl;
}

DecisionProblemEnumerator::DecisionProblemEnumerator
(SoftInstance& ins, Heuristic& h, NTYPE b)
  : instance(ins), selected(ins.n_soft()), h_cost(h), nb_min(b)
{
  // done
}

DecisionProblemEnumerator::~DecisionProblemEnumerator()
{
  // done
}

bool DecisionProblemEnumerator::find_next(bool more)
{
  while (more) {
    selected.current_set(g_sel);
    a_sel.clear();
    v_sel = instance.null_value;
    for (index_type k = 0; k < g_sel.length(); k++) {
      a_sel.insert(instance.soft[g_sel[k]].atoms);
      v_sel += instance.soft[g_sel[k]].weight;
    }
    c_sel = h_cost.eval(a_sel);
    nb_sel = v_sel - c_sel;
    if (nb_sel >= nb_min) return true;
    more = selected.next();
  }
  return false;
}

bool DecisionProblemEnumerator::first()
{
  bool more = selected.first();
  return find_next(more);
}

bool DecisionProblemEnumerator::next()
{
  bool more = selected.next();
  return find_next(more);
}

NTYPE DecisionProblemEnumerator::current_value() const
{
  return v_sel;
}

NTYPE DecisionProblemEnumerator::current_min_cost() const
{
  return c_sel;
}

NTYPE DecisionProblemEnumerator::current_max_cost() const
{
  return v_sel - nb_min;
}

NTYPE DecisionProblemEnumerator::current_max_nb() const
{
  return nb_sel;
}

NTYPE DecisionProblemEnumerator::current_min_nb() const
{
  return nb_min;
}

const index_set& DecisionProblemEnumerator::current_soft_goals() const
{
  return g_sel;
}

const index_set& DecisionProblemEnumerator::current_goal_atoms() const
{
  return a_sel;
}

void  DecisionProblemEnumerator::create_decision_problem(Instance& ins)
{
  instance.create_decision_problem(selected.current_set(), nb_min, ins);
}

index_type MaxValueSearch::print_options_max = 3;

MaxValueSearch::MaxValueSearch
(SoftInstance& i, ACF& c, Statistics& s, Result& r)
  : instance(i), cost(c), stats(s), res(r), search(0), lb(NEG_INF),
    solved_flag(false), trace_level(SearchAlgorithm::default_trace_level)
{
  lb = instance.empty_plan_value();
  HashTable* tt = new HashTable(100007);
  search = new IDA(stats, res, tt);
}

MaxValueSearch::MaxValueSearch
(SoftInstance& i, ACF& c, Statistics& s, Result& r,
 index_type tt_size, bool use_cc)
  : instance(i), cost(c), stats(s), res(r), search(0), lb(NEG_INF),
    solved_flag(false), trace_level(SearchAlgorithm::default_trace_level)
{
  lb = instance.empty_plan_value();
  HashTable* tt = (tt_size > 0 ? new HashTable(tt_size) : 0);
  IDA* i_search = new IDA(stats, res, tt);
  i_search->set_cycle_check(use_cc);
  search = i_search;
}

MaxValueSearch::~MaxValueSearch()
{
  // done - note search alg and tt not deleted!
}

void MaxValueSearch::insert_option_in_list(const option& o)
{
  option_p p = options.begin();
  while (p != options.end()) {
    if (p->est_value() < o.est_value()) {
      options.insert(p, o);
      while (p != options.end()) {
	if (p->goals == o.goals) {
	  option_p q = p++;
	  options.erase(q);
	}
	else {
	  p++;
	}
      }
      return;
    }
    else if (p->goals == o.goals) {
      // option with the same atoms and better or equal value exists
      return;
    }
    p++;
  }
  options.insert(p, o);
}

void MaxValueSearch::print_option_list(std::ostream& s)
{
  option_p p = options.begin();
  count_type n = 0;
  while ((p != options.end()) && (n < print_options_max)) {
    s << *(p->root) << ": " << PRINT_NTYPE(p->est_value())
      << " (" << PRINT_NTYPE(p->goal_value) << " - "
      << PRINT_NTYPE(p->est_cost) << ", "
      << p->n_goals() << " goals)"
      << std::endl;
    p++;
    n++;
  }
}

void MaxValueSearch::init_option_list()
{
  SubsetEnumerator se(instance.n_soft());
  bool more = se.first();
  index_set s;
  count_type n = 0;
  while (more && !stats.break_signal_raised()) {
    se.current_set(s);
    option o;
    init_option(s, o);
    if (o.est_value() > lb) {
      insert_option_in_list(o);
    }
    if ((trace_level > 1) && (options.size() >= (n + 1000))) {
      std::cerr << options.size() << " options created (" << stats << ")"
		<< std::endl;
      n = options.size();
    }
    more = se.next();
  }
}

void MaxValueSearch::make_empty_plan()
{
  ZeroHeuristic hz(instance);
  SeqCRegState s(instance, hz, cost, EMPTYSET);
  res.solution(s, 0);
}

NTYPE MaxNetBenefit::explore_next_option()
{
  option_p p = options.begin();
  if (p == options.end()) {
    std::cerr << "error: explore_next_option called with empty list"
	      << std::endl;
    exit(255);
  }
  option o = *p;
  p++;
  options.pop_front();
  NTYPE next_best_value = lb;
  if (p != options.end()) {
    next_best_value = p->est_value();
  }
  if (o.root == 0) {
    std::cerr << "error: NIL root in option" << std::endl;
    exit(255);
  }
  search->set_cost_limit(o.goal_value - next_best_value);
  if (trace_level > 0) {
    std::cerr << "searching " << *(o.root) << " with cost limit "
	      << PRINT_NTYPE(search->get_cost_limit()) << "..."
	      << std::endl;
  }
  NTYPE new_est_cost = search->resume(*o.root, o.est_cost);
  if ((new_est_cost <= o.est_cost) && !search->solved()) {
    std::cerr << "error: search returned less or equal cost ("
	      << PRINT_NTYPE(new_est_cost) << ") with no solution"
	      << std::endl;
    stats.error();
  }
  o.est_cost = new_est_cost;
  if (search->solved()) {
    std::cerr << "solution found (value "
	      << PRINT_NTYPE(o.est_value())
	      << ", " << stats << ")" << std::endl;
    solved_flag = true;
    return o.est_value();
  }
  else {
    if (trace_level > 0) {
      std::cerr << "no solution (estimated value "
		<< PRINT_NTYPE(o.est_value())
		<< ", " << stats << ")" << std::endl;
    }
    if (o.est_value() > lb) {
      insert_option_in_list(o);
    }
    if (options.empty()) {
      return 0;
    }
    else {
      return options.front().est_value();
    }
  }
}

void MaxNetBenefit::init_option(const index_set& selected, option& o)
{
  index_set g(instance.hard);
  o.goal_value = instance.null_value;
  for (index_type k = 0; k < selected.length(); k++) {
    g.insert(instance.soft[selected[k]].atoms);
    o.goal_value += instance.soft[selected[k]].weight;
  }
  o.goals = g;
  if (root_rs)
    o.root = new SeqCRegState(instance, heuristic, cost, g, root_rs);
  else
    o.root = new SeqCRegState(instance, heuristic, cost, g);
  o.est_cost = o.root->est_cost();
}

NTYPE MaxValueSearch::best_option_estimated_value()
{
  if (options.empty()) {
    return 0;
  }
  else {
    return options.front().est_value();
  }
}

index_type MaxValueSearch::best_option_size()
{
  if (options.empty()) {
    return 0;
  }
  else {
    return options.front().n_goals();
  }
}

void MaxValueSearch::init()
{
  stats.start();
  std::cerr << "initializing option list ("
	    << instance.n_soft() << " soft goals, lb = "
	    << lb << ")..."
	    << std::endl;
  init_option_list();
  solved_flag = false;
  stats.stop();
}

NTYPE MaxValueSearch::main()
{
  NTYPE best_value = best_option_estimated_value();
  stats.start();
  while (!solved_flag && !stats.break_signal_raised() && !options.empty()) {
    if (print_options_max > 0) {
      std::cerr << "top " << print_options_max
		<< " of " << n_options()
		<< " current options:" << std::endl;
      print_option_list(std::cerr);
    }
    best_value = explore_next_option();
    std::cerr << "current best value = " << PRINT_NTYPE(best_value);
    if (print_options_max == 0)
      std::cerr << ", " << n_options() << " options";
    std::cerr << " (solved = " << solved_flag << ", " << stats << ")"
	      << std::endl;
  }
  if (options.empty()) {
    if (instance.empty_plan_valid()) {
      std::cerr << "no more options - empty plan is best" << std::endl;
      make_empty_plan();
      solved_flag = true;
      return 0;
    }
    else {
      std::cerr << "no more options - problem unsolvable" << std::endl;
      return NEG_INF;
    }
  }
  stats.stop();
  return best_value;
}

END_HSPS_NAMESPACE
