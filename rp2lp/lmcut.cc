
#include "lmcut.h"

//#define TRACE_PRINT_LOTS

BEGIN_HSPS_NAMESPACE

LMCutBase::LMCutBase(Instance& i, const ACF& c, Stopwatch& s)
  : Heuristic(i), init_cost(c), stats(s), table(0)
{
  table = new CostTable(instance, s);
}

LMCutBase::~LMCutBase()
{
  delete table;
}

// void LMCutBase::extend_goal_set
// (const index_set& sgset, const AnyACF& costs, bool_vec& ext_goal_set)
// {
//   if (stats.break_signal_raised()) return;
//   // find atom in subgoal set that is not in ext. goal set
//   // with highest non-zero h^1 cost
//   index_type a_max = no_such_index;
//   NTYPE c_max = 0;
//   for (index_type k = 0; k < sgset.size(); k++)
//     if (!ext_goal_set[sgset[k]] && (table->eval(sgset[k]) > c_max))
//       a_max = sgset[k];
//   std::cerr << "extend_goal_set: sgset = ";
//   instance.write_atom_set(std::cerr, sgset);
//   if (a_max != no_such_index)
//     std::cerr << ", max atom = " << instance.atoms[a_max].name;
//   else
//     std::cerr << ", max atom = nil";
//   std::cerr << std::endl;
//   // if such an atom exists, add it to the extended goal set and
//   // recurse on all actions that add it and have zero cost
//   if (a_max != no_such_index) {
//     ext_goal_set[a_max] = true;
//     for (index_type k = 0; k < instance.atoms[a_max].add_by.size(); k++) {
//       index_type i = instance.atoms[a_max].add_by[k];
//       if (IS_ZERO(costs(i))) {
// 	std::cerr << "recursing on " << instance.actions[i].name << std::endl;
// 	extend_goal_set(instance.actions[i].pre, costs, ext_goal_set);
//       }
//     }
//   }
// }

void LMCutBase::extend_goal_set
(const index_set& sgset, const AnyACF& costs, bool_vec& ext_goal_set)
{
  if (stats.break_signal_raised()) return;
  // find all h^1-maximisers in the subgoal set, and check if any of them
  // already belongs to the extended goal set
  index_type a_max = no_such_index;
  NTYPE c_max = 0;
  bool  max_in_set = false;
  for (index_type k = 0; k < sgset.size(); k++) {
    if (table->eval(sgset[k]) > c_max) {
      c_max = table->eval(sgset[k]);
      max_in_set = false;
      if (ext_goal_set[sgset[k]])
	max_in_set = true;
      a_max = sgset[k];
    }
    else if (table->eval(sgset[k]) == c_max) {
      if (ext_goal_set[sgset[k]])
	max_in_set = true;
    }
  }
  // if no maximiser is in the extended goal set, we have to add one
  // (a_max), and recurse on all actions that add it and have zero cost.
  // if a_max == no_such_index, all atoms in sgset have h^1 value zero.
  if (!max_in_set && (a_max != no_such_index)) {
    ext_goal_set[a_max] = true;
    for (index_type k = 0; k < instance.atoms[a_max].add_by.size(); k++) {
      index_type i = instance.atoms[a_max].add_by[k];
      if (IS_ZERO(costs(i))) {
	extend_goal_set(instance.actions[i].pre, costs, ext_goal_set);
      }
    }
  }
}

void LMCutBase::find_cut
(const bool_vec& ext_goal_set, const AnyACF& costs, bool_vec& cut)
{
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (!IS_ZERO(costs(k)))
      if (instance.actions[k].add.have_common_element(ext_goal_set))
	if (FINITE(table->eval(instance.actions[k].pre)))
	  cut[k] = true;
}

NTYPE LMCutBase::compute
(const bool_vec& s, const index_set& g, const ACF& c, const bool_vec* a)
{
  AnyACF costs(instance.n_actions(), c);
  table->compute_H1(costs, s, a);
  if (INFINITE(table->eval(g)))
    return POS_INF;
  NTYPE total_cost = 0;
  while (table->eval(g) > 0) {
#ifdef TRACE_PRINT_LOTS
    std::cerr << "h^1 table:" << std::endl;
    table->write_pddl(std::cerr, (Instance&)instance);
#endif
    if (stats.break_signal_raised()) return 0;
    bool_vec ext_goal_set(false, instance.n_atoms());
    extend_goal_set(g, costs, ext_goal_set);
#ifdef TRACE_PRINT_LOTS
    std::cerr << "extended goal set: ";
    instance.write_atom_set(std::cerr, ext_goal_set);
    std::cerr << std::endl;
#endif
    bool_vec allowed_acts(true, instance.n_actions());
    if (a) {
      for (index_type k = 0; k < instance.n_actions(); k++)
	if (!(*a)[k]) allowed_acts[k] = false;
    }
    for (index_type k = 0; k < instance.n_actions(); k++)
      if (allowed_acts[k] &&
	  instance.actions[k].pre.have_common_element(ext_goal_set))
	allowed_acts[k] = false;
#ifdef TRACE_PRINT_LOTS
    for (index_type k = 0; k < instance.n_actions(); k++)
      if (!allowed_acts[k]) {
	std::cerr << "excluding " << instance.actions[k].name << std::endl;
      }
#endif
    table->compute_H1(costs, s, &allowed_acts);
#ifdef TRACE_PRINT_LOTS
    std::cerr << "revised h^1 table:" << std::endl;
    table->write_pddl(std::cerr, (Instance&)instance);
#endif
    bool_vec cut(false, instance.n_actions());
    find_cut(ext_goal_set, costs, cut);
    cut.intersect(allowed_acts);
    NTYPE c_cut = costs.min_cost(cut);
#ifdef TRACE_PRINT_LOTS
    std::cerr << "cut: ";
    instance.write_action_set(std::cerr, cut);
    std::cerr << ", cost = " << c_cut << std::endl;
#endif
    assert(c_cut > 0);
    total_cost += c_cut;
    costs.decrease(cut, c_cut);
    table->compute_H1(costs, s);
  }
#ifdef TRACE_PRINT_LOTS
  std::cerr << "total cost = " << total_cost << std::endl;
#endif
  return total_cost;
}

NTYPE LMCutBase::compute
(const index_set& s, const index_set& g, const ACF& c, const bool_vec* a)
{
  bool_vec s1(s, instance.n_atoms());
  return compute(s1, g, c, a);
}

NTYPE LMCutBase::compute
(const bool_vec& s, const index_set& g)
{
  return compute(s, g, init_cost, 0);
}

NTYPE LMCutBase::compute
(const index_set& s, const index_set& g)
{
  return compute(s, g, init_cost, 0);
}

void LMCutBase::select_relevant(const index_set& g, bool_vec& rel)
{
  if (g.empty()) return;
  // pick a subgoal with non-zero estimated cost in the set g: this can
  // be done in many ways. here we take the one that adds the fewest new
  // relevant actions, and tie-break on estimated cost (preferring higher)
  index_type g_best = 0;
  index_type r_best = (instance.atoms[g[0]].add_by.size() -
		       instance.atoms[g[0]].add_by.count_common(rel));
  NTYPE c_best = table->eval(g[0]);
  for (index_type k = 1; k < g.size(); k++) {
    NTYPE c_k = table->eval(g[k]);
    if (!IS_ZERO(c_k)) {
      if (IS_ZERO(c_best)) {
	g_best = k;
	r_best = (instance.atoms[g[k]].add_by.size() -
		  instance.atoms[g[k]].add_by.count_common(rel));
	c_best = c_k;
      }
      else {
	index_type r_k = (instance.atoms[g[k]].add_by.size() -
			  instance.atoms[g[k]].add_by.count_common(rel));
	if ((r_k < r_best) || ((r_k == r_best) && (c_k > c_best))) {
	  g_best = k;
	  r_best = r_k;
	  c_best = c_k;
	}
      }
    }
  }
  // all actions that add this subgoal are relevant; for each one that
  // isn't already in the relevant set, we must add it and recurse on
  // its preconditions
  for (index_type k = 0; k < instance.atoms[g[g_best]].add_by.size(); k++)
    if (!rel[instance.atoms[g[g_best]].add_by[k]]) {
      rel[instance.atoms[g[g_best]].add_by[k]] = true;
      select_relevant(instance.actions[instance.atoms[g[g_best]].add_by[k]].pre,
		      rel);
    }
}

NTYPE LMCutBase::compute2(const bool_vec& s, const index_set& g)
{
  AnyACF costs(instance.n_actions(), init_cost);
  table->compute_H1(costs, s);
  if (INFINITE(table->eval(g)))
    return POS_INF;
  NTYPE total_cost = 0;
  NTYPE g_cost = table->eval(g);
  while (g_cost > 0) {
    // std::cerr << "h^1(g) = " << PRINT_NTYPE(g_cost) << std::endl;
    if (stats.break_signal_raised()) return 0;
    // compute a set of "relevant actions"
    bool_vec cut(false, instance.n_actions());
    select_relevant(g, cut);
    // find the subset of relevant actions that have non-zero cost
    // and an estimated precondition cost of zero: this is our cut
    for (index_type k = 0; k < instance.n_actions(); k++)
      if (cut[k]) {
	if (IS_ZERO(costs(k)))
	  cut[k] = false;
	else if (!IS_ZERO(table->eval(instance.actions[k].pre)))
	  cut[k] = false;
      }
    NTYPE c_cut = costs.min_cost(cut);
    // std::cerr << "cut: ";
    // instance.write_action_set(std::cerr, cut);
    // std::cerr << ", cost = " << c_cut << std::endl;
    assert(c_cut > 0);
    total_cost += c_cut;
    costs.decrease(cut, c_cut);
    table->compute_H1(costs, s);
    g_cost = table->eval(g);
  }
  // std::cerr << "total cost = " << total_cost << std::endl;
  return total_cost;
}

NTYPE LMCutBase::get_landmark
(const bool_vec& s, const index_set& g, bool_vec& lm)
{
  table->compute_H1(init_cost, s);
  NTYPE g_cost = table->eval(g);
  if (INFINITE(g_cost) || IS_ZERO(g_cost))
    return g_cost;
  // compute a set of "relevant actions"
  lm.assign_value(false, instance.n_actions());
  select_relevant(g, lm);
  // find the subset of relevant actions that have non-zero cost
  // and an estimated precondition cost of zero: this is our cut,
  // which is also our landmark
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (lm[k]) {
      if (IS_ZERO(init_cost(k)))
	lm[k] = false;
      else if (!IS_ZERO(table->eval(instance.actions[k].pre)))
	lm[k] = false;
    }
  NTYPE c_cut = init_cost.min_cost(lm);
  assert(c_cut > 0);
  return c_cut;
}

ForwardLMCut::ForwardLMCut
(Instance& i, const index_set& g, const ACF& c, Stopwatch& s)
  : LMCutBase(i, c, s), goals(g)
{
  // done
}

NTYPE ForwardLMCut::eval(const index_set& s)
{
  bool_vec s1(s, instance.n_atoms());
  return compute(s1, goals);
}

NTYPE ForwardLMCut::eval(const bool_vec& s)
{
  return compute(s, goals);
}


NTYPE ForwardLMCut2::eval(const index_set& s)
{
  bool_vec s1(s, instance.n_atoms());
  return compute2(s1, goals);
}

NTYPE ForwardLMCut2::eval(const bool_vec& s)
{
  return compute2(s, goals);
}


RegressionLMCut::RegressionLMCut
(Instance& i, const ACF& c, Stopwatch& s)
  : LMCutBase(i, c, s), inits(false, instance.n_atoms())
{
  for (index_type k = 0; k < instance.n_atoms(); k++)
    inits[k] = instance.atoms[k].init;
}

RegressionLMCut::RegressionLMCut
(Instance& i, const index_set& s0, const ACF& c, Stopwatch& s)
  : LMCutBase(i, c, s), inits(s0, instance.n_atoms())
{
  // done
}

NTYPE RegressionLMCut::eval(const index_set& s)
{
  return LMCutBase::compute(inits, s, init_cost);
}

NTYPE RegressionLMCut::eval(const bool_vec& s)
{
  index_set s1(s);
  return LMCutBase::compute(inits, s1, init_cost);
}

NTYPE RegressionLMCut::compute(const index_set& g, const ACF& c)
{
  return LMCutBase::compute(inits, g, c);
}

RegressionLMCutP2::RegressionLMCutP2
(Instance& p2_ins, const ACF& c, Instance& p_ins,
 const s2index& pm, const index_vec& am, Stopwatch& s)
  : LMCutBase(p2_ins, c, s), base_ins(p_ins), atom_map(pm), action_map(am),
    inits(false, p2_ins.n_atoms())
{
  for (index_type i = 0; i < base_ins.n_atoms(); i++) {
    inits[atom_map(i, i)] = base_ins.atoms[i].init;
    // assert(inits[atom_map(i, i)] == p2_ins.atoms[atom_map(i, i)].init);
    for (index_type j = i + 1; j < base_ins.n_atoms(); j++) {
      assert(atom_map(i, j) < inits.length());
      inits[atom_map(i, j)] =
	(base_ins.atoms[i].init && base_ins.atoms[j].init);
      // assert(inits[atom_map(i, j)] == p2_ins.atoms[atom_map(i, j)].init);
    }
  }
}

NTYPE RegressionLMCutP2::compute(const bool_vec& s, const index_set& g)
{
  // std::cerr << "action map: " << action_map << std::endl;
  AnyACF costs(instance.n_actions());
  for (index_type k = 0; k < instance.n_actions(); k++) {
    assert(action_map[k] < base_ins.n_actions());
    costs.set_cost(k, init_cost(action_map[k]));
  }
  table->compute_H1(costs, s);
  if (INFINITE(table->eval(g)))
    return POS_INF;
  NTYPE total_cost = 0;
  while (table->eval(g) > 0) {
    if (stats.break_signal_raised()) return 0;
    bool_vec ext_goal_set(false, instance.n_atoms());
    extend_goal_set(g, costs, ext_goal_set);
    bool_vec allowed_acts(true, instance.n_actions());
    for (index_type k = 0; k < instance.n_actions(); k++)
      if (instance.actions[k].pre.have_common_element(ext_goal_set))
	allowed_acts[k] = false;
    table->compute_H1(costs, s, &allowed_acts);
    bool_vec cut(false, instance.n_actions());
    find_cut(ext_goal_set, costs, cut);
    cut.intersect(allowed_acts);
    NTYPE c_cut = costs.min_cost(cut);
    // std::cerr << "cut: ";
    // instance.write_action_set(std::cerr, cut);
    // std::cerr << ", cost = " << c_cut << std::endl;
    assert(c_cut > 0);
    total_cost += c_cut;
    bool_vec base_cut(false, base_ins.n_actions());
    for (index_type k = 0; k < instance.n_actions(); k++)
      if (cut[k]) base_cut[action_map[k]] = true;
    for (index_type k = 0; k < instance.n_actions(); k++)
      if (base_cut[action_map[k]]) cut[k] = true;
    costs.decrease(cut, c_cut);
    table->compute_H1(costs, s);
  }
  // std::cerr << "total cost = " << total_cost << std::endl;
  return total_cost;
}

NTYPE RegressionLMCutP2::eval(const index_set& s)
{
  index_set s2;
  for (index_type i = 0; i < s.length(); i++) {
    s2.insert(atom_map(s[i], s[i]));
    for (index_type j = i + 1; j < s.length(); j++) {
      s2.insert(atom_map(s[i], s[j]));
    }
  }
  return compute(inits, s2);
}

NTYPE RegressionLMCutP2::eval(const bool_vec& s)
{
  index_set s2;
  for (index_type i = 0; i < base_ins.n_atoms(); i++)
    if (s[i]) {
      s2.insert(atom_map(i, i));
      for (index_type j = i + 1; j < s.length(); j++)
	if (s[j]) {
	  s2.insert(atom_map(i, j));
	}
    }
  return compute(inits, s2);
}

END_HSPS_NAMESPACE
