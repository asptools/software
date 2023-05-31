
#include "ilb.h"
#include "enumerators.h"
#include "plans.h"

#include <fstream>
#include <sstream>
#include <algorithm>

//#define USE_NEWLM_BC
//#define USE_SCORE
//#define USE_LB1
//#define USE_LB2
//#define USE_LB3
//#define USE_LB4
//#define USE_OUT

#define ILA_USE_HS_EXTERN
#define ILA_USE_2ND_APX_HS
#define IMAI_FUKUNAGA

#ifdef IMAI_FUKUNAGA
#define HPLUS_USE_RP_EXTERN
#define USE_CPLEX
#endif

//#define USE_SCIP
//#define RP_SCIP_CG
#define USE_CPLEX
//#define USE_GECODE

// Apply zero cost actions to initial state before computing
// h+. This is different from the zero_cost_fill option, because
// it applies to all h+ computation methods, and before other
// reductions, such as relevance analysis.
// #define INIT_ZERO_FILL

// option HOFFMANN: use FF-style (non-optimal) delete-relaxed plan
// extraction instead of h+ inside hplus_with_ce loop.
//#define HOFFMANN

#define HPLUS_PRINT_STATS
//#define HPLUS_PRINT_STATS_EVERY_ITERATION
#define HPLUSPLUS_PRINT_STATS
//#define MEASURE_DENSITY
//#define MEASURE_STABILITY

#define TRACE_PRINT_HLB
//#define TRACE_PRINT_META_ATOMS_AS_SETS

//#define WSAT_REDUNDANT_CONSTRAINTS

#ifndef HAVE_GLPK
#undef USE_LB3
#endif

#ifdef USE_LB3
#include <glpk.h>
#endif

#ifndef HAVE_SCIP
#undef USE_SCIP
#endif

#ifdef USE_SCIP
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#ifdef HPLUS_USE_RP_EXTERN
#ifndef RP_SCIP_CG
#include "examples/LOP/src/cons_linearordering.h"
#endif
#endif
#endif

#ifndef HAVE_CPLEX
#undef USE_CPLEX
#endif

#ifdef USE_CPLEX
#ifndef CPLEX_INCREMENTAL
#include <ilcplex/ilocplex.h>
#endif
#endif

#ifndef HAVE_GECODE
#undef USE_GECODE
#endif

#ifdef USE_GECODE
#include "gecode/kernel.hh"
#include "gecode/int.hh"
#include "gecode/search.hh"
#include "gecode/bs.h"
#endif

typedef IloArray<IloBoolVarArray> IloBoolVarArray2;



BEGIN_HSPS_NAMESPACE

NTYPE ILB::hlb = 0;

void addedge(int **g, int *gs, int *gn, int s, int t, int check){
  if (check==1)
    for (int i=0; i<gn[s]; i++)
      if (g[s][i]==t)
        return;
      if (gs[s]==gn[s]){
        g[s]= (int*) realloc(g[s],gs[s]*2*sizeof(int) );
        gs[s]=gs[s]*2;
        
      }
      g[s][gn[s]]=t;
      gn[s]++;
}




ILB::ConflictModifier::ConflictModifier()
{
  // dummy constructor
}

ILB::ConflictModifier::~ConflictModifier()
{
  // dummy destructor
}

void ILB::ConflictModifier::apply(ILB& ilb, index_set_vec& cs)
{
  // default modifier does nothing
}

void ILB::new_lb(NTYPE lb)
{
  if (lb > hlb) {
    hlb = lb;
    std::cout << ";; hlb = " << hlb << std::endl;
  }
}

void ILB::reset_hlb()
{
  hlb = 0;
}

bool ILB::set_compare::operator()
(const HSPS::bool_vec& s0, const HSPS::bool_vec& s1) const
{
  for (index_type i = 0; (i < s0.size()) && (i < s1.size()); i++)
    if (s0[i] != s1[i]) return false;
  if (s0.size() < s1.size()) {
    for (index_type i = s0.size(); i < s1.size(); i++)
      if (s1[i]) return false;
  }
  else if (s1.size() < s0.size()) {
    for (index_type i = s1.size(); i < s0.size(); i++)
      if (s0[i]) return false;
  }
  return true;
}

set_hash_function ILB::set_hash::f;

unsigned int ILB::set_hash::operator()
(const HSPS::bool_vec& s0) const
{
  return f(s0);
}

unsigned int ILB::set_hash::operator()
(const HSPS::index_set& s0) const
{
  return f(s0);
}

NTYPE EvalWithMetaAtoms::eval(const index_set& s)
{
  index_set b(s);
  bool done = false;
  while (!done) {
    done = true;
    index_set b1;
    for (index_type k = 0; k < b.size(); k++) {
      index_type i = meta_atom_map.find_rule(b[k]);
      if (i != no_such_index) {
	assert(i < meta_atom_map.size());
	b1.insert(meta_atom_map[i].antecedent);
	done = false;
      }
      else {
	b1.insert(b[k]);
      }
    }
    if (!done)
      b.assign_copy(b1);
  }
  return base_h.eval(b);
}

NTYPE EvalWithMetaAtoms::eval(const bool_vec& s)
{
  index_set b(s);
  return eval(b);
}

ILB::ILB(Instance& i, const ACF& c, Heuristic* n, Stopwatch& s)
  : ins(i), is_sub(false),
    n_real_atoms(i.n_atoms()), n_real_actions(i.n_actions()),
    cost(c), inc(n), stats(s), hs_stats(&s), hss_stats(&hs_stats),
    hsa_stats(&hs_stats), lm_stats(&s), h1_stats(&s), ra_stats(&s),
    rpn_stats(&s), ce_stats(&s), sce_stats(&s),
    calls_to_hplus(0), calls_to_hplus_ce(0), calls_to_newlm(0),
    calls_to_hs_opt(0), relevant_actions_ratio(0), sum_density(0),
    calls_to_hs_apx(0), hs_apx_improve(0), sum_width(0), n_width(0),
    hs_nodes(0), hs_cache_hits(0), hs_cache_miss(0), hs_lb_calls(0),
    hs_lb1_max(0), hs_lb2_max(0), hs_lb3_max(0), hs_lb4_max(0),
    hs_split1(0), hs_splits(0), hs_branch(0),
    hpp_iterations(0), conflict_set_size(0), pc_actions_ratio(0),
    n_ce_compiled(0), max_ce_compiled(0), cce_actions_ratio(0),
    reach(i, false), verbose_level(1), check_rp_strategy(0),
    remove_dominated_conditions(false), prune_relaxed_irrelevant(true),
    prune_relaxed_dominated(true), iterated_preprocessing(true),
    ILA_use_approximate(false), ILA_use_saturation(false),
    ILA_use_lmcut(false), zero_cost_fill(false),
    hitting_set_use_split(false), hitting_set_use_dominance(false),
    conflict_mod(NULL), h1(NULL)
{
#ifdef NTYPE_RATIONAL
  NTYPE d = cost.cost_gcd(ins.n_actions());
  ac_div = d.divisor();
#endif
  set_hash::f.init(5000);
  relevant_atoms.assign_value(true, ins.n_actions());
  relevant_actions.assign_value(true, ins.n_actions());
  h1 = new CostTable(ins, h1_stats);
#ifdef USE_NEWLM_BC
  h1->compute_H1(UnitACF());
#endif
#ifdef RANDOMIZE
  rng.seed_with_pid();
#endif
#ifdef CPLEX_INCREMENTAL
  cplex_model_atoms = 0;
  cplex_model_actions = 0;
  cplex_model_landmarks = 0;
#endif
}

ILB::~ILB()
{
  if (h1) delete h1;
#ifdef CPLEX_INCREMENTAL
  cplex_env.end();
#endif
}

void ILB::copy_options(const ILB& ilb, int vl_adjustment)
{
  verbose_level = ilb.verbose_level + vl_adjustment;
  remove_dominated_conditions = ilb.remove_dominated_conditions;
  prune_relaxed_irrelevant = ilb.prune_relaxed_irrelevant;
  prune_relaxed_dominated = ilb.prune_relaxed_dominated;
  iterated_preprocessing = ilb.iterated_preprocessing;
  ILA_use_approximate = ilb.ILA_use_approximate;
  ILA_use_saturation = ilb.ILA_use_saturation;
  ILA_use_lmcut = ilb.ILA_use_lmcut;
  zero_cost_fill = ilb.zero_cost_fill;
  hitting_set_use_split = ilb.hitting_set_use_split;
  hitting_set_use_dominance = ilb.hitting_set_use_dominance;
}

void ILB::add_hplus_stats(const ILB& ilb)
{
  calls_to_hplus += ilb.calls_to_hplus;
  calls_to_newlm += ilb.calls_to_newlm;
  calls_to_hs_opt += ilb.calls_to_hs_opt;
  relevant_actions_ratio += ilb.relevant_actions_ratio;
  sum_density += ilb.sum_density;
  calls_to_hs_apx += ilb.calls_to_hs_apx;
  hs_apx_improve += ilb.hs_apx_improve;
  sum_width += ilb.sum_width;
  n_width += ilb.n_width;
  hs_nodes += ilb.hs_nodes;
  hs_cache_hits += ilb.hs_cache_hits;
  hs_cache_miss += ilb.hs_cache_miss;
  hs_lb_calls += ilb.hs_lb_calls;
  hs_lb1_max += ilb.hs_lb1_max;
  hs_lb2_max += ilb.hs_lb2_max;
  hs_lb3_max += ilb.hs_lb3_max;
  hs_lb4_max += ilb.hs_lb4_max;
  hs_split1 += ilb.hs_split1;
  hs_splits += ilb.hs_splits;
  hs_branch += ilb.hs_branch;
  hs_stats.add(ilb.hs_stats);
  hss_stats.add(ilb.hss_stats);
  hsa_stats.add(ilb.hsa_stats);
  lm_stats.add(ilb.lm_stats);
  ra_stats.add(ilb.ra_stats);
  rpn_stats.add(ilb.rpn_stats);
}

///
/// OLD VERSION
///
// void ILB::L1C1(const index_set& goal)
// {
//   for (index_type i = 0; i < goal.length(); i++)
//     if (!relevant_atoms[goal[i]]) {
//       relevant_atoms[goal[i]] = true;
//       for (index_type k = 0; k < ins.atoms[goal[i]].add_by.length(); k++)
// 	if (!relevant_actions[ins.atoms[goal[i]].add_by[k]]) {
// 	  Instance::Action& act = ins.actions[ins.atoms[goal[i]].add_by[k]];
// 	  bool ok = true;
// 	  for (index_type j = 0; (j < act.pre.length()) && ok; j++)
// 	    if (atom_lmg.adjacent(goal[i], act.pre[j]))
// 	      ok = false;
// 	  if (ok) {
// 	    relevant_actions[act.index] = true;
// 	    L1C1(act.pre);
// 	  }
// 	}
//     }
// }
//
// void ILB::compute_relevant(const index_set& init, const index_set& goal)
// {
//   ra_stats.start();
//   if (verbose_level > 1) {
//     std::cerr << "computing atom landmark graph..." << std::endl;
//   }
//   reach.recompute(init);
//   bool_vec reach1(reach.reachable());
//   atom_lmg.init(ins.n_atoms());
//   for (index_type i = 0; i < ins.n_atoms(); i++) {
//     bool_vec acts_without_i(true, ins.n_actions());
//     acts_without_i.subtract(ins.atoms[i].req_by);
//     index_set init_without_i(init);
//     init_without_i.subtract(i);
//     reach.recompute(init_without_i, acts_without_i);
//     for (index_type j = 0; j < ins.n_atoms(); j++)
//       if ((j != i) && reach1[j] && reach.unreachable(j))
// 	atom_lmg.add_edge(i, j);
//   }
//   if (verbose_level > 1) {
//     std::cerr << "computing relevance..." << std::endl;
//   }
//   relevant_atoms.assign_value(false, ins.n_atoms());
//   relevant_actions.assign_value(false, ins.n_actions());
//   L1C1(goal);
//   if (verbose_level > 0) {
//     std::cerr << relevant_actions.count(true) << " of " << ins.n_actions()
// 	      << " actions are L1C1 relevant" << std::endl;
//   }
//   ra_stats.stop();
// }
///
/// END OLD
///

void ILB::init_preprocessing()
{
  FA_and_ALM_avail = false;
  relevant_atoms.assign_value(true, ins.n_atoms());
  relevant_actions.assign_value(true, ins.n_actions());
}

bool ILB::compute_FA_and_ALM
(const index_set& init)
{
  ra_stats.start();
  if (verbose_level > 1) {
    std::cerr << "computing first achievers..." << std::endl;
  }
  first_achievers.assign_value(EMPTYSET, ins.n_atoms());
  if (!FA_and_ALM_avail) {
    atom_landmarks.assign_value(EMPTYSET, ins.n_atoms());
  }
  bool changed = false;
  reach.recompute(init, relevant_actions);
  bool_vec reachable_atom(reach.reachable());
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (!init.contains(i)) {
      // compute atoms reachable without actions that add i
      //bool_vec acts(true, ins.n_actions());
      bool_vec acts(relevant_actions);
      acts.subtract(ins.atoms[i].add_by);
      reach.recompute(init, acts);
      for (index_type j = 0; j < ins.atoms[i].add_by.size(); j++)
	if (relevant_actions[ins.atoms[i].add_by[j]])
	  if (!reach.unreachable_action(ins.atoms[i].add_by[j]))
	    first_achievers[i].insert(ins.atoms[i].add_by[j]);
      for (index_type j = 0; j < ins.n_atoms(); j++)
	if ((j != i) && reachable_atom[j] && reach.unreachable(j))
	  if (!atom_landmarks[j].contains(i)) {
	    if (verbose_level > 1) {
	      std::cerr << "new: " << i << "." << ins.atoms[i].name
			<< " lm of " << j << "." << ins.atoms[j].name
			<< std::endl;
	    }
	    atom_landmarks[j].insert(i);
	    changed = true;
	  }
      if (stats.break_signal_raised()) {
	ra_stats.stop();
	return changed;
      }
    }
  FA_and_ALM_avail = true;
  ra_stats.stop();
  if (verbose_level > 2) {
    std::cerr << "first achievers: " << std::endl;
    for (index_type i = 0; i < ins.n_atoms(); i++)
      std::cerr << i << " : " << first_achievers[i] << std::endl;
    std::cerr << "changed = " << changed << std::endl;
  }
  return changed;
}

void ILB::L1C1(const index_set& goal, bool_vec& ratm, bool_vec& ract)
{
  for (index_type i = 0; i < goal.length(); i++)
    if (!ratm[goal[i]]) {
      ratm[goal[i]] = true;
      for (index_type k = 0; k < first_achievers[goal[i]].length(); k++)
	if (!ract[first_achievers[goal[i]][k]]) {
	  ract[first_achievers[goal[i]][k]] = true;
	  L1C1(ins.actions[first_achievers[goal[i]][k]].pre, ratm, ract);
	}
    }
}

void ILB::compute_relevant
(const index_set& init,
 const index_set& goal)
{
  ra_stats.start();
  if (!FA_and_ALM_avail) {
    compute_FA_and_ALM(init);
    if (stats.break_signal_raised()) {
      ra_stats.stop();
      return;
    }
  }
  if (verbose_level > 1) {
    std::cerr << "computing relevance..." << std::endl;
  }
  bool_vec r1atm(false, ins.n_atoms());
  bool_vec r1act(false, ins.n_actions());
  L1C1(goal, r1atm, r1act);
  relevant_atoms.intersect(r1atm);
  relevant_actions.intersect(r1act);
  if (verbose_level > 0) {
    std::cerr << relevant_actions.count(true) << " of " << ins.n_actions()
	      << " actions are L1C1 relevant" << std::endl;
    if (verbose_level > 2) {
      std::cerr << "relevant atoms: " << index_set(relevant_atoms) << std::endl;
      std::cerr << "relevant actions: " << index_set(relevant_actions) << std::endl;
    }
  }
  ra_stats.stop();
}

///
/// OLD VERSION
///
// void ILB::remove_relaxed_redundant_actions(const index_set& init)
// {
//   for (index_type k = 0; k < ins.n_actions(); k++)
//     if (relevant_actions[k]) {
//       bool is_red = false;
//       index_set pre_a(ins.actions[k].pre);
//       pre_a.subtract(init);
//       index_set add_a(ins.actions[k].add);
//       add_a.subtract(init);
//       for (index_type j = 0; (j < ins.n_actions()) && !is_red; j++)
// 	if (relevant_actions[j]) {
// 	  if (cost(j) <= cost(k)) {
// 	    index_set pre_b(ins.actions[j].pre);
// 	    pre_b.subtract(init);
// 	    index_set add_b(ins.actions[j].add);
// 	    add_b.subtract(init);
// 	    if (pre_a.contains(pre_b) && add_b.contains(add_a)) {
// 	      if ((pre_a.size() == pre_b.size()) &&
// 		  (add_a.size() == add_b.size()) &&
// 		  (cost(k) == cost(j))) {
// 		if (j < k)
// 		  is_red = true;
// 	      }
// 	      else {
// 		is_red = true;
// 	      }
// 	    }
// 	  }
// 	}
//       if (is_red)
// 	relevant_actions[k] = false;
//     }
//   if (verbose_level > 0) {
//     std::cerr << relevant_actions.count(true) << " of " << ins.n_actions()
// 	      << " actions are non-redundant" << std::endl;
//   }
// }
///
/// END OLD
///

void ILB::remove_dominated_actions(const index_set& init)
{
  ra_stats.start();
  if (!FA_and_ALM_avail) {
    compute_FA_and_ALM(init);
    if (stats.break_signal_raised()) {
      ra_stats.stop();
      return;
    }
  }
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (relevant_actions[k]) {
      // action k is redundant if there exists a dominating action b;
      // candb is the set of candidates for this action b.
      // candidates must themselves be relevant:
      bool_vec candb(relevant_actions);
      // k cannot dominate itself:
      candb[k] = false;
      // candidates must add everything in add(k) - init:
      // for (index_type i = 0; i < ins.actions[k].add.size(); i++)
      // 	if (!init.contains(ins.actions[k].add[i]))
      // 	  candb.intersect(ins.atoms[ins.actions[k].add[i]].add_by);
      // alt: candidate must be first ach of everything that k is:
      for (index_type i = 0; i < ins.actions[k].add.size(); i++)
	if (relevant_atoms[ins.actions[k].add[i]])
	  if (first_achievers[ins.actions[k].add[i]].contains(k))
	    candb.intersect(first_achievers[ins.actions[k].add[i]]);
      if (candb.count(true) > 0) {
	// candidates cannot have a prec not in pre(k) U init
	bool_vec allowed(init, ins.n_atoms());
	allowed.insert(ins.actions[k].pre);
	for (index_type j = 0; j < ins.actions[k].pre.size(); j++)
	  allowed.insert(atom_landmarks[ins.actions[k].pre[j]]);
	allowed.complement(); // now it's actually "disallowed"
	bool is_red = false;
	for (index_type j = 0; (j < ins.n_actions()) && !is_red; j++)
	  if (candb[j])
	    // candidates must have cost <= cost of action k
	    if (cost(j) <= cost(k)) {
	      // candidates must not have a disallowed precondition
	      if (!ins.actions[j].pre.have_common_element(allowed)) {
		if (verbose_level > 2) {
		  std::cerr << "action " << k << "." << ins.actions[k].name
			    << " is dominated by " << j << "."
			    << ins.actions[j].name << std::endl;
		}
		is_red = true;
	      }
	    }
	// note: don't need to check symmetry; if two actions are
	// equivalent (both dominate the other), one will be tested (and
	// removed) first, which prevents it from dominating the other
	// when that is later tested.
	if (is_red)
	  relevant_actions[k] = false;
      }
      if (stats.break_signal_raised()) {
	ra_stats.stop();
	return;
      }
    }
  if (verbose_level > 0) {
    std::cerr << relevant_actions.count(true) << " of " << ins.n_actions()
	      << " actions are (also) non-dominated" << std::endl;
  }
  ra_stats.stop();
}

NTYPE ILB::weighted_degree(index_type i)
{
  assert(i < landmarks.size());
  NTYPE s = 0;
  for (index_type k = 0; k < lm_conflict_set[i].size(); k++)
    s += lm_min_cost[lm_conflict_set[i][k]];
#ifdef NTYPE_RATIONAL
  NTYPE wd = safemul(s, lm_min_cost[i].invert()).round();
#else
  NTYPE wd = (s / lm_min_cost[i]);
#endif
  return wd;
}

void ILB::update_weighted_degree()
{
  lm_wd.assign_value(0, landmarks.size());
  for (index_type k = 0; k < landmarks.size(); k++)
    lm_wd[k] = weighted_degree(k);
}

NTYPE ILB::hitting_set_lb_1(const bool_vec& rem, const bool_vec& out)
{
  // approximate weighted independent set algorithm
  weighted_vec<index_type,NTYPE> sorted;
  for (index_type k = 0; k < landmarks.size(); k++)
    if (rem[k] && (lm_min_cost[k] > 0))
      sorted.insert_increasing(k, lm_wd[k]);
  bool_vec may_be_in(rem);
  NTYPE w = 0;
  for (index_type k = 0; k < sorted.size(); k++) {
    index_type i = sorted[k].value;
    if (may_be_in[i]) {
      may_be_in.subtract(lm_conflict_set[i]);
      may_be_in[i] = false;
      w += lm_min_cost[i];
    }
  }
  return w;
}

NTYPE ILB::hitting_set_lb_4(const bool_vec& rem, const bool_vec& out)
{
  NTYPE total = 0;
  AnyACF rem_cost(ins.n_actions(), cost);
  for (index_type k = 0; k < landmarks.size(); k++)
    if (rem[k]) {
      index_type a_min = rem_cost.min_cost_action(landmarks[k]);
      assert(a_min != no_such_index);
      assert(!out[a_min]); // no support for action exclusion yet
      NTYPE c_a = rem_cost(a_min);
      if (c_a > 0) {
	total += c_a;
	rem_cost.decrease(landmarks[k], c_a);
      }
    }
  return total;
}

NTYPE ILB::hitting_set_lb_2(const bool_vec& rem, const bool_vec& out)
{
  index_set rx(rem);
  lvector<double> u(0, rx.size());
  lvector<double> a(0, ins.n_actions());
  for (index_type i = 0; i < ins.n_actions(); i++)
    a[i] = N_TO_D(cost(i));
  lvector<double> d(0, rx.size());
  for (index_type i = 0; i < rx.size(); i++) {
    d[i] = a[landmarks[rx[i]][0]];
    for (index_type j = 1; j < landmarks[rx[i]].size(); j++)
      if (a[landmarks[rx[i]][j]] < d[i])
	d[i] = a[landmarks[rx[i]][j]];
  }

  // initialise u
  bool done = false;
  while (!done) {
    index_type i_min = no_such_index;
    index_type s_min = index_type_max;
    for (index_type i = 0; i < rx.size(); i++)
      if (d[i] > 0)
	if (landmarks[rx[i]].size() < s_min)
	  i_min = i;
    if (i_min != no_such_index) {
      double c = d[i_min];
      u[i_min] += c;
      for (index_type j = 0; j < landmarks[rx[i_min]].size(); j++)
	a[landmarks[rx[i_min]][j]] -= c;
      for (index_type i = 0; i < rx.size(); i++)
	for (index_type j = 0; j < landmarks[rx[i]].size(); j++)
	  if (a[landmarks[rx[i]][j]] < d[i])
	    d[i] = a[landmarks[rx[i]][j]];
      //assert(IS_ZERO(d[i_min]));
    }
    else {
      done = true;
    }
  }

  // local improvement loop
  done = false;
  while (!done) {
    done = true;
    for (index_type i1 = 0; (i1 < rx.size()) && done; i1++)
      if (u[i1] > 0) {
	for (index_type i2 = 0; (i2 < rx.size()) && done; i2++)
	  if (i2 != i1) {
	    bool ok = true;
	    for (index_type k = 0; (k < landmarks[rx[i2]].size()) && ok; k++)
	      if ((a[landmarks[rx[i2]][k]] > 0) &&
		  !landmarks[rx[i1]].contains(landmarks[rx[i2]][k]))
		ok = false;
	    if (ok) {
	      for (index_type i3 = 0; (i3 < rx.size()) && done; i3++)
		if ((i3 != i1) && (i3 != i2)) {
		  ok = true;
		  for (index_type k = 0; (k < landmarks[rx[i3]].size()) && ok;
		       k++) {
		    if ((a[landmarks[rx[i3]][k]] > 0) &&
			!landmarks[rx[i1]].contains(landmarks[rx[i3]][k]))
		      ok = false;
		    if ((a[landmarks[rx[i3]][k]] > 0) &&
			landmarks[rx[i2]].contains(landmarks[rx[i3]][k]))
		      ok = false;
		  }
		  if (ok) {
		    double delta = u[i1];
		    index_vec divi(0, ins.n_actions());
		    for (index_type k = 0; k < landmarks[rx[i2]].size(); k++)
		      divi[landmarks[rx[i2]][k]] += 1;
		    for (index_type k = 0; k < landmarks[rx[i3]].size(); k++)
		      divi[landmarks[rx[i3]][k]] += 1;
		    for (index_type k = 0; k < landmarks[rx[i1]].size(); k++)
		      if (divi[landmarks[rx[i1]][k]] > 0)
			divi[landmarks[rx[i1]][k]] -= 1;
		    for (index_type j = 0; j < ins.n_actions(); j++)
		      if (divi[j] > 0) {
			double limj = a[j]/divi[j];
			if (limj < 0.001)
			  delta = 0;
			else if (limj < delta)
			  delta = limj;
		      }
		    if (delta > 0) {
		      u[i1] -= delta;
		      u[i2] += delta;
		      u[i3] += delta;
		      for (index_type k = 0; k < landmarks[rx[i1]].size(); k++)
			a[landmarks[rx[i1]][k]] += delta;
		      for (index_type k = 0; k < landmarks[rx[i2]].size(); k++) {
			a[landmarks[rx[i2]][k]] -= delta;
			assert(a[landmarks[rx[i2]][k]] >= 0);
		      }
		      for (index_type k = 0; k < landmarks[rx[i3]].size(); k++) {
			a[landmarks[rx[i3]][k]] -= delta;
			assert(a[landmarks[rx[i3]][k]] >= 0);
		      }
		      done = false;
		    }
		  }
		}
	    }
	  }
      }
  }

  double s = 0;
  for (index_type i = 0; i < rx.size(); i++)
    s += u[i];
#ifdef NTYPE_RATIONAL
  NTYPE b = rational::dtor(s);
  b = rational::ceil_to(b, ac_div);
  return b;
#else
  return s;
#endif
}

NTYPE ILB::hitting_set_lb_3(const bool_vec& rem, const bool_vec& out)
{
#ifdef USE_LB3
  index_set rx(rem);

  glp_prob* lp = glp_create_prob();
  glp_set_prob_name(lp, "lb");
  glp_set_obj_dir(lp, GLP_MIN);

  // rows:
  // 1 .. #sets
  int n_rows = rx.size();
  glp_add_rows(lp, n_rows);
  for (int i = 1; i <= n_rows; i++) {
    glp_set_row_bnds(lp, i, GLP_LO, 1.0, 0.0);
  }

  // column index/value vector
  int c_ind[n_rows + 2];
  double c_val[n_rows + 2];
  int len;

  // columns:
  // 1 .. #actions (in lm_actions)
  glp_add_cols(lp, lm_actions.size());
  for (HSPS::index_type j = 0; j < lm_actions.size(); j++) {
#ifdef USE_OUT
    if (out[lm_actions[j]])
      glp_set_col_bnds(lp, j + 1, GLP_FX, 0.0, 0.0);
    else
      glp_set_col_bnds(lp, j + 1, GLP_DB, 0.0, 1.0);
#else
    glp_set_col_bnds(lp, j + 1, GLP_DB, 0.0, 1.0);
#endif
    len = 0;
    for (HSPS::index_type i = 0; i < rx.size(); i++)
      if (landmarks[rx[i]].contains(lm_actions[j])) {
	len += 1;
	c_ind[len] = i + 1;
	c_val[len] = 1.0;
      }
    glp_set_mat_col(lp, j + 1, len, c_ind, c_val);
    glp_set_obj_coef(lp, j + 1, N_TO_D(cost(lm_actions[j])));
  }

  // solve
  glp_smcp params;
  glp_init_smcp(&params);
  params.msg_lev = GLP_MSG_ERR;
  glp_simplex(lp, &params);

  int status = glp_get_status(lp);
  if (status != GLP_OPT) {
    std::cerr << "error: lp status is " << status << std::endl;
    exit(1);
  }
  double objective_value = glp_get_obj_val(lp);
  glp_delete_prob(lp);
  NTYPE b = D_TO_N(objective_value);
#ifdef NTYPE_RATIONAL
  b = rational::ceil_to(b, ac_div);
#endif
  return b;
#else
  return ZERO;
#endif
}

bool ILB::hitting_set_fc(const bool_vec& rem, const bool_vec& out)
{
  for (index_type k = 0; k < landmarks.size(); k++)
    if (rem[k]) {
      bool ok = false;
      for (index_type i = 0; (i < lm_branch_order[k].size()) && !ok; i++)
	if (!out[lm_branch_order[k][i]])
	  ok = true;
      if (!ok) return true;
    }
  return false;
}

NTYPE ILB::hitting_set_lb(const bool_vec& rem, const bool_vec& out)
{
  hs_lb_calls += 1;
#ifdef USE_OUT
  if (hitting_set_fc(rem, out))
    return POS_INF;
#endif
  NTYPE max_lb = 0;
#ifdef USE_LB3
  NTYPE lb3 = hitting_set_lb_3(rem, out);
  max_lb = MAX(max_lb, lb3);
#endif
#ifdef USE_LB2
  NTYPE lb2 = hitting_set_lb_2(rem, out);
  max_lb = MAX(max_lb, lb2);
#endif
#ifdef USE_LB1
  NTYPE lb1 = hitting_set_lb_1(rem, out);
  max_lb = MAX(max_lb, lb1);
#endif
#ifdef USE_LB4
  NTYPE lb4 = hitting_set_lb_4(rem, out);
  max_lb = MAX(max_lb, lb4);
#endif

#ifdef USE_LB3
  if (max_lb == lb3) hs_lb3_max += 1;
#endif
#ifdef USE_LB2
  if (max_lb == lb2) hs_lb2_max += 1;
#endif
#ifdef USE_LB1
  if (max_lb == lb1) hs_lb1_max += 1;
#endif
#ifdef USE_LB4
  if (max_lb == lb4) hs_lb4_max += 1;
#endif

  return max_lb;

/// old implement:
// #ifdef USE_LB3
//   NTYPE lb3 = hitting_set_lb_3(rem, out);
// #ifdef USE_LB1
//   NTYPE lb1 = hitting_set_lb_1(rem, out);
//   if (lb3 > lb1) {
//     hs_lb3_max += 1;
//     return lb3;
//   }
//   else {
//     if (lb1 > lb3) hs_lb1_max += 1;
//     return lb1;
//   }
// #else // lb3 only
//   return lb3;
// #endif
// #else // no lb3
// #ifdef USE_LB2
//   NTYPE lb2 = hitting_set_lb_2(rem, out);
// #ifdef USE_LB1
//   NTYPE lb1 = hitting_set_lb_1(rem, out);
//   if (lb2 > lb1) {
//     hs_lb2_max += 1;
//     return lb2;
//   }
//   else {
//     if (lb1 > lb2) hs_lb1_max += 1;
//     return lb1;
//   }
// #else // lb2 only
//   return lb2;
// #endif
// #else // no lb2 (and no lb3)
//   NTYPE lb1 = hitting_set_lb_1(rem, out);
//   hs_lb1_max += 1;
//   return lb1;
// #endif
// #endif
}

NTYPE ILB::hitting_set
(const bool_vec& rem,
 bool_vec& out,
 bool_vec& set,
 NTYPE acc,
 index_set& best,
 NTYPE& ub,
 NTYPE lb)
{
  hs_nodes += 1;
  if (acc >= ub)
    return acc;
  // choose a set that remains to hit
  index_type i_branch = no_such_index;
  NTYPE co_max = -1;
  index_type ca_min = index_type_max;
  index_type n_rem = 0;
  for (index_type i = 0; i < rem.size(); i++)
    if (rem[i]) {
      if (lm_min_cost[i] > co_max) {
	i_branch = i;
	co_max = lm_min_cost[i];
	ca_min = lm_branch_order[i].size();
      }
      else if ((lm_min_cost[i] == co_max) &&
	       (lm_branch_order[i].size() < ca_min)) {
	i_branch = i;
	ca_min = lm_branch_order[i].size();
      }
      n_rem += 1;
    }

  // if no set remains to hit, we have a solution; it should be better
  if (i_branch == no_such_index) {
    assert(n_rem == 0);
    assert(acc < ub);
    set.copy_to(best);
    ub = acc;
    if (verbose_level > 2) {
      std::cerr << "new best rp cost: " << ub << std::endl;
    }
    return acc;
  }
  // lookup or compute lower bound
  NTYPE est = 0;
#ifdef USE_CACHE
  set_map::iterator p = store_lb.find(rem);
  if (p != store_lb.end()) {
    est = p->second;
    hs_cache_hits += 1;
  }
  else {
    est = hitting_set_lb(rem, out);
    hs_cache_miss += 1;
  }
#else
  est = hitting_set_lb(rem, out);
#endif
  if ((acc + est) >= ub) {
    return acc + est;
  }
  // else we have to branch...
  hs_branch += 1;
  bool_vec rem_copy(rem);
#ifdef USE_OUT
  //bool_vec out_reset(false, lm_branch_order[i_branch].size());
  bool_vec out_copy(out);
#endif
  NTYPE v_min = POS_INF;
  for (index_type k = 0; k < lm_branch_order[i_branch].size(); k++)
#ifdef USE_OUT
    if (!out[lm_branch_order[i_branch][k]])
#endif
    {
      assert(!out[lm_branch_order[i_branch][k]]);
      index_type a = lm_branch_order[i_branch][k];
      set[a] = true;
      for (index_type j = 0; j < lm_with[a].size(); j++)
	rem_copy[lm_with[a][j]] = false;
#ifdef USE_OUT
      for (index_type j = 0; j < lm_branch_order[i_branch].size(); j++)
	if (j != k) {
	  //if (!out[lm_branch_order[i_branch][j]])
	  //  out_reset[j] = true;
	  //out[lm_branch_order[i_branch][j]] = true;
	  out[lm_branch_order[i_branch][j]] = true;
	}
#endif
      NTYPE v = hitting_set(rem_copy, out, set, acc + cost(a), best, ub, lb);
      if (stats.break_signal_raised()) return lb;
      set[a] = false;
      // note: this check is redundant, given the one on v_min
      // // if we have found a better solution on this call, ub may have
      // // decreased; if we cannot find an even better solution below this
      // // node, return
      // if ((acc + est) >= ub) {
      //   assert(v == ub);
      //   return acc + est;
      // }
      // if the value returned by last branch equals (global) lower bound,
      // we cannot find any better solution, so return.
      if (v == lb) return v;
      v_min = MIN(v_min, v);
      // if v_min (== best value found below this node) equals acc + est
      // (== best value possible through this node), we have already found
      // the best solution (if any) we can through this node, so return.
      if (v_min == (acc + est)) {
	return v_min;
      }
      if ((k + 1) < lm_branch_order[i_branch].size())
	rem_copy.assign_copy(rem);
#ifdef USE_OUT
      //for (index_type j = 0; j < lm_branch_order[i_branch].size(); j++)
      //if (out_reset[j])
      //  out[lm_branch_order[i_branch][j]] = false;
      //if ((k + 1) < lm_branch_order[i_branch].size())
      //out_reset.assign_value(false, lm_branch_order[i_branch].size());
      out.assign_copy(out_copy);
#endif
    }
#ifdef USE_CACHE
  if ((v_min - acc) > est)
    store_lb[rem] = (v_min - acc);
#endif
  return v_min;
}

// solve the hitting set "problem" optimally for a single set to hit
NTYPE ILB::hitting_set_special_case_1
(index_type rem,
 const bool_vec& out,
 index_type& opt)
{
  assert(rem < landmarks.size());
  opt = cost.min_cost_action(landmarks[rem]);
  return cost(opt);
}

// solve the hitting set problem optimally for two sets to hit
NTYPE ILB::hitting_set_special_case_2
(index_type rem1,
 index_type rem2,
 const bool_vec& out,
 index_set& opt)
{
  assert(rem1 < landmarks.size());
  assert(rem2 < landmarks.size());
  index_type amin1 = cost.min_cost_action(landmarks[rem1]);
  index_type amin2 = cost.min_cost_action(landmarks[rem2]);
  index_type aminc = no_such_index;
  NTYPE cminc = POS_INF;
  index_type i1 = 0;
  index_type i2 = 0;
  while ((i1 < landmarks[rem1].size()) && (i2 < landmarks[rem2].size())) {
    if (landmarks[rem1][i1] == landmarks[rem2][i2]) {
      NTYPE c = cost(landmarks[rem1][i1]);
      if (c < cminc) {
	cminc = c;
	aminc = landmarks[rem1][i1];
      }
      i1 += 1;
      i2 += 1;
    }
    else if (landmarks[rem1][i1] < landmarks[rem2][i2]) {
      i1 += 1;
    }
    else {
      i2 += 1;
    }
  }
  if ((cost(amin1) + cost(amin2)) < cminc) {
    assert(amin1 != amin2);
    opt.set_length(2);
    if (amin1 < amin2) {
      opt[0] = amin1;
      opt[1] = amin2;
    }
    else {
      opt[0] = amin2;
      opt[1] = amin1;
    }
    return cost(amin1) + cost(amin2);
  }
  else {
    assert(aminc < ins.n_actions());
    opt.assign_singleton(aminc);
    return cminc;
  }
}

NTYPE ILB::hitting_set_split
(const bool_vec& rem,
 bool_vec& out,
 bool_vec& set,
 NTYPE acc,
 index_set& best,
 NTYPE& ub,
 NTYPE lb,
 index_type depth)
{
  hs_nodes += 1;

  //bool superverbose = false; //(hs_calls > 138687);
  //DBG/ if (superverbose) {
  //DBG/   std::cerr << "enter: rem = " << rem << ", acc = " << acc
  //DBG/ 	      << ", lb = " << lb << ", ub = " << ub
  //DBG/ 	      << ", solved = " << !best.empty() << std::endl;
  //DBG/ }

  if ((acc > ub) || ((acc >= ub) && !best.empty())) {
    //DBG/ if (superverbose) {
    //DBG/   std::cerr << "return (0): " << acc << std::endl;
    //DBG/ }
    return acc;
  }

  index_type rem1 = no_such_index;
  index_type rem2 = no_such_index;
  index_type n_rem = 0;
  for (index_type i = 0; i < landmarks.size(); i++)
    if (rem[i]) {
      if (n_rem == 0) rem1 = i;
      if (n_rem == 1) rem2 = i;
      n_rem += 1;
    }

  //DBG/ if (superverbose) {
  //DBG/   std::cerr << "n_rem = " << n_rem << std::endl;
  //DBG/ }

  // if no set remains to hit, we have a better solution
  if (n_rem == 0) {
    assert((acc < ub) || ((acc <= ub) && best.empty()));
    set.copy_to(best);
    ub = acc;
    //DBG/ if (superverbose) {
    //DBG/   std::cerr << "return (00): " << acc << std::endl;
    //DBG/ }
    if (verbose_level > 2) {
      std::cerr << "new best rp cost: " << ub << std::endl;
    }
    return acc;
  }

  // check for special cases
  if (n_rem == 1) {
    assert(rem1 < landmarks.size());
    index_type new_a;
    NTYPE new_c = hitting_set_special_case_1(rem1, out, new_a);
    if (((acc + new_c) < ub) ||
	(((acc + new_c) <= ub) && best.empty())) {
      set[new_a] = true;
      set.copy_to(best);
      assert(!best.empty());
      set[new_a] = false;
      ub = (acc + new_c);
      if (verbose_level > 2) {
	std::cerr << "new best rp cost: " << ub << std::endl;
      }
    }
    //DBG/ if (superverbose) {
    //DBG/   std::cerr << "return (1): " << acc + new_c << std::endl;
    //DBG/ }
    return acc + new_c;
  }
  if (rem.size() == 2) {
    assert(rem1 < landmarks.size());
    assert(rem2 < landmarks.size());
    index_set new_s;
    NTYPE new_c = hitting_set_special_case_2(rem1, rem2, out, new_s);
    if (((acc + new_c) < ub) ||
	(((acc + new_c) <= ub) && best.empty())) {
      set.insert(new_s);
      set.copy_to(best);
      assert(!best.empty());
      set.subtract(new_s);
      ub = (acc + new_c);
      if (verbose_level > 2) {
	std::cerr << "new best rp cost: " << ub << std::endl;
      }
    }
    //DBG/ if (superverbose) {
    //DBG/   std::cerr << "return (2): " << acc + new_c << std::endl;
    //DBG/ }
    return acc + new_c;
  }

  // compute components
  index_vec irem(no_such_index, n_rem);
  index_type next = 0;
  for (index_type i = 0; i < landmarks.size(); i++)
    if (rem[i]) {
      assert(next < n_rem);
      irem[next++] = i;
    }

  equivalence e(irem.size());
  for (index_type i = 0; i < irem.size(); i++) {
    index_type ci = e.canonical(i);
    for (index_type j = i + 1; j < irem.size(); j++)
      //assert(lm_conflict(irem[i],irem[j]) ==
      //       landmarks[irem[i]].have_common_element(landmarks[irem[j]]));
      if (lm_conflict(irem[i],irem[j])) {
	//e.merge(i, j);
	index_type cj = e.canonical(j);
	e[cj] = ci;
      }
  }
  index_vec eclass;
  index_vec ecsize;
  index_type n_classes = e.class_number_and_size(eclass, ecsize);

  if (depth == 0) {
    index_type w = index_vec_util::max(ecsize, 0);
    assert(w > 0);
    sum_width += w;
    n_width += 1;
  }

  //DBG/ if (superverbose) {
  //DBG/   std::cerr << "n_classes = " << n_classes
  //DBG/ 	      << ", size = " << ecsize << std::endl;
  //DBG/ }

  // n_classes == # of eq classes
  // eclass[i] == eq class id (in 0 .. n_classes-1) of element i
  // ecsize[i] == # of elems in eq class i
  if (n_classes > 1) {
    hs_split1 += 1;
    // std::cout << "split: " << ecsize << std::endl;
    // for each component...
    bool_vec new_set(set);
    NTYPE new_acc = acc;
    cost_vec class_lb(0, n_classes);
    bool_vec class_solved(false, n_classes);
    index_type n_unsolved = n_classes;
    for (index_type k = 0; k < n_classes; k++) {
      // if it is size one or two, solve it optimally
      if (ecsize[k] == 1) {
	index_type rk = eclass.first(k);
	assert(rk < n_rem);
	index_type a_new;
	class_lb[k] = hitting_set_special_case_1(irem[rk], out, a_new);
	new_set[a_new] = true;
	class_solved[k] = true;
	n_unsolved -= 1;
	new_acc += class_lb[k];
	if ((new_acc > ub) || ((new_acc >= ub) && !best.empty())) {
	  //DBG/ if (superverbose) {
	  //DBG/   std::cerr << "return (3): " << new_acc << std::endl;
	  //DBG/ }
	  return new_acc;
	}
      }
      else if (ecsize[k] == 2) {
	index_type rk1 = eclass.first(k);
	assert(rk1 < n_rem);
	index_type rk2 = eclass.next(k, rk1);
	assert(rk2 < n_rem);
	index_set s_new;
	class_lb[k] =
	  hitting_set_special_case_2(irem[rk1], irem[rk2], out, s_new);
	new_set.insert(s_new);
	class_solved[k] = true;
	n_unsolved -= 1;
	new_acc += class_lb[k];
	if ((new_acc > ub) || ((new_acc >= ub) && !best.empty())) {
	  //DBG/ if (superverbose) {
	  //DBG/   std::cerr << "return (4): " << new_acc << std::endl;
	  //DBG/ }
	  return new_acc;
	}
      }
    }

    // if we have solved every component optimally, and not returned
    // on a bound violation, we must have a new better solution
    if (n_unsolved == 0) {
      assert((new_acc < ub) || best.empty());
      new_set.copy_to(best);
      ub = new_acc;
      //DBG/ if (superverbose) {
      //DBG/ 	std::cerr << "return (000): " << new_acc << std::endl;
      //DBG/ }
      if (verbose_level > 2) {
	std::cerr << "new best rp cost: " << ub << std::endl;
      }
      return new_acc;
    }

    // if there is only one unsolved component (i.e., only one of
    // size > 2), just adjust acc cost/set based on solutions to the
    // small components, and recurse
    if (n_unsolved == 1) {
      bool_vec r(false, landmarks.size());
      for (index_type i = 0; i < irem.size(); i++)
	if (!class_solved[eclass[i]])
	  r[irem[i]] = true;
      NTYPE est = 0;
#ifdef USE_CACHE
      set_map::iterator p = store_lb.find(r);
      if (p != store_lb.end()) {
	est = p->second;
	hs_cache_hits += 1;
      }
      else {
	est = hitting_set_lb(r, out);
	hs_cache_miss += 1;
      }
#else
      est = hitting_set_lb(r, out);
#endif
      if (((new_acc + est) > ub) ||
	  (((new_acc + est) >= ub) && !best.empty())) {
	//DBG/ if (superverbose) {
	//DBG/   std::cerr << "return (5): " << new_acc + est << std::endl;
	//DBG/ }
	return new_acc + est;
      }
      hs_branch += 1;
      NTYPE v = hitting_set_branch(r, out, new_set, new_acc, est, best, ub, lb, depth);
      //DBG/ if (superverbose) {
      //DBG/ 	std::cerr << "return (6): " << v << std::endl;
      //DBG/ }
      return v;
    }

    // otherwise, it's a real split, so we lookup/compute lower
    // bounds for each unsolved component
    NTYPE est = 0;
    for (index_type k = 0; k < n_classes; k++)
      if (!class_solved[k]) {
	bool_vec r(false, landmarks.size());
	for (index_type i = 0; i < irem.size(); i++)
	  if (eclass[i] == k)
	    r[irem[i]] = true;
#ifdef USE_CACHE
	set_map::iterator p = store_lb.find(r);
	if (p != store_lb.end()) {
	  class_lb[k] = p->second;
	  hs_cache_hits += 1;
	}
	else {
	  class_lb[k] = hitting_set_lb(r, out);
	  hs_cache_miss += 1;
	}
#else
	class_lb[k] = hitting_set_lb(r, out);
#endif
	est += class_lb[k];
	if (((new_acc + est) > ub) ||
	    (((new_acc + est) >= ub) && !best.empty())) {
	  //DBG/ if (superverbose) {
	  //DBG/   std::cerr << "return (10): " << new_acc + est << std::endl;
	  //DBG/ }
	  return (new_acc + est);
	}
      }

    hs_splits += 1;

    // for each unsolved component...
    for (index_type k = 0; k < n_classes; k++)
      if (!class_solved[k]) {
	hs_nodes += 1;
	// find the set of sets to hit in this component
	bool_vec r(false, landmarks.size());
	for (index_type i = 0; i < irem.size(); i++)
	  if (eclass[i] == k)
	    r[irem[i]] = true;
	// the upper bound on the solution for this component
	// is the upper bound for the whole collection, minus
	// new_acc and minus the lb's of all other unsolved
	// components
	NTYPE class_ub = (ub - (new_acc + est - class_lb[k]));
	assert(class_ub >= 0);
	// project the current best solution, if it exists, to the
	// component
	index_set class_best;
	if (!best.empty()) {
	  for (index_type i = 0; i < best.size(); i++)
	    if (lm_with[best[i]].have_common_element(r))
	      class_best.insert(best[i]);
	  if (cost.sum(class_best) <= class_ub) {
	    class_ub = cost.sum(class_best);
	  }
	  else {
	    // if the projected solution is not within the upper
	    // bound, it is not a solution (but that doesn't mean
	    // one doesn't exist).
	    class_best.clear();
	  }
	}
	//count_type marker = hs_nodes;
	bool_vec class_set(false, ins.n_actions());
	hs_branch += 1;
	NTYPE v = hitting_set_branch(r, out, class_set, 0, class_lb[k],
				     class_best, class_ub, class_lb[k],
				     depth);
	if (stats.break_signal_raised()) return (new_acc + est);
	// if there was no solution within the bound of this
	// component, there is no solution to the whole problem
	if (class_best.empty()) {
	  //DBG/ if (!(v > class_ub)) {
	  //DBG/   std::cerr << "pre-call mark: " << marker << std::endl;
	  //DBG/   std::cerr << "post-call mark: " << hs_calls << std::endl;
	  //DBG/   std::cerr << "r = " << r << std::endl;
	  //DBG/   std::cerr << "acc = " << new_acc << ", lbs = " << class_lb
	  //DBG/ 	      << ", solved = " << class_solved
	  //DBG/ 	      << ", est = " << est << ", lb = " << lb
	  //DBG/ 	      << ", ub = " << ub << std::endl;
	  //DBG/   std::cerr << "class #" << k << ": lb = " << class_lb[k]
	  //DBG/ 	      << ", ub = " << class_ub << std::endl;
	  //DBG/   std::cerr << "best = " << best << std::endl;
	  //DBG/   if (!best.empty()) {
	  //DBG/     for (index_type i = 0; i < best.size(); i++)
	  //DBG/ 	if (lm_with[best[i]].have_common_element(r))
	  //DBG/ 	  class_best.insert(best[i]);
	  //DBG/     std::cerr << "class initial best = " << class_best
	  //DBG/ 		<< ", cost = " << cost.sum(class_best)
	  //DBG/ 		<< std::endl;
	  //DBG/   }
	  //DBG/ }
	  assert(v > class_ub);
	  assert((new_acc + (est - class_lb[k]) + v) > ub);
	  return new_acc + (est - class_lb[k]) + v;
	}
	// else, this component is now solved
	assert(!class_best.empty());
	class_solved[k] = true;
	n_unsolved -= 1;
	est = (est - class_lb[k]);
	assert(class_ub == cost.sum(class_best));
	new_acc += class_ub;
	new_set.insert(class_best);
      }
    // if every component is solved, we have a solution to the
    // whole problem, but not necessarily a better one.
    if ((new_acc < ub) || ((new_acc <= ub) && best.empty())) {
      new_set.copy_to(best);
      assert(!best.empty());
      ub = new_acc;
      if (verbose_level > 2) {
	std::cerr << "new best rp cost: " << ub << std::endl;
      }
    }
    //DBG/ if (superverbose) {
    //DBG/   std::cerr << "return (7): " << new_acc << std::endl;
    //DBG/ }
    return new_acc;
  }

  // problem is a single component: normal lb check and branch
  NTYPE est = 0;
#ifdef USE_CACHE
  set_map::iterator p = store_lb.find(rem);
  if (p != store_lb.end()) {
    est = p->second;
    hs_cache_hits += 1;
  }
  else {
    est = hitting_set_lb(rem, out);
    hs_cache_miss += 1;
  }
#else
  est = hitting_set_lb(rem, out);
#endif
  if (((acc + est) > ub) || (((acc + est) >= ub) && !best.empty())) {
    //DBG/ if (superverbose) {
    //DBG/   std::cerr << "return (8): " << acc + est << std::endl;
    //DBG/ }
    return acc + est;
  }
  hs_branch += 1;
  NTYPE v = hitting_set_branch(rem, out, set, acc, est, best, ub, lb, depth);
  //DBG/ if (superverbose) {
  //DBG/   std::cerr << "return (9): " << v << std::endl;
  //DBG/ }
  return v;  
}

// assumption: hitting_set_branch is called with a non-empty
// collection of sets to hit.
NTYPE ILB::hitting_set_branch
(const bool_vec& rem,
 bool_vec& out,
 bool_vec& set,
 NTYPE acc,
 NTYPE est,
 index_set& best,
 NTYPE& ub,
 NTYPE lb,
 index_type depth)
{
  // choose a set to branch on: take the set with the most expensive
  // cheapest element first, tie-breaking on min cardinality; in the
  // future, one could perhaps do something smarter, aimed at creating
  // early, balanced splits...
  index_type i_branch = no_such_index;
  NTYPE co_max = -1;
  index_type ca_min = index_type_max;
  for (index_type i = 0; i < rem.size(); i++)
    if (rem[i]) {
      if (lm_min_cost[i] > co_max) {
	i_branch = i;
	co_max = lm_min_cost[i];
	ca_min = lm_branch_order[i].size();
      }
      else if ((lm_min_cost[i] == co_max) &&
	       (lm_branch_order[i].size() < ca_min)) {
	i_branch = i;
	ca_min = lm_branch_order[i].size();
      }
    }
  assert(i_branch < landmarks.size());

  bool_vec rem_copy(rem);
  NTYPE v_min = POS_INF;
  for (index_type k = 0; k < lm_branch_order[i_branch].size(); k++) {
    index_type a = lm_branch_order[i_branch][k];
    set[a] = true;
    for (index_type j = 0; j < lm_with[a].size(); j++)
      rem_copy[lm_with[a][j]] = false;
    NTYPE v = hitting_set_split(rem_copy, out, set, acc + cost(a), best, ub, lb, depth + 1);
    if (stats.break_signal_raised()) return (acc + est);
    set[a] = false;
    if ((v == lb) && !best.empty()) return v;
    v_min = MIN(v_min, v);
    if ((v_min == (acc + est)) && !best.empty()) {
      return v_min;
    }
    if ((k + 1) < lm_branch_order[i_branch].size())
      rem_copy.assign_copy(rem);
  }
#ifdef USE_CACHE
  if ((v_min - acc) > est)
    store_lb[rem] = (v_min - acc);
#endif
  return v_min;
}

NTYPE ILB::apx_hitting_set(index_set& hs, NTYPE bound)
{
  bool_vec rem_to_hit(true, landmarks.size());
  // if input hs is non-empty, clear lms already hit.
  for (index_type i = 0; i < hs.size(); i++)
    for (index_type k = 0; k < lm_with[hs[i]].size(); k++)
      rem_to_hit[lm_with[hs[i]][k]] = false;
  index_type n_rem = rem_to_hit.count(true);
  index_vec n_hit(0, ins.n_actions());
  NTYPE sum = cost.sum(hs);
  while (n_rem > 0) {
    n_hit.assign_value(0, ins.n_actions());
    for (index_type k = 0; k < landmarks.size(); k++)
      if (rem_to_hit[k])
    	for (index_type i = 0; i < landmarks[k].size(); i++)
    	  n_hit[landmarks[k][i]] += 1;
    // for (index_type i = 0; i < ins.n_actions(); i++)
    //   if (relevant_actions[i]) {
    // 	n_hit[i] = lm_with[i].count_common(rem_to_hit);
    //   }
    index_type i_best = no_such_index;
    NTYPE w_best = POS_INF;
    index_type n_best = 0;
    for (index_type i = 0; i < ins.n_actions(); i++)
      if (relevant_actions[i]) {
	NTYPE w = cost(i) / n_hit[i];
	if ((w < w_best) || (IS_ZERO(w) && (n_hit[i] > n_best))) {
	  i_best = i;
	  w_best = w;
	  n_best = n_hit[i];
	}
      }
    assert((i_best != no_such_index) && (i_best < ins.n_actions()));
    hs.insert(i_best);
    sum += cost(i_best);
    if (sum >= bound) return POS_INF;
    for (index_type k = 0; k < lm_with[i_best].size(); k++)
      rem_to_hit[lm_with[i_best][k]] = false;
    index_type r = rem_to_hit.count(true);
    // if (r >= n_rem) {
    //   std::cerr << "to hit = " << rem_to_hit << std::endl;
    //   std::cerr << "n_rem = " << n_rem << std::endl;
    //   std::cerr << "i_best = " << i_best << std::endl;
    //   std::cerr << "lm_with[i_best] = " << lm_with[i_best] << std::endl;
    //   std::cerr << "w_best = " << w_best << std::endl;
    //   std::cerr << "n_best = " << n_best << std::endl;
    // }
    assert(r < n_rem);
    n_rem = r;
  }
  assert(sum == cost.sum(hs));
  return sum;
}

#ifdef USE_SCIP
NTYPE ILB::hitting_set_extern
(index_set& best,
 NTYPE& ub,
 NTYPE lb)
{
  // set of active actions (vars in hs problem) and mapping from
  // instance action id to active action number.
  index_vec aa;
  index_vec amap(no_such_index, ins.n_actions());
  for (index_type k = 0; k < lm_actions.size(); k++)
    if (!dominated.contains(lm_actions[k])) {
      amap[lm_actions[k]] = aa.size();
      aa.append(lm_actions[k]);
    }
  SCIP* scip;
  SCIPcreate(&scip);
  if (verbose_level < 2) {
    SCIPsetMessagehdlrQuiet(scip, TRUE);
  }
  SCIPincludeDefaultPlugins(scip);
  SCIP_RETCODE res;
  // create the problem
  res = SCIPcreateProbBasic(scip, "hitting set");
  assert(res == SCIP_OKAY);
  SCIP_VAR* vars[aa.size() + 1];
  for (index_type k = 0; k < aa.size(); k++) {
    double w = N_TO_D(ins.actions[aa[k]].cost);
    SCIP_VAR* v;
    res = SCIPcreateVarBasic(scip, &v, NULL, 0, 1, w, SCIP_VARTYPE_BINARY);
    assert(res == SCIP_OKAY);
    res = SCIPaddVar(scip, v);
    assert(res == SCIP_OKAY);
    vars[k] = v;
  }
  //std::cerr << "done creating variables..." << std::endl;
  SCIP_VAR* vs[aa.size() + 1];
  for (index_type k = 0; k < landmarks.size(); k++) {
    int n = 0;
    for (index_type i = 0; i < landmarks[k].size(); i++) {
      index_type p = amap[landmarks[k][i]];
      // redundant check:
      //index_type q = aa.first(landmarks[k][i]);
      //assert(p == q);
      if (p != no_such_index) {
	vs[n] = vars[p];
	n += 1;
      }
    }
    assert(n > 0);
    SCIP_CONS* c;
    res = SCIPcreateConsBasicLogicor(scip, &c, "lm", n, vs);
    assert(res == SCIP_OKAY);
    res = SCIPaddCons(scip, c);
    assert(res == SCIP_OKAY);
    SCIPreleaseCons(scip, &c);
  }
  // done creating the problem; tell the solver that
  res = SCIPtransformProb(scip);
  assert(res == SCIP_OKAY);
  // solve
  res = SCIPsolve(scip);
  assert(res == SCIP_OKAY);
  SCIP_Status status = SCIPgetStatus(scip);
  if (status != SCIP_STATUS_OPTIMAL) {
    std::cerr << "error: not solved to optimality!" << std::endl;
    exit(1);
  }
  // extract solution
  SCIP_SOL* sol = SCIPgetBestSol(scip);
  assert(sol != NULL);
  //SCIPprintSol(scip, sol, NULL, FALSE);
  double obj_val = SCIPgetSolOrigObj(scip, sol);
  NTYPE val = D_TO_N(obj_val);
  if (verbose_level > 0) {
    std::cerr << "obj val = " << obj_val << " (R), " << val << " (Q)"
	      << std::endl;
  }
#ifdef NTYPE_RATIONAL
  val = rational::ceil_to(val, ac_div);
#endif
  //assert(SAFE_EQ(val, lb) || SAFE_GT(val, lb));
  //assert(SAFE_EQ(ub, val) || SAFE_GT(ub, val));
  NTYPE sum = 0;
  best.clear();
  for (index_type k = 0; k < aa.size(); k++) {
    double v = SCIPgetSolVal(scip, sol, vars[k]);
    if (v > 0.5) {
      sum += ins.actions[aa[k]].cost;
      best.insert(aa[k]);
    }
  }
  //assert(SAFE_EQ(sum, val));
  // free memory
  for (index_type k = 0; k < aa.size(); k++)
    SCIPreleaseVar(scip, &(vars[k]));
  SCIPfreeTransform(scip);
  SCIPfree(&scip);
  // return optimal cost
  return sum;
}
#endif

#ifdef USE_CPLEX
#ifdef CPLEX_INCREMENTAL

/// hitting-set solver based on incremental CPLEX

NTYPE ILB::hitting_set_extern
(index_set& best,
 NTYPE& ub,
 NTYPE lb)
{
  assert(ins.n_actions() > 0);

  // first call, set env parameters
  if (cplex_model_actions == 0) {
    if (verbose_level < 2) {
      cplex_env.setOut(cplex_env.getNullStream());
      if (verbose_level < 1) {
	cplex_env.setWarning(cplex_env.getNullStream());
      }
    }
  }

  // if instance has changed since last call, reset the model
  if ((cplex_model_actions > 0) &&
      ((cplex_model_actions != ins.n_actions()) ||
       (cplex_model_atoms != ins.n_atoms()))) {
    cplex_solver.end();
    cplex_model.end();
    cplex_vars.end();
    cplex_model_actions = 0;
    cplex_model_atoms = 0;
  }

  // create the model if necessary
  if (cplex_model_actions == 0) {
    cplex_model = IloModel(cplex_env);
    cplex_vars = IloBoolVarArray(cplex_env, ins.n_actions());
    IloNumArray costs(cplex_env);
    for (index_type k = 0; k < ins.n_actions(); k++) {
      costs.add(N_TO_D(ins.actions[k].cost));
    }
    cplex_model.add(cplex_vars);
    cplex_model.add(IloObjective(cplex_env, IloScalProd(costs, cplex_vars)));
    cplex_model_actions = ins.n_actions();
    cplex_model_atoms = ins.n_atoms();
    cplex_model_landmarks = 0;
    cplex_solver = IloCplex(cplex_model);
    cplex_solver.setParam(IloCplex::Threads, 1);
  }

  // add new landmark constraints
  for (index_type k = cplex_model_landmarks; k < landmarks.size(); k++) {
    IloExpr expr(cplex_env);
    for (index_type i = 0; i < landmarks[k].size(); i++) {
      IloBoolVar &var = cplex_vars[landmarks[k][i]];
      expr += var;
    }
    IloRange con(cplex_env, 1, expr);
    cplex_model.add(con);
    cplex_model_landmarks = landmarks.size();
  }

  // solve
  NTYPE sum = 0;
  try {
    cplex_solver.solve();
    if (cplex_solver.getStatus() != IloAlgorithm::Optimal) {
      std::cerr << "error: not solved to optimality!" << std::endl;
      std::cerr << "cplex status = " << cplex_solver.getStatus() << std::endl;
      exit(1);
    }

    // double obj_val = cplex_solver.getObjValue();
    // NTYPE val = D_TO_N(obj_val);
    // if (verbose_level > 0) {
    //   std::cerr << "obj val = " << obj_val << " (R), " << val << " (Q)"
    // 		<< std::endl;
    // }
    // #ifdef NTYPE_RATIONAL
    // val = rational::ceil_to(val, ac_div);
    // #endif
    // assert(SAFE_EQ(val, lb) || SAFE_GT(val, lb));
    // assert(SAFE_EQ(ub, val) || SAFE_GT(ub, val));
    best.clear();
    for (index_type k = 0; k < ins.n_actions(); k++) {
      double v = cplex_solver.getValue(cplex_vars[k]);
      if (v > 0.5) {
	sum += ins.actions[k].cost;
	best.insert(k);
      }
    }
    // assert(SAFE_EQ(sum, val));
  }
  catch (IloException& e) {
    std::cerr << "cplex exception: " << e << std::endl;
    exit(1);
  }
  catch (...) {
    std::cerr << "unknown exception" << std::endl;
    exit(1);
  }

  if (verbose_level > 0) {
    std::cerr << "val = " << sum << std::endl;
  }

  // return optimal cost
  return sum;
}

#else // !CPLEX_INCREMENTAL

/// hitting-set solver based on non-incremental CPLEX

NTYPE ILB::hitting_set_extern
(index_set& best,
 NTYPE& ub,
 NTYPE lb)
{
  // set of active actions (vars in hs problem) and mapping from
  // instance action id to active action number.
  exit(0);
  index_vec aa;
  index_vec amap(no_such_index, ins.n_actions());
  for (index_type k = 0; k < lm_actions.size(); k++)
    if (!dominated.contains(lm_actions[k])) {
      amap[lm_actions[k]] = aa.size();
      aa.append(lm_actions[k]);
    }

  // have to declare sum outside the try-catch scope
  NTYPE sum = 0;

  IloEnv env;
  try {
    if (verbose_level < 2) {
      env.setOut(env.getNullStream());
      if (verbose_level < 1) {
	env.setWarning(env.getNullStream());
      }
    }
    IloModel model(env);

    IloBoolVarArray vars(env, aa.size());
    IloNumArray costs(env);
    for (index_type k = 0; k < aa.size(); k++) {
      costs.add(N_TO_D(ins.actions[aa[k]].cost));
    }
    model.add(vars);
    model.add(IloObjective(env, IloScalProd(costs, vars))); 

    for (index_type k = 0; k < landmarks.size(); k++) {
      IloExpr expr(env);
      for (index_type i = 0; i < landmarks[k].size(); i++)
	if (amap[landmarks[k][i]] != no_such_index) {
	  assert(amap[landmarks[k][i]] < aa.size());
	  IloBoolVar &var = vars[amap[landmarks[k][i]]];
	  expr += var;
	}
      IloRange con(env, 1, expr);
      model.add(con);
    }

    IloCplex cplex(model);
    cplex.setParam(IloCplex::Threads, 1);
    cplex.solve();
    if (cplex.getStatus() != IloAlgorithm::Optimal) {
      std::cerr << "error: not solved to optimality!" << std::endl;
      std::cerr << "cplex status = " << cplex.getStatus() << std::endl;
      exit(1);
    }

    // double obj_val = cplex.getObjValue();
    // NTYPE val = D_TO_N(obj_val);
    // if (verbose_level > 0) {
    //   std::cerr << "obj val = " << obj_val << " (R), " << val << " (Q)"
    // 		<< std::endl;
    // }
    // #ifdef NTYPE_RATIONAL
    // val = rational::ceil_to(val, ac_div);
    // #endif
    // assert(SAFE_EQ(val, lb) || SAFE_GT(val, lb));
    // assert(SAFE_EQ(ub, val) || SAFE_GT(ub, val));
    best.clear();
    for (index_type k = 0; k < aa.size(); k++) {
      double v = cplex.getValue(vars[k]);
      if (v > 0.5) {
	sum += ins.actions[aa[k]].cost;
	best.insert(aa[k]);
      }
    }
    // assert(SAFE_EQ(sum, val));
  }
  catch (IloException& e) {
    std::cerr << "cplex exception: " << e << std::endl;
    exit(1);
  }
  catch (...) {
    std::cerr << "unknown exception" << std::endl;
    exit(1);
  }
  env.end();

  if (verbose_level > 0) {
    std::cerr << "val = " << sum << std::endl;
  }

  // return optimal cost
  return sum;
}
#endif // #else
#endif // USE_CPLEX

void ILB::dump_hitting_set_problem_wcnf(NTYPE lb, NTYPE ub, std::ostream& to)
{
  to << "c optimal hitting set problem #" << calls_to_hs_opt << std::endl;
  to << "c lower bound = " << PRINT_NTYPE(lb)
     << ", upper bound = " << PRINT_NTYPE(ub) << std::endl;
  // construct an index to name the actions used in the hsp 1,2,...
  index_set els(lm_actions);
  els.subtract(dominated);
  index_vec rn(no_such_index, ins.n_actions());
  NTYPE sum_el_cost = 0;
  for (index_type i = 0; i < els.size(); i++) {
    assert(els[i] < ins.n_actions());
    rn[els[i]] = i + 1;
    assert(INTEGRAL(ins.actions[els[i]].cost));
    sum_el_cost += ins.actions[els[i]].cost;
  }
  NTYPE hard_cost = sum_el_cost + 1;
  to << "p wcnf " << els.size()
     << " " << landmarks.size() + els.size()
     << " " << hard_cost << std::endl;
  for (index_type k = 0; k < landmarks.size(); k++) {
    to << hard_cost << " ";
    for (index_type i = 0; i < lm_branch_order[k].size(); i++) {
      index_type v = rn[lm_branch_order[k][i]];
      assert(v <= els.size());
      to << v << " ";
    }
    to << "0" << std::endl;
  }
  for (index_type k = 0; k < els.size(); k++) {
    to << ins.actions[els[k]].cost << " -" << (k + 1) << " 0" << std::endl;
  }
}

void ILB::dump_hitting_set_problem_lp(NTYPE lb, NTYPE ub, std::ostream& to)
{
  to << "\\ optimal hitting set problem #" << calls_to_hs_opt << std::endl;
  to << "\\ lower bound = " << PRINT_NTYPE(lb)
     << ", upper bound = " << PRINT_NTYPE(ub) << std::endl;
  to << "Minimize" << std::endl << " cost :";
  for (index_type k = 0; k < lm_actions.size(); k++)
    if (!dominated.contains(lm_actions[k])) {
      to << " + " << PRINT_NTYPE(ins.actions[lm_actions[k]].cost)
	 << " a" << lm_actions[k];
    }
  to << std::endl << "subject to" << std::endl;
  for (index_type k = 0; k < landmarks.size(); k++)
    if (!dominated.contains(landmarks[k])) {
      to << " lm" << k << " :";
      for (index_type i = 0; i < landmarks[k].size(); i++)
	if (!dominated.contains(landmarks[k][i])) {
	  to << " + a" << landmarks[k][i];
	}
      to << " >= 1" << std::endl;
    }
  to << "binary" << std::endl;
  for (index_type k = 0; k < lm_actions.size(); k++)
    if (!dominated.contains(lm_actions[k])) {
      to << " a" << lm_actions[k] << std::endl;
    }
}

// check if action a is dominated by any action in set s w.r.t.
// the current landmark collection; action a is dominated by action
// b iff every landmark that contains a also contains b, and cost(b)
// <= cost(a); in case of a tie (equal set of landmarks, equal cost),
// the lower-indexed action dominates the higher-indexed one.
bool ILB::is_dominated(index_type a, const index_set& s) const
{
  for (index_type j = 0; j < s.size(); j++)
    if ((s[j] != a) && (cost(s[j]) <= cost(a)) &&
	lm_with[s[j]].contains(lm_with[a])) {
      assert(lm_with[s[j]].size() >= lm_with[a].size());
      if (lm_with[s[j]].size() < lm_with[a].size())
	return true;
      else if (cost(s[j]) < cost(a))
	return true;
      else if (s[j] < a)
	return true;
    }
  return false;
}

// add a new landmark to the collection, and update some stuff
index_type ILB::add_new_landmark(const index_set& new_lm)
{
  landmarks.append(new_lm);
  lm_min_cost.append(cost.min_cost(new_lm));
  assert(lm_min_cost.size() == landmarks.size());
  lm_actions.insert(new_lm);
  index_type new_i = landmarks.size() - 1;
  for (index_type k = 0; k < new_lm.size(); k++)
    lm_with[new_lm[k]].insert(new_i);
#ifndef ILA_USE_HS_EXTERN
  lm_conflict.extend_to(landmarks.size());
  lm_conflict_set.append(EMPTYSET);
  assert(lm_conflict_set.size() == landmarks.size());
  for (index_type k = 0; k < new_i; k++)
    if (landmarks[k].have_common_element(landmarks[new_i])) {
      lm_conflict.set(new_i, k);
      lm_conflict_set[k].insert(new_i);
      lm_conflict_set[new_i].insert(k);
    }
#endif
#ifdef USE_LB1
  update_weighted_degree();
#endif
  action_cost_increasing o(cost);
  // check if the new landmark breaks some dominance relationship
  index_type i = 0;
  while (i < dominated.size()) {
    if (new_lm.contains(dominated[i]))
      if (!is_dominated(dominated[i], new_lm)) {
	for (index_type k = 0; k < lm_with[dominated[i]].size(); k++) {
	  index_type j = lm_with[dominated[i]][k];
	  if (j < new_i) {
	    assert(landmarks[j].contains(dominated[i]));
#ifndef ILA_USE_HS_EXTERN
	    assert(!lm_branch_order[j].contains(dominated[i]));
	    lm_branch_order[j].insert_ordered(dominated[i], o);
#endif
	  }
	}
	dominated.remove(i);
      }
      else i += 1;
    else i += 1;
  }
  index_vec b;
  for (index_type i = 0; i < new_lm.size(); i++) {
    // if new_lm[i] is dominated, a dominating action must be in new_lm.
    if (hitting_set_use_dominance)
      if (is_dominated(new_lm[i], new_lm))
	dominated.insert(new_lm[i]);
#ifndef ILA_USE_HS_EXTERN
      else
	b.insert_ordered(new_lm[i], o);
    else
      b.insert_ordered(new_lm[i], o);
#endif
  }
#ifndef ILA_USE_HS_EXTERN
  lm_branch_order.append(b);
  assert(lm_branch_order.size() == landmarks.size());
#endif
  return new_i;
}

// b is a superset of a
void ILB::add_redundant_actions_to_sets(bool_vec& a, bool_vec& b)
{
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (!b[k])
      if (!reach.unreachable(ins.actions[k].add)) {
	a[k] = true;
	b[k] = true;
      }
}

index_type ILB::make_new_landmark_ST
(const index_set& init, const index_set& goal, const index_set& rp_set)
{
  calls_to_newlm += 1;

  // update scores
#ifdef USE_SCORE
  for (index_type i = 0; i < ins.n_actions(); i++)
    score[i] = (score[i] / 2);
  for (index_type i = 0; i < lm_actions.size(); i++)
    score[i] += 1;
  for (index_type i = 0; i < rp_set.size(); i++)
    score[i] += 8;
#endif

  // build the vector of actions to test
  bool_vec in(rp_set, ins.n_actions());
  bool_vec out(relevant_actions);
  out.complement();
  out.insert(rp_set);

  // intialise reach computation (before testing redudancy)
  reach.recompute(init, in);

  // // restrict search to executable landmarks
  // for (index_type i = 0; i < ins.n_actions(); i++)
  //   if (reach.unreachable_action(i)) {
  //     in[i] = true;
  //     out[i] = true;
  //   }

  add_redundant_actions_to_sets(in, out);

  index_vec acts_to_test;
  acts_to_test.reserve(ins.n_actions());
  for (index_type i = 0; i < ins.n_actions(); i++)
    if (!out[i])
      acts_to_test.append(i);
#ifdef USE_SCORE
  index_vec::decreasing_by_value o(score);
  std::sort(acts_to_test.begin(), acts_to_test.end(), o);
#endif
#ifdef RANDOMIZE
  rng.permute(acts_to_test);
#endif

  index_set new_lm;
  index_type n = 0;
  Reachability::reachability_state s;

  while (n < acts_to_test.size()) {
    if (stats.break_signal_raised())
      return no_such_index;
    if (!out[acts_to_test[n]]) {
      out[acts_to_test[n]] = true;
      // if the actions preconditions are unreachable, then addding
      // this action will not change reachability of any atom
      if (reach.unreachable_action(acts_to_test[n])) {
	in[acts_to_test[n]] = true;
      }
      // otherwise, we have to check goal reachability
      else {
	in[acts_to_test[n]] = true;
	reach.save_state(s);
	reach.update(in, acts_to_test[n]);
	if (!reach.unreachable(goal)) {
	  in[acts_to_test[n]] = false;
	  new_lm.insert(acts_to_test[n]);
	  // if the test succeeds, we have to restore the reachability
	  // state after removing the action from set 'in'
	  reach.restore_state(s);
	  //reach.recompute(init, in);
	}
	else {
	  add_redundant_actions_to_sets(in, out);
	}
      }
    }
    n += 1;
  }
  assert(!new_lm.empty());
  return add_new_landmark(new_lm);
}

///
// simplifed version of landmark generation, with extra debug checks
///
//
// index_type ILB::make_new_landmark_ST
// (const index_set& init, const index_set& goal, const index_set& rp_set)
// {
//   calls_to_newlm += 1;
//   bool_vec in(rp_set, ins.n_actions());
//   index_set new_lm;
//   bool_vec s0(init, ins.n_atoms());
//   if (!h1) {
//     h1 = new CostTable(ins, h1_stats);
//   }
//   for (index_type k = 0; k < ins.n_actions(); k++)
//     if (!in[k]) {
//       in[k] = true;
//       reach.recompute(init, in);
//       h1->compute_H1(UnitACF(), s0, &in);
//       if (!reach.unreachable(goal)) {
// 	assert(FINITE(h1->eval(goal)));
// 	in[k] = false;
// 	new_lm.insert(k);
//       }
//       else {
// 	if (FINITE(h1->eval(goal))) {
// 	  std::cerr << "init = " << init << std::endl;
// 	  std::cerr << "in = " << index_set(in) << std::endl;
// 	  std::cerr << "reach = " << reach.reachable() << std::endl;
// 	  for (index_type i = 0; i < ins.n_atoms(); i++) {
// 	    bool x = reach.reachable()[i];
// 	    bool y = !reach.unreachable(i);
// 	    NTYPE z1 = h1->eval(i);
// 	    bool z = z1.finite();
// 	    std::cerr << i << " : " << x << " " << y
// 		      << " " << z1 << " " << z << std::endl;
// 	  }
// 	  std::cerr << "reach = " << reach.reachable() << std::endl;
// 	  assert(0);
// 	}
//       }
//     }
//   assert(!new_lm.empty());
//   return add_new_landmark(new_lm);
// }

void ILB::select_pcf(const index_set& g, bool_vec& rel)
{
  if (g.empty()) return;
  if (!reach.unreachable(g)) {
    std::cerr << "entered select with g = " << g << " which is reachable"
	      << std::endl;
  }
  // pick an unreachable subgoal in set g: this can be done in many ways.
  // here we take the one that has the fewest adders in common with existing
  // landmarks, adds the fewest new relevant actions, and tie-break on
  // estimated h^1 (unit) cost (preferring smaller)
  index_type g_best = no_such_index;
  index_type o_best = 0;
  index_type r_best = 0;
  NTYPE c_best = 0;
  for (index_type k = 0; k < g.size(); k++)
    if (reach.unreachable(g[k])) {
      if (g_best == no_such_index) {
	g_best = k;
	o_best = ins.atoms[g[k]].add_by.count_common(lm_actions);
	r_best = (ins.atoms[g[k]].add_by.size() -
		  ins.atoms[g[k]].add_by.count_common(rel));
	c_best = h1->eval(g[k]);
      }
      else {
	index_type o_k = ins.atoms[g[k]].add_by.count_common(lm_actions);
	index_type r_k = (ins.atoms[g[k]].add_by.size() -
			  ins.atoms[g[k]].add_by.count_common(rel));
	NTYPE c_k = h1->eval(g[k]);
	if ((o_k < o_best) || ((o_k == o_best) && (r_k < r_best)) ||
	    ((o_k == o_best) && (r_k == r_best) && (c_k < c_best))) {
	  g_best = k;
	  r_best = r_k;
	  o_best = o_k;
	  c_best = c_k;
	}
      }
    }
  assert(g_best != no_such_index);
  assert(g_best < g.size());
  // all actions that add this subgoal are relevant; for each one that
  // isn't already in the relevant set, we must add it, and if some of
  // it's preconditions are not (all) reachable, recurse
  // to minimise the relevant set, we go through all adders and
  // (greedily) cover their unreached preconditions with a smaller
  // number of sets; in particular, if there's an unreached precond
  // that is common to all adders, we should only recurse on that
  index_set_vec precs;
  for (index_type k = 0; k < ins.atoms[g[g_best]].add_by.size(); k++)
    if (!rel[ins.atoms[g[g_best]].add_by[k]]) {
      rel[ins.atoms[g[g_best]].add_by[k]] = true;
      index_set upre;
      reach.unreachable_subset
	(ins.actions[ins.atoms[g[g_best]].add_by[k]].pre, upre);
      if (!upre.empty()) {
	bool placed = false;
	for (index_type i = 0; (i < precs.size()) && !placed; i++)
	  if (precs[i].have_common_element(upre)) {
	    precs[i].intersect(upre);
	    placed = true;
	  }
	if (!placed)
	  precs.append(upre);
      }
    }
  for (index_type i = 0; i < precs.size(); i++)
    select_pcf(precs[i], rel);
}

index_type ILB::make_new_landmark_BC
(const index_set& init, const index_set& goal, const index_set& rp_set)
{
  calls_to_newlm += 1;
  // // saturation idea: first, check reachability with all
  // // current landmark actions; if the goal remains unreachable,
  // // we'll find a better landmark using this
  // bool_vec in(lm_actions, ins.n_actions());
  // reach.recompute(init, in);
  // // however, if the goal is reachable using all lm actions,
  // // we have to try with only the current rp set.
  // if (!reach.unreachable(goal)) {
  //   //std::cerr << "not good!" << std::endl;
  //   bool_vec in(rp_set, ins.n_actions());
  //   reach.recompute(init, in);
  // }

  // backchain from the (unsatisfied) goals to produce a set
  // of relevant actions
  bool_vec cut(false, ins.n_actions());
  select_pcf(goal, cut);

  // the landmark is the set of relevant actions whose preconds
  // are reachable.
  index_set new_lm;
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (cut[k] && !reach.unreachable(ins.actions[k].pre))
      new_lm.insert(k);

  assert(!new_lm.empty());
  return add_new_landmark(new_lm);
}

NTYPE ILB::initialise_with_lmcut
(const index_set& init, const index_set& goal)
{
  bool_vec s0(init, ins.n_atoms());
  AnyACF costs(ins.n_actions(), cost);
  h1->compute_H1(costs, s0, &relevant_actions);
  if (INFINITE(h1->eval(goal)))
    return POS_INF;
  NTYPE total_cost = 0;
  while (h1->eval(goal) > 0) {
    if (stats.break_signal_raised()) return 0;
    bool_vec ext_goal_set(false, ins.n_atoms());
    h1->extend_goal_set(goal, costs, ext_goal_set);
    bool_vec allowed_acts(relevant_actions);
    for (index_type k = 0; k < ins.n_actions(); k++)
      if (allowed_acts[k] &&
	  ins.actions[k].pre.have_common_element(ext_goal_set))
	allowed_acts[k] = false;
    h1->compute_H1(costs, s0, &allowed_acts);
    bool_vec cut(false, ins.n_actions());
    h1->find_cut(ext_goal_set, costs, cut);
    cut.intersect(allowed_acts);
    NTYPE c_cut = costs.min_cost(cut);
    assert(c_cut > 0);
    total_cost += c_cut;
    costs.decrease(cut, c_cut);
    index_set new_lm(cut);
    add_new_landmark(new_lm);
    h1->compute_H1(costs, s0);
  }
  return total_cost;
}

NTYPE ILB::compute_relaxed_plan_ILA
(const index_set& init, const index_set& goal, NTYPE bound)
{
#ifdef CPLEX_INCREMENTAL
  cplex_model_atoms = 0;
  cplex_model_landmarks = 0;
#endif
  // initialisation: check if goal is reachable with all (relevant) actions
  reach.recompute(init, relevant_actions);
  if (reach.unreachable(goal)) {
    return POS_INF;
  }
  // initialise landmark collection (empty) and hitting set (empty)
  landmarks.clear();
  lm_actions.clear();
  lm_min_cost.clear();
  lm_with.assign_value(EMPTYSET, ins.n_actions());
  lm_conflict.clear();
  lm_conflict_set.clear();
  lm_branch_order.clear();
  dominated.clear();
#ifdef USE_CACHE
  store_lb.clear();
#endif
#ifdef USE_SCORE
  score.assign_value(0, ins.n_actions());
#endif
  index_set all_rp_actions;
  rp_actions.clear();
  index_set zero_cost_actions;
  if (zero_cost_fill) {
    for (index_type k = 0; k < ins.n_actions(); k++)
      if (IS_ZERO(cost(k)) && relevant_actions[k])
	zero_cost_actions.insert(k);
    rp_actions.insert(zero_cost_actions);
    all_rp_actions.insert(zero_cost_actions);
  }
  rp_cost = 0;
  // allocate some variables used by the optimal hitting set algorithm
  bool_vec lm_to_hit(true, landmarks.size());
  bool_vec out(false, ins.n_actions());
  bool_vec set(false, ins.n_actions());
  NTYPE v_max = 0;
  bool init_phase =
    (ILA_use_saturation && ILA_use_approximate && !ILA_use_lmcut);
  if (ILA_use_lmcut) {
    NTYPE v = initialise_with_lmcut(init, goal);
    lm_to_hit.assign_value(true, landmarks.size());
    v_max = hitting_set_lb(lm_to_hit, out);
    assert(v <= v_max);
#ifdef TRACE_PRINT_HLB
    new_lb(v_max);
#endif
    hsa_stats.start();
    calls_to_hs_apx += 1;
    hs_apx_improve += 1;
    NTYPE c = apx_hitting_set(rp_actions);
    rp_cost = cost.sum(rp_actions);
    hsa_stats.stop();
    // need to compute an initial rp over the landmark set...
  }
  // recompute reachability using only hitting set actions
  reach.recompute(init, rp_actions);
  // as long as goal is not reachable...
  bool reached_goal = !reach.unreachable(goal);
  while (!reached_goal) {
    if (stats.break_signal_raised()) return v_max;
    if (v_max >= bound) return v_max;
    if (verbose_level > 1) {
      std::cerr << landmarks.size() << " landmarks (d: "
		<< dominated.size() << "/" << lm_actions.size()
		<< "), rp size = "
		<< (ILA_use_saturation ?
		    all_rp_actions.size() : rp_actions.size())
		<< ", cost = [" << PRINT_NTYPE(v_max)
		<< ", " << PRINT_NTYPE(rp_cost) << "]"
		<< std::endl;
      if (verbose_level > 2) {
	std::cerr << "partial rp = ";
	std::cerr << rp_actions;
	// ins.write_action_set(std::cerr << " = ", rp_actions);
	std::cerr << std::endl;
      }
    }
    // find a new landmark
    lm_stats.start();
#ifdef USE_NEWLM_BC
    index_type new_lm =
      make_new_landmark_BC
      (init, goal, (ILA_use_saturation ? all_rp_actions : rp_actions));
#else
    index_type new_lm =
      make_new_landmark_ST
      (init, goal, (ILA_use_saturation ? all_rp_actions : rp_actions));
#endif
    lm_stats.stop();
    if (stats.break_signal_raised()) return v_max;
    assert(new_lm < landmarks.size());
    if (verbose_level > 2) {
      std::cerr << "new landmark: ";
      std::cerr << landmarks[new_lm];
      if (verbose_level > 3)
	ins.write_action_set(std::cerr << " = ", landmarks[new_lm]);
      std::cerr << std::endl;
    }
    lm_to_hit.append(true);
    out.assign_value(false, ins.n_actions());
    assert(lm_to_hit.size() == landmarks.size());
    if (stats.break_signal_raised()) return v_max;
    // in "initialisation phase", use all actions appearing in any
    // landmark as the hitting set
    if (init_phase) {
      all_rp_actions.insert(landmarks[new_lm]);
      reach.recompute(init, all_rp_actions);
      if (!reach.unreachable(goal)) {
	if (verbose_level > 0) {
	  std::cerr << "finished init phase with " << landmarks.size()
		    << " landmarks" << std::endl;
	}
	rp_cost = cost.sum(all_rp_actions);
	hss_stats.start();
	calls_to_hs_opt += 1;
#ifdef ILA_USE_HS_EXTERN
	NTYPE new_v = hitting_set_extern(all_rp_actions, rp_cost, v_max);
	if (zero_cost_fill)
	  all_rp_actions.insert(zero_cost_actions);
#else
	if (zero_cost_fill)
	  set.insert(zero_cost_actions);
	NTYPE new_v =
	  (hitting_set_use_split ?
	   hitting_set_split(lm_to_hit, out, set, 0, all_rp_actions, rp_cost, v_max, 0)
	   : hitting_set(lm_to_hit, out, set, 0, all_rp_actions, rp_cost, v_max));
#endif
	hss_stats.stop();
	if (stats.break_signal_raised()) return v_max;
	v_max = MAX(v_max, new_v);
#ifdef TRACE_PRINT_HLB
	new_lb(v_max);
#endif
	rp_actions.assign_copy(all_rp_actions);
#ifdef MEASURE_STABILITY
	last_rp.assign_copy(rp_actions);
#endif
	reach.recompute(init, all_rp_actions);
	if (!reach.unreachable(goal))
	  reached_goal = true;
	init_phase = false;
      }
    }
    else {
      // find a new hitting set:
      // as an initial ub, pick the action in the new lm with the
      // smallest cost, and add it to the current hitting set
      index_type new_a = cost.min_cost_action(landmarks[new_lm]);
      if (zero_cost_fill) {
	assert(cost(new_a) > 0);
      }
      rp_actions.insert(new_a);
      rp_cost = cost.sum(rp_actions);
      if (ILA_use_approximate) {
#ifdef ILA_USE_2ND_APX_HS
	hsa_stats.start();
	calls_to_hs_apx += 1;
	index_set rp2;
	if (zero_cost_fill) {
	  rp2.assign_copy(zero_cost_actions);
	}
	NTYPE c2 = apx_hitting_set(rp2, rp_cost);
	if (c2 < rp_cost) {
	  hs_apx_improve += 1;
	  // std::cerr << "rp_cost: " << rp_cost << " -> " << c2
	  // 	    << std::endl << "rp1 = " << rp_actions
	  // 	    << std::endl << "rp2 = " << rp2
	  // 	    << std::endl;
	  rp_actions.assign_copy(rp2);
	  rp_cost = c2;
	}
	hsa_stats.stop();
#endif
	NTYPE new_v = hitting_set_lb(lm_to_hit, out);
	assert(new_v <= rp_cost);
	v_max = MAX(v_max, new_v);
#ifdef TRACE_PRINT_HLB
	new_lb(v_max);
#endif
      }
      else {
	// find an optimal (min cost) hitting set:
	if (verbose_level > 2) {
	  std::cerr << "initial rp cost: " << rp_cost
		    << ", lb = " << v_max << std::endl;
	}
#ifdef MEASURE_DENSITY
	long n = 0;
	for (index_type k = 0; k < landmarks.size(); k++)
	  n += lm_branch_order[k].size();
	double d = n/(double)((lm_actions.size() - dominated.size()) *
			      landmarks.size());
	sum_density += d;
#endif
	calls_to_hs_opt += 1;
	hss_stats.start();
#ifdef ILA_USE_HS_EXTERN
	NTYPE new_v = hitting_set_extern(rp_actions, rp_cost, v_max);
	if (zero_cost_fill)
	  rp_actions.insert(zero_cost_actions);
#else
	set.assign_value(false, ins.n_actions());
	if (zero_cost_fill)
	  set.insert(zero_cost_actions);
	NTYPE new_v =
	  (hitting_set_use_split ?
	   hitting_set_split(lm_to_hit, out, set, 0, rp_actions, rp_cost, v_max, 0)
	   : hitting_set(lm_to_hit, out, set, 0, rp_actions, rp_cost, v_max));
#endif
	hss_stats.stop();
#ifdef MEASURE_STABILITY
	index_type nc = rp_actions.count_common(last_rp);
	double s = (2*nc)/((double)(rp_actions.size() + last_rp.size()));
	sum_stability += s;
	n_stability += 1;
	last_rp.assign_copy(rp_actions);
#endif
	if (stats.break_signal_raised()) return v_max;
	v_max = MAX(v_max, new_v);
#ifdef TRACE_PRINT_HLB
	new_lb(v_max);
#endif
      }
      if (ILA_use_saturation) {
	all_rp_actions.insert(rp_actions);
	reach.recompute(init, all_rp_actions);
      }
      else {
	reach.recompute(init, rp_actions);
      }
      // if the goal is now reachable...
      if (!reach.unreachable(goal)) {
	// if our current hitting set is approximate, we now have to
	// find a min cost set, and check if the goal is still reachable
	// using only those actions:
	if (ILA_use_approximate) {
	  if (verbose_level > 1) {
	    std::cerr << "checking if rp cost " << PRINT_NTYPE(rp_cost)
		      << " is min..." << std::endl;
	    // std::cerr << "rp_actions = " << rp_actions << std::endl;
	  }
	  calls_to_hs_opt += 1;
	  //dump_hitting_set_problem_wcnf(v_max, rp_cost, std::cout);
	  //dump_hitting_set_problem_lp(v_max, rp_cost, std::cout);
#ifdef MEASURE_DENSITY
	  long n = 0;
	  for (index_type k = 0; k < landmarks.size(); k++)
	    n += lm_branch_order[k].size();
	  double d = n/(double)((lm_actions.size() - dominated.size()) *
				landmarks.size());
	  sum_density += d;
#endif
	  hss_stats.start();
#ifdef ILA_USE_HS_EXTERN
	  NTYPE new_v = hitting_set_extern(rp_actions, rp_cost, v_max);
	  if (zero_cost_fill)
	    rp_actions.insert(zero_cost_actions);
#else
	  set.assign_value(false, ins.n_actions());
	  if (zero_cost_fill)
	    set.insert(zero_cost_actions);
	  NTYPE new_v =
	    (hitting_set_use_split ?
	     hitting_set_split(lm_to_hit, out, set, 0, rp_actions, rp_cost, v_max, 0)
	     : hitting_set(lm_to_hit, out, set, 0, rp_actions, rp_cost, v_max));
#endif
	  // if (!(new_v <= rp_cost)) {
	  //   std::cerr << "new_v = " << new_v << std::endl;
	  //   std::cerr << "rp_cost = " << rp_cost << std::endl;
	  // }
	  assert(new_v <= rp_cost);
	  hss_stats.stop();
#ifdef MEASURE_STABILITY
	  index_type nc = rp_actions.count_common(last_rp);
	  double s = (2*nc)/((double)(rp_actions.size() + last_rp.size()));
	  sum_stability += s;
	  n_stability += 1;
	  last_rp.assign_copy(rp_actions);
#endif
	  if (stats.break_signal_raised()) return v_max;
	  assert(new_v == cost.sum(rp_actions));
	  v_max = MAX(v_max, new_v);
#ifdef TRACE_PRINT_HLB
	  new_lb(v_max);
#endif
	  reach.recompute(init, rp_actions);
	  all_rp_actions.assign_copy(rp_actions);
	  if (!reach.unreachable(goal))
	    reached_goal = true;
	}
	else if (ILA_use_saturation) {
	  reach.recompute(init, rp_actions);
	  all_rp_actions.assign_copy(rp_actions);
	  if (!reach.unreachable(goal))
	    reached_goal = true;
	}
	else {
	  reached_goal = true;
	}
      }
    }
  }
  if (zero_cost_fill) {
    bool_vec ua(false, rp_actions.size());
    for (index_type k = 0; k < rp_actions.size(); k++)
      if (reach.unreachable_action(rp_actions[k]))
	ua[k] = true;
    rp_actions.remove(ua);
  }
  if (verbose_level > 0) {
    std::cerr << "valid rp cost = " << PRINT_NTYPE(v_max)
	      << " (" << stats.time() << " seconds)" << std::endl;
    if (verbose_level > 2) {
      std::cerr << "valid rp:" << std::endl;
      print_rp(std::cerr, rp_actions);
    }
  }
  return v_max;
}

#ifdef USE_SCIP
#ifdef HPLUS_USE_RP_EXTERN
NTYPE ILB::compute_relaxed_plan_extern
(const index_set& init, const index_set& goal)
{
  NTYPE ub = compute_relaxed_plan_FF(init, goal);
  std::cerr << "initial rp cost = " << PRINT_NTYPE(ub) << std::endl;
  if (INFINITE(ub)) return ub;

  // atoms in the encoding are relevant, non-initial atoms
  index_set atoms(ins.n_atoms(), init);
  atoms.intersect(relevant_atoms);
  index_vec atom_map;
  mapping::invert_map(atoms, atom_map);
  // actions in the encoding are all relevant actions
  index_set acts(relevant_actions);
  index_vec act_map;
  mapping::invert_map(acts, act_map);

  /// debug printing
  // for (index_type i = 0; i < atoms.size(); i++)
  //   std::cout << "atom" << i << " = " << ins.atoms[atoms[i]].name
  // 	      << std::endl;
  // for (index_type k = 0; k < acts.size(); k++)
  //   std::cout << "act" << k << " = " << ins.actions[acts[k]].name
  // 	      << std::endl;

  // standard init stuff
  SCIP* scip;
  SCIPcreate(&scip);
  if (verbose_level < 2) {
    SCIPsetMessagehdlrQuiet(scip, TRUE);
  }
  SCIPincludeDefaultPlugins(scip);
  SCIP_RETCODE res = SCIPincludeConshdlrLinearOrdering(scip);
  assert(res == SCIP_OKAY);
  res = SCIPcreateProbBasic(scip, "h+");
  assert(res == SCIP_OKAY);

  // variables of the encoding
  SCIP_VAR* p_vars[atoms.size()];
  SCIP_VAR* a_vars[acts.size()];
  // add indexes into acts:
  index_set_vec add(EMPTYSET, atoms.size());
  SCIP_VAR** add_vars[atoms.size()];
  SCIP_VAR** cl_vars[atoms.size()];
	
  // create variables
  std::cerr << "creating variables..." << std::endl;
  for (index_type i = 0; i < atoms.size(); i++) {
    // force goal vars to be 1 by setting lower bound.
    double l = (goal.contains(atoms[i]) ? 1 : 0);
    SCIP_VAR* v;
    std::ostringstream vname;
    vname << "atom" << i;
    res = SCIPcreateVarBasic(scip, &v, vname.str().c_str(),
			     l, 1, 0, SCIP_VARTYPE_BINARY);
    assert(res == SCIP_OKAY);
    res = SCIPaddVar(scip, v);
    assert(res == SCIP_OKAY);
    p_vars[i] = v;
  }
  for (index_type k = 0; k < acts.size(); k++) {
    // set action costs in objective.
    double w = N_TO_D(ins.actions[acts[k]].cost);
    std::ostringstream vname;
    vname << "act" << k;
    res = SCIPcreateVarBasic(scip, &(a_vars[k]), vname.str().c_str(),
			     0, 1, w, SCIP_VARTYPE_BINARY);
    assert(res == SCIP_OKAY);
    res = SCIPaddVar(scip, a_vars[k]);
    assert(res == SCIP_OKAY);
  }
  for (index_type i = 0; i < atoms.size(); i++) {
    for (index_type k = 0; k < acts.size(); k++)
      if (ins.actions[acts[k]].add.contains(atoms[i]))
	add[i].insert(k);
    // for (index_type k = 0; k < ins.atoms[atoms[i]].add_by.size(); k++)
    //   if (act_map[ins.atoms[atoms[i]].add_by[k]] != no_such_index)
    // 	add[i].insert(act_map[ins.atoms[atoms[i]].add_by[k]]);
    add_vars[i] = new SCIP_VAR*[add[i].size() + 1];
    for (index_type k = 0; k < add[i].size(); k++) {
      std::ostringstream vname;
      vname << "add(" << i << "," << add[i][k] << ")";
      res = SCIPcreateVarBasic(scip, &(add_vars[i][k]), vname.str().c_str(),
			       0, 1, 0, SCIP_VARTYPE_BINARY);
      assert(res == SCIP_OKAY);
      res = SCIPaddVar(scip, add_vars[i][k]);
      assert(res == SCIP_OKAY);
    }
  }
  for (index_type i = 0; i < atoms.size(); i++) {
    cl_vars[i] = new SCIP_VAR*[atoms.size() + 1];
    for (index_type j = 0; j < atoms.size(); j++) { 
      std::ostringstream vname;
      vname << "cl(" << i << "," << j << ")";
      res = SCIPcreateVarBasic(scip, &(cl_vars[i][j]), vname.str().c_str(),
			       0, 1, 0, SCIP_VARTYPE_BINARY);
      assert(res == SCIP_OKAY);
      res = SCIPaddVar(scip, cl_vars[i][j]);
      assert(res == SCIP_OKAY);
    }
  }

  // create constraints

  // (1) support: p -> \/ add(p, a)
  std::cerr << "creating type 1 constraints..." << std::endl;
  for (index_type i = 0; i < atoms.size(); i++) {
    double cs[add[i].size() + 2];
    SCIP_VAR* vs[add[i].size() + 2];
    cs[0] = -1;
    vs[0] = p_vars[i];
    for (index_type k = 0; k < add[i].size(); k++) {
      cs[k + 1] = 1;
      vs[k + 1] = add_vars[i][k];
    }
    SCIP_CONS* c;
    res = SCIPcreateConsBasicLinear(scip, &c, "type1", add[i].size() + 1,
				    vs, cs, 0, SCIPinfinity(scip));
    assert(res == SCIP_OKAY);
    res = SCIPsetUpgradeConsLinear(scip, c, TRUE);
    assert(res == SCIP_OKAY);
    res = SCIPaddCons(scip, c);
    assert(res == SCIP_OKAY);
    SCIPreleaseCons(scip, &c);
  }

  // (2) support implies action: add(p, a) -> a
  std::cerr << "creating type 2 constraints..." << std::endl;
  for (index_type i = 0; i < atoms.size(); i++)
    for (index_type k = 0; k < add[i].size(); k++) {
      double cs[2];
      cs[0] = -1;
      cs[1] = 1;
      SCIP_VAR* vs[2];
      vs[0] = add_vars[i][k];
      assert(add[i][k] < acts.size());
      vs[1] = a_vars[add[i][k]];
      SCIP_CONS* c;
      res = SCIPcreateConsBasicLinear(scip, &c, "type2", 2,
				      vs, cs, 0, SCIPinfinity(scip));
      assert(res == SCIP_OKAY);
      res = SCIPsetUpgradeConsLinear(scip, c, TRUE);
      assert(res == SCIP_OKAY);
      res = SCIPaddCons(scip, c);
      assert(res == SCIP_OKAY);
      SCIPreleaseCons(scip, &c);
    }

  // (3) action implies prec: a -> q forall q in pre(a).
  std::cerr << "creating type 3 constraints..." << std::endl;
  for (index_type k = 0; k < acts.size(); k++)
    for (index_type i = 0; i < ins.actions[acts[k]].pre.size(); i++)
      if (atom_map[ins.actions[acts[k]].pre[i]] != no_such_index) {
	double cs[2];
	cs[0] = -1;
	cs[1] = 1;
	SCIP_VAR* vs[2];
	vs[0] = a_vars[k];
	vs[1] = p_vars[atom_map[ins.actions[acts[k]].pre[i]]];
	SCIP_CONS* c;
	res = SCIPcreateConsBasicLinear(scip, &c, "type3", 2,
					vs, cs, 0, SCIPinfinity(scip));
	assert(res == SCIP_OKAY);
	res = SCIPsetUpgradeConsLinear(scip, c, TRUE);
	assert(res == SCIP_OKAY);
	res = SCIPaddCons(scip, c);
	assert(res == SCIP_OKAY);
	SCIPreleaseCons(scip, &c);
      }

  // (4) support implies cl: add(p, a) -> cl(q,p) forall q in pre(a).
  std::cerr << "creating type 4 constraints..." << std::endl;
  for (index_type i = 0; i < atoms.size(); i++)
    for (index_type k = 0; k < add[i].size(); k++) {
      index_type a = acts[add[i][k]];
      for (index_type j = 0; j < ins.actions[a].pre.size(); j++)
	if (atom_map[ins.actions[a].pre[j]] != no_such_index) {
	  double cs[2];
	  cs[0] = -1;
	  cs[1] = 1;
	  SCIP_VAR* vs[2];
	  vs[0] = add_vars[i][k];
	  vs[1] = cl_vars[atom_map[ins.actions[a].pre[j]]][i];
	  SCIP_CONS* c;
	  res = SCIPcreateConsBasicLinear(scip, &c, "type4", 2,
					  vs, cs, 0, SCIPinfinity(scip));
	  assert(res == SCIP_OKAY);
	  res = SCIPsetUpgradeConsLinear(scip, c, TRUE);
	  assert(res == SCIP_OKAY);
	  res = SCIPaddCons(scip, c);
	  assert(res == SCIP_OKAY);
	  SCIPreleaseCons(scip, &c);
	  assert(res == SCIP_OKAY);
	}
    }

#ifndef RP_SCIP_CG
  // (5) LO(cl).
  std::cerr << "creating LO constraint..." << std::endl;
  SCIP_CONS* clo;
  res =  SCIPcreateConsLinearOrdering(scip, &clo, "type5",
  				      atoms.size(), cl_vars,
  				      TRUE, // initial
  				      TRUE, // separate
  				      TRUE, // enforce
  				      TRUE, // check
  				      TRUE, // propagate
  				      FALSE, // local
  				      FALSE, // modifiable
  				      FALSE, // dynamic
  				      FALSE, // removable
  				      FALSE  // stickingatnode
  				      );
  assert(res == SCIP_OKAY);
  res = SCIPaddCons(scip, clo);
  assert(res == SCIP_OKAY);
  SCIPreleaseCons(scip, &clo);
  assert(res == SCIP_OKAY);

  // std::cerr << "creating type 5 constraints..." << std::endl;
  // double cs[3] = {1,1,1};
  // SCIP_VAR* vs[3];
  // for (index_type i = 0; i < atoms.size(); i++)
  //   for (index_type j = 0; j < atoms.size(); j++)
  //     if (j != i) {
  // 	// 2-cycle: cl[i][j] + cl[j][i] <= 1.
  // 	vs[0] = cl_vars[i][j];
  // 	vs[1] = cl_vars[j][i];
  // 	cs[0] = 1; cs[1] = 1;
  // 	SCIP_CONS* c;
  // 	res = SCIPcreateConsBasicLinear(scip, &c, "type5:2", 2,
  // 					vs, cs, -SCIPinfinity(scip), 1);
  // 	assert(res == SCIP_OKAY);
  // 	res = SCIPsetUpgradeConsLinear(scip, c, TRUE);
  // 	assert(res == SCIP_OKAY);
  // 	res = SCIPaddCons(scip, c);
  // 	assert(res == SCIP_OKAY);
  // 	SCIPreleaseCons(scip, &c);
  // 	assert(res == SCIP_OKAY);
  // 	for (index_type k = 0; k < atoms.size(); k++)
  // 	  if ((k != j) && (k != i)) {
  // 	    // transitivity: -cl[i][j] + -cl[j][k] + cl[i][k] >= -1.
  // 	    vs[0] = cl_vars[i][j];
  // 	    vs[1] = cl_vars[j][k];
  // 	    vs[2] = cl_vars[i][k];
  // 	    cs[0] = -1; cs[1] = -1; cs[2] = 1;
  // 	    SCIP_CONS* c;
  // 	    res = SCIPcreateConsBasicLinear(scip, &c, "type5:3", 3,
  // 					    vs, cs, -1, SCIPinfinity(scip));
  // 	    assert(res == SCIP_OKAY);
  // 	    res = SCIPsetUpgradeConsLinear(scip, c, TRUE);
  // 	    assert(res == SCIP_OKAY);
  // 	    res = SCIPaddCons(scip, c);
  // 	    assert(res == SCIP_OKAY);
  // 	    SCIPreleaseCons(scip, &c);
  // 	    assert(res == SCIP_OKAY);
  // 	  }
  //     }
#else
  std::cerr << "creating type 5:2 constraints..." << std::endl;
  double cs[2];
  SCIP_VAR* vs[2];
  for (index_type i = 0; i < atoms.size(); i++)
    for (index_type j = 0; j < atoms.size(); j++)
      if (j != i) {
  	// 2-cycle: cl[i][j] + cl[j][i] <= 1.
  	vs[0] = cl_vars[i][j];
  	vs[1] = cl_vars[j][i];
  	cs[0] = 1; cs[1] = 1;
  	SCIP_CONS* c;
  	res = SCIPcreateConsBasicLinear(scip, &c, "type5:2", 2,
  					vs, cs, -SCIPinfinity(scip), 1);
  	assert(res == SCIP_OKAY);
  	res = SCIPsetUpgradeConsLinear(scip, c, TRUE);
  	assert(res == SCIP_OKAY);
  	res = SCIPaddCons(scip, c);
  	assert(res == SCIP_OKAY);
  	SCIPreleaseCons(scip, &c);
  	assert(res == SCIP_OKAY);
      }
#endif

  // debug printing:
  // SCIPprintOrigProblem(scip, NULL, NULL, FALSE);
  
#ifndef RP_SCIP_CG
  // create an initial (non-optimal) solution
  index_vec rp_seq;
  reach.order_relaxed_plan(init, goal, rp_actions, rp_seq);

  SCIP_SOL* insol;
  res = SCIPcreateSol(scip, &insol, NULL);
  assert(res == SCIP_OKAY);

  bool_vec s(false, ins.n_atoms());
  index_vec atom_seq;
  for (index_type i = 0; i < init.size(); i++) {
    atom_seq.append(init[i]);
    s[init[i]] = true;
  }
  for (index_type k = 0; k < rp_seq.size(); k++) {
    index_type a = acts.first(rp_seq[k]);
    if (a != no_such_index) {
      for (index_type i = 0; i < ins.actions[rp_seq[k]].add.size(); i++)
	if (!s[ins.actions[rp_seq[k]].add[i]]) {
	  atom_seq.append(ins.actions[rp_seq[k]].add[i]);
	  s[ins.actions[rp_seq[k]].add[i]] = true;
	  index_type p = atoms.first(ins.actions[rp_seq[k]].add[i]);
	  if (p != no_such_index) {
	    index_type ai = add[p].first(a);
	    assert(ai != no_such_index);
	    //std::cerr << "action " << a << " (#" << ai << ") adds "
	    // 	      << p << std::endl;
	    res = SCIPsetSolVal(scip, insol, add_vars[p][ai], 1);
	    assert(res == SCIP_OKAY);
	    //std::cerr << "atom " << p << " holds" << std::endl;
	    res = SCIPsetSolVal(scip, insol, p_vars[p], 1);
	    assert(res == SCIP_OKAY);
	  }
	}
      //std::cerr << "action " << a << " included" << std::endl;
      res = SCIPsetSolVal(scip, insol, a_vars[a], 1);
      assert(res == SCIP_OKAY);
    }
    else {
      std::cerr << "error: irrelevant action " << rp_seq[k]
		<< "." << ins.actions[rp_seq[k]].name
		<< " in initial rp" << std::endl;
    }
  }
  // add atoms not made true in rp to end of sequence
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (!s[i]) atom_seq.append(i);
  for (index_type i = 0; i < atom_seq.size(); i++) {
    index_type pi = atoms.first(atom_seq[i]);
    if (pi != no_such_index) {
      for (index_type j = i + 1; j < atom_seq.size(); j++) {
  	index_type pj = atoms.first(atom_seq[j]);
  	if (pj != no_such_index) {
  	  //std::cerr << "atom " << pi << " before " << pj << std::endl;
	  res = SCIPsetSolVal(scip, insol, cl_vars[pi][pj], 1);
	  assert(res == SCIP_OKAY);
	  res = SCIPsetSolVal(scip, insol, cl_vars[pj][pi], 0);
	  assert(res == SCIP_OKAY);
  	}
      }
    }
  }

  // SCIPprintSol(scip, insol, NULL, FALSE);

  std::cerr << "testing solution..." << std::endl;
  // must call transformProb before storing a solution!
  res = SCIPtransformProb(scip);
  assert(res == SCIP_OKAY);
  SCIP_Bool insolok;
  res = SCIPtrySolFree(scip, &insol, TRUE, TRUE, TRUE, TRUE, &insolok);
  assert(res == SCIP_OKAY);
  if (insolok == FALSE) {
    std::cerr << "error: input solution not accepted" << std::endl;
  }
#endif // !RP_SCIP_CG

  // solve
  std::cerr << "solving..." << std::endl;
#ifdef RP_SCIP_CG
  graph prec;
  bool ok = false;
  while (!ok) {
    ok = true;
    res = SCIPsolve(scip);
    assert(res == SCIP_OKAY);
    SCIP_Status status = SCIPgetStatus(scip);
    if (status != SCIP_STATUS_OPTIMAL) {
      std::cerr << "error: not solved to optimality!" << std::endl;
      exit(1);
    }
    SCIP_SOL* sol = SCIPgetBestSol(scip);
    double obj_val = SCIPgetSolOrigObj(scip, sol);
    if (verbose_level > 0) {
      std::cerr << "obj val = " << obj_val << std::endl;
    }
    prec.init(atoms.size());
    for (index_type i = 0; i < atoms.size(); i++)
      for (index_type j = 0; j < atoms.size(); j++)
	if (i != j) {
	  double v = SCIPgetSolVal(scip, sol, cl_vars[i][j]);
	  if (v > 0.5)
	    prec.add_edge(i, j);
	}
    res = SCIPfreeTransform(scip);
    assert(res == SCIP_OKAY);
    std::cerr << "extracted " << prec.n_edges() << " precedences" << std::endl;
    index_vec c_el;
    index_type c_len = prec.shortest_cycle(c_el);
    // prec.strongly_connected_components();
    // index_type cy_min = no_such_index;
    // index_type cy_min_size = 0;
    // std::cerr << "component sizes:";
    // for (index_type i = 0; i < prec.n_components(); i++) {
    //   std::cerr << " " << prec.component_size(i);
    //   if (prec.component_size(i) > 1) {
    // 	if ((cy_min == no_such_index) ||
    // 	    (prec.component_size(i) < cy_min_size)) {
    // 	  cy_min = i;
    // 	  cy_min_size = prec.component_size(i);
    // 	}
    //   }
    // }
    // std::cerr << std::endl;
    if (c_len != no_such_index) {
      assert(c_len == c_el.size());
      std::cerr << "found cycle of length " << c_len
		<< ": " << c_el << std::endl;
      assert(c_len > 1);
      double cs[c_len];
      SCIP_VAR* vs[c_len];
      cs[0] = 1;
      vs[0] = cl_vars[c_el[c_len - 1]][c_el[0]];
      for (index_type i = 1; i < c_len; i++) {
	cs[i] = 1;
	vs[i] = cl_vars[c_el[i - 1]][c_el[i]];
      }
      SCIP_CONS* c;
      res = SCIPcreateConsBasicLinear
	(scip, &c, "type5:c", c_len, vs, cs, 0, c_len - 1);
      assert(res == SCIP_OKAY);
      res = SCIPsetUpgradeConsLinear(scip, c, TRUE);
      assert(res == SCIP_OKAY);
      res = SCIPaddCons(scip, c);
      assert(res == SCIP_OKAY);
      SCIPreleaseCons(scip, &c);
      assert(res == SCIP_OKAY);
      ok = false;
    }
    else {
      assert(prec.acyclic());
    }
  }
#else
  res = SCIPsolve(scip);
  assert(res == SCIP_OKAY);
  SCIP_Status status = SCIPgetStatus(scip);
  if (status != SCIP_STATUS_OPTIMAL) {
    std::cerr << "error: not solved to optimality!" << std::endl;
    exit(1);
  }
#endif

  // extract solution
  SCIP_SOL* sol = SCIPgetBestSol(scip);
  assert(sol != NULL);
  // debug printing:
  // SCIPprintSol(scip, sol, NULL, FALSE);
  double obj_val = SCIPgetSolOrigObj(scip, sol);
  NTYPE val = D_TO_N(obj_val);
  if (verbose_level > 0) {
    std::cerr << "obj val = " << obj_val << " (R), " << val << " (Q)"
	      << std::endl;
  }
#ifdef NTYPE_RATIONAL
  val = rational::ceil_to(val, ac_div);
#endif
  NTYPE sum = 0;
  rp_actions.clear();
  for (index_type k = 0; k < acts.size(); k++) {
    double v = SCIPgetSolVal(scip, sol, a_vars[k]);
    if (v > 0.5) {
      sum += ins.actions[acts[k]].cost;
      rp_actions.insert(acts[k]);
    }
  }
  assert(sum == val);

  // free memory
  for (index_type i = 0; i < atoms.size(); i++)
    SCIPreleaseVar(scip, &(p_vars[i]));
  for (index_type k = 0; k < acts.size(); k++)
    SCIPreleaseVar(scip, &(a_vars[k]));
  for (index_type i = 0; i < atoms.size(); i++) {
    for (index_type k = 0; k < add[i].size(); k++)
      SCIPreleaseVar(scip, &(add_vars[i][k]));
    delete[] add_vars[i];
  }
  for (index_type i = 0; i < atoms.size(); i++) {
    for (index_type j = 0; j < atoms.size(); j++)
      SCIPreleaseVar(scip, &(cl_vars[i][j]));
    delete[] cl_vars[i];
  }

  SCIPfree(&scip);

  // return optimal cost
  return sum;
}
#endif
#endif

#ifdef IMAI_FUKUNAGA
// The Imai-Fukunaga (ECAI 2014) formulation of h+ as a MIP, solved
// with CPLEX.
NTYPE ILB::compute_relaxed_plan_extern
(const index_set& init, const index_set& goal)
{
  
  int smodelvars=1;
  int smfalse=1;
  std::cerr << "applying Imai-Fukunaga encoding of h+" << std::endl;
  NTYPE ub = compute_relaxed_plan_FF(init, goal);
  std::cerr << "initial rp cost = " << PRINT_NTYPE(ub) << std::endl;
  if (INFINITE(ub)) return ub;

  if (!FA_and_ALM_avail) {
    compute_FA_and_ALM(init);
    if (stats.break_signal_raised()) {
      return 0;
    }
  }

  NTYPE sum = 0;
	FILE *Flp;
	char fname[1000];
	char lpname[1000];
	sprintf(fname,"hsplp");
	
	
	if (probnum < 10)
		sprintf(lpname,"%s0%d.lp",fname, probnum);
		
	else
		sprintf(lpname,"%s%d.lp",fname, probnum);
	Flp = fopen(lpname,"w");
  // atoms in the encoding are relevant, non-initial atoms
  index_set atoms(ins.n_atoms(), init); // complement constructor!!
  atoms.intersect(relevant_atoms);
  // atom_map maps from orignal (full) set of atoms to indices
  // in the reduced set of relevant atoms (index_set atoms)
  index_vec atom_map;
  mapping::invert_map(atoms, atom_map, ins.n_atoms());
  // actions in the encoding are all relevant actions
  index_set acts(relevant_actions);
	//index_vec ac_map;
	//mapping::invert_map(acts, ac_map, ins.n_actions());
  index_set glms; // atom landmarks of the goal
  //std::cerr << "atom landmarks = " << atom_landmarks << std::endl;
  //std::cerr << "goal = " << goal << std::endl;
  index_set alms; // singleton action landmarks of the goal
  for (index_type i = 0; i < goal.size(); i++)
    for (index_type j = 0; j < atom_landmarks[goal[i]].size(); j++)
      if (relevant_atoms[atom_landmarks[goal[i]][j]]) {
	glms.insert(atom_landmarks[goal[i]][j]);
	if (first_achievers[atom_landmarks[goal[i]][j]].size() == 1)
	  alms.insert(first_achievers[atom_landmarks[goal[i]][j]]);
		
      }
  std::cerr << "found " << glms.size() << " atom and " << alms.size()
	    << " (singleton) action landmarks" << std::endl;


  // standard init stuff
  IloEnv env;
  try {


    int **SstepOrder =(int **)malloc(atoms.size() * sizeof(int *));
    for (int prop1=0;prop1<atoms.size();prop1++){
      SstepOrder[prop1]=(int *)malloc(atoms.size() * sizeof(int ));
      for (int prop2=0;prop2<atoms.size();prop2++)
        SstepOrder[prop1][prop2]=-1;
    }

	//fprintf(Flp,"landmarksize(%d).\n",landmarks.size());
    //IloBoolVar **graphEdge = (IloBoolVar **)malloc(atoms.size() * sizeof(IloBoolVar *));
	int **graphEdge = (int **)malloc(atoms.size() * sizeof(int *));
    for (int prop1=0;prop1<atoms.size();prop1++){
      graphEdge[prop1] =(int *)malloc(atoms.size() * sizeof(int ));
	  for (int prop2=0;prop2<atoms.size();prop2++){
		  graphEdge[prop1][prop2]=-1;
	  }
    }
    std::cerr << "creating variables..." << std::endl;
    // variables of the encoding
    int p_vars[atoms.size()];
    for (index_type i = 0; i < atoms.size(); i++){
		smodelvars++;
		p_vars[i]= smodelvars;
	}
		


    int a_vars[acts.size()];
	
    for (index_type k = 0; k < acts.size(); k++){
		
		smodelvars++;
		a_vars[k]= smodelvars;
	}
      //a_vars[k].setName(ins.actions[acts[k]].name->to_cstring());
    //model.add(a_vars);
    
    //model.add(at_vars);
    // add indexes into acts:
    index_set_vec add(EMPTYSET, atoms.size());
    int * add_vars[atoms.size()];
    //exit(0);
    for (index_type i = 0; i < atoms.size(); i++) {
      for (index_type k = 0; k < acts.size(); k++)
        if (first_achievers[atoms[i]].contains(acts[k])) {
          assert(ins.actions[acts[k]].add.contains(atoms[i]));
          add[i].insert(k);
      }
      IloBoolVarArray tmp(env, add[i].size());
        //add_vars[i] = tmp;
	  add_vars[i]=(int *)malloc(add[i].size() * sizeof(int ));
	  for (int k=0; k<add[i].size(); k++){
		  smodelvars++;
		  add_vars[i][k]=smodelvars;
	  }
		  
    }

        
        //model.add(graphEdge[i][j]);

    for (int a=0; a<acts.size(); a++){
      for (index_type i = 0; i < ins.actions[acts[a]].pre.size(); i++){
        if (atom_map[ins.actions[acts[a]].pre[i]] == no_such_index)
          continue;
        for (index_type j = 0; j < ins.actions[acts[a]].add.size(); j++){
          if (atom_map[ins.actions[acts[a]].add[j]] == no_such_index)
            continue;
          if (ins.actions[acts[a]].pre[i]!=ins.actions[acts[a]].add[j]){
            SstepOrder[atom_map[ins.actions[acts[a]].pre[i]]][atom_map[ins.actions[acts[a]].add[j]]]=1;
			if (graphEdge[atom_map[ins.actions[acts[a]].pre[i]]][atom_map[ins.actions[acts[a]].add[j]]]==-1){
				smodelvars++;
				graphEdge[atom_map[ins.actions[acts[a]].pre[i]]][atom_map[ins.actions[acts[a]].add[j]]]=smodelvars;
			}
            //if(0)
            if (SstepOrder[atom_map[ins.actions[acts[a]].add[j]]][atom_map[ins.actions[acts[a]].pre[i]]]==1){
              //model.add(graphEdge[atom_map[ins.actions[acts[a]].pre[i]]][atom_map[ins.actions[acts[a]].add[j]]]+graphEdge[atom_map[ins.actions[acts[a]].add[j]]][atom_map[ins.actions[acts[a]].pre[i]]] < 2);
			  //fprintf(Flp,":- p%dp%d, p%dp%d.\n",atom_map[ins.actions[acts[a]].pre[i]],atom_map[ins.actions[acts[a]].add[j]],atom_map[ins.actions[acts[a]].add[j]],atom_map[ins.actions[acts[a]].pre[i]]);	
            //model.add(a_vars[a] <= graphEdge[atom_map[ins.actions[acts[a]].pre[i]]][atom_map[ins.actions[acts[a]].add[j]]]);
			}
          }
        }

      }

    }
    
    
    
    

    for (index_type i = 0; i < atoms.size(); i++)
      for (index_type k = 0; k < add[i].size(); k++) {
        int a = add[i][k];
        for (index_type j = 0; j < ins.actions[acts[a]].pre.size(); j++){
          if (atom_map[ins.actions[acts[a]].pre[j]] == no_such_index)
            continue;
          if (atom_map[ins.actions[acts[a]].pre[j]]==i)
			  //exit(0);
            continue;
          //model.add(add_vars[i][k] <= graphEdge[atom_map[ins.actions[acts[a]].pre[j]]][i]);
		  
		  //fprintf(Flp,"p%dp%d :- a%dp%d.\n",atom_map[ins.actions[acts[a]].pre[j]],i,a,i);
		  if (graphEdge[atom_map[ins.actions[acts[a]].pre[j]]][i] == -1){
			  smodelvars++;
			  graphEdge[atom_map[ins.actions[acts[a]].pre[j]]][i]=smodelvars;
		  }
		  fprintf(Flp,"1 %d 1 0 %d\n",graphEdge[atom_map[ins.actions[acts[a]].pre[j]]][i], add_vars[i][k]);

          //int x = (acts.size() + 1);
          //model.add(pt_vars[atom_map[ins.actions[acts[a]].pre[j]]]+1 <= (pt_vars[i] + x - (x * add_vars[i][k])));

        }
        
      }
    
    
    
    
    
    
    
    

 
    int **adjorddes =(int **)malloc(atoms.size() * sizeof(int *));
    int *nadjorddes = (int *)malloc(atoms.size() * sizeof(int ));
    int *sadjorddes = (int *)malloc(atoms.size() * sizeof(int ));
    int **adjorddesrev =(int **)malloc(atoms.size() * sizeof(int *));;
    int *nadjorddesrev = (int *)malloc(atoms.size() * sizeof(int ));
    int *sadjorddesrev = (int *)malloc(atoms.size() * sizeof(int ));
    for (int i=0; i<atoms.size(); i++){
      adjorddes[i] = (int *)malloc(100*sizeof(int ));
      adjorddesrev[i] = (int *)malloc(100*sizeof(int ));
      nadjorddes[i]=0;
      sadjorddes[i]=100;
      nadjorddesrev[i]=0;
      sadjorddesrev[i]=100;            
    }
    
    for (int i=0; i<atoms.size(); i++)
      for (int j=0; j<atoms.size(); j++){
		  //fprintf(Flp,":- p%dp%d, p%dp%d.\n",i,j,j,i);
        if (i==j)
          continue;
        if (SstepOrder[i][j]==1) {
          //exit(0);
						if (SstepOrder[j][i]==1)
							if (j>i)
								//fprintf(Flp,":- p%dp%d, p%dp%d.\n",i,j,j,i);
								fprintf(Flp,"1 1 2 0 %d %d\n",graphEdge[i][j],graphEdge[j][i]);

          addedge(adjorddes,sadjorddes,nadjorddes,i,j,1);
          addedge(adjorddesrev,sadjorddesrev,nadjorddesrev,j,i,1);
        }
        
        
        
      }
    int degreesOut[atoms.size()];
    int degreesIn[atoms.size()];
    for (int i=0; i<atoms.size(); i++){
      degreesOut[i]=0;
      degreesIn[i]=0;
    }
    for (int i=0; i<atoms.size(); i++){
      
      degreesOut[i]=nadjorddes[i];
      degreesIn[i]=nadjorddesrev[i];
      
    }  
    int delta=-1;
    int ts=0;
    //if(0)
    for (int i=0; i<atoms.size(); i++){
      
      long int mindegree=1000000000;
      int minv =-1;
      for (int j=0; j<atoms.size(); j++){
        if (degreesOut[j]==-1)
          continue;
        if (degreesIn[j]*degreesOut[j]<mindegree){
          mindegree=degreesIn[j]*degreesOut[j];
          minv=j;
        }
        //minv=i;
      }
      ts+=mindegree;
      printf("\n %d",  degreesOut[minv]);
      if (degreesOut[minv]>delta)
        delta=degreesOut[minv];



      for (int j=0; j<nadjorddesrev[minv]; j++){
        if ((SstepOrder[adjorddesrev[minv][j]][minv]!=1)) 
          continue;
        if (adjorddesrev[minv][j]==minv)
          continue;
        
        for (int k=0; k<nadjorddes[minv]; k++){
          if (adjorddes[minv][k]==adjorddesrev[minv][j])
            continue;
          if ((SstepOrder[minv][adjorddes[minv][k]]!=1) )
            continue;
          if (adjorddes[minv][k]==minv)
            continue;
          
          
          
          //exit(0);
          //if (k==j)
          //continue;
	  
          if ((SstepOrder[adjorddesrev[minv][j]][minv]==1) && (SstepOrder[minv][adjorddes[minv][k]]==1)){
            //exit(0);
            

            //addformula(transition,Fimpl(graphEdge(adjorddesrev[minv][j],minv),Fimpl(graphEdge(minv,adjorddes[minv][k]),graphEdge(adjorddesrev[minv][j],adjorddes[minv][k])) ));
            //model.add(a_vars[k] <= expr);
            //model.add(IloIfThen(env, graphEdge[adjorddesrev[minv][j]][minv]==IloTrue && graphEdge[minv][adjorddes[minv][k]]==IloTrue, graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]]==IloTrue ));
            //model.add(pt_vars[atom_map[ins.actions[acts[a]].pre[j]]]+1 <= (pt_vars[i] + x - (x * add_vars[i][k])));
            int x = (acts.size() + 1);
            //model.add(1 <= (graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]] + x - (x * (graphEdge[adjorddesrev[minv][j]][minv] + graphEdge[minv][adjorddes[minv][k]] -1))));
			//fprintf(Flp,"p%dp%d :- p%dp%d, p%dp%d.\n",adjorddesrev[minv][j],adjorddes[minv][k],adjorddesrev[minv][j],minv,minv,adjorddes[minv][k]);
                            if (adjorddesrev[minv][j]!=adjorddes[minv][k]){
								  if (graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]]==-1){
									  smodelvars++;
									  graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]]=smodelvars;
								  }
								  if (graphEdge[adjorddesrev[minv][j]][minv]==-1){
									  smodelvars++;
									  graphEdge[adjorddesrev[minv][j]][minv]=smodelvars;
								  }
								  if (graphEdge[minv][adjorddes[minv][k]]==-1){
									  smodelvars++;
									  graphEdge[minv][adjorddes[minv][k]]=smodelvars;
								  }								  
								//fprintf(Flp,"p%dp%d :- p%dp%d, p%dp%d.\n",adjorddesrev[minv][j],adjorddes[minv][k],adjorddesrev[minv][j],minv,minv,adjorddes[minv][k]); //done
								  fprintf(Flp,"1 %d 2 0 %d %d\n",graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]],graphEdge[adjorddesrev[minv][j]][minv],graphEdge[minv][adjorddes[minv][k]]); //done
								
							//else
								  if (graphEdge[adjorddes[minv][k]][adjorddesrev[minv][j]]==-1){
									  smodelvars++;
									  graphEdge[adjorddes[minv][k]][adjorddesrev[minv][j]]=smodelvars;
								  }
								  fprintf(Flp,"1 1 2 0 %d %d\n",graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]],graphEdge[adjorddes[minv][k]][adjorddesrev[minv][j]]);	
								//fprintf(Flp,":- p%dp%d, p%dp%d.\n",adjorddesrev[minv][j],adjorddes[minv][k],adjorddes[minv][k],adjorddesrev[minv][j]);	
							}
            //model.add(IloIfThen(env, graphEdge[adjorddesrev[minv][j]][minv] + graphEdge[minv][adjorddes[minv][k]] -1 >= 1, graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]] >= 1 ));
            if ((SstepOrder[adjorddesrev[minv][j]][adjorddes[minv][k]]!=1)){
              //if (SstepOrder[adjorddes[minv][k]][adjorddesrev[minv][j]]==-1){
				//model.add(graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]]+graphEdge[adjorddes[minv][k]][adjorddesrev[minv][j]] < 2);
				
				
				//fprintf(Flp,":- p%dp%d, p%dp%d.\n",adjorddesrev[minv][j],minv,minv,adjorddes[minv][k]);
			  //}
              SstepOrder[adjorddesrev[minv][j]][adjorddes[minv][k]]=1;
			  if (graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]]==-1){
				  smodelvars++;
				  graphEdge[adjorddesrev[minv][j]][adjorddes[minv][k]]=smodelvars;
			  }
              addedge(adjorddesrev,sadjorddesrev,nadjorddesrev,adjorddes[minv][k],adjorddesrev[minv][j],1);
              addedge(adjorddes,sadjorddes,nadjorddes,adjorddesrev[minv][j],adjorddes[minv][k],1);
              degreesOut[adjorddesrev[minv][j]]++;
              degreesIn[adjorddes[minv][k]]++;
            }
            continue;
          } 
          
        }
      }                      
      
      
      
      
      
      for (int k=0; k<atoms.size(); k++){
        // if ((SstepOrder[minv][k]!=1) && (SdestroyedTemp[minv][k]!=1) )
        // continue;
        if (k==minv)
          continue;
        
        if ((SstepOrder[minv][k]==1) ){
          SstepOrder[minv][k]=-1;
          degreesOut[minv]--;
          degreesIn[k]--;
        }
        
        if ((SstepOrder[k][minv]==1)){
          SstepOrder[k][minv]=-1;
          degreesOut[k]--;
          degreesIn[minv]--;
        }
        
      }
      if ((degreesOut[minv]!=0) || (degreesIn[minv]!=0))
      {
        printf("\n %d %d %d DEGREES ERROR2\n",degreesOut[minv], degreesIn[minv],i);
        exit(0);
      }          
      degreesOut[minv]=-1;
      degreesIn[minv]=-1;  
      
    }
    
    
    printf("\n delta = %d\n number of atoms = %d\n",delta,atoms.size());
    // create variables


    // add the objective: minimise action costs
    IloNumArray costs(env, acts.size());
    for (index_type k = 0; k < acts.size(); k++)
      costs[k] = N_TO_D(ins.actions[acts[k]].cost);
    //model.add(IloObjective(env, IloScalProd(costs, a_vars))); 

    // create constraints

    // (1) goals: p for all p in G.
    std::cerr << "creating type 1 constraints..." << std::endl;
    for (index_type i = 0; i < goal.size(); i++)
      if (atom_map[goal[i]] != no_such_index) {
		  fprintf(Flp,"1 %d 0 0\n",p_vars[atom_map[goal[i]]]);
		  //fprintf(Flp,"p%d.\n",atom_map[goal[i]]);
	//model.add(p_vars[atom_map[goal[i]]] == 1);
      }
    // (1') p for all p in landmarks(G).
    std::cerr << "creating type 1' constraints..." << std::endl;
    for (index_type i = 0; i < glms.size(); i++) {
      assert(atom_map[glms[i]] != no_such_index);
      //model.add(p_vars[atom_map[glms[i]]] == 1);
	  fprintf(Flp,"1 %d 0 0\n",p_vars[atom_map[glms[i]]]);
	  //fprintf(Flp,"p%d.\n",atom_map[glms[i]]);
    }
	for (index_type i = 0; i < atoms.size(); i++)
		fprintf(Flp,"3 1 %d 0 0\n",p_vars[i]);

    // (2) action implies precondition
    // ... and precondition not established by an inverse action
    std::cerr << "creating type 2 constraints..." << std::endl;
    for (index_type k = 0; k < acts.size(); k++) {
      //std::cerr << "atoms = " << atoms << std::endl;
      //std::cerr << "acts = " << acts << std::endl;
      //std::cerr << "k = " << k << std::endl;
      // find the inverses of acts[k]; note, invs will be a set
      // of action indices in the reduced set.
      
      //std::cerr << "add[k] = " << ins.actions[acts[k]].add << std::endl;
      //std::cerr << "rel-add[k] = " << radd << std::endl;
      //std::cerr << "pre[k] = " << ins.actions[acts[k]].pre << std::endl;

      // then construct a constraint for each precondition of acts[k]
      for (index_type i = 0; i < ins.actions[acts[k]].pre.size(); i++)
		if (atom_map[ins.actions[acts[k]].pre[i]] != no_such_index) {
		  index_type q = atom_map[ins.actions[acts[k]].pre[i]];
		  //model.add(a_vars[k] <= p_vars[q]);
		  //fprintf(Flp,"p%d :- a%d.\n",q,k);
		  fprintf(Flp,"1 %d 1 0 %d\n",p_vars[q],a_vars[k]);
		  
		  //expr = expr && p_vars[q];

		  //model.add(a_vars[k] <= expr);
		}
    }

    // (3) support implies action: add(p, a) -> a
    std::cerr << "creating type 3 constraints..." << std::endl;
    for (index_type i = 0; i < atoms.size(); i++)
      for (index_type k = 0; k < add[i].size(); k++) {
	      //model.add(add_vars[i][k] <= a_vars[add[i][k]]);
		  //for (index_type j = 0; j < ins.actions[acts[add[i][k]]].pre.size(); j++)
	  //if (atom_map[ins.actions[acts[add[i][k]]].pre[j]] != no_such_index) {
	  //index_type q = atom_map[ins.actions[acts[add[i][k]]].pre[j]];
		//if ((int)q== (int)i)
			//exit(0);
	  //}
		  //fprintf(Flp,"a%d :- a%dp%d.\n",add[i][k],add[i][k],i); 
		  fprintf(Flp,"1 %d 1 0 %d\n",a_vars[add[i][k]],add_vars[i][k]); 
        //model.add(add_vars[i][k] <= p_vars[i]);
      }

    // (4) support: p -> \/ add(p, a)
    std::cerr << "creating type 4 constraints..." << std::endl;
	
    for (index_type i = 0; i < atoms.size(); i++) {
		//if (add[i].size()==0)
			//exit(0);

		//fprintf(Flp,"1{");
      IloExpr expr(env);
      for (index_type k = 0; k < add[i].size(); k++) {

		//if (k==0)
			//fprintf(Flp,"{a%dp%d} :- p%d.\n", add[i][k],i, i );
		fprintf(Flp,"3 1 %d 1 0 %d\n", add_vars[i][k], p_vars[i] );
		//else
			//fprintf(Flp,";a%dp%d", add[i][k],i );
      }
	  //fprintf(Flp,"}1 :- p%d.\n",i);
      //model.add(p_vars[i] == expr); // note: == instead of <= !!
    }
	
	for (index_type i = 0; i < atoms.size(); i++) {
		//if (add[i].size()==0)
			//exit(0);

		//fprintf(Flp,":- p%d",i);
		fprintf(Flp,"1 1 %d %d",add[i].size()+1,add[i].size());
      IloExpr expr(env);
      for (index_type k = 0; k < add[i].size(); k++) {

		//if (k==0)
			//fprintf(Flp,", not a%dp%d", add[i][k],i );
		fprintf(Flp," %d", add_vars[i][k]);
		//else
			//fprintf(Flp,";a%dp%d", add[i][k],i );
      }
	  fprintf(Flp," %d", p_vars[i]);
	  fprintf(Flp,"\n");
      //model.add(p_vars[i] == expr); // note: == instead of <= !!
    }
	
	
/*
    // (5) T(pre(a) <= T(a)
    std::cerr << "creating type 5 constraints..." << std::endl;
    for (index_type k = 0; k < acts.size(); k++)
      for (index_type i = 0; i < ins.actions[acts[k]].pre.size(); i++)
	if (atom_map[ins.actions[acts[k]].pre[i]] != no_such_index) {
	  index_type q = atom_map[ins.actions[acts[k]].pre[i]];
	  model.add(pt_vars[q] <= at_vars[k]);
	}

    // (6) if add(p, a), then T(a) + 1 <= T(p)
    std::cerr << "creating type 6 constraints..." << std::endl;
    for (index_type i = 0; i < atoms.size(); i++)
      for (index_type k = 0; k < add[i].size(); k++) {
	int x = (acts.size() + 1);
	model.add((at_vars[add[i][k]] + 1) <=
		  (pt_vars[i] + x - (x * add_vars[i][k])));
      }
*/



							//fprintf(Flp,"#minimize {");
							fprintf(Flp,"6 0 %d 0", acts.size());

							//fprintf(Flp,"%d,a%d: a%d", (int)N_TO_D(ins.actions[acts[0]].cost), 0, 0);
							for (int i=0; i<acts.size(); i++){
								//if (costs[i]!=1)
									//exit(0);
								//if (subacts[i]==1)
									//continue;
								//if (lmactions[i]==1)
									//fprintf(Flp,"; %d,a%d: not a%d", 1000, i, i);
									//continue;
								//fprintf(Flp,"; %d,a%d: a%d", (int)N_TO_D(ins.actions[acts[i]].cost), i, i);
								fprintf(Flp," %d",  a_vars[i]);
							}
							for (int i=0; i<acts.size(); i++){
								//if (costs[i]!=1)
									//exit(0);
								//if (subacts[i]==1)
									//continue;
								//if (lmactions[i]==1)
									//fprintf(Flp,"; %d,a%d: not a%d", 1000, i, i);
									//continue;
								//fprintf(Flp,"; %d,a%d: a%d", (int)N_TO_D(ins.actions[acts[i]].cost), i, i);
								fprintf(Flp," %d",  (int)N_TO_D(ins.actions[acts[i]].cost));
							}
							fprintf(Flp,"\n");
							
							
							
							
							
							
							
							
							fprintf(Flp,"0\n0\nB+\n0\nB-\n1\n0\n1");
							
							fclose(Flp);
							//exit(0);



 








    std::cerr << "solving..." << std::endl;
    //IloCplex cplex(model);
    //cplex.setParam(IloCplex::Threads, 1);
    //cplex.writeParam("myparam.prm");
    //cplex.exportModel("hplus.lp");
    //cplex.solve();
    //if (cplex.getStatus() != IloAlgorithm::Optimal) {
      //std::cerr << "error: not solved to optimality!" << std::endl;
      ////std::cerr << "cplex status = " << cplex.getStatus() << std::endl;
      //exit(1);
    //}

    /// time-sorted vector of actions for debug solution print:
    //std::vector< std::pair<double, unsigned int> > timed_rp_acts;
    //rp_actions.clear();
    //printf("\n ******************************************************** \n", v);
    //for (index_type k = 0; k < acts.size(); k++) {
      ////double v = cplex.getValue(a_vars[k]);
      //if (v > 0.5) {
        //printf("\n %f", v);
	//sum += ins.actions[acts[k]].cost;
	//rp_actions.insert(acts[k]);
	//double w = cplex.getValue(at_vars[k]);
	//timed_rp_acts.push_back(std::pair<double, unsigned int>(w, k));
      //}
    //}
    //std::sort(timed_rp_acts.begin(), timed_rp_acts.end());
    //
    // for (index_type k = 0; k < timed_rp_acts.size(); k++) {
    //   index_type a = timed_rp_acts[k].second;
    //   std::cerr << "at " << timed_rp_acts[k].first << ": "
    // 		<< ins.actions[acts[a]].name
    // 		<< std::endl;
    //   for (index_type i = 0; i < ins.actions[acts[a]].pre.size(); i++)
    // 	if (atom_map[ins.actions[acts[a]].pre[i]] != no_such_index) {
    // 	  index_type q = atom_map[ins.actions[acts[a]].pre[i]];
    // 	  std::cerr << " - pre " << ins.atoms[ins.actions[acts[a]].pre[i]].name
    // 		    << " (active = " << cplex.getValue(p_vars[q])
    // 		    << ", at = " << cplex.getValue(pt_vars[q]) << ")"
    // 		    << std::endl;
    // 	}
    //   for (index_type i = 0; i < atoms.size(); i++)
    // 	for (index_type j = 0; j < add[i].size(); j++)
    // 	  if (add[i][j] == a)
    // 	    if (cplex.getValue(add_vars[i][j]) > 0.5)
    // 	      std::cerr << " - add " << ins.atoms[atoms[i]].name
    // 			<< " at " << cplex.getValue(pt_vars[i])
    // 			<< std::endl;
    // }
    // for (index_type i = 0; i < atoms.size(); i++)
    //   if (cplex.getValue(p_vars[i]) > 0.5) {
    // 	std::cerr << "atom " << ins.atoms[atoms[i]].name
    // 		  << " at " << cplex.getValue(pt_vars[i]) << std::endl;
    // 	for (index_type j = 0; j < add[i].size(); j++)
    // 	  if (cplex.getValue(add_vars[i][j]) > 0.5) {
    // 	    std::cerr << " - add by " << ins.actions[acts[add[i][j]]].name
    // 		      << std::endl;
    // 	  }
    //   }
  }
  catch (IloException& e) {
    std::cerr << "cplex exception: " << e << std::endl;
    exit(1);
  }
  catch (...) {
    std::cerr << "unknown exception" << std::endl;
    exit(1);
  }
  env.end();

  if (verbose_level > 0) {
    std::cerr << "val = " << sum << std::endl;
  }

  // return optimal cost
  return sum;
}

////
// old version
////
// NTYPE ILB::compute_relaxed_plan_extern
// (const index_set& init, const index_set& goal)
// {
//   NTYPE ub = compute_relaxed_plan_FF(init, goal);
//   std::cerr << "initial rp cost = " << PRINT_NTYPE(ub) << std::endl;
//   if (INFINITE(ub)) return ub;
// 
//   NTYPE sum = 0;
// 
//   // atoms in the encoding are relevant, non-initial atoms
//   index_set atoms(ins.n_atoms(), init);
//   atoms.intersect(relevant_atoms);
//   index_vec atom_map;
//   mapping::invert_map(atoms, atom_map);
//   // actions in the encoding are all relevant actions
//   index_set acts(relevant_actions);
//   index_vec act_map;
//   mapping::invert_map(acts, act_map);
// 
//   /// debug printing
//   // for (index_type i = 0; i < atoms.size(); i++)
//   //   std::cout << "atom" << i << " = " << ins.atoms[atoms[i]].name
//   // 	      << std::endl;
//   // for (index_type k = 0; k < acts.size(); k++)
//   //   std::cout << "act" << k << " = " << ins.actions[acts[k]].name
//   // 	      << std::endl;
// 
//   // standard init stuff
//   IloEnv env;
//   try {
//     IloModel model(env);
// 
//     // create variables
//     std::cerr << "creating variables..." << std::endl;
//     // variables of the encoding
//     IloBoolVarArray p_vars(env, atoms.size());
//     model.add(p_vars);
//     IloBoolVarArray a_vars(env, acts.size());
//     model.add(a_vars);
//     // add indexes into acts:
//     index_set_vec add(EMPTYSET, atoms.size());
//     IloBoolVarArray add_vars[atoms.size()];
//     IloBoolVarArray cl_vars[atoms.size()];
// 
//     for (index_type i = 0; i < atoms.size(); i++) {
//       for (index_type k = 0; k < acts.size(); k++)
// 	if (ins.actions[acts[k]].add.contains(atoms[i]))
// 	  add[i].insert(k);
//       IloBoolVarArray tmp(env, add[i].size());
//       add_vars[i] = tmp;
//       model.add(add_vars[i]);
//     }
//     for (index_type i = 0; i < atoms.size(); i++) {
//       IloBoolVarArray tmp(env, atoms.size());
//       cl_vars[i] = tmp;
//       model.add(cl_vars[i]);
//     }
// 
//     // add the objective: minimise action costs
//     IloNumArray costs(env, acts.size());
//     for (index_type k = 0; k < acts.size(); k++)
//       costs[k] = N_TO_D(ins.actions[acts[k]].cost);
//     model.add(IloObjective(env, IloScalProd(costs, a_vars))); 
// 
//     // create constraints
// 
//     // (0) goals: p for all p in G.
//     for (index_type i = 0; i < goal.size(); i++)
//       if (atom_map[goal[i]] != no_such_index) {
// 	model.add(p_vars[atom_map[goal[i]]] == 1);
//       }
// 
//     // (1) support: p -> \/ add(p, a)
//     std::cerr << "creating type 1 constraints..." << std::endl;
//     for (index_type i = 0; i < atoms.size(); i++) {
//       IloExpr expr(env);
//       for (index_type k = 0; k < add[i].size(); k++) {
// 	IloBoolVar &var = add_vars[i][k];
// 	expr += var;
//       }
//       model.add(p_vars[i] <= expr);
//     }
// 
//     // (2) support implies action: add(p, a) -> a
//     std::cerr << "creating type 2 constraints..." << std::endl;
//     for (index_type i = 0; i < atoms.size(); i++)
//       for (index_type k = 0; k < add[i].size(); k++) {
// 	model.add(add_vars[i][k] <= a_vars[add[i][k]]);
//       }
// 
//     // (3) action implies prec: a -> q forall q in pre(a).
//     std::cerr << "creating type 3 constraints..." << std::endl;
//     for (index_type k = 0; k < acts.size(); k++)
//       for (index_type i = 0; i < ins.actions[acts[k]].pre.size(); i++)
// 	if (atom_map[ins.actions[acts[k]].pre[i]] != no_such_index) {
// 	  index_type q = atom_map[ins.actions[acts[k]].pre[i]];
// 	  model.add(a_vars[k] <= p_vars[q]);
// 	}
// 
//     // (4) support implies cl: add(p, a) -> cl(q,p) forall q in pre(a).
//     std::cerr << "creating type 4 constraints..." << std::endl;
//     for (index_type i = 0; i < atoms.size(); i++)
//       for (index_type k = 0; k < add[i].size(); k++) {
// 	index_type a = acts[add[i][k]];
// 	for (index_type j = 0; j < ins.actions[a].pre.size(); j++)
// 	  if (atom_map[ins.actions[a].pre[j]] != no_such_index) {
// 	    index_type q = atom_map[ins.actions[a].pre[j]];
// 	    model.add(add_vars[i][k] <= cl_vars[q][i]);
// 	  }
//       }
// 
//     // (5) TC(CL)
//     std::cerr << "creating type 5 constraints..." << std::endl;
//     for (index_type i = 0; i < atoms.size(); i++)
//       for (index_type j = 0; j < atoms.size(); j++)
// 	if (j != i) {
// 	  // 2-cycle: cl[i][j] + cl[j][i] <= 1.
// 	  model.add(cl_vars[i][j] + cl_vars[j][i] <= 1);
// 	  for (index_type k = 0; k < atoms.size(); k++)
// 	    if ((i != k) && (j != k)) {
// 	      // transitivity: cl[i][k] + cl[k][j] <= cl[i][j] + 1
// 	      model.add(cl_vars[i][k] + cl_vars[k][j] <= cl_vars[i][j] + 1);
// 	    }
// 	}
// 
//     std::cerr << "solving..." << std::endl;
//     IloCplex cplex(model);
//     cplex.setParam(IloCplex::Threads, 1);
//     cplex.solve();
//     if (cplex.getStatus() != IloAlgorithm::Optimal) {
//       std::cerr << "error: not solved to optimality!" << std::endl;
//       std::cerr << "cplex status = " << cplex.getStatus() << std::endl;
//       exit(1);
//     }
// 
//     rp_actions.clear();
//     for (index_type k = 0; k < acts.size(); k++) {
//       double v = cplex.getValue(a_vars[k]);
//       if (v > 0.5) {
// 	sum += ins.actions[acts[k]].cost;
// 	rp_actions.insert(acts[k]);
//       }
//     }
//   }
//   catch (IloException& e) {
//     std::cerr << "cplex exception: " << e << std::endl;
//     exit(1);
//   }
//   catch (...) {
//     std::cerr << "unknown exception" << std::endl;
//     exit(1);
//   }
//   env.end();
// 
//   if (verbose_level > 0) {
//     std::cerr << "val = " << sum << std::endl;
//   }
// 
//   // return optimal cost
//   return sum;
// }
#endif // IMAI_FUKUNAGA


////
// OLD VERSION
////
// // standard back-chaining rp extraction, not guaranteed to be optimal.
// // h1c is h^1 computed with action costs (which may be zero), h1u is
// // h^1 computed with unit costs.
// void ILB::extract_rp
// (Heuristic* h1c, Heuristic* h1u,
//  const index_set& g, bool_vec& holds, index_vec& stack)
// {
//   weighted_vec<index_type, NTYPE> g_sort;
//   for (index_type i = 0; i < g.size(); i++)
//     if (!holds[g[i]]) {
//       g_sort.insert_increasing(g[i], h1c->eval(g[i]));
//     }
//   for (index_type i = 0; i < g_sort.size(); i++) {
//     index_set mca;
//     NTYPE c_min = POS_INF;
//     NTYPE l_min = POS_INF;
//     for (index_type k = 0; k < ins.atoms[g_sort[i].value].add_by.size(); k++) {
//       index_type a = ins.atoms[g_sort[i].value].add_by[k];
//       if (relevant_actions[a] &&
// 	  (ins.actions[a].pre.first_common_element(stack) == no_such_index)) {
// 	NTYPE c_a = h1c->eval(ins.actions[a].pre) + cost(a);
// 	if (c_a < c_min) {
// 	  c_min = c_a;
// 	  l_min = h1u->eval(ins.actions[a].pre);
// 	  mca.assign_singleton(a);
// 	}
// 	else if (c_a == c_min) {
// 	  NTYPE l_a = h1u->eval(ins.actions[a].pre);
// 	  if (l_a < l_min) {
// 	    l_min = l_a;
// 	    mca.assign_singleton(a);
// 	  }
// 	  else if (l_a == l_min) {
// 	    mca.insert(a);
// 	  }
// 	}
//       }
//     }
//     if (mca.empty()) {
//       std::cerr << "goal = " << g_sort[i]
// 		<< "." << ins.atoms[g_sort[i].value].name
// 		<< std::endl;
//       index_set reladd(ins.atoms[g_sort[i].value].add_by);
//       reladd.intersect(relevant_actions);
//       std::cerr << "relevant adds = " << reladd << std::endl;
//       std::cerr << "stack = " << stack << std::endl;
//       ((CostTable*)h1c)->write(std::cerr);
//     }
//     assert(!mca.empty());
//     index_type a = mca[0];
//     stack.append(g_sort[i].value);
//     extract_rp(h1c, h1u, ins.actions[a].pre, holds, stack);
//     stack.dec_length();
//     if (!holds[g_sort[i].value]) {
//       rp_actions.insert(a);
//       holds.insert(ins.actions[a].add);
//     }
//   }
// }
//
// NTYPE ILB::compute_relaxed_plan_FF
// (const index_set& init, const index_set& goal)
// {
//   // initialise H1 cost, and check for unsolvability
//   bool_vec s0(init, ins.n_atoms());
//   CostTable* h1c = new CostTable(ins, h1_stats);
//   h1c->compute_H1(cost, s0, &relevant_actions);
//   if (INFINITE(h1c->eval(goal))) {
//     delete h1c;
//     return POS_INF;
//   }
//   CostTable* h1u = new CostTable(ins, h1_stats);
//   h1c->compute_H1(UnitACF(), s0, &relevant_actions);
//   // call rp extraction
//   rp_actions.clear();
//   index_vec stack;
//   extract_rp(h1c, h1u, goal, s0, stack);
//   delete h1c;
//   delete h1u;
//   return cost.sum(rp_actions);
// }

NTYPE ILB::compute_relaxed_plan_FF
(const index_set& init, const index_set& goal)
{
  // initialise H1 cost, and check for unsolvability; use action cost +1
  // to ensure strictly well-founded min cost achievers.
  bool_vec s0(init, ins.n_atoms());
  CostTable* h1p = new CostTable(ins, h1_stats);
  PlusOne p1cost(cost);
  h1p->compute_H1(p1cost, s0, &relevant_actions);
  if (INFINITE(h1p->eval(goal))) {
    delete h1p;
    return POS_INF;
  }
  rp_actions.clear();
  index_vec open(goal);
  index_type next = 0;
  index_vec tail;
  while (next < open.size()) {
    assert(next < open.size());
    index_type p = open[next++];
    if (!s0[p]) {
      index_type mca = no_such_index;
      index_type n_pre = 0;
      for (index_type k = 0; k < ins.atoms[p].add_by.size(); k++) {
	index_type a = ins.atoms[p].add_by[k];
	if (relevant_actions[a])
	  if (h1p->eval(ins.actions[a].pre) + p1cost(a) == h1p->eval(p)) {
	    index_type n =
	      (ins.actions[a].pre.size() - ins.actions[a].pre.count_common(s0));
	    if (mca == no_such_index) {
	      mca = a;
	      n_pre = n;
	    }
	    else if (n < n_pre) {
	      mca = a;
	      n_pre = n;
	    }
	  }
      }
      assert(mca != no_such_index);
      tail.append(mca);
      for (index_type i = 0; i < ins.actions[mca].pre.size(); i++) {
	index_type q = ins.actions[mca].pre[i];
	if (!s0[q])
	  open.append(q);
      }
    }
    bool done = (tail.size() == 0);
    while (!done) {
      index_type a = tail[tail.size() - 1];
      if (s0.contains(ins.actions[a].pre)) {
	rp_actions.insert(a);
	s0.insert(ins.actions[a].add);
	tail.dec_length();
	if (tail.size() == 0) done = true;
      }
      else
	done = true;
    }
  }
  assert(tail.size() == 0);
  delete h1p;
  return cost.sum(rp_actions);
}

void ILB::rp_bb
(const index_set& goal,
 bool_vec& holds,
 index_type next,
 NTYPE lb,
 const index_vec& action_choice_order,
 const index_set& zero_cost_actions,
 bool_vec& rem_acts,
 bool_vec& plan,
 NTYPE acc,
 index_set& best,
 NTYPE& ub)
{
  if (next >= action_choice_order.size()) return;
  // std::cerr << " " << next;
  // assumption: goal test and bounds checking have been done,
  // so here we go straight to branching:
  index_type a = action_choice_order[next];
  rem_acts[a] = false;
  // if the action is redundant at this point, we don't include it
  // and we don't need to update anything; the action is redundant
  // if it does not add any relevant atom that is not already true.
  bool redundant = true;
  for (index_type i = 0; (i < ins.actions[a].add.size()) && redundant; i++)
    if (relevant_atoms[ins.actions[a].add[i]] && !holds[ins.actions[a].add[i]])
      redundant = false;
  if (redundant) {
    rp_bb(goal, holds, next + 1, lb, action_choice_order,
	  zero_cost_actions, rem_acts, plan, acc, best, ub);
  }
  else {
    // first choice is always to not include next action:
    NTYPE new_lb = h1->compute_lmcut(cost, holds, goal, &rem_acts);
    if ((acc + new_lb) < ub) {
      rp_bb(goal, holds, next + 1, new_lb, action_choice_order,
	    zero_cost_actions, rem_acts, plan, acc, best, ub);
    }
    // after first branch, ub may have lowered, so check if this node
    // is still worth exploring
    if ((acc + lb) >= ub) return;
    // if the action is applicable, the second branch is to include it:
    if (holds.contains(ins.actions[a].pre) && ((acc + cost(a)) < ub)) {
      plan[a] = true;
      bool_vec new_holds(holds);
      new_holds.insert(ins.actions[a].add);
      // update with zero-cost actions
      bool done = false;
      index_vec reset_zca;
      while (!done) {
	done = true;
	for (index_type k = 0; k < zero_cost_actions.size(); k++) {
	  index_type b = zero_cost_actions[k];
	  if (rem_acts[b]) {
	    if (new_holds.contains(ins.actions[b].add)) {
	      rem_acts[b] = false;
	      reset_zca.append(b);
	    }
	    else if (new_holds.contains(ins.actions[b].pre)) {
	      new_holds.insert(ins.actions[b].add);
	      plan[b] = true;
	      rem_acts[b] = false;
	      reset_zca.append(b);
	      done = false;
	    }
	  }
	}
      }
      // check if we've reached the goal
      if (new_holds.contains(goal)) {
	plan.copy_to(best);
	ub = acc + cost(a);
	std::cerr << "new rp cost = " << PRINT_NTYPE(ub) << std::endl;
      }
      // if not, check if we're still within bound, and if so recurse
      else {
	NTYPE new_lb = h1->compute_lmcut(cost, new_holds, goal, &rem_acts);
	if ((acc + cost(a) + new_lb) < ub) {
	  rp_bb(goal, new_holds, next + 1, new_lb, action_choice_order,
		zero_cost_actions, rem_acts, plan, acc + cost(a), best, ub);
	}
      }
      // reset
      for (index_type k = 0; k < reset_zca.size(); k++) {
	rem_acts[reset_zca[k]] = true;
	plan[reset_zca[k]] = false;
      }
      plan[a] = false;
    }
  }
  rem_acts[a] = true;
}

class precondition_cost_order : public lvector<index_type>::order
{
  Instance&  ins;
  Heuristic& pc;
public:
  precondition_cost_order(Instance& i, Heuristic& h) : ins(i), pc(h) { };
  virtual bool operator()(const index_type& v0, const index_type& v1) const;
};

bool precondition_cost_order::operator()
(const index_type& v0, const index_type& v1) const
{
  assert(v0 < ins.n_actions());
  assert(v1 < ins.n_actions());
  return (pc.eval_precondition(ins.actions[v0]) <
	  pc.eval_precondition(ins.actions[v1]));
}

NTYPE ILB::compute_relaxed_plan_BB
(const index_set& init, const index_set& goal)
{
  h1->compute_H1(UnitACF());
  if (INFINITE(h1->eval(goal))) {
    return POS_INF;
  }
  rp_actions.clear();
  rp_cost = POS_INF;

  if (prune_relaxed_irrelevant)
    compute_relevant(init, goal);
  if (prune_relaxed_dominated)
    remove_dominated_actions(init);
  if (stats.break_signal_raised())
    return rp_cost;

  index_vec action_choice_order;
  index_set zero_cost_actions;
  precondition_cost_order o(ins, *h1);
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (relevant_actions[k]) {
      if (IS_ZERO(cost(k)))
	zero_cost_actions.insert(k);
      else
	action_choice_order.insert_ordered(k, o);
    }

//   for (index_type k = 0; k < action_choice_order.size(); k++) {
//     std::cerr << action_choice_order[k] << ". "
// 	      << ins.actions[action_choice_order[k]].name << ", "
// 	      << PRINT_NTYPE(h1->eval_precondition(ins.actions[action_choice_order[k]]))
// 	      << std::endl;
//   }

  //h1->compute_H1(cost);
  bool_vec holds(init, ins.n_atoms());
  //index_vec stack;
  //extract_relaxed_plan(goal, holds, stack);
  //rp_cost = cost.sum(rp_actions);

  //holds.assign_value(false, ins.n_atoms());
  //holds.insert(init);
  bool_vec rem_acts(relevant_actions);
  NTYPE lb = h1->compute_lmcut(cost, holds, goal, &rem_acts);
  std::cerr << "initial rp cost = " << PRINT_NTYPE(rp_cost)
	    << ", lb = " << PRINT_NTYPE(lb) << std::endl;
  if (lb == rp_cost) {
    return rp_cost;
  }
  bool_vec plan(false, ins.n_actions());
  rp_bb(goal, holds, 0, lb, action_choice_order, zero_cost_actions,
	rem_acts, plan, 0, rp_actions, rp_cost);
  std::cerr << "final rp cost = " << PRINT_NTYPE(rp_cost) << std::endl;
  return rp_cost;
}

NTYPE ILB::rp_regress
(bool_vec& goals,
 bool_vec& supported,
 graph& prec,
 bool_vec& plan,
 NTYPE acc,
 NTYPE& ub,
 index_set& best)
{
  // select next open goal
  index_type next_goal = no_such_index;
  index_type min_domain = ins.n_actions();
  index_type min_degree = prec.size();
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (goals[i] && !supported[i]) {
      index_type n = 0;
      for (index_type k = 0; k < ins.atoms[i].add_by.size(); k++)
	if (relevant_actions[ins.atoms[i].add_by[k]] &&
	    !prec.adjacent(i, ins.n_atoms() + ins.atoms[i].add_by[k]))
	  n += 1;
      if (n < min_domain) {
	next_goal = i;
	min_domain = i;
	min_degree = prec.out_degree(i);
      }
      // use out-degree in prec as a proxy for distance from top-level goal
      else if ((n == min_domain) && (prec.out_degree(i) < min_degree)) {
	next_goal = i;
	min_domain = i;
	min_degree = prec.out_degree(i);
      }
    }
  // if no open goal, plan is done
  if (next_goal == no_such_index) {
    best = index_set(plan);
    std::cerr << "plan found: " << best << ", cost = " << acc << std::endl;
    assert(cost.sum(best) == acc);
    ub = acc;
    return acc;
  }
  //std::cerr << "next goal = " << next_goal << std::endl;
  // remaining possible establishers for the goal:
  bool_vec rem_est(true, ins.atoms[next_goal].add_by.size());
  NTYPE c_min = POS_INF;
  pair_set e0;
  while (true) {
    // select next establisher to try
    index_type next_action = no_such_index;
    bool in_plan = false;
    min_degree = prec.size();
    for (index_type k = 0; k < ins.atoms[next_goal].add_by.size(); k++)
      if (rem_est[k] && relevant_actions[ins.atoms[next_goal].add_by[k]] &&
	  !prec.adjacent(next_goal,
			 ins.n_atoms() + ins.atoms[next_goal].add_by[k]) &&
	  ((acc + cost(ins.atoms[next_goal].add_by[k])) < ub))
	if (plan[ins.atoms[next_goal].add_by[k]] && !in_plan) {
	  next_action = k;
	  in_plan = true;
	  min_degree = prec.in_degree(ins.n_atoms() + ins.atoms[next_goal].add_by[k]);
	}
	else if ((plan[ins.atoms[next_goal].add_by[k]] == in_plan) &&
		 (prec.in_degree(ins.n_atoms() + ins.atoms[next_goal].add_by[k])
		  < min_degree)) {
	  next_action = k;
	  min_degree = prec.in_degree(ins.n_atoms() + ins.atoms[next_goal].add_by[k]);
	}
    // if no more possible establishers, fail
    if (next_action == no_such_index) {
      prec.remove_edges(e0);
      return c_min;
    }
    //std::cerr << "next_action = " << next_action << ", rem = " << rem_est
    //	      << std::endl;
    // else, try next action:
    rem_est[next_action] = false;
    index_type a = ins.atoms[next_goal].add_by[next_action];
    supported[next_goal] = true;
    bool_vec new_goals(goals);
    new_goals.insert(ins.actions[a].pre);
    pair_set e1;
    prec.add_edge_to_transitive_closure(ins.n_atoms() + a, next_goal, e1);
    in_plan = plan[a];
    NTYPE new_acc = acc;
    if (!plan[a]) {
      plan[a] = true;
      new_acc += cost(a);
    }
    NTYPE c_a = rp_regress(new_goals, supported, prec, plan, new_acc, ub, best);
    c_min = MIN(c_min, c_a);
    // roll-back
    prec.remove_edges(e1);
    supported[next_goal] = false;
    plan[a] = in_plan;
    // add C1 constraint:
    prec.add_edge_to_transitive_closure(next_goal, ins.n_atoms() + a, e0);
  }
}

NTYPE ILB::compute_relaxed_plan_BDGBT
(const index_set& init, const index_set& goal)
{
  bool_vec goals(goal, ins.n_atoms());
  bool_vec supported(init, ins.n_atoms());
  bool_vec plan(false, ins.n_actions());
  rp_cost = POS_INF;
  rp_actions.clear();
  graph prec(ins.n_atoms() + ins.n_actions() + 1);
  for (index_type k = 0; k < ins.n_actions(); k++)
    for (index_type i = 0; i < ins.actions[k].pre.size(); i++) {
      prec.add_edge(ins.actions[k].pre[i], ins.n_atoms() + k);
    }
  rp_regress(goals, supported, prec, plan, 0, rp_cost, rp_actions);
  return rp_cost;
}

#ifdef USE_GECODE

class RPCP : public Gecode::Space
{
  const Instance& ins;
  const bool_vec& relevant_actions;
  const index_set& init;
  const index_set& goal;
  index_type lmax;
  Gecode::IntVarArray  atom_level;
  Gecode::BoolVarArray sub_goal;
  Gecode::IntVarArray  establisher;
  Gecode::IntVar       goal_level;
  Gecode::IntVarArray  action_level;
  Gecode::BoolVarArray in_plan;
  Gecode::IntVar       total_cost;

  struct branch_info {
    index_type next_goal;
    index_type next_action;
  };

 public:
  RPCP(const Instance& in,
       const bool_vec& ra,
       const index_set& is,
       const index_set& gs,
       const ACF& cost);
  RPCP(bool share, RPCP& s);
  virtual ~RPCP();
  virtual Gecode::Space* copy(bool share);
  virtual void constrain(const Space& best);
  unsigned int cb_status(branch_info& b) const;
  Gecode::ExecStatus cb_commit(const branch_info& b, unsigned int a);
  virtual void print(std::ostream& s);
  int plan_cost() const;
};

RPCP::RPCP
(const Instance& in,
 const bool_vec& ra,
 const index_set& is,
 const index_set& gs,
 const ACF& cost)
  : ins(in),
    relevant_actions(ra),
    init(is),
    goal(gs),
    atom_level(*this, in.n_atoms(), 0, in.n_actions()),
    sub_goal(*this, in.n_atoms(), 0, 1),
    establisher(*this, in.n_atoms(), 0, in.n_actions()),
    goal_level(*this, 0, in.n_actions()),
    action_level(*this, in.n_actions(), 0, in.n_actions()),
    in_plan(*this, in.n_actions(), 0, 1),
    total_cost(*this, 0, INT_MAX - 2)
{
  index_type m = relevant_actions.count(true);
  lmax = (ins.n_atoms() < m ? ins.n_atoms() : m);
  for (index_type i = 0; i < ins.n_atoms(); i++) {
    //Gecode::dom(*this, atom_level[i], 0, lmax);
    Gecode::dom(*this, establisher[i], 0, ins.n_actions());
    if (init.contains(i)) {
      // if initial, the atom does not need an establisher
      // and its level is zero
      Gecode::dom(*this, atom_level[i], 0, 0);
    }
    else {
      // atom_level[i] = action_level[establisher[i]]
      Gecode::element(*this, action_level, establisher[i], atom_level[i]);
      // establisher[i] in { relevant adders of i }
      int es[ins.atoms[i].add_by.size() + 1];
      index_type ne = 0;
      for (index_type k = 0; k < ins.atoms[i].add_by.size(); k++)
	if (relevant_actions[ins.atoms[i].add_by[k]])
	  es[ne++] = ins.atoms[i].add_by[k];
      Gecode::IntSet s(es, ne);
      Gecode::dom(*this, establisher[i], s);
    }
  }
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (relevant_actions[k]) {
      //Gecode::dom(*this, action_level, 0, lmax);
      // action_level[k] > atom_level[i] for all i in pre(k)
      for (index_type i = 0; i < ins.actions[k].pre.size(); i++) {
	Gecode::rel(*this, action_level[k], Gecode::IRT_GR,
		    atom_level[ins.actions[k].pre[i]]);
      }
      // action is in plan iff action_level < goal_level
      Gecode::rel(*this, action_level[k], Gecode::IRT_LQ, goal_level,
		  in_plan[k]);
      // in_plan[k] <= sub_goal[i] for all i in pre(k)
      for (index_type i = 0; i < ins.actions[k].pre.size(); i++) {
	Gecode::rel(*this, in_plan[k], Gecode::IRT_LQ,
		    sub_goal[ins.actions[k].pre[i]]);
      }
    }
  // goal level is the highest level where a goal is achieved;
  // all goals are subgoals.
  for (index_type i = 0; i < goal.size(); i++) {
    Gecode::rel(*this, goal_level, Gecode::IRT_GQ, atom_level[goal[i]]);
    Gecode::rel(*this, sub_goal[goal[i]], Gecode::IRT_EQ, 1);
  }
  // plan cost is the sum of costs of actions in the plan
  Gecode::IntArgs costs(ins.n_actions());
  for (index_type k = 0; k < ins.n_actions(); k++) {
    costs[k] = FLOOR_TO_INT(cost(k));
  }
  Gecode::linear(*this, costs, in_plan, Gecode::IRT_EQ, total_cost);

  // first, custom branching
  Gecode::BSB<RPCP,RPCP::branch_info>::post
    (*this, &RPCP::cb_status, &RPCP::cb_commit);
}

RPCP::RPCP(bool share, RPCP& s)
  : Gecode::Space(share, s),
    ins(s.ins),
    relevant_actions(s.relevant_actions),
    init(s.init),
    goal(s.goal)
{
  atom_level.update(*this, share, s.atom_level);
  sub_goal.update(*this, share, s.sub_goal);
  establisher.update(*this, share, s.establisher);
  goal_level.update(*this, share, s.goal_level);
  action_level.update(*this, share, s.action_level);
  in_plan.update(*this, share, s.in_plan);
  total_cost.update(*this, share, s.total_cost);
}

RPCP::~RPCP()
{
  // done
}

Gecode::Space* RPCP::copy(bool share)
{
  return new RPCP(share, *this);
}

unsigned int RPCP::cb_status(branch_info& b) const
{
  index_type next_goal = no_such_index;
  index_type min_domain = ins.n_actions();
  index_type max_level = 0;
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (sub_goal[i].one() && !init.contains(i)) {
      if (establisher[i].size() > 1) {
	if (establisher[i].size() < min_domain) {
	  next_goal = i;
	  min_domain = establisher[i].size();
	  max_level = atom_level[i].min();
	}
	else if ((establisher[i].size() == min_domain) &&
		 (atom_level[i].min() > max_level)) {
	  next_goal = i;
	  min_domain = establisher[i].size();
	  max_level = atom_level[i].min();
	}
      }
    }
  if (next_goal != no_such_index) {
    // std::cerr << "status: " << next_goal << " / " << b.next_action
    // 	      << std::endl;
    b.next_goal = next_goal;
    b.next_action = establisher[next_goal].min();
    return 2;
  }
  if (goal_level.min() < goal_level.max()) {
    b.next_goal = no_such_index;
    b.next_action = goal_level.min();
    return 1;
  }
  return 0;
}

Gecode::ExecStatus RPCP::cb_commit(const branch_info& b, unsigned int a)
{
  if (b.next_goal == no_such_index) {
    assert(a == 0);
    return SET_INT_EQ(goal_level, static_cast<int>(b.next_action));
  }
  assert(b.next_goal < ins.n_atoms());
  if (a == 0) {
    // std::cerr << "commit: establisher[" << b.next_goal << "] = "
    // 	      << b.next_action << std::endl;
    SET_INT_EQ(establisher[b.next_goal], static_cast<int>(b.next_action));
    for (index_type i = 0; i < ins.actions[b.next_action].pre.size(); i++)
      SET_BOOL_EQ(sub_goal[ins.actions[b.next_action].pre[i]], 1);
  }
  else {
    // std::cerr << "commit: " << b.next_action << " will not establish "
    // 	      << b.next_goal << std::endl;
    SET_INT_NEQ(establisher[b.next_goal], static_cast<int>(b.next_action));
    Gecode::rel(*this, atom_level[b.next_goal], Gecode::IRT_LE,
		action_level[b.next_action]);
  }
  return Gecode::ES_OK;
}

void RPCP::constrain(const Space& best)
{
  const RPCP& best_rp = static_cast<const RPCP&>(best);
  Gecode::rel(*this, total_cost, Gecode::IRT_LE, best_rp.total_cost.min());
}

void RPCP::print(std::ostream& s)
{
  s << "sub_goal: " << sub_goal << std::endl;
  s << "atom_level: " << atom_level << std::endl;
  s << "establisher: " << establisher << std::endl;
  s << "goal_level: " << goal_level << std::endl;
  s << "action_level: " << action_level << std::endl;
  s << "in_plan: " << in_plan << std::endl;
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (in_plan[k].one()) {
      s << k << ". " << ins.actions[k].name << std::endl;
    }
  s << "plan cost = " << total_cost << std::endl;
}

int RPCP::plan_cost() const
{
  return total_cost.min();
}

#ifdef HPLUS_USE_RP_EXTERN
NTYPE ILB::compute_relaxed_plan_extern
(const index_set& init, const index_set& goal)
{
  RPCP* root = new RPCP(ins, relevant_actions, init, goal, cost);
  Gecode::BAB<RPCP> engine(root);
  RPCP* best = NULL;
  while (RPCP* sol = engine.next()) {
    std::cerr << "found solution:" << std::endl;
    sol->print(std::cerr);
    if (best) delete best;
    best = sol;
  }
  if (best == NULL) return POS_INF;
  int val = best->plan_cost();
  delete best;
  return R_TO_N(val, 1);
}
#endif // HPLUS_USE_RP_EXTERN
#endif // USE_GECODE

bool ILB::write_hplus_wcnf
(const index_set& init, const index_set& goal, std::ostream& to,
 count_type climit)
{
  // an EncodingIndex provides methods to index variables
  // in the SAT encoding (and some other utils)
  struct EncodingIndex {
    const ILB& con; // the ILB object that is the context
    const index_set& init;
    const index_set& goal;
    index_set atoms; // relevant atoms
    index_set acts; // relevant actions
    index_vec atm_map; // back-map for atoms
    index_vec act_map; // back-map for actions
    index_set_vec adds; // rel. atoms added by rel. actions
    index_vec first_add;
    index_type n_add; // number of add-vars (total)
    index_type n_pre_ax; // number of precondition axioms
    index_type n_nz; // number of non-zero cost actions
    NTYPE sum_cost;

    EncodingIndex(const ILB& c, const index_set& i, const index_set& g)
      : con(c), init(i), goal(g),
	atoms(con.ins.n_atoms(), init), acts(con.relevant_actions)
    {
      // atoms in the encoding are relevant, non-initial atoms
      atoms.intersect(con.relevant_atoms);
      bool ok = mapping::invert_map(atoms, atm_map, con.ins.n_atoms());
      assert(ok);
      // actions in the encoding are all relevant actions
      assert(acts.size() == con.relevant_actions.count(true));
      ok = mapping::invert_map(acts, act_map, con.ins.n_actions());
      assert(ok);
      //std::cerr << "acts = " << acts << std::endl;
      //std::cerr << "map = " << act_map << std::endl;
      adds.assign_value(EMPTYSET, acts.size());
      first_add.assign_value(0, acts.size());
      // initialise the add-vars index, and compute sum-cost
      n_add = 0;
      sum_cost = 0;
      n_pre_ax = 0;
      n_nz = 0;
      for (index_type k = 0; k < acts.size(); k++) {
	assert(INTEGRAL(con.ins.actions[acts[k]].cost));
	if (!IS_ZERO(con.ins.actions[acts[k]].cost)) {
	  sum_cost += con.ins.actions[acts[k]].cost;
	  n_nz += 1;
	}
	if (k > 0)
	  first_add[k] = (first_add[k - 1] + adds[k - 1].size());
	for (index_type i = 0; i < con.ins.actions[acts[k]].add.size(); i++)
	  if (atm_map[con.ins.actions[acts[k]].add[i]] != no_such_index)
	    adds[k].insert(atm_map[con.ins.actions[acts[k]].add[i]]);
	assert(adds[k].size() ==
	       con.ins.actions[acts[k]].add.count_common(atoms));
	index_type q = con.ins.actions[acts[k]].pre.count_common(atoms);
	n_add += adds[k].size();
	n_pre_ax += (adds[k].size() * q);
      }
    }

    bool atom_in(index_type p) {
      assert(p < con.ins.n_atoms());
      return (atm_map[p] != no_such_index);
    }

    bool action_in(index_type a) {
      assert(a < con.ins.n_actions());
      return (act_map[a] != no_such_index);
    }

    NTYPE hard_cost() {
      return (sum_cost + 1);
    }

    // variable indexing functions: args are action/atom indices
    // in the original problem.

    index_type prop_var(index_type p) {
      assert(p < con.ins.n_atoms());
      assert(atm_map[p] != no_such_index);
      return (acts.size() + atm_map[p] + 1);
    }

    index_type act_var(index_type a) {
      assert(a < con.ins.n_actions());
      assert(act_map[a] != no_such_index);
      return (act_map[a] + 1);
    }

    index_type add_var(index_type a, index_type p) {
      assert(a < con.ins.n_actions());
      assert(act_map[a] != no_such_index);
      assert(act_map[a] < adds.size());
      assert(p < con.ins.n_atoms());
      assert(atm_map[p] != no_such_index);
      index_type i = adds[act_map[a]].first(atm_map[p]);
      assert(i != no_such_index);
      return (atoms.size() + acts.size() + first_add[act_map[a]] + i + 1);
    }

    index_type cl_var(index_type p, index_type q) {
      assert(p < con.ins.n_atoms());
      assert(atm_map[p] != no_such_index);
      assert(q < con.ins.n_atoms());
      assert(atm_map[q] != no_such_index);
      index_type i = (atoms.size() * atm_map[p]) + atm_map[q];
      return (atoms.size() + acts.size() + n_add + i + 1);
    }

    index_type n_vars() {
      return (acts.size() + atoms.size() + n_add +
	      (atoms.size() * atoms.size()));
    }

    count_type n_type_1() { return atoms.size(); }
    count_type n_type_2() { return n_pre_ax; }
    count_type n_type_3() { return n_add; }
    count_type n_type_4() { return (atoms.size() * atoms.size()); }
    count_type n_type_5() {
      return (atoms.size() * (atoms.size() - 1) * (atoms.size() - 1));
    }
    count_type n_type_6() { return atoms.size(); }
    count_type n_type_7() { return n_nz; }
    count_type n_type_8() {
      count_type n = 0;
      for (index_type i = 0; i < goal.size(); i++)
	if (atom_in(goal[i]))
	  n += 1;
      return n;
    }

    count_type n_type_10() {
      count_type n = 0;
      for (index_type i = 0; i < atoms.size(); i++) {
	index_type a = con.ins.atoms[atoms[i]].add_by.count_common(acts);
	n += ((a * (a - 1)) / 2);
      }
      return n;
    }

    count_type n_clauses() {
      return (n_type_1() + n_type_2() + n_type_3() + n_type_4()
	      + n_type_5() + n_type_6() + n_type_7() + n_type_8()
	      // redundant constraints:
#ifdef WSAT_REDUNDANT_CONSTRAINTS
	      + n_type_4() + n_type_10()
#endif
	      );
    }
  };

  EncodingIndex enc(*this, init, goal);

  to << "c hplus problem #" << calls_to_hplus << std::endl;
  to << "c " << enc.atoms.size() << " atoms, "
     << enc.acts.size() << " actions" << std::endl;
  // for (index_type i = 0; i < ins.n_atoms(); i++)
  //   if (enc.atom_in(i)) {
  //     to << "c " << enc.prop_var(i) << " = atom #" << i << std::endl;
  //   }
  // for (index_type k = 0; k < ins.n_actions(); k++)
  //   if (enc.action_in(k)) {
  //     to << "c " << enc.act_var(k) << " = action #" << k << std::endl;
  //   }
  to << "p wcnf " << enc.n_vars() << " " << enc.n_clauses()
     << " " << enc.hard_cost() << std::endl;
  if (climit > 0)
    if (enc.n_clauses() > climit) {
      to << "error: too many clauses!!" << std::endl;
      return false;
    }

  // (1) p -> \/ add(a,p)
  count_type n = 0;
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (enc.atom_in(i)) {
      to << enc.hard_cost() << " -" << enc.prop_var(i);
      for (index_type k = 0; k < ins.atoms[i].add_by.size(); k++)
	if (enc.action_in(ins.atoms[i].add_by[k]))
	  to << " " << enc.add_var(ins.atoms[i].add_by[k], i);
      to << " 0" << std::endl;
      n += 1;
    }
  assert(n == enc.n_type_1());

  // (2) add(a,p) -> CL(q,p) for all q in pre(a)
  n = 0;
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (enc.action_in(k))
      for (index_type i = 0; i < ins.actions[k].add.size(); i++)
	if (enc.atom_in(ins.actions[k].add[i]))
	  for (index_type j = 0; j < ins.actions[k].pre.size(); j++)
	    if (enc.atom_in(ins.actions[k].pre[j])) {
	      to << enc.hard_cost()
		 << " -" << enc.add_var(k, ins.actions[k].add[i])
		 << " " << enc.cl_var(ins.actions[k].pre[j],
				      ins.actions[k].add[i])
		 << " 0" << std::endl;
	      n += 1;
	    }
  assert(n == enc.n_type_2());

  // (3) add(a,p) -> a
  n = 0;
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (enc.action_in(k))
      for (index_type i = 0; i < ins.actions[k].add.size(); i++)
	if (enc.atom_in(ins.actions[k].add[i])) {
	  to << enc.hard_cost()
	     << " -" << enc.add_var(k, ins.actions[k].add[i])
	     << " " << enc.act_var(k)
	     << " 0" << std::endl;
	  n += 1;
	}
  assert(n == enc.n_type_3());

  // (4) CL(p,q) -> p
  n = 0;
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (enc.atom_in(i))
      for (index_type j = 0; j < ins.n_atoms(); j++)
	if (enc.atom_in(j)) {
	  to << enc.hard_cost()
	     << " -" << enc.cl_var(i, j)
	     << " " << enc.prop_var(i) << " 0" << std::endl;
	  n += 1;
	}
  assert(n == enc.n_type_4());

  // (5) TC(CL)
  n = 0;
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (enc.atom_in(i))
      for (index_type j = 0; j < ins.n_atoms(); j++)
	if (enc.atom_in(j) && (j != i))
	  for (index_type k = 0; k < ins.n_atoms(); k++)
	    if (enc.atom_in(k) && (k != j)) {
	      to << enc.hard_cost()
		 << " -" << enc.cl_var(i, j)
		 << " -" << enc.cl_var(j, k)
		 << " " << enc.cl_var(i, k) << " 0" << std::endl;
	      n += 1;
	    }
  assert(n == enc.n_type_5());

  // (6) -CL(p,p)
  n = 0;
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (enc.atom_in(i)) {
      to << enc.hard_cost()
	 << " -" << enc.cl_var(i, i) << " 0" << std::endl;
      n += 1;
    }
  assert(n == enc.n_type_6());

  // (7) soft: -a for all (relevant) a with cost(a) > 0.
  n = 0;
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (enc.action_in(k))
      if (!IS_ZERO(ins.actions[k].cost)) {
	to << ins.actions[k].cost
	   << " -" << enc.act_var(k) << " 0" << std::endl;
	n += 1;
      }
  assert(n == enc.n_type_7());

  // (8) p for all (relevant) p in goal
  n = 0;
  for (index_type i = 0; i < goal.size(); i++)
    if (enc.atom_in(goal[i])) {
      to << enc.hard_cost()
	 << " " << enc.prop_var(goal[i]) << " 0" << std::endl;
      n += 1;
    }
  assert(n == enc.n_type_8());

  // redundant constraints:

#ifdef WSAT_REDUNDANT_CONSTRAINTS
  // (9) CL(p,q) -> q
  n = 0;
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (enc.atom_in(i))
      for (index_type j = 0; j < ins.n_atoms(); j++)
	if (enc.atom_in(j)) {
	  to << enc.hard_cost()
	     << " -" << enc.cl_var(i, j)
	     << " " << enc.prop_var(j) << " 0" << std::endl;
	  n += 1;
	}
  assert(n == enc.n_type_4());

  // (10) add(a,p) -> -add(a',p)
  n = 0;
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (enc.atom_in(i))
      for (index_type k1 = 0; k1 < ins.atoms[i].add_by.size(); k1++)
	if (enc.action_in(ins.atoms[i].add_by[k1]))
	  for (index_type k2 = k1 + 1; k2 < ins.atoms[i].add_by.size(); k2++)
	    if (enc.action_in(ins.atoms[i].add_by[k2])) {
	      to << enc.hard_cost()
		 << " -" << enc.add_var(ins.atoms[i].add_by[k1], i)
		 << " -" << enc.add_var(ins.atoms[i].add_by[k2], i)
		 << " 0" << std::endl;
	      n += 1;
	    }
  assert(n == enc.n_type_10());
#endif
  return true;
}

void ILB::save_hplus_wcnf
(const index_set& init, const index_set& goal, const char* pname,
 count_type climit)
{
  if (prune_relaxed_irrelevant)
    compute_relevant(init, goal);
  if (prune_relaxed_dominated)
    remove_dominated_actions(init);
  if (stats.break_signal_raised())
    return;
  if (prune_relaxed_irrelevant || prune_relaxed_dominated) {
    count_type n_rel = relevant_actions.count(true);
    relevant_actions_ratio += (n_rel/(double)ins.n_actions());
    if (verbose_level > 0)
      std::cerr << n_rel << " of " << ins.n_actions()
		<< " actions considered relevant" << std::endl;
  }
  else {
    relevant_actions_ratio += 1;
  }
  std::string wcnf_file_name =
    std::string("hplus-") + std::string(pname) + std::string(".wcnf");
  std::ofstream wcnf_file(wcnf_file_name.c_str());
  bool ok = write_hplus_wcnf(init, goal, wcnf_file, climit);
  wcnf_file.close();
  if (!ok) return;
  std::string index_file_name =
    std::string("hplus-") + std::string(pname) + std::string(".acts");
  std::ofstream index_file(index_file_name.c_str());
  index_type n = 0;
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (relevant_actions[k]) {
      index_file << "/" << n << "/" << ins.actions[k].name << "/" << std::endl;
      n += 1;
    }
  index_file.close();
}

NTYPE ILB::hplus
(const index_set& init, const index_set& goal, NTYPE bound, Plan** sol)
{
  calls_to_hplus += 1;
  stats.start();

  // check that the goal is actually reachable, because a lot of
  // stuff might fail if it's not.
  reach.recompute(init);
  if (reach.unreachable(goal))
    return POS_INF;

#ifdef INIT_ZERO_FILL
  bool_vec zero_cost_actions(false, ins.n_actions());
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (IS_ZERO(cost(k)))
      zero_cost_actions[k] = true;
  reach.recompute(init, zero_cost_actions);
  // remove zero cost actions that have not been applied
  for (index_type k = 0; k < ins.n_actions(); k++)
    if (zero_cost_actions[k])
      if (reach.unreachable(ins.actions[k].pre))
	zero_cost_actions[k] = false;
  if (verbose_level > 0)
    std::cerr << zero_cost_actions.count(true)
	      << " zero cost actions applied to init state"
	      << std::endl;
  index_set einit(reach.reachable());
  const index_set& _init = einit;
#else
  const index_set& _init = init;  
#endif

  init_preprocessing();
  // first pass
  if (prune_relaxed_irrelevant) {
    compute_relevant(_init, goal);
    if (stats.break_signal_raised()) return 0;
  }
  index_type nrel1 = relevant_actions.count(true);
  if (prune_relaxed_dominated) {
    remove_dominated_actions(_init);
    if (stats.break_signal_raised()) return 0;
  }
  index_type nrel2 = relevant_actions.count(true);
  // verify goal is still reachable
  reach.recompute(_init, relevant_actions);
  assert(!reach.unreachable(goal));
  bool done = (nrel2 == nrel1) | !iterated_preprocessing;
  while (!done) {
    bool lmch = compute_FA_and_ALM(init);
    index_type nrp_before = relevant_atoms.count(true);
    compute_relevant(_init, goal);
    if (stats.break_signal_raised()) return 0;
    index_type nrp_after = relevant_atoms.count(true);
    nrel1 = relevant_actions.count(true);
    lmch = compute_FA_and_ALM(init) ? true : lmch;
    if (lmch | (nrp_after < nrp_before))
      remove_dominated_actions(_init);
    if (stats.break_signal_raised()) return 0;
    nrel2 = relevant_actions.count(true);
    done = (nrel2 == nrel1);
  }
  if (prune_relaxed_irrelevant || prune_relaxed_dominated) {
    relevant_actions_ratio += (nrel2/(double)ins.n_actions());
    if (verbose_level > 0) {
      std::cerr << nrel2 << " of " << ins.n_actions()
		<< " actions considered relevant" << std::endl;
      if (verbose_level > 2) {
	for (index_type k = 0; k < ins.n_actions(); k++)
	  if (relevant_actions[k])
	    std::cerr << " " << ins.actions[k].name << std::endl;
      }
    }
  }
  else {
    relevant_actions_ratio += 1;
  }

  rp_actions.clear();
  // save_hplus_wcnf(_init, goal);
#ifdef HPLUS_USE_RP_EXTERN
  NTYPE v = compute_relaxed_plan_extern(_init, goal);
#else
  NTYPE v = compute_relaxed_plan_ILA(_init, goal, bound);
#endif
  stats.stop();
#ifdef HPLUS_PRINT_STATS
  print_hplus_stats(std::cerr);
#endif
#ifdef INIT_ZERO_FILL
  // must add zero cost actions that were applied before computing
  // the relaxed plan:
  if (!stats.break_signal_raised() && FINITE(v)) {
    for (index_type k = 0; k < ins.n_actions(); k++)
      if (zero_cost_actions[k])
	rp_actions.insert(k);
  }
#endif
  if (!stats.break_signal_raised() && FINITE(v) &&
      ((sol != NULL) || (verbose_level > 1))) {
    if (verbose_level > 0)
      std::cerr << "computing non-redundant rp..." << std::endl;
    compute_non_redundant_rp(_init, goal, true);
    /// Annoying: we want to do an rp validation at verbose_level 2+
    /// if h+ was called from test_ilb, but not if it was called from
    /// h++ or h+/ce...
    // if (verbose_level > 1) {
    //   validate_rp(rp_actions);
    // }
    if (sol != NULL) {
      if (verbose_level > 0)
	std::cerr << "ordering rp..." << std::endl;
      index_vec seq;
      bool ok = reach.order_relaxed_plan(_init, goal, rp_actions, seq);
      assert(ok);
      ActionSequence* plan = new ActionSequence();
      for (index_type k = 0; k < seq.length(); k++)
	plan->insert(seq[k]);
      plan->end();
      *sol = plan;
    }
  }
  return v;
}

NTYPE ILB::hplus_with_selected_actions
(const index_set& init,
 const index_set& goal,
 const index_set& acts,
 NTYPE bound,
 Plan** sol)
{
  calls_to_hplus += 1;
  stats.start();

  // check that the goal is actually reachable, because a lot of
  // stuff might fail if it's not.
  reach.recompute(init, acts);
  if (reach.unreachable(goal))
    return POS_INF;

  init_preprocessing();
  relevant_actions.assign_value(false);
  relevant_atoms.assign_value(false);
  for (index_type i = 0; i < acts.size(); i++) {
    assert(acts[i] < ins.n_actions());
    relevant_actions[acts[i]] = true;
    relevant_atoms.insert(ins.actions[acts[i]].pre);
    relevant_atoms.insert(ins.actions[acts[i]].add);
  }
  relevant_actions_ratio += (acts.size()/(double)ins.n_actions());
  std::cerr << relevant_actions.count(true) << " of " << ins.n_actions()
	    << " actions and " << relevant_atoms.count(true) << " of "
	    << ins.n_atoms() << " atoms relevant" << std::endl;

  rp_actions.clear();
#ifdef HPLUS_USE_RP_EXTERN
  NTYPE v = compute_relaxed_plan_extern(init, goal);
#else
  NTYPE v = compute_relaxed_plan_ILA(init, goal, bound);
#endif
  stats.stop();
#ifdef HPLUS_PRINT_STATS
  print_hplus_stats(std::cerr);
#endif
  if (!stats.break_signal_raised() && FINITE(v) &&
      ((sol != NULL) || (verbose_level > 1))) {
    if (verbose_level > 0)
      std::cerr << "computing non-redundant rp..." << std::endl;
    compute_non_redundant_rp(init, goal, true);
    if (sol != NULL) {
      if (verbose_level > 0)
	std::cerr << "ordering rp..." << std::endl;
      index_vec seq;
      bool ok = reach.order_relaxed_plan(init, goal, rp_actions, seq);
      assert(ok);
      ActionSequence* plan = new ActionSequence();
      for (index_type k = 0; k < seq.length(); k++)
	plan->insert(seq[k]);
      plan->end();
      *sol = plan;
    }
  }
  return v;
}

///
// hplus_with_ce
///

index_type ILB::get_ce_atom
(const Instance& irce, index_type first_rce_atom, index_type act)
{
  assert(act < irce.n_actions());
  index_type atom = no_such_index;
  for (index_type k = 0; k < irce.actions[act].pre.size(); k++)
    if (irce.actions[act].pre[k] >= first_rce_atom) {
      assert(atom == no_such_index);
      atom = irce.actions[act].pre[k];
    }
  assert(atom != no_such_index);
  return atom;
}

bool ILB::sequence_dg
(const Instance& irce, index_type first_rce_atom,
 index_set_graph& dg, index_vec& seq)
{
  bool_vec rem(true, dg.size());
  bool_vec s(irce.init_atoms, irce.n_atoms());
  // rce-flag atoms should not be considered, so set them all true
  for (index_type i = first_rce_atom; i < irce.n_atoms(); i++)
    s[i] = true;
  bool done = false;
  seq.clear();
  while (!done) {
    done = true;
    // search for an applicable remaining node
    for (index_type k = 0; (k < dg.size()) && done; k++)
      if (rem[k]) {
	// node with empty label is the goal, handle separately
	if (dg.node_label(k).empty()) {
	  if (s.contains(irce.goal_atoms)) {
	    seq.append(k);
	    rem[k] = false;
	    done = false;
	  }
	}
	// normal node with set of actions
	else {
	  bool app = true;
	  for (index_type i = 0; (i < dg.node_label(k).size()) && app; i++)
	    if (!s.contains(irce.actions[dg.node_label(k)[i]].pre))
	      app = false;
	  if (app) {
	    for (index_type i = 0; i < dg.node_label(k).size(); i++)
	      s.insert(irce.actions[dg.node_label(k)[i]].add);
	    seq.append(k);
	    rem[k] = false;
	    done = false;
	  }
	}
      }
  }
  // exited loop; if any nodes remain, sequencing has failed:
  return (rem.count(true) == 0);
}

bool ILB::schedule_ce_rec
(const Instance& irce, index_type first_rce_atom, index_type first_rce_action,
 index_set_graph& dg, index_set& failed)
{
  // find floating ce's, and their candidate action nodes:
  index_vec ce_nodes;
  index_set_vec ca_nodes;
  for (index_type i = 0; i < dg.size(); i++)
    if ((dg.node_label(i).size() == 1) &&
	(dg.node_label(i)[0] >= first_rce_action)) {
      index_type ce_act = dg.node_label(i)[0];
      index_type ce_atom = get_ce_atom(irce, first_rce_atom, ce_act);
      ce_nodes.append(i);
      index_set cands;
      for (index_type j = 0; j < dg.size(); j++)
	if (!dg.node_label(j).empty()) {
	  index_type a = dg.node_label(j)[0];
	  if (a < first_rce_action)
	    if (irce.actions[a].add.contains(ce_atom))
	      cands.insert(j);
	}
      assert(!cands.empty());
      ca_nodes.append(cands);
    }
  assert(ce_nodes.size() == ca_nodes.size());
  // if there are no ce's, we're done
  if (ce_nodes.size() == 0) return true;
  bool ok = true;
  index_type ca_min = 0;
  index_type ca_min_size = ca_nodes[0].size();
  // first, prune impossible candidates:
  for (index_type k = 0; k < ce_nodes.size(); k++) {
    index_set poss_ca;
    index_type ce_node = ce_nodes[k];
    for (index_type i = 0; i < ca_nodes[k].size(); i++) {
      // if the ce is unordered wrt the ca, it's possible
      if (!dg.reachable(ca_nodes[k][i], ce_node) &&
	  !dg.reachable(ce_node, ca_nodes[k][i]))
	poss_ca.insert(ca_nodes[k][i]);
    }
    ca_nodes[k] = poss_ca;
    // if ce has no possible candidates, it's failed
    if (ca_nodes[k].empty()) {
      failed.insert(dg.node_label(ce_node)[0]);
      ok = false;
    }
    else if (ca_nodes[k].size() < ca_min_size) {
      ca_min = k;
      ca_min_size = ca_nodes[k].size();
    }
  }
  if (!ok) {
    assert(!failed.empty());
    return false;
  }
  // now, check ce with minimal candidate set:
  assert(ca_min < ce_nodes.size());
  // if there is only one candidate:
  if (ca_nodes[ca_min].size() == 1) {
    index_type ce_node = ce_nodes[ca_min];
    index_type ce_act = dg.node_label(ce_node)[0];
    index_type ca_node = ca_nodes[ca_min][0];
    // merge the nodes and try sequencing the resulting graph:
    index_set_graph new_dg(dg);
    new_dg.contract(ca_node, ce_node);
    assert(new_dg.acyclic());
    index_vec seq;
    ok = sequence_dg(irce, first_rce_atom, new_dg, seq);
    // if this works, we recurse with the new graph:
    if (ok) {
      ok = schedule_ce_rec(irce, first_rce_atom, first_rce_action,
			   new_dg, failed);
      // if that worked solution is now in new_dg
      if (ok) {
	dg = new_dg;
      }
      return ok;
    }
    // else, fail:
    else {
      failed.insert(ce_act);
      return false;
    }
  }
  // if there is more than one candidate, we have to try them all:
  else {
    index_type ce_node = ce_nodes[ca_min];
    index_type ce_act = dg.node_label(ce_node)[0];
    for (index_type k = 0; k < ca_nodes[ca_min].size(); k++) {
      index_type ca_node = ca_nodes[ca_min][k];
      // merge the nodes and try sequencing the resulting graph:
      index_set_graph new_dg(dg);
      new_dg.contract(ca_node, ce_node);
      assert(new_dg.acyclic());
      index_vec seq;
      ok = sequence_dg(irce, first_rce_atom, new_dg, seq);
      // if this works, we recurse with the new graph:
      if (ok) {
	ok = schedule_ce_rec(irce, first_rce_atom, first_rce_action,
			     new_dg, failed);
	failed.clear();
	// if that worked, return success; note that solution is now
	// in new_dg, so we have to copy it.
	if (ok) {
	  dg = new_dg;
	  return true;
	}
      }
    }
    // if no candidate worked, fail:
    failed.insert(ce_act);
    return false;
  }
}

bool ILB::schedule_ce
(const Instance& irce, index_type first_rce_atom, index_type first_rce_action,
 index_set_graph& dg, index_set& failed)
{
  sce_stats.start();
  //std::cerr << "1. dg = " << dg << std::endl;
  //std::cerr << "first rce atom = " << first_rce_atom << std::endl;
  //std::cerr << "first rce action = " << first_rce_action << std::endl;
  // remove rce atoms from dg edges
  for (index_type i = 0; i < dg.size(); i++)
    for (index_type j = 0; j < dg.size(); j++)
      if (dg.adjacent(i, j)) {
	index_set l = dg.edge_label(i, j);
	assert(!l.empty());
	dg.edge_label(i, j).clear();
	for (index_type k = 0; k < l.size(); k++)
	  if (l[k] < first_rce_atom)
	    dg.edge_label(i, j).insert(l[k]);
	// if an edge label becomes empty, remove the edge
	if (dg.edge_label(i, j).empty())
	  dg.remove_edge(i, j);
      }
  //std::cerr << "2. dg = " << dg << std::endl;
  bool ok =
    schedule_ce_rec(irce, first_rce_atom, first_rce_action, dg, failed);
  sce_stats.stop();
  return ok;
}

void ILB::schedule_ce_nonopt
(const Instance& irce,
 index_type first_rce_atom,
 index_type first_rce_action,
 index_set_graph& dg,
 index_set& failed)
{
  // try to schedule the current dg
  failed.clear();
  bool ok = schedule_ce(irce, first_rce_atom, first_rce_action, dg, failed);
  while (!ok) {
    // if this does not work, find the failed ce node...
    assert(!failed.empty());
    index_type ce_act = failed[0];
    assert(ce_act >= first_rce_action);
    index_type ce_node = no_such_index;
    for (index_type i = 0; (i < dg.size()) && (ce_node == no_such_index); i++)
      if (dg.node_label(i).size() == 1)
	if (dg.node_label(i)[0] == ce_act)
	  ce_node = i;
    // (there must be one)
    //std::cerr << "ce_act = " << ce_act << std::endl;
    //std::cerr << "dg = " << dg << std::endl;      
    assert(ce_node != no_such_index);
    assert(ce_node < dg.size());
    index_type ce_atom = get_ce_atom(irce, first_rce_atom, ce_act);
    // ...find an action that the ce is associated with...
    index_type assoc_act = no_such_index;
    for (index_type k = 0;
	 (k < first_rce_action) && (assoc_act == no_such_index); k++)
      if (irce.actions[k].add.contains(ce_atom))
	assoc_act = k;
    // (there must be one)
    assert((assoc_act != no_such_index) && (assoc_act < first_rce_action));
    // ...and add that action to the failed dg node:
    dg.node_label(ce_node).insert(assoc_act);
    // then try scheduling again:
    failed.clear();
    ok = schedule_ce(irce, first_rce_atom, first_rce_action, dg, failed);
  }
  assert(failed.empty());
}

NTYPE ILB::hplus_with_ce
(const index_set& init, const index_set& goal, NTYPE bound,
 index_vec& rp, index_set_vec& rp_ce)
{
  calls_to_hplus_ce += 1;
  count_type iterations = 0;
  count_type n_compiled = 0;
  stats.start();
  Instance irce;
  irce.copy(ins);
  // mapping to original action indices; this is needed because
  // compiling out ce's creates new actions:
  index_vec original(no_such_index, irce.n_actions());
  mapping::identity_map(irce.n_actions(), original);
  // keep track of which ce's have been compiled into new actions
  // (sets are indices into the cadd array below):
  index_set_vec compiled_ce(EMPTYSET, irce.n_actions());
  // all conditional adds, the original action they belong to, their
  // position in that actions list of cadds, and their corresponding
  // "a_done" atom:
  lvector<rule> cadd;
  index_vec cadd_action;
  index_vec cadd_index;
  index_vec rce_atom;
  index_type first_rce_atom = irce.n_atoms();
  index_type first_rce_action = irce.n_actions();
  for (index_type k = 0; k < irce.n_actions(); k++)
    if (!irce.actions[k].cadd.empty()) {
      Instance::Atom& a_done = irce.new_atom(irce.actions[k].name);
      irce.actions[k].add.insert(a_done.index);
      for (index_type i = 0; i < irce.actions[k].cadd.size(); i++) {
	cadd.append(irce.actions[k].cadd[i]);
	cadd_action.append(k);
	cadd_index.append(i);
	rce_atom.append(a_done.index);
      }
      assert(cadd.size() == cadd_action.size());
    }
  // remaining (not compiled) ce's
  bool_vec cadd_rem(true, cadd.size());
  // track how many ce's are compiled per action (for stats
  // purposes only)
  index_vec act_ce_compiled(0, ins.n_actions());
  NTYPE v_max = 0;
  bool solved = false;
  while (!solved) {
    iterations += 1;
    if (verbose_level > 0)
      std::cerr << "constructing relaxed conditional effect problem..."
		<< std::endl;
    // construct rce actions for remaining conditional effects
    first_rce_action = irce.n_actions();
    // local index, mapping rce actions back to cadd array
    index_vec rce_index;
    for (index_type k = 0; k < cadd.size(); k++)
      if (cadd_rem[k]) {
	assert(rce_atom[k] != no_such_index);
	Instance::Action& a_ce =
	  irce.new_action(irce.actions[cadd_action[k]].name);
	a_ce.pre = cadd[k].antecedent;
	a_ce.pre.insert(rce_atom[k]);
	a_ce.add.assign_singleton(cadd[k].consequent);
	a_ce.cost = 0;
	rce_index.append(k);
	assert((rce_index.size() + first_rce_action) == irce.n_actions());
      }
    irce.clear_cross_reference();
    irce.cross_reference();
    if (verbose_level > 0) {
      std::cerr << first_rce_action << " actions and "
		<< irce.n_actions() - first_rce_action
		<< " remaining ce's" << std::endl;
      if (verbose_level > 4) {
	irce.print(std::cerr);
	std::cerr << "first rce action = " << first_rce_action << std::endl;
	std::cerr << "original = " << original << std::endl;
	std::cerr << "compiled = " << compiled_ce << std::endl;
      }
    }
    // compute h+ on rce instance:
    AnyACF rce_cost(ins.n_actions(), cost);
    rce_cost.extend_to(irce.n_actions(), 0);
    for (index_type i = ins.n_actions(); i < first_rce_action; i++)
      rce_cost.set_cost(i, cost(original[i]));
    if (verbose_level > 0)
      std::cerr << "computing relaxed plan..." << std::endl;
    ILB hp_rce(irce, rce_cost, NULL, stats);
    hp_rce.copy_options(*this, -1);
    hp_rce.is_sub = true;
#ifdef HOFFMANN
    v_max = hp_rce.compute_relaxed_plan_FF(init, goal);
    calls_to_hplus += 1;
#else
    NTYPE v_new = hp_rce.hplus(init, goal, bound, NULL);
    assert((v_new >= v_max) || stats.break_signal_raised());
    v_max = MAX(v_max, v_new);
#ifdef TRACE_PRINT_HLB
    new_lb(v_max);
#endif
#endif
    add_hplus_stats(hp_rce);
#ifdef HPLUS_PRINT_STATS_EVERY_ITERATION
    print_hplus_stats(std::cerr);
#endif
    if (stats.break_signal_raised() || (v_max >= bound) || INFINITE(v_max)) {
      stats.stop();
      double size_increase = first_rce_action/(double)ins.n_actions();
      cce_actions_ratio += size_increase;
      return v_max;
    }
    if (verbose_level == 1)
      std::cerr << "rp cost = " << PRINT_NTYPE(v_max) << std::endl;
    // post-process and try to schedule the rp:
    if (verbose_level > 0)
      std::cerr << "computing non-redundant relaxed plan..." << std::endl;
#ifdef HOFFMANN
    hp_rce.compute_non_redundant_rp(init, goal, false);
#else
    hp_rce.compute_non_redundant_rp(init, goal, true);
#endif
    if (verbose_level > 2) {
      for (index_type i = 0; i < hp_rce.rp_actions.size(); i++)
	std::cerr << hp_rce.rp_actions[i] << ". "
		  << irce.actions[hp_rce.rp_actions[i]].name << std::endl;
    }
    if (stats.break_signal_raised()) {
      stats.stop();
      return v_max;
    }
    rp_actions.assign_copy(hp_rce.relaxed_plan());
    if (verbose_level > 0)
      std::cerr << "computing dependency graph..." << std::endl;
    hp_rce.compute_dependency_graph(init, goal, rp_dg);
    if (verbose_level > 0)
      std::cerr << "scheduling conditional effects..." << std::endl;
    index_set failed;
#ifdef HOFFMANN
    schedule_ce_nonopt(irce, first_rce_atom, first_rce_action, rp_dg, failed);
    bool ok = true;
#else
    bool ok =
      schedule_ce(irce, first_rce_atom, first_rce_action, rp_dg, failed);
#endif
    if (verbose_level > 2) {
      std::cerr << "dependency graph = " << rp_dg << std::endl;
      std::cerr << "failed = " << failed << std::endl;
    }
    if (stats.break_signal_raised()) {
      stats.stop();
      return v_max;
    }
    // if scheduling worked, we're done:
    if (ok) {
      assert(failed.empty());
      std::cerr << "valid rp cost = " << PRINT_NTYPE(v_max)
		<< " (" << stats.time() << " seconds)" << std::endl;
      //hp_rce.validate_rp(hp_rce.rp_actions);
      solved = true;
      index_vec seq;
      bool ok = sequence_dg(irce, first_rce_atom, rp_dg, seq);
      assert(ok);
      // printing the relaxed plan
      if (verbose_level > 1) {
	std::cerr << "sequenced relaxed plan:" << std::endl;
	for (index_type k = 0; k < seq.length(); k++) {
	  if (rp_dg.node_label(seq[k]).empty()) {
	    std::cerr << " <Goal>" << std::endl;
	  }
	  else {
	    index_type a = rp_dg.node_label(seq[k])[0];
	    assert(a < first_rce_action);
	    std::cerr << " " << irce.actions[a].name << std::endl;
	    for (index_type i = 0; i < compiled_ce[a].size(); i++) {
	      index_type oa = original[a];
	      assert(oa < ins.n_actions());
	      index_type e = compiled_ce[a][i];
	      assert(e < cadd.size());
	      assert(cadd_action[e] == oa);
	      assert(!cadd_rem[e]);
	      std::cerr << " + ";
	      ins.write_atom_set(std::cerr, cadd[e].antecedent);
	      std::cerr << " -> " << ins.atoms[cadd[e].consequent].name
			<< " (compiled)" << std::endl;
	    }
	    for (index_type i = 1; i < rp_dg.node_label(seq[k]).size(); i++) {
	      index_type a = rp_dg.node_label(seq[k])[i];
	      assert(a >= first_rce_action);
	      index_type e = a - first_rce_action;
	      assert(e < rce_index.size());
	      assert(cadd_rem[rce_index[e]]);
	      std::cerr << " + ";
	      irce.write_atom_set(std::cerr, irce.actions[a].pre);
	      std::cerr << " -> ";
	      irce.write_atom_set(std::cerr, irce.actions[a].add);
	      std::cerr << std::endl;
	    }
	  }
	}
      }
      // extracting the relaxed plan into rp/rp_ce
      for (index_type k = 0; k < seq.length(); k++)
	if (!rp_dg.node_label(seq[k]).empty()) {
	  index_type a = rp_dg.node_label(seq[k])[0];
	  assert(a < first_rce_action);
	  index_type oa = original[a];
	  assert(oa < ins.n_actions());
	  index_set oa_ces;
	  // std::cerr << k << "." << seq[k] << " = "
	  // 	    << rp_dg.node_label(seq[k])
	  // 	    << ", oa = " << oa << ", cce = " << compiled_ce[a]
	  // 	    << std::endl;
	  for (index_type i = 0; i < compiled_ce[a].size(); i++) {
	    index_type e = compiled_ce[a][i];
	    assert(e < cadd.size());
	    assert(cadd_action[e] == oa);
	    assert(!cadd_rem[e]);
	    index_type oe = cadd_index[e];
	    assert(oe < ins.actions[oa].cadd.size());
	    oa_ces.insert(oe);
	  }
	  for (index_type i = 1; i < rp_dg.node_label(seq[k]).size(); i++) {
	    index_type a = rp_dg.node_label(seq[k])[i];
	    assert(a >= first_rce_action);
	    index_type e = a - first_rce_action;
	    assert(e < rce_index.size());
	    assert(cadd_rem[rce_index[e]]);
	    assert(cadd_action[rce_index[e]] == oa);
	    index_type oe = cadd_index[rce_index[e]];
	    assert(oe < ins.actions[oa].cadd.size());
	    oa_ces.insert(oe);
	  }
	  rp.append(oa);
	  rp_ce.append(oa_ces);
	}
    }
    // else, compile out failed ce's
    else {
      assert(!failed.empty());
      if (verbose_level > 0)
	std::cerr << "scheduling failed! " << failed.size()
		  << " unscheduled effect(s)" << std::endl;
      index_vec failed_act;
      index_set_vec failed_ce;
      for (index_type k = 0; k < failed.size(); k++) {
	assert(failed[k] >= first_rce_action);
	assert((failed[k] - first_rce_action) < rce_index.size());
	index_type e = rce_index[failed[k] - first_rce_action];
	assert(e < cadd.size());
	assert(cadd_rem[e]);
	if (verbose_level > 1) {
	  std::cerr << " " << e << ".";
	  ins.print_rule(std::cerr, cadd[e]);
	  std::cerr << " of " << cadd_action[e] << "."
		    << irce.actions[cadd_action[e]].name
		    << std::endl;
	}
	index_type i = failed_act.first(cadd_action[e]);
	if (i != no_such_index) {
	  failed_ce[i].insert(e);
	}
	else {
	  failed_act.append(cadd_action[e]);
	  index_set s;
	  s.assign_singleton(e);
	  failed_ce.append(s);
	  assert(failed_act.size() == failed_ce.size());
	}
      }
      // remove rce actions
      irce.actions.set_length(first_rce_action);
      // deal with actions one at a time
      for (index_type k = 0; k < failed_act.size(); k++) {
	assert(!failed_ce[k].empty());
	// strategy: compile out one ce for each action that has failed ce's
	if (failed_ce[k].size() > 1) {
	  if (verbose_level > 0)
	    std::cerr << "limiting compilation to 1 effect" << std::endl;
	  failed_ce[k].set_length(1);
	}
	// enumerate subsets of (to-be-compiled) ce's for this action
	// and create a copy for each non-empty subset:
	SubsetEnumerator se(failed_ce[k].size());
	bool more = se.first();
	while (more) {
	  if (se.current_set_size() > 0) {
	    index_set sel;
	    se.current_set(failed_ce[k], sel);
	    index_set ce_pre;
	    index_set ce_add;
	    for (index_type j = 0; j < sel.size(); j++) {
	      ce_pre.insert(cadd[sel[j]].antecedent);
	      ce_add.insert(cadd[sel[j]].consequent);
	    }
	    // for each action that is a compiled version of the action
	    // with the failed ce's, create a copy with these ce's:
	    index_type n = irce.n_actions();
	    for (index_type i = 0; i < n; i++)
	      if (original[i] == failed_act[k]) {
		if (irce.actions[i].pre.contains(ce_pre)) {
		  irce.actions[i].add.insert(ce_add);
		  compiled_ce[i].insert(sel);
		  if (verbose_level > 3) {
		    std::cerr << "added ce's " << sel << " into ";
		    irce.print_action(std::cerr, irce.actions[i]);
		  }
		}
		else {
		  Instance::Action& a = irce.copy_action(i);
		  a.pre.insert(ce_pre);
		  a.add.insert(ce_add);
		  if (verbose_level > 3) {
		    std::cerr << "created ";
		    irce.print_action(std::cerr, a);
		    std::cerr << "by compiling ce's " << sel << " into ";
		    irce.print_action(std::cerr, irce.actions[i]);
		  }
		  original.append(failed_act[k]); // == original[i]
		  assert(original.size() == irce.n_actions());
		  compiled_ce.append(compiled_ce[i]);
		  assert(compiled_ce.size() == irce.n_actions());
		  compiled_ce[compiled_ce.size() - 1].insert(sel);
		}
		if (stats.break_signal_raised()) {
		  stats.stop();
		  double size_increase =
		    irce.n_actions()/(double)ins.n_actions();
		  cce_actions_ratio += size_increase;
		  return v_max;
		}
	      }
	  }
	  more = se.next();
	}
	// remove compiled ce's:
	for (index_type j = 0; j < failed_ce[k].size(); j++)
	  cadd_rem[failed_ce[k][j]] = false;
	n_compiled += failed_ce[k].size();
	n_ce_compiled += failed_ce[k].size();
	act_ce_compiled[failed_act[k]] += failed_ce[k].size();
	if (act_ce_compiled[failed_act[k]] > max_ce_compiled)
	  max_ce_compiled = act_ce_compiled[failed_act[k]];
      }
    }
  }
  stats.stop();
#ifdef HPLUS_PRINT_STATS
  std::cerr << "h+/ce calls: " << calls_to_hplus_ce << std::endl;
  std::cerr << "h+/ce: " << iterations << " iterations ("
	    << calls_to_hplus << " total, avg "
	    << PRINT_NTYPE(R_TO_N(calls_to_hplus, calls_to_hplus_ce))
	    << " its/call)" << std::endl;
  std::cerr << "h+/ce: time (scheduling rp) = "
	    << sce_stats.total_time() << std::endl;
  double size_increase = first_rce_action/(double)ins.n_actions();
  cce_actions_ratio += size_increase;
  std::cerr << "h+/ce: " << n_compiled
	    << " ce's compiled (total), " << max_ce_compiled
	    << " (max), size increase = " << size_increase << std::endl;
  print_hplus_stats(std::cerr);
#endif
  return v_max;
}

NTYPE ILB::hplus_with_ce
(const index_set& init, const index_set& goal, NTYPE bound, Plan** sol)
{
  index_vec rp;
  index_set_vec rp_ce;
  NTYPE val = hplus_with_ce(init, goal, bound, rp, rp_ce);
  if (FINITE(val) && !stats.break_signal_raised() && (sol != NULL)) {
    ActionSequence* plan = new ActionSequence(rp);
    *sol = plan;
  }
  return val;
}

NTYPE ILB::hplus_with_ce_time_indexed
(const index_set& init, const index_set& goal, NTYPE bound,
 index_vec& rp, index_set_vec& rp_ce)
{
  calls_to_hplus_ce += 1;
  stats.start();
  // if (prune_relaxed_irrelevant)
  //   compute_relevant(init, goal);
  // if (prune_relaxed_dominated)
  //   remove_dominated_actions(init);

  // construct time-indexed problem
  if (verbose_level > 0)
    std::cerr << "constructing time-indexed instance..." << std::endl;
  index_type tmax = (ins.n_atoms() - init.size());

  Instance i2;
  index_vec ce_act_map(no_such_index, ins.n_actions());
  index_type n_ce_act = 0;
  pair_vec dc_act;
  index_set init2;
  index_set goal2;

  for (index_type k = 0; k < ins.n_actions(); k++)
    if (relevant_actions[k] && !ins.actions[k].cadd.empty())
      ce_act_map[k] = n_ce_act++;

  #define atom_at_t(p, t) (((ins.n_atoms() + n_ce_act) * t) + p)
  #define name_at_t(n, t) n
  //#define name_at_t(n, t) (new NameAtIndex(n, t))

  // make atoms
  for (index_type t = 0; t <= tmax; t++) {
    for (index_type k = 0; k < ins.n_atoms(); k++)
      i2.new_atom(name_at_t(ins.atoms[k].name, t));
    for (index_type i = 0; i < init.size(); i++) {
      i2.atoms[atom_at_t(init[i], t)].init = true;
      init2.insert(atom_at_t(init[i], t));
    }
    if (t < tmax)
      for (index_type k = 0; k < ins.n_actions(); k++)
	if (ce_act_map[k] != no_such_index)
	  i2.new_atom(name_at_t(ins.actions[k].name, t));
  }
  for (index_type i = 0; i < goal.size(); i++) {
    i2.atoms[atom_at_t(goal[i], tmax)].goal = true;
    goal2.insert(atom_at_t(goal[i], tmax));
  }
  for (index_type t = 0; t < tmax; t++) {
    // make standard actions
    for (index_type k = 0; k < ins.n_actions(); k++)
      if (relevant_actions[k]) {
	Instance::Action& a = i2.new_action(name_at_t(ins.actions[k].name, t));
	for (index_type i = 0; i < ins.actions[k].pre.size(); i++)
	  a.pre.insert(atom_at_t(ins.actions[k].pre[i], t));
	for (index_type i = 0; i < ins.actions[k].add.size(); i++)
	  for (index_type t2 = t + 1; t2 <= tmax; t2++)
	    a.add.insert(atom_at_t(ins.actions[k].add[i], t2));
	if (ce_act_map[k] != no_such_index)
	  a.add.insert(atom_at_t(ins.n_atoms() + ce_act_map[k], t));
	a.cost = cost(k);
      }
    // make de-conditionalising actions
    for (index_type k = 0; k < ins.n_actions(); k++)
      if (ce_act_map[k] != no_such_index)
	for (index_type j = 0; j < ins.actions[k].cadd.size(); j++) {
	  Instance::Action& a =
	    i2.new_action(name_at_t(ins.actions[k].name, t));
	  const index_set& c = ins.actions[k].cadd[j].antecedent;
	  for (index_type i = 0; i < c.size(); i++)
	    a.pre.insert(atom_at_t(c[i], t));
	  a.pre.insert(atom_at_t(ins.n_atoms() + ce_act_map[k], t));
	  for (index_type t2 = t + 1; t2 <= tmax; t2++)
	    a.add.insert(atom_at_t(ins.actions[k].cadd[j].consequent, t2));
	  a.cost = ZERO;
	  if (t == 0) dc_act.append(index_pair(k, j));
	}
  }
  if (verbose_level > 0)
    std::cerr << "cross-referencing..." << std::endl;
  i2.cross_reference();
  // i2.print(std::cerr);

  // bool_vec s0(init2, i2.n_atoms());
  // UnitACF c0;
  // CostTable h0(i2, h1_stats);
  // h0.compute_H1(c0, s0);
  // index_set_graph g;
  // h0.compute_hm_graph(c0, goal2, true, g);
  // ((HSPS::labeled_graph<HSPS::index_set, HSPS::index_set>&)g).
  //   write_digraph(std::cerr, false, true, true, false, "h-m Graph");

  // compute h+ on time-indexed instance:
  std::cerr << "computing relaxed plan..." << std::endl;
  CostACF c2(i2);
  ILB hp2(i2, c2, NULL, stats);
  hp2.copy_options(*this, 0);
  hp2.is_sub = true;
  NTYPE v_max = hp2.hplus(init2, goal2, bound, NULL);
  add_hplus_stats(hp2);
#ifdef HPLUS_PRINT_STATS
  print_hplus_stats(std::cerr);
#endif
  if (stats.break_signal_raised() || (v_max >= bound) || INFINITE(v_max)) {
    stats.stop();
    return v_max;
  }

  // print the plan
  if (verbose_level > 1) {
    for (index_type k = 0; k < hp2.rp_actions.size(); k++) {
      index_type ca = hp2.rp_actions[k];
      index_type t = ca / (ins.n_actions() + dc_act.size());
      assert(t < tmax);
      index_type act_at_t = ca - (t * (ins.n_actions() + dc_act.size()));
      assert(act_at_t < ins.n_actions() + dc_act.size());
      if (act_at_t < ins.n_actions()) {
	std::cerr << t << " : " << act_at_t << "."
		  << ins.actions[act_at_t].name << std::endl;
      }
      else {
	index_type dc_at_t = act_at_t - ins.n_actions();
	assert(dc_at_t < dc_act.size());
	index_type ce_act = dc_act[dc_at_t].first;
	index_type ce_num = dc_act[dc_at_t].second;
	assert(ce_num < ins.actions[ce_act].cadd.size());
	std::cerr << t << " : ce " << ce_num << " ";
	ins.print_rule(std::cerr, ins.actions[ce_act].cadd[ce_num]);
	std::cerr << " of " << ce_act << "." << ins.actions[ce_act].name
		  << std::endl;
      }
    }
  }
  stats.stop();
  return v_max;
}

///
// hplusplus
///

void ILB::compute_non_redundant_rp
(const index_set& init, const index_set& goal, bool is_opt)
{
  rpn_stats.start();
  index_set non_red_rp;
  bool_vec red_acts(false, ins.n_actions());
  for (index_type k = 0; k < rp_actions.size(); k++) {
    if (!is_opt || IS_ZERO(cost(rp_actions[k]))) {
      bool_vec in(rp_actions, ins.n_actions());
      in.subtract(red_acts);
      in[rp_actions[k]] = false;
      reach.recompute(init, in);
      if (reach.unreachable(goal))
	non_red_rp.insert(k);
      else
	red_acts[rp_actions[k]] = true;
    }
    else {
      non_red_rp.insert(k);
    }
    if (stats.break_signal_raised()) return;
  }
  if (non_red_rp.size() < rp_actions.size()) {
    if (verbose_level > 0)
      std::cerr << "note: only " << non_red_rp.size()
		<< " of " << rp_actions.size()
		<< " actions in rp are necessary"
		<< std::endl;
    rp_actions = index_set(rp_actions, non_red_rp);
    reach.recompute(init, rp_actions);
    assert(!reach.unreachable(goal));
    // index_type rp_size = rp_actions.size();
    // compute_non_redundant_rp(init, goal, is_opt);
    // assert(rp_actions.size() == rp_size);
  }
  rpn_stats.stop();
}

void ILB::compute_dependency_graph
(const index_set& init, const index_set& goal, index_set_graph& dg)
{
  reach.recompute(init);
  assert(!reach.unreachable(goal));
  dg.init(rp_actions.size() + 1);
  for (index_type k = 0; k < rp_actions.size(); k++) {
    dg.node_label(k).assign_singleton(rp_actions[k]);
    bool_vec in(rp_actions, ins.n_actions());
    in[rp_actions[k]] = false;
    reach.recompute(init, in);
    index_set l;
    reach.unreachable_subset(goal, l);
    if (l.empty()) {
      std::cerr << "redundant action:" << std::endl;
      ins.print_action(std::cerr, ins.actions[rp_actions[k]]);
    }
    assert(!l.empty());
    dg.add_edge(k, rp_actions.size(), l);
    for (index_type i = 0; i < rp_actions.size(); i++)
      if (i != k) {
	reach.unreachable_subset(ins.actions[rp_actions[i]].pre, l);
	if (!l.empty()) {
	  dg.add_edge(k, i, l);
	}
      }
  }
  bool is_acyc = dg.acyclic();
  if (!is_acyc) {
    dg.write_digraph(std::cerr, "RPDG");
  }
  assert(is_acyc);
  /////
  /// note: compute_dependency_graph should NOT perform transitive
  /// reduction on the graph! this must be done by the caller on
  /// return, if desired.
  /// dg.transitive_reduction();
  ////
  // debug checking and trace printing:
  // index_set_vec prod(EMPTYSET, ins.n_atoms());
  // for (index_type k = 0; k < rp_actions.size(); k++) {
  //   std::cerr << k << ". " << ins.actions[rp_actions[k]].name << std::endl;
  //   for (index_type i = 0; i < rp_dg.successors(k).size(); i++) {
  //     index_type n = rp_dg.successors(k)[i];
  //     std::cerr << " -> " << n << ". ";
  //     if (n < rp_actions.size())
  // 	std::cerr << ins.actions[rp_actions[n]].name << std::endl;
  //     else
  // 	std::cerr << "goal" << std::endl;
  //     assert(rp_dg.edge_has_label(k, n));
  //     assert(!rp_dg.edge_label(k, n).empty());
  //     for (index_type j = 0; j < rp_dg.edge_label(k, n).size(); j++) {
  // 	index_type p = rp_dg.edge_label(k, n)[j];
  // 	std::cerr << "    + " << ins.atoms[p].name
  // 		  << std::endl;
  // 	assert(ins.actions[rp_actions[k]].add.contains(p));
  // 	prod[p].insert(k);
  //     }
  //   }
  // }
  // for (index_type i = 0; i < ins.n_atoms(); i++) {
  //   std::cerr << "producers of " << ins.atoms[i].name
  // 	      << ": " << prod[i] << std::endl;
  //   assert(prod[i].size() <= 1);
  // }
}

void ILB::partition_applicable
(const index_vec& rp,
 const bool_vec& rem_rp,
 const bool_vec& rstate,
 index_set_vec& apps)
{
  for (index_type k = 0; k < rp.size(); k++)
    if (rem_rp[k] && rstate.contains(ins.actions[rp[k]].pre)) {
      bool placed = false;
      for (index_type i = 0; (i < apps.size()) && !placed; i++) {
	bool can_go = true;
	for (index_type j = 0; (j < apps[i].size()) && can_go; j++)
	  if (!ins.commutative(rp[apps[i][j]], rp[k]))
	    can_go = false;
	if (can_go) {
	  apps[i].insert(k);
	  placed = true;
	}
      }
      if (!placed) {
	apps.append(EMPTYSET);
	apps[apps.length() - 1].insert(k);
      }
    }
}

void ILB::conflicts2
(index_type n_fail,
 const index_set& c_fail,
 lvector<bool_vec>& state,
 lvector<bool_vec>& rstate,
 const index_vec& rp,
 index_type s_fail,
 const index_set_vec& seq,
 bool check_ce,
 node_conflict_set& cs)
{
  // const index_set& c_fail =
  // (n_fail < rp_actions.size() ? ins.actions[rp_actions[n_fail]].pre : goal);
  if (verbose_level > 2) {
    std::cerr << "extracting conflicts: " << n_fail << ", c = ";
    ins.write_atom_set(std::cerr, c_fail);
    std::cerr << ", n_fail at step " << s_fail << std::endl;
  }
  assert(s_fail > 0);
  // for every atom (p) in the failed condition that is not true in
  // the current state (at step s_fail)...
  for (index_type i = 0; i < c_fail.size(); i++) {
    index_type p = c_fail[i];
    if (!state[s_fail][p]) {
      if (verbose_level > 2) {
	std::cerr << "p = " << p << ": " << ins.atoms[p].name << std::endl;
      }
      // find the last step where the atom was true (s_true)
      index_type s_true = s_fail;
      while ((s_true > 0) && !state[s_true - 1][p])
	s_true -= 1;
      assert(s_true > 0);
      s_true -= 1;
      assert(s_true < s_fail);
      assert(state[s_true][p]);
      if (verbose_level > 2) {
	std::cerr << "s_true = " << s_true << std::endl;
      }
      // rem stores "future" relaxed plan actions
      bool_vec rem(true, rp.size());
      for (index_type s = 0; s < s_true; s++)
	for (index_type k = 0; k < seq[s].size(); k++)
	  rem[seq[s][k]] = false;
      // for every action taking place between s_true and current step...
      for (index_type s = s_true; s < s_fail; s++) {
	if (verbose_level > 3) {
	  std::cerr << "step " << s << ":";
	  for (index_type k = 0; k < seq[s].size(); k++)
	    std::cerr << " " << ins.actions[rp[seq[s][k]]].name;
	  std::cerr << std::endl;
	}
	for (index_type k = 0; k < seq[s].size(); k++)
	  rem[seq[s][k]] = false;
	for (index_type k = 0; k < seq[s].size(); k++) {
	  index_type a = rp[seq[s][k]];
	  // if this action deletes p...
	  if (ins.actions[a].del.contains(p)) {
	    if (verbose_level > 2) {
	      std::cerr << seq[s][k] << " = " << ins.actions[a].name
			<< " deletes p at " << s << std::endl;
	      std::cerr << "generated node conflict: ("
			<< index_pair(seq[s][k], n_fail) << "," << p << ")"
			<< std::endl;
	    }
	    cs.insert(node_conflict(index_pair(seq[s][k], n_fail), p));
	  }
	  // else, if the action conditionally deletes p...
	  else if (check_ce) {
	    for (index_type i = 0; i < ins.actions[a].cdel.size(); i++)
	      if (ins.actions[a].cdel[i].consequent == p)
		if (state[s].contains(ins.actions[a].cdel[i].antecedent)) {
		  if (verbose_level > 2) {
		    std::cerr << seq[s][k] << " = " << ins.actions[a].name
			      << " conditionally deletes p at " << s
			      << std::endl;
		    std::cerr << "generated node conflict: ("
			      << index_pair(seq[s][k], n_fail)
			      << "," << p << ")" << std::endl;
		  }
		  cs.insert(node_conflict(index_pair(seq[s][k], n_fail), p));
		}
	  }
	  // end of if/else
	}
      }
    }
  }
  assert(!cs.empty());
}

bool ILB::exec_rp2
(bool_vec& rem_rp,
 lvector<bool_vec>& state,
 lvector<bool_vec>& rstate,
 const index_set& goal,
 index_type step,
 index_set_vec& seq,
 lvector<node_conflict_set>& cs)
{
  //std::cerr << " -> " << seq << std::endl;
  if (rem_rp.count(true) == 0) {
    assert(rstate[step].contains(goal));
    if (state[step].contains(goal)) {
      return true;
    }
    else {
      node_conflict_set goal_cs;
      conflicts2(rp_actions.size(), goal, state, rstate, rp_actions, step,
		 seq, false, goal_cs);
      assert(cs.size() > 0);
      cs[cs.size() - 1].insert(goal_cs);
      return false;
    }
  }
  index_set_vec pnext;
  partition_applicable(rp_actions, rem_rp, rstate[step], pnext);
  assert(pnext.length() > 0);
  assert(cs.size() > 0);
  node_conflict_set cur_cs(cs[cs.size() - 1]);
  for (index_type k = 0; k < pnext.length(); k++) {
    state[step + 1] = state[step];
    rstate[step + 1] = rstate[step];
    bool ok_now = true;
    node_conflict_set step_cs;
    for (index_type i = 0; i < pnext[k].size(); i++) {
      if (!state[step].contains(ins.actions[rp_actions[pnext[k][i]]].pre)) {
	node_conflict_set tmp_cs;
	conflicts2(pnext[k][i], ins.actions[rp_actions[pnext[k][i]]].pre,
		   state, rstate, rp_actions, step, seq, false, tmp_cs);
	step_cs.insert(tmp_cs);
	ok_now = false;
      }
      state[step + 1].subtract(ins.actions[rp_actions[pnext[k][i]]].del);
      state[step + 1].insert(ins.actions[rp_actions[pnext[k][i]]].add);
      rstate[step + 1].insert(ins.actions[rp_actions[pnext[k][i]]].add);
      rem_rp[pnext[k][i]] = false;
    }
    //std::cerr << step << "/" << k << ": " << ok_now << std::endl;
    seq[step] = pnext[k];
    cs[cs.size() - 1].insert(step_cs);
    bool ok_rest =
      exec_rp2(rem_rp, state, rstate, goal, step + 1, seq, cs);
    if (cs[cs.size() - 1].empty()) return true;
    if (stats.break_signal_raised()) return false;
    for (index_type i = 0; i < pnext[k].size(); i++)
      rem_rp[pnext[k][i]] = true;
    if (k + 1 < pnext.length())
      cs.append(cur_cs);
    // HACK! HACK! HACK!
    //assert(cs.size() < 10000);
  }
  return false;
}

index_type ILB::estimate_conflict_weight
(const index_set_graph& dg, const index_pair& c, index_type& ncd)
{
  index_type deler = c.first;
  index_type failed = c.second;
  index_vec d_deler;
  dg.distance(deler, d_deler);
  // case 1:
  if (d_deler[failed] != no_such_index) {
    assert(dg.reachable(deler, failed));
    ncd = no_such_index;
    return d_deler[failed];
  }
  // case 2:
  bool_vec ncds;
  dg.nearest_common_descendants(deler, failed, ncds);
  assert(ncds.first(true) != no_such_index);
  // remove ncd's where paths from deler and failed have a common atom:
  for (index_type i = 0; i < ncds.size(); i++)
    if (ncds[i]) {
      index_set_vec ps;
      atom_sets_on_path(dg, deler, i, ps);
      index_set u1;
      ps.union_set(u1);
      ps.clear();
      atom_sets_on_path(dg, failed, i, ps);
      index_set u2;
      ps.union_set(u2);
      if (u1.have_common_element(u2))
	ncds[i] = false;
    }
  // now find cheapest ncd among remaining
  index_vec d_failed;
  dg.distance(failed, d_failed);
  index_type i_min = ncds.first(true);
  assert(i_min != no_such_index);
  index_type d_min = d_deler[i_min] * d_failed[i_min];
  for (index_type i = 0; i < ncds.size(); i++)
    if (ncds[i]) {
      assert(d_deler[i] != no_such_index);
      assert(d_failed[i] != no_such_index);
      index_type d = d_deler[i] * d_failed[i];
      if (d < d_min) {
	i_min = i;
	d_min = d;
      }
    }
  ncd = i_min;
  return d_min;
}

void ILB::choose_atoms_on_path
(const index_set_graph& dg,
 const index_type n_from,
 const index_type n_to,
 index_set& as)
{
  index_vec p;
  index_type d = dg.shortest_path(n_from, n_to, p);
  assert(d != no_such_index);
  for (index_type i = 0; (i + 1) < p.length(); i++) {
    const index_set& e = dg.edge_label(p[i], p[i + 1]);
    assert(!e.empty());
    // here we could try to do something clever... but for now,
    // let's just pick one atom from the set:
    as.insert(e[0]);
  }
}

void ILB::atom_sets_on_path
(const index_set_graph& dg,
 const index_type n_from,
 const index_type n_to,
 index_set_vec& ps)
{
  index_vec p;
  index_type d = dg.shortest_path(n_from, n_to, p);
  assert(d != no_such_index);
  for (index_type i = 0; (i + 1) < p.length(); i++) {
    const index_set& e = dg.edge_label(p[i], p[i + 1]);
    assert(!e.empty());
    ps.append(e);
  }
}

// a "smarter" variant of node_conflict_to_atom_conflicts: here, we try
// to minimise the size of the conflict set, by always checking if
// there is a choice that will result in an atom set that is already
// in it.
void ILB::node_conflict_to_atom_conflicts
(const index_set_graph& dg,
 const node_conflict& nc,
 index_set_vec& acs)
{
  index_type deler = nc.first.first;
  index_type failed = nc.first.second;
  index_type p = nc.second;
  if (verbose_level > 2) {
    std::cerr << "deler = " << deler << ", failed = " << failed
	      << ", p = " << p << std::endl;
  }
  if (dg.reachable(deler, failed)) {
    index_set_vec ps;
    atom_sets_on_path(dg, deler, failed, ps);
    if (verbose_level > 2) {
      std::cerr << "case 1: ps = " << ps << std::endl;
    }
    for (index_type j = 0; j < ps.size(); j++) {
      bool found = false;
      for (index_type k = 0; (k < ps[j].size()) && !found; k++)
	if (acs.contains(index_set(p, ps[j][k])))
	  found = true;
      if (!found)
	acs.append(index_set(p, ps[j][0]));
    }
  }
  else {
    index_type ncd;
    estimate_conflict_weight(dg, index_pair(deler, failed), ncd);
    assert((ncd != no_such_index) && (ncd < dg.size()));
    // std::cerr << "deler = " << deler << ", failed = " << failed
    //	<< ", ncd = " << ncd << std::endl;
    index_set_vec ps1;
    atom_sets_on_path(dg, deler, ncd, ps1);
    index_set_vec ps2;
    atom_sets_on_path(dg, failed, ncd, ps2);
    if (verbose_level > 2) {
      std::cerr << "case 2: ncd = " << ncd
		<< std::endl << "ps1 = " << ps1
		<< std::endl << "ps2 = " << ps2
		<< std::endl;
    }
    ps2.append(index_set(p));
    for (index_type j1 = 0; j1 < ps1.size(); j1++)
      for (index_type j2 = 0; j2 < ps2.size(); j2++) {
	bool found = false;
	for (index_type k1 = 0; (k1 < ps1[j1].size()) && !found; k1++)
	  for (index_type k2 = 0; (k2 < ps2[j2].size()) && !found; k2++) {
	    if (ps1[j1][k1] == ps2[j2][k2]) {
	      std::cerr << "deler = " << deler
			<< ", failed = " << failed
			<< ", ncd = " << ncd
			<< std::endl;
	      std::cerr << "ps1 = " << ps1 << std::endl;
	      std::cerr << "ps2 = " << ps2 << std::endl;
	    }
	    assert(ps1[j1][k1] != ps2[j2][k2]);
	    if (acs.contains(index_set(ps1[j1][k1], ps2[j2][k2])))
	      found = true;
	  }
	if (!found)
	  acs.append(index_set(ps1[j1][0], ps2[j2][0]));
      }
  }
}

void ILB::compute_dependency_closure
(const index_set& rp,
 const index_set_graph& dg,
 index_type n_from,
 index_type n_to,
 index_set& dc_labels)
{
  bool_vec dc(false, dg.size());
  dc[n_from] = true;
  index_vec p;
  index_type d = dg.shortest_path(n_from, n_to, p);
  assert(d != no_such_index);
  for (index_type i = 0; (i + 1) < p.length(); i++) {
    dc[p[i + 1]] = true;
    const index_set& e = dg.edge_label(p[i], p[i + 1]);
    assert(!e.empty());
    // just pick an arbitrary atom from the label
    dc_labels.insert(e[0]);
  }
  bool done = false;
  while (!done) {
    done = true;
    for (index_type i = 0; i < dc_labels.size(); i++)
      for (index_type k = 0; k < rp.size(); k++)
	if (ins.actions[rp[k]].add.contains(dc_labels[i]) && !dc[k]) {
	  p.clear();
	  d = dg.shortest_path(n_from, k, p);
	  assert(d != no_such_index);
	  for (index_type j = 0; (j + 1) < p.length(); j++) {
	    dc[p[j + 1]] = true;
	    const index_set& e = dg.edge_label(p[j], p[j + 1]);
	    assert(!e.empty());
	    dc_labels.insert(e[0]);
	  }
	  done = false;
	}
  }
}

void ILB::node_conflict_to_atom_conflicts_with_dc
(const index_set& rp,
 const index_set_graph& dg,
 const node_conflict& nc,
 index_set_vec& acs)
{
  index_type deler = nc.first.first;
  index_type failed = nc.first.second;
  index_type p = nc.second;
  if (verbose_level > 2) {
    std::cerr << "(full dc) deler = " << deler << ", failed = " << failed
	      << ", p = " << p << std::endl;
  }
  if (dg.reachable(deler, failed)) {
    index_set ps;
    compute_dependency_closure(rp, dg, deler, failed, ps);
    if (verbose_level > 2) {
      std::cerr << "(full dc) case 1: ps = " << ps << std::endl;
    }
    for (index_type j = 0; j < ps.size(); j++)
      acs.append(index_set(p, ps[j]));
  }
  else {
    index_type ncd;
    estimate_conflict_weight(dg, index_pair(deler, failed), ncd);
    assert((ncd != no_such_index) && (ncd < dg.size()));
    // std::cerr << "deler = " << deler << ", failed = " << failed
    //	<< ", ncd = " << ncd << std::endl;
    index_set ps1;
    compute_dependency_closure(rp, dg, deler, ncd, ps1);
    index_set ps2;
    compute_dependency_closure(rp, dg, failed, ncd, ps2);
    if (verbose_level > 2) {
      std::cerr << "(full dc) case 2: ncd = " << ncd
		<< std::endl << "ps1 = " << ps1
		<< std::endl << "ps2 = " << ps2
		<< std::endl;
    }
    assert(!ps1.contains(p));
    ps2.insert(p);
    for (index_type j1 = 0; j1 < ps1.size(); j1++)
      for (index_type j2 = 0; j2 < ps2.size(); j2++) {
	assert(ps1[j1] != ps2[j2]);
	acs.append(index_set(ps1[j1], ps2[j2]));
      }
  }
}

void ILB::node_conflicts_to_atom_conflicts
(const index_set_graph& dg,
 const node_conflict_set& ncs,
 index_set_vec& acs)
{
  for (index_type k = 0; k < ncs.size(); k++) {
    node_conflict_to_atom_conflicts(dg, ncs[k], acs);
  }
}

bool ILB::check_rp2
(const index_set& init,
 const index_set& goal,
 index_set_vec& seq,
 index_set_vec& cs)
{
  if (verbose_level > 0)
    std::cerr << "reducing relaxed plan..." << std::endl;
#ifdef HOFFMANN
  compute_non_redundant_rp(init, goal, false);
#else
  compute_non_redundant_rp(init, goal, true);
#endif
  if (stats.break_signal_raised()) return false;
  if (verbose_level > 0)
    std::cerr << "computing RPDG..." << std::endl;
  compute_dependency_graph(init, goal, rp_dg);
  if (stats.break_signal_raised()) return false;
  rp_dg.transitive_reduction();
  if (verbose_level > 2) {
    std::cerr << "dependency graph: " << rp_dg << std::endl;
    print_rpdg(std::cerr, rp_actions, rp_dg);
  }
  if (verbose_level > 0)
    std::cerr << "checking validity & extracting conflicts..." << std::endl;
  bool_vec rem_rp(true, rp_actions.size());
  lvector<bool_vec>
    state(bool_vec(init, ins.n_atoms()), rp_actions.size() + 1);
  lvector<bool_vec>
    rstate(bool_vec(init, ins.n_atoms()), rp_actions.size() + 1);
  seq.assign_value(EMPTYSET, rp_actions.size());
  lvector<node_conflict_set> exec_cs(node_conflict_set(), 1);
  bool solved = exec_rp2(rem_rp, state, rstate, goal, 0, seq, exec_cs);
  if (verbose_level > 0) {
    std::cerr << "examined " << exec_cs.size() << " executions" << std::endl;
  }
  if (stats.break_signal_raised()) return false;
  if (!solved) {
    if (verbose_level > 1) {
      std::cerr << "conflict sets:" << std::endl;
      for (index_type k = 0; k < exec_cs.size(); k++) {
	std::cerr << k << ":";
	for (index_type i = 0; i < exec_cs[k].size(); i++) {
	  index_type ncd;
	  index_type w =
	    estimate_conflict_weight(rp_dg, exec_cs[k][i].first, ncd);
	  std::cerr << " (" << exec_cs[k][i].first << ","
		    << exec_cs[k][i].second << ")/" << w;
	}
	std::cerr << std::endl;
      }
    }
    // and now we have to solve a hitting set problem again...
    node_conflict_set sel_cs;
    bool_vec rem_to_hit(true, exec_cs.size());
    for (index_type k = 0; k < exec_cs.size(); k++)
      if (rem_to_hit[k]) {
	assert(!exec_cs[k].empty());
	index_type i_best = no_such_index;
	NTYPE w_best = POS_INF;
	for (index_type i = 0; i < exec_cs[k].size(); i++) {
	  index_type n = 1;
	  for (index_type j = k + 1; j < exec_cs.size(); j++)
	    if (rem_to_hit[j] && exec_cs[j].contains(exec_cs[k][i])) n += 1;
	  index_type ncd;
	  NTYPE w = R_TO_N(estimate_conflict_weight(rp_dg, exec_cs[k][i].first, ncd), n);
	  if (w < w_best) {
	    i_best = i;
	    w_best = w;
	  }
	}
	assert((i_best != no_such_index) && (i_best < exec_cs[k].size()));
	sel_cs.insert(exec_cs[k][i_best]);
	rem_to_hit[k] = false;
	for (index_type j = k + 1; j < exec_cs.size(); j++)
	  if (rem_to_hit[j] && exec_cs[j].contains(exec_cs[k][i_best]))
	    rem_to_hit[j] = false;
      }
    if (verbose_level > 1) {
      std::cerr << "selected conflict set:";
      for (index_type i = 0; i < sel_cs.size(); i++)
	std::cerr << " (" << sel_cs[i].first << "," << sel_cs[i].second << ")";
      std::cerr << std::endl;
    }
    assert(!sel_cs.empty());
    node_conflicts_to_atom_conflicts(rp_dg, sel_cs, cs);
    assert(!cs.empty());
  }
  return solved;
}

// note: this method is not used.
void ILB::meta_atom_to_atom_set
(index_type p,
 const rule_set& map,
 index_set& set)
{
  index_type i = map.find_rule(p);
  if (i != no_such_index) {
    assert(i < map.size());
    set.insert(map[i].antecedent);
  }
  else {
    set.insert(p);
  }
}

bool ILB::simple_exec_rp
(bool_vec& rem_rp,
 index_type step,
 lvector<bool_vec>& state,
 const index_set& goal,
 const index_vec& rp,
 index_set_vec& seq)
{
  if (rem_rp.count(true) == 0) {
    if (state[step].contains(goal))
      return true;
    else
      return false;
  }
  assert((step + 1) < state.size());
  assert(step < seq.size());
  index_set_vec pnext;
  partition_applicable(rp, rem_rp, state[step], pnext);
  // std::cerr << "step = " << step << ", seq = " << seq
  // 	    << ", pnext = " << pnext << std::endl;
  for (index_type k = 0; k < pnext.length(); k++) {
    state[step + 1] = state[step];
    for (index_type i = 0; i < pnext[k].size(); i++) {
      index_type a = rp[pnext[k][i]];
      state[step + 1].subtract(ins.actions[a].del);
      for (index_type j = 0; j < ins.actions[a].cdel.size(); j++)
	if (state[step].contains(ins.actions[a].cdel[j].antecedent))
	  state[step + 1][ins.actions[a].cdel[j].consequent] = false;
      state[step + 1].insert(ins.actions[a].add);
      for (index_type j = 0; j < ins.actions[a].cadd.size(); j++)
	if (state[step].contains(ins.actions[a].cadd[j].antecedent))
	  state[step + 1][ins.actions[a].cadd[j].consequent] = true;
      rem_rp[pnext[k][i]] = false;
    }
    seq[step] = pnext[k];
    bool ok = simple_exec_rp(rem_rp, step + 1, state, goal, rp, seq);
    if (ok) return true;
    for (index_type i = 0; i < pnext[k].size(); i++)
      rem_rp[pnext[k][i]] = true;
  }
  return false;
}

bool ILB::rp_is_plan
(const index_set& init,
 const index_set& goal,
 const index_vec& rp,
 index_set_vec& seq)
{
  bool_vec rem(true, rp.size());
  lvector<bool_vec> state(bool_vec(init, ins.n_atoms()), rp.size() + 1);
  seq.assign_value(EMPTYSET, rp.size());
  return simple_exec_rp(rem, 0, state, goal, rp, seq);
}

bool ILB::PNCWeightIncreasing::operator()
(const PotentialNodeConflict& v0, const PotentialNodeConflict& v1) const
{
  return (v0.weight < v1.weight);
}

void ILB::remove_self_pairs(index_set_vec& cs)
{
  index_type k = 0;
  while (k < cs.size()) {
    assert(cs[k].size() == 2);
    bool remove = false;
    index_type m0 = meta_atom_map.find_rule(cs[k][0]);
    index_type m1 = meta_atom_map.find_rule(cs[k][1]);
    if (m0 != no_such_index) {
      if (meta_atom_map[m0].antecedent.contains(cs[k][1])) {
	remove = true;
	std::cerr << "removing pair " << cs[k] << " with "
		  << cs[k][0] << " = " << meta_atom_map[m0].antecedent
		  << std::endl;
      }
    }
    if (m1 != no_such_index) {
      if (meta_atom_map[m1].antecedent.contains(cs[k][0])) {
	remove = true;
	std::cerr << "removing pair " << cs[k] << " with "
		  << cs[k][1] << " = " << meta_atom_map[m1].antecedent
		  << std::endl;
      }
    }
    if (remove)
      cs.remove(k);
    else
      k += 1;
  }
}

bool ILB::contains_new_meta_atom(const index_set_vec& cs) const
{
  for (index_type k = 0; k < cs.size(); k++) {
    assert(cs[k].size() >= 2);
    index_type m = meta_atom_map.find_rule(cs[k]);
    if (m == no_such_index)
      return true;
  }
  return false;
}

void ILB::conflicts3
(const index_set& init,
 const index_set& goal,
 index_set_vec& cs)
{
  graph prec(rp_dg);
  prec.transitive_closure();
  if (verbose_level > 1) {
    std::cerr << "initial precedence graph: " << prec << std::endl;
  }
  lvector<PotentialNodeConflict> pc;
  PNCWeightIncreasing o;
  index_type g_node = rp_actions.size();
  index_type ncd;
  // initialise sorted list of potential conflicts
  for (index_type i = 0; i < rp_actions.size(); i++) {
    for (index_type j = 0; j < rp_actions.size(); j++)
      if ((i != j) && !prec.adjacent(j, i)) {
	index_set dps(ins.actions[rp_actions[j]].pre);
	dps.intersect(ins.actions[rp_actions[i]].del);
	for (index_type l = 0; l < dps.size(); l++) {
	  bool nec_rescue = false;
	  PotentialNodeConflict nc(i, j, dps[l]);
	  for (index_type k = 0; (k < rp_actions.size()) && !nec_rescue; k++)
	    if ((k != i) && (k != j) &&
		!prec.adjacent(k, i) && !prec.adjacent(j, k))
	      if (ins.actions[rp_actions[k]].add.contains(dps[l])) {
		nc.rescue.insert(k);
		if (prec.adjacent(i, k) && prec.adjacent(k, j))
		  nec_rescue = true;
	      }
	  if (!nec_rescue) {
	    nc.weight = estimate_conflict_weight(rp_dg, nc.pair, ncd);
	    pc.insert_ordered(nc, o);
	  }
	}
      }
    index_set dps(goal);
    dps.intersect(ins.actions[rp_actions[i]].del);
    for (index_type l = 0; l < dps.size(); l++) {
      bool nec_rescue = false;
      PotentialNodeConflict nc(i, g_node, dps[l]);
      for (index_type k = 0; (k < rp_actions.size()) && !nec_rescue; k++)
	if ((k != i) && !prec.adjacent(k, i))
	  if (ins.actions[rp_actions[k]].add.contains(dps[l])) {
	    nc.rescue.insert(k);
	    if (prec.adjacent(i, k))
	      nec_rescue = true;
	  }
      if (!nec_rescue) {
	nc.weight = estimate_conflict_weight(rp_dg, nc.pair, ncd);
	pc.insert_ordered(nc, o);
      }
    }
  }
  bool_vec rem(true, pc.size());
  index_set  best;
  best.fill(pc.size());
  // initialise best solution to be the set of all potential conflicts
  index_type w_best = 0;
  for (index_type k = 0; k < pc.size(); k++)
    w_best += pc[k].weight;
  index_set current;
  index_type w_current = 0;
  // find a (hopefully) better solution
  choose_next_conflict(pc, rem, prec, current, w_current, best, w_best);
  // convert the chosen set of conflicts to atom conflict sets
  for (index_type k = 0; k < best.size(); k++) {
    node_conflict nc(pc[best[k]].pair, pc[best[k]].dp);
    node_conflict_to_atom_conflicts(rp_dg, nc, cs);
  }
  // remove self-pairs
  remove_self_pairs(cs);
  bool new_c = contains_new_meta_atom(cs);
  if ((verbose_level > 0) && !new_c) {
    std::cerr << "no new conflict found, retrying with full dc..." << std::endl;
  }
  // try again, with full-dc option on
  for (index_type k = 0; (k < best.size()) && !new_c; k++) {
    node_conflict nc(pc[best[k]].pair, pc[best[k]].dp);
    node_conflict_to_atom_conflicts_with_dc(rp_actions, rp_dg, nc, cs);
    remove_self_pairs(cs);
    new_c = contains_new_meta_atom(cs);
  }
}

void ILB::choose_next_conflict
(const lvector<PotentialNodeConflict>& pc,
 bool_vec& rem,
 graph& prec,
 index_set& current,
 index_type w_current,
 index_set& best,
 index_type w_best)
{
  // update rem to reflect changes in prec
  for (index_type k = 0; k < pc.size(); k++)
    if (rem[k]) {
      if (prec.adjacent(pc[k].pair.second, pc[k].pair.first))
	rem[k] = false;
      for (index_type i = 0; (i < pc[k].rescue.size()) && rem[k]; i++)
	if (prec.adjacent(pc[k].pair.first, pc[k].rescue[i]) &&
	    prec.adjacent(pc[k].rescue[i], pc[k].pair.second))
	  rem[k] = false;
    }
  // for each remaining conflict, check how many ways we can avoid it
  // at the same time, find the minimum and second minimum weight
  // among remaining conflicts
  index_vec na(0, pc.size());
  index_vec nr(0, pc.size());
  index_type w_min = w_best + 1;
  index_type w_2nd_min = w_best + 1;
  for (index_type k = 0; k < pc.size(); k++)
    if (rem[k]) {
      if (!prec.adjacent(pc[k].pair.first, pc[k].pair.second))
	na[k] += 1;
      for (index_type i = 0; i < pc[k].rescue.size(); i++)
	if (!prec.adjacent(pc[k].rescue[i], pc[k].pair.first) &&
	    !prec.adjacent(pc[k].pair.second, pc[k].rescue[i]))
	  nr[k] += 1;
      na[k] += nr[k];
      if (pc[k].weight < w_min) {
	w_2nd_min = w_min;
	w_min = pc[k].weight;
      }
      else if (pc[k].weight < w_2nd_min) {
	w_2nd_min = pc[k].weight;
      }
    }
  // check if we can get a better solution by adding one unavoidable
  // conflict to the current set
  for (index_type k = 0; k < pc.size(); k++)
    if (rem[k] && (na[k] == 0))
      if ((w_current + pc[k].weight) < w_best) {
	best = current;
	best.insert(k);
	w_best = w_current + pc[k].weight;
      }
  // check if we can get a better solution by adding two conflicts that
  // are together unavoidable to the current set
  for (index_type k1 = 0; k1 < pc.size(); k1++)
    if (rem[k1])
      for (index_type k2 = k1 + 1; k2 < pc.size(); k2++)
	if (rem[k2])
	  if ((pc[k1].pair.first == pc[k2].pair.second) &&
	      (pc[k1].pair.second == pc[k2].pair.first) &&
	      (nr[k1] == 0) && (nr[k2] == 0) &&
	      (na[k1] > 0) && (na[k2] > 0) &&
	      ((w_current + pc[k1].weight + pc[k2].weight) < w_best)) {
	    best = current;
	    best.insert(k1);
	    best.insert(k2);
	    w_best = w_current + pc[k1].weight + pc[k2].weight;
	  }
  // any solution not considered by the two options above must add
  // at least two more conflicts to the current set; thus, if the
  // weight of the best solution is no more than the current plus
  // the minimum and second minimum remaining, then we can't do better
  if (w_best <= (w_current + w_min + w_2nd_min))
    return;
  // now, we explore further options...
  index_type w_best_choice = index_type_max;
  index_type na_best_choice = index_type_max;
  index_type best_choice = no_such_index;
  for (index_type k = 0; k < pc.size(); k++)
    if (rem[k] && (na[k] > 0)) {
      // the estimated additional cost of a conflict is its weight
      // + the weights of the next N remaining conflicts, where N
      // is the number of ways it can be avoided.
      index_type w = pc[k].weight;
      index_type n = na[k];
      for (index_type k2 = k + 1; (k2 < pc.size()) && (n > 0); k2++)
	if (rem[k2]) {
	  w += pc[k2].weight;
	  n -= 1;
	}
      if ((w < w_best_choice) ||
	  ((w == w_best_choice) && (na[k] < na_best_choice))) {
	w_best_choice = w;
	na_best_choice = na[k];
	best_choice = k;
      }
    }
  if (best_choice != no_such_index) {
    current.insert(best_choice);
    w_current += pc[best_choice].weight;
    assert(na[best_choice] > 0);
    // if there is only one way to avoid the chosen conflict, we won't
    // have to backtrack
    if (na[best_choice] == 1) {
      // can we put the failed node before the deleter?
      if (!prec.adjacent(pc[best_choice].pair.first,
			 pc[best_choice].pair.second)) {
	pair_set e;
	prec.add_edge_to_transitive_closure
	  (pc[best_choice].pair.second, pc[best_choice].pair.first, e);
	choose_next_conflict(pc, rem, prec, current, w_current, best, w_best);
	return;
      }
      // can we place a rescuer between them?
      for (index_type i = 0; i < pc[best_choice].rescue.size(); i++)
	if (!prec.adjacent(pc[best_choice].rescue[i],
			   pc[best_choice].pair.first) &&
	    !prec.adjacent(pc[best_choice].pair.second,
			   pc[best_choice].rescue[i])) {
	  pair_set e;
	  prec.add_edge_to_transitive_closure
	    (pc[best_choice].pair.first, pc[best_choice].rescue[i], e);
	  e.clear();
	  prec.add_edge_to_transitive_closure
	    (pc[best_choice].rescue[i], pc[best_choice].pair.second, e);
	  choose_next_conflict(pc, rem, prec, current, w_current, best, w_best);
	  return;
	}
      assert(0);
    }
    else {
      // if there are multiple resolvers, we have to find a best solution
      // for each of them, and merge the results.
      index_set merged_best;
      // can we put the failed node before the deleter?
      if (!prec.adjacent(pc[best_choice].pair.first,
			 pc[best_choice].pair.second)) {
	graph new_prec(prec);
	index_set new_current(best);
	index_set new_best(best);
	bool_vec new_rem(rem);
	pair_set e;
	new_prec.add_edge_to_transitive_closure
	  (pc[best_choice].pair.second, pc[best_choice].pair.first, e);
	choose_next_conflict(pc, new_rem, new_prec, new_current, w_current,
			     new_best, w_best);
	// if there is at least one branch whose best solution is
	// the current best, the merged best cannot be better, so
	// we can stop.
	if (new_best == best)
	  return;
	merged_best.insert(new_best);
      }
      // can we place a rescuer between them?
      for (index_type i = 0; i < pc[best_choice].rescue.size(); i++)
	if (!prec.adjacent(pc[best_choice].rescue[i],
			   pc[best_choice].pair.first) &&
	    !prec.adjacent(pc[best_choice].pair.second,
			   pc[best_choice].rescue[i])) {
	  graph new_prec(prec);
	  index_set new_current(best);
	  index_set new_best(best);
	  bool_vec new_rem(rem);
	  pair_set e;
	  prec.add_edge_to_transitive_closure
	    (pc[best_choice].pair.first, pc[best_choice].rescue[i], e);
	  e.clear();
	  prec.add_edge_to_transitive_closure
	    (pc[best_choice].rescue[i], pc[best_choice].pair.second, e);
	  choose_next_conflict(pc, new_rem, new_prec, new_current, w_current,
			       new_best, w_best);
	  if (new_best == best)
	    return;
	  merged_best.insert(new_best);
	}
      index_type w_merged = 0;
      for (index_type k = 0; k < merged_best.size(); k++)
	w_merged += pc[merged_best[k]].weight;
      if (w_merged < w_best) {
	best = merged_best;
	w_best = w_merged;
      }
    }
  }
}

bool ILB::check_rp3
(const index_set& init,
 const index_set& goal,
 index_set_vec& seq,
 index_set_vec& cs)
{
  if (verbose_level > 0)
    std::cerr << "reducing relaxed plan..." << std::endl;
#ifdef HOFFMANN
  compute_non_redundant_rp(init, goal, false);
#else
  compute_non_redundant_rp(init, goal, true);
#endif
  if (stats.break_signal_raised()) return false;
  if (verbose_level > 2) {
    std::cerr << "non-redundant rp:" << std::endl;
    print_rp(std::cerr, rp_actions);
    if (verbose_level > 3) {
      validate_rp(rp_actions);
    }
  }
  if (verbose_level > 0)
    std::cerr << "checking real validity..." << std::endl;
  if (rp_is_plan(init, goal, rp_actions, seq)) {
    if (stats.break_signal_raised()) return false;
    std::cerr << "rp is valid!" << std::endl;
    return true;
  }
  if (stats.break_signal_raised()) return false;
  if (verbose_level > 0)
    std::cerr << "computing RPDG..." << std::endl;
  compute_dependency_graph(init, goal, rp_dg);
  if (stats.break_signal_raised()) return false;
  rp_dg.transitive_reduction();
  if (verbose_level > 2) {
    std::cerr << "dependency graph: " << rp_dg << std::endl;
    print_rpdg(std::cerr, rp_actions, rp_dg);
  }
  if (verbose_level > 0)
    std::cerr << "extracting conflicts..." << std::endl;
  conflicts3(init, goal, cs);
  if (stats.break_signal_raised()) return false;
  assert(!cs.empty());
  return false;
}

// bool ILB::check_rp_both
// (const index_set& init,
//  const index_set& goal,
//  index_set_vec& seq,
//  index_set_vec& cs)
// {
//   index_set_vec seq2;
//   index_set_vec cs2;
//   bool ok2 = check_rp2(init, goal, seq2, cs2);
//   bool ok3 = rp_is_plan(init, goal, rp_actions, seq);
//   if (ok3) {
//     if (!ok2) {
//       std::cerr << "error: check_rp2 = " << ok2 << ", rp_is_plan = " << ok3
// 		<< std::endl;
//       exit(1);
//     }
//     std::cerr << "rp is valid!" << std::endl;
//     return true;
//   }
//   else if (ok2) {
//     std::cerr << "error: check_rp2 = " << ok2 << ", rp_is_plan = " << ok3
// 	      << std::endl;
//     std::cerr << "extracting conflicts..." << std::endl;
//     conflicts3(init, goal, cs);
//     exit(1);
//   }
//   std::cerr << "extracting conflicts..." << std::endl;
//   conflicts3(init, goal, cs);
//   assert(!cs.empty());
//   return false;
// }

bool ILB::check_rp
(const index_set& init,
 const index_set& goal,
 index_set_vec& seq,
 index_set_vec& cs)
{
  switch (check_rp_strategy) {
  case 3:
    return check_rp3(init, goal, seq, cs);
  default:
    return check_rp2(init, goal, seq, cs);
  }
}

NTYPE ILB::hplusplus(NTYPE bound, Plan** sol)
{
  stats.start();
  reset_hlb();
  assert(meta_atom_map.empty());
  // to be able to distinguish new meta-actions from actions in the
  // original instance, we need to clear the src field before starting
  for (index_type k = 0; k < ins.n_actions(); k++)
    ins.actions[k].src = 0;
  std::cerr << "computing relaxed plan..." << std::endl;
#ifdef HOFFMANN
  NTYPE v_max = compute_relaxed_plan_FF(ins.init_atoms, ins.goal_atoms);
  calls_to_hplus += 1;
#else
  NTYPE v_max = hplus(ins.init_atoms, ins.goal_atoms, bound, NULL);
#endif
  if (!stats.break_signal_raised()) {
    std::cerr << "rp cost = " << PRINT_NTYPE(v_max) << std::endl;
#ifndef HOFFMANN
#ifdef TRACE_PRINT_HLB
    new_lb(v_max);
#endif
#endif
  }
  //if (stats.break_signal_raised() || (v_max >= bound)) {
    stats.stop();
    return v_max;
  //}
  index_set_vec seq;
  index_set_vec new_cs;
  ce_stats.start();
  bool solved = check_rp(ins.init_atoms, ins.goal_atoms, seq, new_cs);
  if ((conflict_mod != NULL) && !solved) {
    conflict_mod->apply(*this, new_cs);
  }
  ce_stats.stop();
  while (!solved) {
    if (stats.break_signal_raised()) {
      stats.stop();
      return v_max;
    }
    if (verbose_level > 0) {
      std::cerr << "plan failed, conflicts: " << std::endl;
      for (index_type k = 0; k < new_cs.size(); k++) {
	// std::cerr << " " << ins.atoms[new_cs[k].first].name
	// 	  << " & " << ins.atoms[new_cs[k].second].name;
	// if (!cs.contains(new_cs[k]))
	//   std::cerr << " (new)";
	// std::cerr << std::endl;
	std::cerr << " ";
	ins.write_atom_set(std::cerr, new_cs[k]);
	if (meta_atom_map.find_rule(new_cs[k]) == no_such_index)
	  std::cerr << " (new)";
	std::cerr << std::endl;
      }
    }
    conflict_set_size += new_cs.size();
    index_type pre_mod_n_actions = ins.n_actions();
    EvalWithMetaAtoms* mx =
      ((inc != 0) ? new EvalWithMetaAtoms(*inc, meta_atom_map) : 0);
    std::cerr << "modifying problem..." << std::endl;
    // don't clear the meta_atom_map!!
    // meta_atom_map.clear();
#ifdef TRACE_PRINT_META_ATOMS_AS_SETS
    std::cout << ";; iterations = " << hpp_iterations << std::endl;
#endif
    bool new_conflict_found = false;
    for (index_type k = 0; k < new_cs.size(); k++) {
      assert(new_cs[k].size() >= 2);
      index_type m = meta_atom_map.find_rule(new_cs[k]);
      if (m == no_such_index) {
#ifdef TRACE_PRINT_META_ATOMS_AS_SETS
	index_set s(new_cs[k]);
	meta_atom_map.backchain_to_fixpoint(s);
	std::cout << "(:set";
	for (HSPS::index_type i = 0; i < s.size(); i++) {
	  assert(s[i] < ins.n_atoms());
	  std::cout << " " << ins.atoms[s[i]].name;
	}
	std::cout << ")" << std::endl;
#endif
	if (verbose_level > 1) {
	  std::cerr << "creating meta atom ";
	  ins.write_atom_set(std::cerr, new_cs[k]);
	  std::cerr << "..." << std::endl;
	}
	ins.create_meta_atom(new_cs[k], meta_atom_map, mx);
	if (stats.break_signal_raised()) {
	  stats.stop();
	  return v_max;
	}
	if (verbose_level > 1) {
	  std::cerr << "now " << ins.n_atoms() << " atoms and "
		    << ins.n_actions() << " actions..." << std::endl;
	}
	// std::cerr << "now " << ins.n_atoms() << " atoms and "
	// 		<< ins.n_actions() << " actions" << std::endl;
	// ins.cross_reference();
	// reach.init_structs();
	// reach.recompute(init);
	// bool_vec dead_acts(false, ins.n_actions());
	// for (index_type i = 0; i < ins.n_actions(); i++)
	// 	if (reach.unreachable(ins.actions[i].pre))
	// 	  dead_acts[i] = true;
	// if (dead_acts.count(true) > 0) {
	// 	std::cerr << "removing " << dead_acts.count(true)
	// 		  << " actions with unreachable preconditions..."
	// 		  << std::endl;
	// 	index_vec dummy_map;
	// 	ins.remove_actions(dead_acts, dummy_map);
	// 	ins.cross_reference();
	// }
	new_conflict_found = true;
      }
    }
    if (!new_conflict_found) {
      std::cerr << "NO NEW CONFLICT FOUND!!!" << std::endl;
      std::cerr << "meta-atom map:" << std::endl;
      for (index_type i = 0; i < meta_atom_map.size(); i++)
	std::cerr << " " << meta_atom_map[i] << std::endl;
      verbose_level = 4;
      solved = check_rp(ins.init_atoms, ins.goal_atoms, seq, new_cs);
      assert(!solved);
    }
    assert(new_conflict_found);
    if (remove_dominated_conditions) {
      ins.remove_dominated_conditions(meta_atom_map);
    }
    ins.cross_reference();
    pc_actions_ratio += (ins.n_actions() / (double)pre_mod_n_actions);
    hpp_iterations += 1;
    std::cerr << "new problem has " << ins.n_atoms() << " atoms and "
    	      << ins.n_actions() << " actions" << std::endl;
    // exit(0);
#ifdef HPLUSPLUS_PRINT_STATS
    print_hplusplus_stats(std::cerr);
#endif
    if (mx) delete mx;
    if (verbose_level > 2) {
      std::cerr << "meta-atom map:" << std::endl;
      for (index_type i = 0; i < meta_atom_map.size(); i++)
	std::cerr << " " << meta_atom_map[i] << std::endl;
      if (verbose_level > 4) {
	std::cerr << "new domain:" << std::endl;
	ins.print(std::cerr);
      }
    }
    reach.init_structs();
#ifdef USE_NEWLM_BC
    delete h1;
    h1 = new CostTable(ins, stats);
    h1->compute_H1(UnitACF());
#else
    if (ILA_use_lmcut) {
      delete h1;
      h1 = new CostTable(ins, stats);
    }
#endif
    std::cerr << "computing relaxed plan..." << std::endl;
#ifdef HOFFMANN
    init_preprocessing();
    NTYPE v = compute_relaxed_plan_FF(ins.init_atoms, ins.goal_atoms);
    calls_to_hplus += 1;
#else
    NTYPE v = hplus(ins.init_atoms, ins.goal_atoms, bound, NULL);
    assert((v >= v_max) || stats.break_signal_raised());
#endif
    if (!stats.break_signal_raised()) {
      std::cerr << "rp cost = " << PRINT_NTYPE(v) << std::endl;
#ifndef HOFFMANN
#ifdef TRACE_PRINT_HLB
      new_lb(v_max);
#endif
#endif
    }
    v_max = MAX(v_max, v);
    if (stats.break_signal_raised() || (v_max >= bound)) {
      stats.stop();
      return v_max;
    }
    new_cs.clear();
    ce_stats.start();
    solved = check_rp(ins.init_atoms, ins.goal_atoms, seq, new_cs);
    if (((conflict_mod != NULL) && !solved)) {
      conflict_mod->apply(*this, new_cs);
    }
    ce_stats.stop();
  }

  if (solved) {
    ActionSequence* plan = (sol != NULL ? new ActionSequence() : NULL);
    std::cerr << "plan succeeded!" << std::endl;
    for (index_type k = 0; k < seq.length(); k++)
      for (index_type i = 0; i < seq[k].size(); i++) {
	index_type a = rp_actions[seq[k][i]];
	std::cerr << ins.actions[a].name << std::endl;
	if (plan != NULL) {
	  if (ins.actions[a].src)
	    a = *((index_type*)ins.actions[a].src);
	  assert(a < ins.n_actions());
	  plan->insert(a);
	}
      }
    if (plan != NULL) {
      plan->end();
      *sol = plan;
    }
  }

  stats.stop();
  return v_max;
}

NTYPE ILB::hplusplus1
(const index_set_vec& cs_in,
 index_set_vec& cs_out,
 Plan** sol,
 bool& completed,
 bool initialise)
{
  cs_out.clear();
  completed = false;
  stats.start();
  // if this is the first call, initialise
  if (initialise) {
    reset_hlb();
    assert(meta_atom_map.empty());
    for (index_type k = 0; k < ins.n_actions(); k++)
      ins.actions[k].src = 0;
  }
  // if input conflict set is non-empty, modify the problem...
  if (cs_in.size() > 0) {
    index_type pre_mod_n_actions = ins.n_actions();
    EvalWithMetaAtoms* mx =
      ((inc != 0) ? new EvalWithMetaAtoms(*inc, meta_atom_map) : 0);
    std::cerr << "modifying problem..." << std::endl;
    for (index_type k = 0; k < cs_in.size(); k++) {
      assert(cs_in[k].size() >= 2);
      index_type m = meta_atom_map.find_rule(cs_in[k]);
      if (m == no_such_index) {
	if (verbose_level > 1) {
	  std::cerr << "creating meta atom ";
	  ins.write_atom_set(std::cerr, cs_in[k]);
	  std::cerr << "..." << std::endl;
	}
	ins.create_meta_atom(cs_in[k], meta_atom_map, mx);
	if (stats.break_signal_raised()) {
	  stats.stop();
	  return 0;
	}
	if (verbose_level > 1) {
	  std::cerr << "now " << ins.n_atoms() << " atoms and "
		    << ins.n_actions() << " actions..." << std::endl;
	}
      }
    }
    // maintainance required after problem modification
    ins.cross_reference();
    pc_actions_ratio += (ins.n_actions() / (double)pre_mod_n_actions);
    hpp_iterations += 1;
    std::cerr << "new problem has " << ins.n_atoms() << " atoms and "
    	      << ins.n_actions() << " actions" << std::endl;
    if (mx) delete mx;
    if (verbose_level > 2) {
      std::cerr << "meta-atom map:" << std::endl;
      for (index_type i = 0; i < meta_atom_map.size(); i++)
	std::cerr << " " << meta_atom_map[i] << std::endl;
      if (verbose_level > 4) {
	std::cerr << "new domain:" << std::endl;
	ins.print(std::cerr);
      }
    }
    reach.init_structs();
#ifdef USE_NEWLM_BC
    delete h1;
    h1 = new CostTable(ins, stats);
    h1->compute_H1(UnitACF());
#else
    if (ILA_use_lmcut) {
      delete h1;
      h1 = new CostTable(ins, stats);
    }
#endif
  }
  // then, compute a new relaxed plan...
  std::cerr << "computing relaxed plan..." << std::endl;
  NTYPE v_max = hplus(ins.init_atoms, ins.goal_atoms, POS_INF, NULL);
  if (!stats.break_signal_raised()) {
    std::cerr << "rp cost = " << PRINT_NTYPE(v_max) << std::endl;
#ifdef TRACE_PRINT_HLB
    new_lb(v_max);
#endif
  }
  if (INFINITE(v_max)) {
    stats.stop();
    completed = true;
    return v_max;
  }
  if (stats.break_signal_raised()) {
    stats.stop();
    return v_max;
  }
  // and check relaxed plan...
  index_set_vec seq;
  index_set_vec check_cs;
  ce_stats.start();
  bool solved = check_rp(ins.init_atoms, ins.goal_atoms, seq, check_cs);
  ce_stats.stop();
  // if not a solution, verify that we have found a new conflict
  if (!solved) {
    bool new_conflict_found = false;
    for (index_type k = 0; k < check_cs.size(); k++) {
      assert(check_cs[k].size() >= 2);
      index_type m = meta_atom_map.find_rule(check_cs[k]);
      if (m == no_such_index) {
	new_conflict_found = true;
	cs_out.append(check_cs[k]);
      }
    }
    if (!new_conflict_found) {
      std::cerr << "NO NEW CONFLICT FOUND!!!" << std::endl;
      std::cerr << "meta-atom map:" << std::endl;
      for (index_type i = 0; i < meta_atom_map.size(); i++)
	std::cerr << " " << meta_atom_map[i] << std::endl;
      verbose_level = 4;
      solved = check_rp(ins.init_atoms, ins.goal_atoms, seq, cs_out);
      assert(!solved);
    }
    assert(new_conflict_found);
  }
  // if it is a solution, extract plan
  else {
    std::cerr << "plan succeeded!" << std::endl;
    ActionSequence* plan = (sol != NULL ? new ActionSequence() : NULL);
    for (index_type k = 0; k < seq.length(); k++)
      for (index_type i = 0; i < seq[k].size(); i++) {
	index_type a = rp_actions[seq[k][i]];
	std::cerr << ins.actions[a].name << std::endl;
	if (plan != NULL) {
	  if (ins.actions[a].src)
	    a = *((index_type*)ins.actions[a].src);
	  assert(a < ins.n_actions());
	  plan->insert(a);
	}
      }
    if (plan != NULL) {
      plan->end();
      *sol = plan;
    }
  }
#ifdef HPLUSPLUS_PRINT_STATS
    print_hplusplus_stats(std::cerr);
#endif
  stats.stop();
  completed = true;
  return v_max;
}


void ILB::print_hplus_stats(std::ostream& s) const
{
  if (is_sub) return;
  s << "h+ calls: " << calls_to_hplus << std::endl;
  s << "ILA: lms generated = " << landmarks.size()
    << " (total = " << calls_to_newlm << ", avg "
    << PRINT_NTYPE(R_TO_N(calls_to_newlm, calls_to_hplus))
    << " lms/h+ call)" << std::endl;
  s << "ILA: hs(opt) calls = " << calls_to_hs_opt
    << ", nodes = " << hs_nodes << " (avg "
    << PRINT_NTYPE(R_TO_N(hs_nodes, calls_to_hs_opt))
    << " nodes/call)" << std::endl;
#ifdef MEASURE_STABILITY
  s << "average stability = " << sum_stability / n_stability << std::endl;
  std::cout << ";; average stability = " << sum_stability / n_stability
	    << std::endl;
#endif
  s << "ILA: hs(apx) calls = " << calls_to_hs_apx << " ("
    << hs_apx_improve << " improved)" << std::endl;
#ifdef MEASURE_DENSITY
  s << "ILA: average density = " << (sum_density/calls_to_hs_opt)*100
    << "%" << std::endl;
#endif
#ifndef ILA_USE_HS_EXTERN
  s << "ILA: hs lbs:";
#ifdef USE_LB1
  s << " lb1 = " << PRINT_NTYPE(R_TO_N(hs_lb1_max, hs_lb_calls)*100);
#endif
#ifdef USE_LB2
  s << " lb2 = " << PRINT_NTYPE(R_TO_N(hs_lb2_max, hs_lb_calls)*100);
#endif
#ifdef USE_LB3
  s << " lb3 = " << PRINT_NTYPE(R_TO_N(hs_lb3_max, hs_lb_calls)*100);
#endif
#ifdef USE_LB4
  s << " lb4 = " << PRINT_NTYPE(R_TO_N(hs_lb4_max, hs_lb_calls)*100);
#endif
  s << std::endl;
  if (hitting_set_use_split) {
    s << "ILA: " << PRINT_NTYPE(R_TO_N(hs_split1, hs_nodes)*100)
      << "% hs nodes split" << std::endl;
    s << "ILA: " << PRINT_NTYPE(R_TO_N(hs_splits, hs_nodes)*100)
      << "% hs nodes properly split" << std::endl;
  }
  s << "ILA: " << PRINT_NTYPE(R_TO_N(hs_branch, hs_nodes)*100)
    << "% hs nodes branched" << std::endl;
  s << "ILA: "
    << PRINT_NTYPE(R_TO_N(hs_cache_hits, hs_cache_hits + hs_cache_miss)*100)
    << "% hs cache hits" << std::endl;
#endif // !ILA_USE_HS_EXTERN
  s << "time (ILB) = " << stats.total_time() << std::endl;
  s << "time (relevance analysis) = "
    << ra_stats.total_time() << std::endl;
  s << "time (approx hitting set) = " << hsa_stats.total_time()
    << " (" << hsa_stats.total_time() / calls_to_hs_apx
    << " sec/call)" << std::endl;
  s << "time (optimal hitting set) = " << hss_stats.total_time()
    << " (" << hss_stats.total_time() / calls_to_hs_opt
    << " sec/call)" << std::endl;
  s << "time (lm generation) = " << lm_stats.total_time()
    << " (" << lm_stats.total_time() / calls_to_newlm
    << " sec/call)" << std::endl;
}

void ILB::print_hplusplus_stats(std::ostream& s) const
{
  s << "average conflict set size = "
    << PRINT_NTYPE(R_TO_N(conflict_set_size, hpp_iterations))
    << std::endl;
  s << "average domain size increase = x"
    << pc_actions_ratio/hpp_iterations
    << std::endl;
  s << "time (rp normalisation) = "
    << rpn_stats.total_time() << std::endl;
  s << "time (conflict extraction) = "
    << ce_stats.total_time() << std::endl;
}

///
// hplusplus_with_ce
///


// A simplified "first fail" variant of exec_rp2, adapted to condeffs

void ILB::partition_applicable_with_ce
(const index_vec& rp,
 const index_set_vec& rp_ce,
 const bool_vec& rem_rp,
 const bool_vec& state,
 index_set_vec& apps)
{
  for (index_type k = 0; k < rp.size(); k++)
    if (rem_rp[k]) {
      bool app = state.contains(ins.actions[rp[k]].pre);
      for (index_type i = 0; (i < rp_ce[k].size()) && app; i++)
	if (!state.contains(ins.actions[rp[k]].cadd[rp_ce[k][i]].antecedent))
	  app = false;
      if (app) {
	bool placed = false;
	for (index_type i = 0; (i < apps.size()) && !placed; i++) {
	  bool can_go = true;
	  for (index_type j = 0; (j < apps[i].size()) && can_go; j++)
	    if (!ins.commutative(rp[apps[i][j]], rp[k]))
	      can_go = false;
	  if (can_go) {
	    apps[i].insert(k);
	    placed = true;
	  }
	}
	if (!placed) {
	  apps.append(EMPTYSET);
	  apps[apps.length() - 1].insert(k);
	}
      }
    }
}

void ILB::choose_from_node_conflict_set
(const index_set_graph& dg,
 const node_conflict_set& ncs,
 index_set_vec& cs)
{
  assert(ncs.size() > 0);
  bool_vec case1(false, ncs.size());
  bool all_are_case_2 = true;
  for (index_type k = 0; k < ncs.size(); k++)
    if (dg.reachable(ncs[k].first.first, ncs[k].first.second)) {
      case1[k] = true;
      all_are_case_2 = false;
    }
  // std::cerr << "case1 = " << case1
  // 	    << ", all2 = " << all_are_case_2
  // 	    << std::endl;
  bool first = true;
  index_set_vec best;
  index_set_vec tmp(cs);
  for (index_type k = 0; k < ncs.size(); k++)
    if (case1[k] || all_are_case_2) {
      if (verbose_level > 3) {
	std::cerr << "set in = " << tmp << std::endl;
      }
      node_conflict_to_atom_conflicts(dg, ncs[k], tmp);
      if (verbose_level > 3) {
	std::cerr << "set out = " << tmp << std::endl;
      }
      index_set_vec new_cs;
      for (index_type i = 0; i < tmp.size(); i++)
	if (cs.first(tmp[i]) == no_such_index)
	  new_cs.append(tmp[i]);
      if (new_cs.size() == 0) return;
      if (first || (new_cs.size() < best.size()))
	best.assign_copy(new_cs);
      first = false;
    }
  assert(!first);
  cs.append(best);
}

void ILB::compute_dependency_closure_in_split_dg
(const index_vec& rp,
 const index_set_vec& rp_ce,
 const pair_vec& split_index,
 const index_set_graph& dg,
 index_type n_from,
 index_type n_to,
 index_set& dc_labels)
{
  bool_vec dc(false, dg.size());
  dc[n_from] = true;
  index_vec p;
  index_type d = dg.shortest_path(n_from, n_to, p);
  assert(d != no_such_index);
  for (index_type i = 0; (i + 1) < p.length(); i++) {
    dc[p[i + 1]] = true;
    if (dg.edge_has_label(p[i], p[i + 1])) {
      const index_set& e = dg.edge_label(p[i], p[i + 1]);
      assert(!e.empty());
      dc_labels.insert(e[0]);
    }
    else {
      // an edge with an empty label should be a back-edge from a
      // ce to its action; verify that
      assert((split_index[p[i]].first == split_index[p[i + 1]].first) &&
	     (split_index[p[i + 1]].second == no_such_index));
    }
  }
  if (verbose_level > 3) {
    std::cerr << "dc from " << n_from << " to " << n_to << " in split dg:"
	      << std::endl << "path #0 = " << p << std::endl
	      << " nodes = " << index_set(dc) << ", labels = " << dc_labels
	      << std::endl;
  }
  bool done = false;
  while (!done) {
    done = true;
    for (index_type i = 0; i < dc_labels.size(); i++)
      for (index_type k = 0; k < split_index.size(); k++)
	if (!dc[k]) {
	  index_type action_in_rp = split_index[k].first;
	  assert(action_in_rp < rp.size());
	  index_type a = rp[action_in_rp];
	  bool check = false;
	  if (split_index[k].second == no_such_index) {
	    if (ins.actions[a].add.contains(dc_labels[i]))
	      check = true;
	  }
	  else {
	    assert(split_index[k].second < rp_ce[action_in_rp].size());
	    index_type j = rp_ce[action_in_rp][split_index[k].second];
	    assert(j < ins.actions[a].cadd.size());
	    if (ins.actions[a].cadd[j].consequent == dc_labels[i])
	      check = true;
	  }
	  if (check) {
	    p.clear();
	    d = dg.shortest_path(n_from, k, p);
	    assert(d != no_such_index);
	    for (index_type j = 0; (j + 1) < p.length(); j++) {
	      dc[p[j + 1]] = true;
	      if (dg.edge_has_label(p[j], p[j + 1])) {
		const index_set& e = dg.edge_label(p[j], p[j + 1]);
		assert(!e.empty());
		dc_labels.insert(e[0]);
	      }
	      else {
		assert((split_index[p[j]].first
			== split_index[p[j + 1]].first) &&
		       (split_index[p[j + 1]].second == no_such_index));
	      }
	    }
	    if (verbose_level > 3) {
	      std::cerr << "dc from " << n_from << " to " << n_to
			<< " in split dg:" << std::endl << "new path= " << p
			<< std::endl << " nodes = " << index_set(dc)
			<< ", labels = " << dc_labels << std::endl;
	    }
	    done = false;
	  }
	}
  }
}

void ILB::node_conflict_to_atom_conflicts_with_dc_in_split_dg
(const index_vec& rp,
 const index_set_vec& rp_ce,
 const pair_vec& split_index,
 const index_set_graph& dg1,
 const index_set_graph& dg2, // dg2 has the back-edges
 const node_conflict& nc,
 index_set_vec& acs)
{
  index_type deler = nc.first.first;
  index_type failed = nc.first.second;
  index_type p = nc.second;
  if (verbose_level > 2) {
    std::cerr << "(full dc) deler = " << deler << ", failed = " << failed
	      << ", p = " << p << std::endl;
  }
  if (dg1.reachable(deler, failed)) {
    index_set ps;
    compute_dependency_closure_in_split_dg
      (rp, rp_ce, split_index, dg2, deler, failed, ps);
    if (verbose_level > 2) {
      std::cerr << "(full dc) case 1: ps = " << ps << std::endl;
    }
    for (index_type j = 0; j < ps.size(); j++)
      acs.append(index_set(p, ps[j]));
  }
  else {
    index_type ncd;
    estimate_conflict_weight(dg1, index_pair(deler, failed), ncd);
    assert((ncd != no_such_index) && (ncd < dg1.size()));
    // std::cerr << "deler = " << deler << ", failed = " << failed
    //	<< ", ncd = " << ncd << std::endl;
    index_set ps1;
    compute_dependency_closure_in_split_dg
      (rp, rp_ce, split_index, dg2, deler, ncd, ps1);
    index_set ps2;
    compute_dependency_closure_in_split_dg
      (rp, rp_ce, split_index, dg2, failed, ncd, ps2);
    if (verbose_level > 2) {
      std::cerr << "(full dc) case 2: ncd = " << ncd
		<< std::endl << "ps1 = " << ps1
		<< std::endl << "ps2 = " << ps2
		<< std::endl;
    }
    assert(!ps1.contains(p));
    ps2.insert(p);
    for (index_type j1 = 0; j1 < ps1.size(); j1++)
      for (index_type j2 = 0; j2 < ps2.size(); j2++) {
	assert(ps1[j1] != ps2[j2]);
	acs.append(index_set(ps1[j1], ps2[j2]));
      }
  }
}

void ILB::choose_from_node_conflict_set_with_dc
(const index_vec& rp,
 const index_set_vec& rp_ce,
 const pair_vec& split_index,
 const index_set_graph& dg1,
 const index_set_graph& dg2,
 const node_conflict_set& ncs,
 index_set_vec& cs)
{
  assert(ncs.size() > 0);
  bool_vec case1(false, ncs.size());
  bool all_are_case_2 = true;
  for (index_type k = 0; k < ncs.size(); k++)
    if (dg1.reachable(ncs[k].first.first, ncs[k].first.second)) {
      case1[k] = true;
      all_are_case_2 = false;
    }
  // std::cerr << "case1 = " << case1
  // 	    << ", all2 = " << all_are_case_2
  // 	    << std::endl;
  bool first = true;
  index_set_vec best;
  index_set_vec tmp(cs);
  for (index_type k = 0; k < ncs.size(); k++)
    if (case1[k] || all_are_case_2) {
      if (verbose_level > 3) {
	std::cerr << "set in = " << tmp << std::endl;
      }
      node_conflict_to_atom_conflicts_with_dc_in_split_dg
	(rp, rp_ce, split_index, dg1, dg2, ncs[k], tmp);
      if (verbose_level > 3) {
	std::cerr << "set out = " << tmp << std::endl;
      }
      index_set_vec new_cs;
      for (index_type i = 0; i < tmp.size(); i++)
	if (cs.first(tmp[i]) == no_such_index)
	  new_cs.append(tmp[i]);
      if (new_cs.size() == 0) return;
      if (first || (new_cs.size() < best.size()))
	best.assign_copy(new_cs);
      first = false;
    }
  assert(!first);
  cs.append(best);
}

void ILB::split_deleters
(const pair_vec& split_index,
 const index_set_graph& split_dg,
 node_conflict_set& ncs)
{
  if (verbose_level > 2) {
    std::cerr << "splitting deleters:" << std::endl << " in =";
    for (index_type k = 0; k < ncs.size(); k++)
      std::cerr << " (" << ncs[k].first << "," << ncs[k].second << ")";
    std::cerr << std::endl;
  }
  node_conflict_set ncs1(ncs);
  ncs.clear();
  for (index_type k = 0; k < ncs1.size(); k++) {
    index_type n1 = ncs1[k].first.first; // original deleter node
    for (index_type i = 0; i < split_index.size(); i++)
      if (split_index[i].first == n1)
	if (split_dg.out_degree(i) > 0)
	  ncs.append(node_conflict(index_pair(i, ncs1[k].first.second),
				   ncs1[k].second));
  }
  if (verbose_level > 2) {
    std::cerr << " out =";
    for (index_type k = 0; k < ncs.size(); k++)
      std::cerr << " (" << ncs[k].first << "," << ncs[k].second << ")";
    std::cerr << std::endl;
  }
}

bool ILB::exec_rp_with_ce
(bool_vec& rem_rp,
 lvector<bool_vec>& state,
 lvector<bool_vec>& rstate,
 const index_set& goal,
 const index_vec& rp,
 const index_set_vec& rp_ce,
 const pair_vec& split_index,
 const index_set_graph& split_dg,
 const index_set_graph& dg2, // split dg with back-edges
 bool with_full_dc,
 index_type step,
 index_set_vec& seq,
 index_set_vec& cs)
{
  if (verbose_level > 2) {
    for (index_type i = 0; i < step; i++) std::cerr << "  ";
    std::cerr << "-> " << seq << std::endl;
    std::cerr << rstate[step] << std::endl;
    std::cerr << state[step] << std::endl;
  }
  if (rem_rp.count(true) == 0) {
    assert(rstate[step].contains(goal));
    if (state[step].contains(goal)) {
      return true;
    }
    else {
      node_conflict_set ncs;
      assert((split_index.size() + 1) == split_dg.size());
      conflicts2(split_index.size(), goal, state, rstate, rp,
		 step, seq, true, ncs);
      split_deleters(split_index, split_dg, ncs);
      if (with_full_dc) {
	choose_from_node_conflict_set_with_dc
	  (rp, rp_ce, split_index, split_dg, dg2, ncs, cs);
      }
      else {
	choose_from_node_conflict_set(split_dg, ncs, cs);
      }
      return false;
    }
  }
  index_set_vec pnext;
  partition_applicable_with_ce(rp, rp_ce, rem_rp, rstate[step], pnext);
  if (verbose_level > 3) {
    std::cerr << "partitioned applicable: " << pnext << std::endl;
  }
  assert(pnext.length() > 0);
  for (index_type k = 0; k < pnext.length(); k++) {
    state[step + 1] = state[step];
    rstate[step + 1] = rstate[step];
    bool ok = true;
    for (index_type i = 0; (i < pnext[k].size()) && ok; i++) {
      index_type rpi = pnext[k][i]; // relaxed plan index
      index_type a = rp[rpi]; // action at rpi
      assert(rstate[step].contains(ins.actions[a].pre));
      // if the actions prec are not true, the action itself fails,
      // as do all of it's condeffs:
      if (!state[step].contains(ins.actions[a].pre)) {
	bool non_red_node_found = false;
	index_type nf = split_index.first(index_pair(rpi, no_such_index));
	assert(nf != no_such_index);
	node_conflict_set ncs;
	conflicts2(nf, ins.actions[a].pre, state, rstate, rp, step,
		   seq, true, ncs);
	split_deleters(split_index, split_dg, ncs);
	// we can only generate flaws from the action itself if it
	// has essential effects:
	if (split_dg.out_degree(nf) > 0) {
	  if (with_full_dc) {
	    choose_from_node_conflict_set_with_dc
	      (rp, rp_ce, split_index, split_dg, dg2, ncs, cs);
	  }
	  else {
	    choose_from_node_conflict_set(split_dg, ncs, cs);
	  }
	  non_red_node_found = true;
	}
	for (index_type j = 0; j < rp_ce[rpi].size(); j++) {
	  nf = split_index.first(index_pair(rpi, j));
	  assert(nf != no_such_index);
	  // some ce's can be redundant, but not all:
	  if (split_dg.out_degree(nf) > 0) {
	    for (index_type l = 0; l < ncs.size(); l++)
	      ncs[l].first.second = nf;
	    if (with_full_dc) {
	      choose_from_node_conflict_set_with_dc
		(rp, rp_ce, split_index, split_dg, dg2, ncs, cs);
	    }
	    else {
	      choose_from_node_conflict_set(split_dg, ncs, cs);
	    }
	    non_red_node_found = true;
	  }
	}
	assert(non_red_node_found);
	ok = false;
      }
      if (ok) {
	bool ce_ok = true;
	bool non_red_node_found = false;
	for (index_type j = 0; j < rp_ce[rpi].size(); j++) {
	  index_type cj = rp_ce[rpi][j];
	  assert(rstate[step].contains(ins.actions[a].cadd[cj].antecedent));
	  // if one of its (included) effects conditions is not true,
	  // the condeff fails. note that there can be several failed
	  // ce's, and we need to choose the "right" one; here, we take
	  // the easy approach of generating separate conflicts for each
	  // failed ce.
	  if (!state[step].contains(ins.actions[a].cadd[cj].antecedent)) {
	    index_type nf = split_index.first(index_pair(rpi, j));
	    assert(nf != no_such_index);
	    if (split_dg.out_degree(nf) > 0) {
	      node_conflict_set ncs;
	      conflicts2(nf, ins.actions[a].cadd[cj].antecedent,
			 state, rstate, rp, step, seq, true, ncs);
	      split_deleters(split_index, split_dg, ncs);
	      if (with_full_dc) {
		choose_from_node_conflict_set_with_dc
		  (rp, rp_ce, split_index, split_dg, dg2, ncs, cs);
	      }
	      else {
		choose_from_node_conflict_set(split_dg, ncs, cs);
	      }
	      non_red_node_found = true;
	    }
	    ce_ok = false;
	  }
	}
	// if some ce failed, then some non-redundant ce should fail:
	if (!ce_ok) {
	  assert(non_red_node_found);
	  ok = false;
	}
      }
    }
    if (ok) {
      // if we get to here, every action in pnext[k] is applicable
      // (in real state); only now can we go on and apply them.
      for (index_type i = 0; i < pnext[k].size(); i++) {
	index_type rpi = pnext[k][i]; // relaxed plan index
	index_type a = rp[rpi]; // action at rpi
	state[step + 1].subtract(ins.actions[a].del);
	for (index_type j = 0; j < ins.actions[a].cdel.size(); j++)
	  if (state[step].contains(ins.actions[a].cdel[j].antecedent))
	    state[step + 1][ins.actions[a].cdel[j].consequent] = false;
	state[step + 1].insert(ins.actions[a].add);
	for (index_type j = 0; j < ins.actions[a].cadd.size(); j++)
	  if (state[step].contains(ins.actions[a].cadd[j].antecedent))
	    state[step + 1][ins.actions[a].cadd[j].consequent] = true;
	rstate[step + 1].insert(ins.actions[a].add);
	for (index_type j = 0; j < rp_ce[rpi].size(); j++) {
	  index_type cj = rp_ce[rpi][j];
	  rstate[step + 1][ins.actions[a].cadd[cj].consequent] = true;
	}
	// for (index_type j = 0; j < ins.actions[a].cadd.size(); j++)
	// 	if (rstate[step].contains(ins.actions[a].cadd[j].antecedent))
	// 	  rstate[step + 1][ins.actions[a].cadd[j].consequent] = true;
	rem_rp[rpi] = false;
      }
      seq[step] = pnext[k];
      ok = exec_rp_with_ce(rem_rp, state, rstate, goal, rp, rp_ce,
			   split_index, split_dg, dg2, with_full_dc,
			   step + 1, seq, cs);
      if (ok) return true;
      if (stats.break_signal_raised()) return false;
      for (index_type i = 0; i < pnext[k].size(); i++)
	rem_rp[pnext[k][i]] = true;
    }
  }
  return false;
}

void ILB::make_rp_instance
(const index_vec& rp,
 const index_set_vec& rp_ce,
 Instance& rp_ins,
 pair_vec& split_index)
{
  split_index.clear();
  rp_ins.copy_atoms(ins);
  for (index_type k = 0; k < rp.length(); k++) {
    split_index.append(index_pair(k, no_such_index));
    Instance::Action& a = rp_ins.new_action(ins.actions[rp[k]].name);
    a.pre = ins.actions[rp[k]].pre;
    a.add = ins.actions[rp[k]].add;
    for (index_type i = 0; i < rp_ce[k].size(); i++) {
      split_index.append(index_pair(k, i));
      a.pre.insert(ins.actions[rp[k]].cadd[rp_ce[k][i]].antecedent);
      a.add.insert(ins.actions[rp[k]].cadd[rp_ce[k][i]].consequent);
    }
    a.cost = ins.actions[rp[k]].cost;
  }
  rp_ins.cross_reference();
}

void ILB::normalise_rp_with_ce
(index_vec& rp, index_set_vec& rp_ce)
{
  rpn_stats.start();
  index_type rp_in = rp.size();
  index_type rp_ce_in = 0;
  index_type rp_ce_out = 0;
  if (verbose_level > 0) {
    for (index_type k = 0; k < rp.size(); k++)
      rp_ce_in += rp_ce[k].size();
  }
  // first, add implied ce's
  for (index_type k = 0; k < rp.size(); k++) {
    index_set c(ins.actions[rp[k]].pre);
    for (index_type i = 0; i < rp_ce[k].size(); i++)
      c.insert(ins.actions[rp[k]].cadd[rp_ce[k][i]].antecedent);
    for (index_type i = 0; i < ins.actions[rp[k]].cadd.size(); i++)
      if (c.contains(ins.actions[rp[k]].cadd[i].antecedent)) {
	rp_ce[k].insert(i);
	if (verbose_level > 2) {
	  std::cerr << "adding implied ce ";
	  ins.print_rule(std::cerr, ins.actions[rp[k]].cadd[i]);
	  std::cerr << " to action " << k << "."
		    << ins.actions[rp[k]].name << std::endl;
	}
      }
  }
  // then, remove redundant ce's (and actions)
  pair_vec split_index;
  Instance rp_ins;
  //make_rp_instance(rp, rp_ce, rp_ins, split_index);
  bool done = false;
  while (!done) {
    done = true;
    //Reachability r(rp_ins);
    //assert(!r.unreachable(rp_ins.goal_atoms));
    index_type n_red = 0;
    for (index_type k = 0; (k < rp.size()) && done; k++) {
      // assert(rp_ins.actions[k].name == ins.actions[rp[k]].name);
      // if action k has ce's, chec if any of them are redundant:
      if (!rp_ce[k].empty()) {
	// only non-implied ce's are candidates for removal:
	// build a graph of implications between ce antecedents...
	graph ce_imp(rp_ce[k].size());
	for (index_type i = 0; i < rp_ce[k].size(); i++)
	  for (index_type j = i + 1; j < rp_ce[k].size(); j++)
	    if (ins.actions[rp[k]].cadd[rp_ce[k][i]].antecedent.
		contains(ins.actions[rp[k]].cadd[rp_ce[k][j]].antecedent))
	      ce_imp.add_edge(i, j);
	    else if (ins.actions[rp[k]].cadd[rp_ce[k][j]].antecedent.
		     contains(ins.actions[rp[k]].cadd[rp_ce[k][i]].antecedent))
	      ce_imp.add_edge(j, i);
	// ...and compute the scc tree of this graph...
	ce_imp.strongly_connected_components();
	graph ce_imp_tree;
	ce_imp.component_tree(ce_imp_tree);
	for (index_type i = 0; (i < ce_imp_tree.size()) && done; i++)
	  // each leaf in this scc tree is a set of mutually implying
	  // ce's (i.e., ce's with identical antecedents) that are not
	  // implied by other ce's; now we test this set for redundancy
	  if (ce_imp_tree.in_degree(i) == 0) {
	    index_set ce_set;
	    for (index_type j = 0; j < rp_ce[k].size(); j++)
	      if (ce_imp.component(j) == i)
		ce_set.insert(rp_ce[k][j]);
	    rp_ce[k].subtract(ce_set);
	    make_rp_instance(rp, rp_ce, rp_ins, split_index);
	    Reachability r(rp_ins);
	    if (!r.unreachable(rp_ins.goal_atoms)) {
	      if (verbose_level > 2) {
		std::cerr << "removing redundant ce's";
		for (index_type j = 0; j < ce_set.size(); j++)
		  ins.print_rule(std::cerr << " ", ins.actions[rp[k]].cadd[j]);
		std::cerr << " from action " << k << "."
			  << ins.actions[rp[k]].name << std::endl;
	      }
	      done = false;
	    }
	    else {
	      rp_ce[k].insert(ce_set);
	    }
	  }
      }
      // if action k has no ce's, check if the action itself is redundant:
      else {
	index_vec new_rp(rp);
	new_rp.remove(k);
	index_set_vec new_rp_ce(rp_ce);
	new_rp_ce.remove(k);
	make_rp_instance(new_rp, new_rp_ce, rp_ins, split_index);
	Reachability r(rp_ins);
	if (!r.unreachable(rp_ins.goal_atoms)) {
	  if (verbose_level > 2) {
	    std::cerr << "removing redundant action "
		      << ins.actions[rp[k]].name << std::endl;
	  }
	  rp.remove(k);
	  rp_ce.remove(k);
	  done = false;
	}
      }
    }
  }
  if (verbose_level > 0) {
    for (index_type k = 0; k < rp.size(); k++)
      rp_ce_out += rp_ce[k].size();
    std::cerr << "normalise rp: (" << rp_in << "," << rp_ce_in
	      << ") -> (" << rp.size() << "," << rp_ce_out << ")"
	      << std::endl;
  }
  rpn_stats.stop();
}

void ILB::compute_split_dependency_graph
(const index_vec& rp,
 index_set_vec& rp_ce,
 pair_vec& split_index,
 index_set_graph& dg)
{
  Instance rp_ins;
  make_rp_instance(rp, rp_ce, rp_ins, split_index);
  if (verbose_level > 3) {
    std::cerr << "splint index:" << std::endl;
    for (index_type k = 0; k < split_index.size(); k++)
      std::cerr << k << " = " << split_index[k] << std::endl;
    std::cerr << "rp instance:" << std::endl;
    rp_ins.print(std::cerr);
  }
  // compute the merged RPDG
  CostACF i_cost(rp_ins);
  ILB ilb2(rp_ins, i_cost, NULL, stats);
  ilb2.rp_actions.fill(rp_ins.n_actions());
  ilb2.reach.recompute();
  ////
  // extended debug printout for assertion below
  ////
  // if (ilb2.reach.unreachable(rp_ins.goal_atoms)) {
  //   rp_ins.print(std::cerr);
  //   std::cerr << "goal unreachability:" << std::endl;
  //   for (index_type i = 0; i < rp_ins.goal_atoms.size(); i++) {
  //     bool u = ilb2.reach.unreachable(rp_ins.goal_atoms[i]);
  //     std::cerr << rp_ins.goal_atoms[i] << " = "
  // 		<< rp_ins.atoms[rp_ins.goal_atoms[i]].name
  // 		<< ": " << u << std::endl;
  //   }
  //   CostTable* h1 = new CostTable(rp_ins, stats);
  //   UnitACF unit_cost;
  //   h1->compute_H1(unit_cost);
  //   index_set_graph hmg;
  //   // h1->compute_hm_graph(unit_cost, rp_ins.goal_atoms, true, hmg);
  //   h1->compute_complete_hm_graph(unit_cost, 1, rp_ins.goal_atoms, true, hmg);
  //   ((labeled_graph<index_set, index_set>&)hmg).
  //     write_digraph(std::cerr, false, true, true, false, "h-m Graph");
  // }
  assert(!ilb2.reach.unreachable(rp_ins.goal_atoms));
#ifdef HOFFMANN
  ilb2.compute_non_redundant_rp(rp_ins.init_atoms, rp_ins.goal_atoms, false);
#else
  ilb2.compute_non_redundant_rp(rp_ins.init_atoms, rp_ins.goal_atoms, true);
#endif
  if (stats.break_signal_raised()) return;
  assert(ilb2.rp_actions.size() == rp_ins.n_actions());
  if (verbose_level > 0)
    std::cerr << "computing RPDG..." << std::endl;
  index_set_graph dg1;
  ilb2.compute_dependency_graph(rp_ins.init_atoms, rp_ins.goal_atoms, dg1);
  if (verbose_level > 2) {
    ((labeled_graph<index_set, index_set>)dg1).
      write_digraph(std::cerr, true, false, true, false,
		    "merged RPDG (before t.r.)");
  }
  // here, we need a "label-sensitive" transitive reduction of dg1
  dg1.edge_label_preserving_transitive_reduction();
  if (verbose_level > 2) {
    ((labeled_graph<index_set, index_set>)dg1).
      write_digraph(std::cerr, true, false, true, false, "merged RPDG");
    index_set atoms;
    for (index_type i = 0; i < dg1.size(); i++)
      for (index_type j = 0; j < dg1.size(); j++)
  	if (dg1.adjacent(i,j))
  	  if (dg1.edge_has_label(i,j))
  	    atoms.insert(dg1.edge_label(i,j));
    for (index_type i = 0; i < atoms.size(); i++)
      std::cerr << atoms[i] << " = " << ins.atoms[atoms[i]].name << std::endl;
  }
  assert(dg1.size() == (rp.size() + 1));
  // now, to compute the split RPDG:
  dg.init(split_index.size() + 1);
  // consider two nodes, i and j, such that the "merged" RPDG (dg1)
  // has an edge i -> j with label L; this edge will be broken up
  // into edges from i and i's ce-nodes to j and j's ce-node as
  // follows:
  // 1(a). an atom p added by a ce of action[i] will label outgoing
  //  edges from the corresponding ce node (this works because ce's
  //  are single-effect, so the atom added by a ce is always directly
  //  needed by some other action/ce, or the ce would be redundant);
  // 1(a). all other atoms will label outgoing edges from the
  //  actions own node;
  // 2(a). in the label on the edge to node i, keep only atoms that
  //  are in the actions precondition;
  // 2(b). in the label on the edge to the j:th ce node, keep only
  //  atoms that are in the actions precondition or the ce's antecedent;
  // finally, don't add any edge with empty label.
  for (index_type i = 0; i < rp.size(); i++) {
    for (index_type j = 0; j < rp.size(); j++)
      if (dg1.adjacent(i, j)) {
	assert(dg1.edge_has_label(i, j));
	assert(!dg1.edge_label(i, j).empty());
	const index_set& label(dg1.edge_label(i, j));
	for (index_type k = 0; k < label.size(); k++) {
	  // find which ce (if any) of i adds label[k]:
	  index_type ci = no_such_index;
	  for (index_type l = 0; l < rp_ce[i].size(); l++)
	    if (ins.actions[rp[i]].cadd[rp_ce[i][l]].consequent == label[k])
	      ci = l;
	  index_type src_node = split_index.first(index_pair(i, ci));
	  assert(src_node != no_such_index);
	  // if label[k] belongs to action[i]'s pre, add it to all
	  // edges from src_node to action i and it's ce nodes:
	  if (ins.actions[rp[j]].pre.contains(label[k])) {
	    index_type dst_node =
	      split_index.first(index_pair(j, no_such_index));
	    assert(dst_node != no_such_index);
	    dg.add_edge(src_node, dst_node, label[k]);
	    for (index_type l = 0; l < rp_ce[j].size(); l++) {
	      index_type dst_node = split_index.first(index_pair(j, l));
	      assert(dst_node != no_such_index);
	      dg.add_edge(src_node, dst_node, label[k]);
	    }
	  }
	  // else, add it only to those ce nodes where label[k]
	  // belongs to the antecendent:
	  else {
	    for (index_type l = 0; l < rp_ce[j].size(); l++)
	      if (ins.actions[rp[j]].cadd[rp_ce[j][l]].antecedent.
		  contains(label[k])) {
		index_type dst_node = split_index.first(index_pair(j, l));
		assert(dst_node != no_such_index);
		dg.add_edge(src_node, dst_node, label[k]);
	      }
	  }
	}
      }
    // now, we also need to take care of edges to the goal node:
    index_type ng1 = rp.size(); // goal node in dg1
    index_type dst_node = split_index.size(); // goal node in split dg
    if (dg1.adjacent(i, rp.size())) {
      const index_set& label(dg1.edge_label(i, ng1));
      for (index_type k = 0; k < label.size(); k++) {
	// find which ce (if any) of i adds label[k]:
	index_type ci = no_such_index;
	for (index_type l = 0; l < rp_ce[i].size(); l++)
	  if (ins.actions[rp[i]].cadd[rp_ce[i][l]].consequent == label[k])
	    ci = l;
	index_type src_node = split_index.first(index_pair(i, ci));
	assert(src_node != no_such_index);
	dg.add_edge(src_node, dst_node, label[k]);
      }
    }
  }
  if (verbose_level > 2) {
    ((labeled_graph<index_set, index_set>)dg).
      write_digraph(std::cerr, true, false, true, false, "split RPDG");
    for (index_type i = 0; i < split_index.size(); i++) {
      if (split_index[i].second == no_such_index) {
	std::cerr << i << " = " << split_index[i] << " = "
		  << ins.actions[rp[split_index[i].first]].name
		  << std::endl;
      }
      else {
	std::cerr << i << " = " << split_index[i] << " = ";
	index_type k = split_index[i].first;
	index_type j = split_index[i].second;
	ins.print_rule(std::cerr, ins.actions[rp[k]].cadd[rp_ce[k][j]]);
	std::cerr << std::endl;
      }
    }
  }
  index_set_vec node_add(EMPTYSET, dg.size());
  for (index_type i = 0; i < split_index.size(); i++) {
    index_type k = split_index[i].first;
    index_type j = split_index[i].second;
    if (j == no_such_index)
      node_add[i] = ins.actions[rp[k]].add;
    else
      node_add[i].assign_singleton(ins.actions[rp[k]].cadd[rp_ce[k][j]].consequent);
  }
  // for (index_type i = 0; i < dg.size(); i++)
  //   for (index_type j = 0; j < dg.size(); j++)
  //     if (dg.adjacent(i, j)) {
  // 	assert(dg.edge_has_label(i,j));
  // 	for (index_type k = 0; k < dg.edge_label(i,j).size(); k++) {
  // 	  index_type p = dg.edge_label(i,j)[k];
  // 	  for (index_type n = 0; n < dg.size(); n++)
  // 	    if ((n != i) && node_add[n].contains(p)) {
  // 	      std::cerr << "p = " << p << " added by " << n << std::endl;
  // 	      assert(dg.reachable(i,n));
  // 	    }
  // 	}
  //     }
}

// add edges, with empty label, from each ce node back to it's action.
void ILB::add_back_edges
(const pair_vec& split_index,
 index_set_graph& split_dg)
{
  assert(split_dg.size() == (split_index.size() + 1));
  for (index_type k = 0; k < split_index.size(); k++)
    if (split_index[k].second != no_such_index) {
      index_type l = split_index.first(index_pair(split_index[k].first,
						  no_such_index));
      assert((l != no_such_index) && (l < split_index.size()));
      split_dg.add_edge(k, l);
    }
}

void ILB::filter_flaws(index_set_vec& cs)
{
  bool_vec cs_to_remove(false, cs.size());
  for (index_type k = 0; k < cs.size(); k++) {
    index_set_vec u(EMPTYSET, cs[k].size());
    for (index_type i = 0; i < cs[k].size(); i++) {
      u[i].assign_singleton(cs[k][i]);
      meta_atom_map.backchain_to_fixpoint(u[i]);
    }
    bool_vec atoms_to_remove(false, cs[k].size());
    for (index_type i = 0; i < cs[k].size(); i++) {
      bool found = false;
      for (index_type j = 0; (j < cs[k].size()) && !found; j++)
	if ((i != j) && !atoms_to_remove[j] && u[j].contains(u[i]))
	  found = true;
      if (found)
	atoms_to_remove[i] = true;
    }
    if ((verbose_level > 1) && (atoms_to_remove.count(true) > 0)) {
      std::cerr << "filter_conflicts: ";
      ins.print_atom_set(std::cerr, cs[k]);
    }
    cs[k].remove(atoms_to_remove);
    if ((verbose_level > 1) && (atoms_to_remove.count(true) > 0)) {
      std::cerr << " => ";
      ins.print_atom_set(std::cerr, cs[k]);
      std::cerr << std::endl;
    }
    if (cs[k].size() < 2)
      cs_to_remove[k] = true;
  }
  cs.remove(cs_to_remove);
}

bool ILB::check_rp_with_ce
(const index_set& init,
 const index_set& goal,
 const index_vec& rp,
 const index_set_vec& rp_ce,
 bool with_full_dc,
 index_set_vec& seq,
 index_set_vec& cs)
{
  if (verbose_level > 1) {
    std::cerr << "checking relaxed plan:" << std::endl;
    print_rp_with_ce(std::cerr, rp, rp_ce);
    if (verbose_level > 4) {
      ins.print(std::cerr);
    }
  }
  if (verbose_level > 0)
    std::cerr << "normalising rp (with ce)..." << std::endl;
  index_vec new_rp(rp);
  index_set_vec new_rp_ce(rp_ce);
  normalise_rp_with_ce(new_rp, new_rp_ce);
  if (verbose_level > 1) {
    std::cerr << "rp after normalisation:" << std::endl;
    print_rp_with_ce(std::cerr, new_rp, new_rp_ce);
    //bool ok = validate_rp_with_ce(new_rp, new_rp_ce);
    //assert(ok);
  }
  if (stats.break_signal_raised()) return false;
  // compute split RPDG, index and non-redundant set of ce's
  pair_vec split_index;
  index_set_graph split_dg;
  if (verbose_level > 0)
    std::cerr << "computing split RPDG..." << std::endl;
  compute_split_dependency_graph(new_rp, new_rp_ce, split_index, split_dg);
  if (stats.break_signal_raised()) return false;
  index_set_graph dg2(split_dg);
  add_back_edges(split_index, dg2);
  // then call exec_rp_with_ce
  if (verbose_level > 0)
    std::cerr << "checking rp valididty..." << std::endl;
  bool_vec rem_rp(true, new_rp.size());
  lvector<bool_vec> state(bool_vec(init, ins.n_atoms()), new_rp.size() + 1);
  lvector<bool_vec> rstate(bool_vec(init, ins.n_atoms()), new_rp.size() + 1);
  seq.assign_value(EMPTYSET, new_rp.size());
  bool solved = exec_rp_with_ce(rem_rp, state, rstate, goal, new_rp, new_rp_ce,
				split_index, split_dg, dg2, with_full_dc,
				0, seq, cs);
  if (solved) {
    if (stats.break_signal_raised()) return false;
    std::cerr << "relaxed plan is valid!" << std::endl;
    return true;
  }
  if (verbose_level > 0)
    std::cerr << "rp not valid, " << cs.size() << " conflicts" << std::endl;
  filter_flaws(cs);
  if (verbose_level > 0)
    std::cerr << cs.size() << " conflicts after filtering" << std::endl;
  assert(!cs.empty());
  return false;
}

NTYPE ILB::hplusplus_with_ce(NTYPE bound, Plan** sol)
{
  stats.start();
  reset_hlb();
  assert(meta_atom_map.empty());
  std::cerr << "computing relaxed plan..." << std::endl;
  index_vec rp;
  index_set_vec rp_ce;
  NTYPE v_max =
    hplus_with_ce(ins.init_atoms, ins.goal_atoms, bound, rp, rp_ce);
  if (stats.break_signal_raised() || INFINITE(v_max) || (v_max >= bound)) {
    stats.stop();
    return v_max;
  }
  bool ok = validate_rp_with_ce(rp, rp_ce);
  assert(ok);
  index_set_vec seq;
  index_set_vec new_cs;
  ce_stats.start();
  bool solved = check_rp_with_ce(ins.init_atoms, ins.goal_atoms, rp, rp_ce,
				 false, seq, new_cs);
  if ((conflict_mod != NULL) && !solved) {
    conflict_mod->apply(*this, new_cs);
  }
  ce_stats.stop();

  while (!solved && FINITE(v_max)) {
    if (stats.break_signal_raised()) {
      stats.stop();
      return v_max;
    }

    if (verbose_level > 0) {
      std::cerr << "plan failed, conflicts: " << std::endl;
      for (index_type k = 0; k < new_cs.size(); k++) {
	std::cerr << " ";
	ins.print_atom_set(std::cerr, new_cs[k]);
	if (meta_atom_map.find_rule(new_cs[k]) == no_such_index)
	  std::cerr << " (new)";
	std::cerr << std::endl;
      }
    }

    if (!contains_new_meta_atom(new_cs)) {
      if (verbose_level > 0) {
	std::cerr << "no new conflicts found, re-trying with full dc..."
		  << std::endl;
      }
      ce_stats.start();
      solved = check_rp_with_ce(ins.init_atoms, ins.goal_atoms, rp, rp_ce,
				true, seq, new_cs);
      ce_stats.stop();
      assert(!solved);
      if (verbose_level > 0) {
	std::cerr << "new conflict set: " << std::endl;
	for (index_type k = 0; k < new_cs.size(); k++) {
	  std::cerr << " ";
	  ins.print_atom_set(std::cerr, new_cs[k]);
	  if (meta_atom_map.find_rule(new_cs[k]) == no_such_index)
	    std::cerr << " (new)";
	  std::cerr << std::endl;
	}
      }
    }

    conflict_set_size += new_cs.size();
    index_type pre_mod_n_ce = ins.n_conditional_effects();
    EvalWithMetaAtoms* mx =
      ((inc != 0) ? new EvalWithMetaAtoms(*inc, meta_atom_map) : 0);
    std::cerr << "modifying problem..." << std::endl;
    bool new_conflict_found = false;
    for (index_type k = 0; k < new_cs.size(); k++) {
      assert(new_cs[k].size() == 2);
      index_type m = meta_atom_map.find_rule(new_cs[k]);
      if (m == no_such_index) {
	if (verbose_level > 1) {
	  std::cerr << "creating meta atom ";
	  ins.write_atom_set(std::cerr, new_cs[k]);
	  std::cerr << "..." << std::endl;
	  //ins.trace_level = 1;
	}
	ins.create_meta_atom_with_ce(new_cs[k], meta_atom_map, mx);
	if (stats.break_signal_raised()) {
	  stats.stop();
	  return v_max;
	}
	new_conflict_found = true;
      }
    }
    assert(new_conflict_found);
    if (remove_dominated_conditions) {
      ins.remove_dominated_conditions(meta_atom_map);
    }
    ins.cross_reference();
    pc_actions_ratio +=
      ((ins.n_actions() + ins.n_conditional_effects()) /
       (ins.n_actions() + (double)pre_mod_n_ce));
    hpp_iterations += 1;
#ifdef HPLUSPLUS_PRINT_STATS
    print_hplusplus_stats(std::cerr);
#endif
    if (mx) delete mx;
    if (verbose_level > 2) {
      std::cerr << "meta-atom map:" << std::endl;
      for (index_type i = 0; i < meta_atom_map.size(); i++)
	std::cerr << " " << meta_atom_map[i] << std::endl;
      if (verbose_level > 4) {
	std::cerr << "new domain:" << std::endl;
	ins.print(std::cerr);
      }
    }
    reach.init_structs();
#ifdef USE_NEWLM_BC
    delete h1;
    h1 = new CostTable(ins, stats);
    h1->compute_H1(UnitACF());
#else
    if (ILA_use_lmcut) {
      delete h1;
      h1 = new CostTable(ins, stats);
    }
#endif
    std::cerr << "computing relaxed plan..." << std::endl;
    rp.clear();
    rp_ce.clear();
    NTYPE v = hplus_with_ce(ins.init_atoms, ins.goal_atoms, bound, rp, rp_ce);
#ifndef HOFFMANN
    assert((v >= v_max) || stats.break_signal_raised());
#endif
    v_max = MAX(v_max, v);
    if (stats.break_signal_raised() || (v_max >= bound)) {
      stats.stop();
      return v_max;
    }
    ok = validate_rp_with_ce(rp, rp_ce);
    assert(ok);
    new_cs.clear();
    ce_stats.start();
    solved = check_rp_with_ce(ins.init_atoms, ins.goal_atoms, rp, rp_ce,
			      false, seq, new_cs);
    if ((conflict_mod != NULL) && !solved) {
      conflict_mod->apply(*this, new_cs);
    }
    ce_stats.stop();
  }

  if (solved) {
    ActionSequence* plan = (sol != NULL ? new ActionSequence() : NULL);
    std::cerr << "plan succeeded!" << std::endl;
    for (index_type k = 0; k < seq.length(); k++)
      for (index_type i = 0; i < seq[k].size(); i++) {
	index_type a = rp[seq[k][i]];
	std::cerr << ins.actions[a].name << std::endl;
	if (plan != NULL) {
	  plan->insert(a);
	}
      }
    if (plan != NULL) {
      plan->end();
      *sol = plan;
    }
  }

  stats.stop();
  return v_max;
}

bool ILB::validate_rp(const index_set& rp)
{
  index_vec seq;
  bool ok = reach.order_relaxed_plan(ins.init_atoms, ins.goal_atoms, rp, seq);
  bool_vec state(ins.init_atoms, ins.n_atoms());
  std::cerr << "init: ";
  ins.print_atom_set(std::cerr, ins.init_atoms);
  std::cerr << std::endl;
  for (index_type k = 0; k < seq.length(); k++) {
    // std::cerr << seq[k] << " = " << ins.actions[seq[k]].name << std::endl;
    std::cerr << seq[k] << " = ";
    ins.print_action(std::cerr, ins.actions[seq[k]]);
    //std::cerr << std::endl;
    assert(state.contains(ins.actions[seq[k]].pre));
    // std::cerr << " - add: " << ins.actions[seq[k]].add << " = ";
    // ins.write_atom_set(std::cerr, ins.actions[seq[k]].add);
    // std::cerr << std::endl;
    state.insert(ins.actions[seq[k]].add);
  }
  if (ok) {
    std::cerr << "goal: ";
    ins.print_atom_set(std::cerr, ins.goal_atoms);
    std::cerr << std::endl;
    assert(seq.length() == rp.size());
    assert(state.contains(ins.goal_atoms));
  }
  return ok;
}

bool ILB::validate_rp_with_ce
(const index_vec& rp, const index_set_vec& rp_ce)
{
  bool_vec state(ins.init_atoms, ins.n_atoms());
  if (verbose_level > 2) {
    std::cerr << "init: " << ins.init_atoms << " = ";
    ins.write_atom_set(std::cerr, ins.init_atoms);
    std::cerr << std::endl;
  }
  for (index_type k = 0; k < rp.size(); k++) {
    if (verbose_level > 2) {
      std::cerr << k << ": " << ins.actions[rp[k]].name << std::endl;
      std::cerr << "pre: ";
      ins.print_atom_set(std::cerr, ins.actions[rp[k]].pre);
      std::cerr << std::endl;
    }
    for (index_type i = 0; i < rp_ce[k].size(); i++) {
      if (verbose_level > 2) {
	std::cerr << " + ";
	ins.print_rule(std::cerr, ins.actions[rp[k]].cadd[rp_ce[k][i]]);
	std::cerr << std::endl;
      }
    }
    if (!state.contains(ins.actions[rp[k]].pre)) {
      if (verbose_level > 2) {
	std::cerr << " - precondition is false" << std::endl;
	for (index_type i = 0; i < ins.actions[rp[k]].pre.size(); i++)
	  if (!state[ins.actions[rp[k]].pre[i]])
	    std::cerr << " - " << ins.atoms[ins.actions[rp[k]].pre[i]].name
		      << std::endl;
      }
      return false;
    }
    for (index_type j = 0; j < rp_ce[k].size(); j++) {
      const index_set& ra = ins.actions[rp[k]].cadd[rp_ce[k][j]].antecedent;
      if (!state.contains(ra)) {
	if (verbose_level > 2) {
	  std::cerr << " - effect condition ";
	  ins.write_atom_set(std::cerr, ra);
	  std::cerr << " is false" << std::endl;
	  for (index_type i = 0; i < ra.size(); i++)
	    if (!state[ra[i]])
	      std::cerr << " - " << ins.atoms[ra[i]].name << std::endl;
	}
	return false;
      }
    }
    if (verbose_level > 2) {
      std::cerr << " - add: " << ins.actions[rp[k]].add << " = ";
      ins.write_atom_set(std::cerr, ins.actions[rp[k]].add);
      std::cerr << std::endl;
    }
    state.insert(ins.actions[rp[k]].add);
    for (index_type j = 0; j < rp_ce[k].size(); j++) {
      index_type rc = ins.actions[rp[k]].cadd[rp_ce[k][j]].consequent;
      if (verbose_level > 2) {
	std::cerr << " - ce add: " << rc << " = " << ins.atoms[rc].name
		  << std::endl;
      }
      state[rc] = true;
    }
  }
  if (verbose_level > 2) {
    std::cerr << "goal: " << ins.goal_atoms << " = ";
    ins.write_atom_set(std::cerr, ins.goal_atoms);
    std::cerr << std::endl;
  }
  if (!state.contains(ins.goal_atoms)) {
    if (verbose_level > 2) {
      std::cerr << " - goal is false" << std::endl;
      for (index_type i = 0; i < ins.goal_atoms.size(); i++)
	if (!state[ins.goal_atoms[i]])
	  std::cerr << " - " << ins.atoms[ins.goal_atoms[i]].name
		    << std::endl;
    }
    return false;
  }
  if (verbose_level > 2) {
    std::cerr << " - goal is true" << std::endl;
  }
  return true;
}

void ILB::print_rp(std::ostream& s, const index_set& rp)
{
  for (index_type i = 0; i < rp.size(); i++)
    s << i << ": " << rp[i] << ". " << ins.actions[rp[i]].name
      << " (" << PRINT_NTYPE(cost(rp[i])) << ")"
      << std::endl;
}

void ILB::print_rpdg
(std::ostream& s, const index_set& rp, const index_set_graph& dg)
{
  dg.write_digraph(s, "RPDG");
  index_set dg_atoms;
  for (index_type i = 0; i < rp.size(); i++)
    for (index_type j = 0; j < rp.size() + 1; j++)
      if (dg.adjacent(i, j))
	dg_atoms.insert(dg.edge_label(i, j));
  for (index_type i = 0; i < dg_atoms.size(); i++)
    s << dg_atoms[i] << " = " << ins.atoms[dg_atoms[i]].name << std::endl;
}

void ILB::print_rp_with_ce
(std::ostream& s, const index_vec& rp, const index_set_vec& rp_ce)
{
  for (index_type k = 0; k < rp.size(); k++) {
    s << k << "." << ins.actions[rp[k]].name << std::endl;
    for (index_type i = 0; i < rp_ce[k].size(); i++) {
      s << " + ";
      ins.print_rule(s, ins.actions[rp[k]].cadd[rp_ce[k][i]]);
      s << std::endl;
    }
  }
  s << rp.size() << ". <goal> = ";
  ins.write_atom_set(s, ins.goal_atoms);
  s << std::endl;
}


///
// encapsulation by-pass
///

bool ILB::is_relaxed_plan(const bool_vec& set)
{
  reach.recompute(ins.init_atoms, set);
  return !reach.unreachable(ins.goal_atoms);
}

const index_set& ILB::make_new_landmark(const index_set& set)
{
  if (lm_with.size() < ins.n_actions())
    lm_with.assign_value(EMPTYSET, ins.n_actions());
  //std::cerr << "set in: " << set << std::endl;
  index_type i = make_new_landmark_ST(ins.init_atoms, ins.goal_atoms, set);
  return landmarks[i];
}

void ILB::compute_relevant_actions()
{
  compute_relevant(ins.init_atoms, ins.goal_atoms);
}

///
// Junk
///

// // construct RPDG with ce's as separate nodes...
// std::cerr << "computing split RPDG..." << std::endl;
// Instance split_ins;
// split_ins.copy_atoms(ins);
// pair_vec split_acts;
// for (index_type k = 0; k < rp.length(); k++) {
//   if (split_acts.first(index_pair(rp[k], no_such_index)) == no_such_index)
//     split_acts.append(index_pair(rp[k], no_such_index));
//   for (index_type i = 0; i < rp_ce[k].size(); i++) {
//   assert(split_acts.first(index_pair(rp[k], rp_ce[k][i])) == no_such_index);
//     split_acts.append(index_pair(rp[k], rp_ce[k][i]));
//   }
// }
// std::cerr << split_acts << std::endl;
// for (index_type k = 0; k < split_acts.length(); k++) {
//   index_type ai = split_acts[k].first;
//   assert(ai < ins.n_actions());
//   if (split_acts[k].second == no_such_index) {
//     Instance::Action& a = split_ins.new_action(ins.actions[ai].name);
//     a.pre = ins.actions[ai].pre;
//     a.add = ins.actions[ai].add;
//     a.cost = ins.actions[ai].cost;
//   }
//   else {
//     index_type ci = split_acts[k].second;
//     assert(ci < ins.actions[ai].cadd.size());
//     Instance::Action& a = split_ins.new_action(ins.actions[ai].name);
//     a.pre = ins.actions[ai].cadd[ci].antecedent;
//     a.add.assign_singleton(ins.actions[ai].cadd[ci].consequent);
//     a.cost = ZERO;
//   }
// }
// split_ins.cross_reference();
// //split_ins.print(std::cerr);
// CostACF split_ins_cost(split_ins);
// ILB ilb2(split_ins, split_ins_cost, NULL, stats);
// ilb2.rp_actions.fill(split_ins.n_actions());
// std::cerr << "computing RPDG..." << std::endl;
// ilb2.compute_dependency_graph(init, goal, ilb2.rp_dg);
// ilb2.rp_dg.transitive_reduction();
// if (verbose_level > 2) {
//   ((labeled_graph<index_set, index_set>)ilb2.rp_dg).
//     write_digraph(std::cerr, true, false, true, false, "RPDG");
//   std::cerr << "nodes:" << std::endl;
//   for (index_type k = 0; k < split_acts.size(); k++) {
//     if (split_acts[k].second == no_such_index) {
// 	std::cerr << k << " = " << ins.actions[split_acts[k].first].name
// 		  << std::endl;
//     }
//     else {
// 	std::cerr << k << " = ";
// 	ins.print_rule(std::cerr, ins.actions[split_acts[k].first].
// 		       cadd[split_acts[k].second]);
// 	std::cerr << std::endl;
//     }
//   }
//   index_set atoms;
//   for (index_type i = 0; i < ilb2.rp_dg.size(); i++)
//     for (index_type j = 0; j < ilb2.rp_dg.size(); j++)
// 	if (ilb2.rp_dg.adjacent(i,j))
// 	  if (ilb2.rp_dg.edge_has_label(i,j))
// 	    atoms.insert(ilb2.rp_dg.edge_label(i,j));
//   std::cerr << "atoms:" << std::endl;
//   for (index_type i = 0; i < atoms.size(); i++)
//   std::cerr << atoms[i] << " = " << ins.atoms[atoms[i]].name << std::endl;
// }

void ILB::test(const index_vec& rp)
{
  // index_vec rp;
  // index_set_vec rp_ce;
  // NTYPE v_max =
  //   hplus_with_ce(ins.init_atoms, ins.goal_atoms, POS_INF, rp, rp_ce);
  // std::cerr << "found relaxed plan with cost " << v_max << std::endl;
  // print_rp_with_ce(std::cerr, rp, rp_ce);

  index_set_vec rp_ce(EMPTYSET, rp.size());
  print_rp_with_ce(std::cerr, rp, rp_ce);

  // add all pairs:

  assert(meta_atom_map.empty());
  index_type n = ins.n_atoms();
  for (index_type i = 0; i < n; i++)
    for (index_type j = i + 1; j < n; j++) {
      index_set s;
      s.insert(i);
      s.insert(j);
      ins.create_meta_atom_with_ce(s, meta_atom_map, NULL);
    }
  ins.cross_reference();

  // check if rp still works:

  bool_vec rstate(ins.init_atoms, ins.n_atoms());
  for (index_type k = 0; k < rp.size(); k++) {
    std::cerr << rstate << std::endl;
    std::cerr << "step " << k << " = ";
    // << ins.actions[rp[k]].name << std::endl;
    ins.print_action(std::cerr, ins.actions[rp[k]]);
    if (rstate.contains(ins.actions[rp[k]].pre)) {
      bool_vec triggered(false, ins.actions[rp[k]].cadd.size());
      for (index_type i = 0; i < ins.actions[rp[k]].cadd.size(); i++) {
	if (rstate.contains(ins.actions[rp[k]].cadd[i].antecedent)) {
	  std::cerr << " + ";
	  ins.print_rule(std::cerr, ins.actions[rp[k]].cadd[i]);
	  std::cerr << " triggered" << std::endl;
	  triggered[i] = true;
	}
	else {
	  std::cerr << " + ";
	  ins.print_rule(std::cerr, ins.actions[rp[k]].cadd[i]);
	  std::cerr << " NOT triggered" << std::endl;
	}
      }
      rstate.insert(ins.actions[rp[k]].add);
      for (index_type i = 0; i < ins.actions[rp[k]].cadd.size(); i++)
	if (triggered[i])
	  rstate[ins.actions[rp[k]].cadd[i].consequent] = true;
    }
    else {
      std::cerr << " - precondition not true!" << std::endl;
      return;
    }
  }
  std::cerr << rstate << std::endl;
  if (rstate.contains(ins.goal_atoms)) {
    std::cerr << "goal achieved!" << std::endl;
  }
  else {
    std::cerr << " - goal = " << ins.goal_atoms << " not true!" << std::endl;
  }
}

///
/// things that didn't work (or, in some cases, that did work but
/// are now obsolete):
///

// NTYPE ILB::hitting_set_rc(index_set& hs)
// {
//   weighted_graph g(landmarks.size());
//   for (index_type i = 0; i < landmarks.size(); i++)
//     for (index_type j = i + 1; j < landmarks.size(); j++)
//       if (landmarks[i].have_common_element(landmarks[j]))
// 	g.add_undirected_edge(i, j);
//   for (index_type i = 0; i < landmarks.size(); i++)
//     g.set_weight(i, ilog(landmarks[i].size()));
//   //index_set_graph dtree;
//   //g.recursive_tree_decomposition(1, 0, dtree);
//   //dtree.write_DOT(std::cerr, false, "dtree");
//   index_set s;
//   index_set_vec p;
//   g.min_vs(s, p, 1, 0);
//   std::cerr << "s = " << s << ", |p| = " << p.size() << std::endl;
//   assert(0);
// }

// void ILB::node_conflict_to_atom_conflicts
// (const node_conflict& nc, index_set_vec& acs)
// {
//   index_type deler = nc.first.first;
//   index_type failed = nc.first.second;
//   index_type p = nc.second;
//   if (rp_dg.reachable(deler, failed)) {
//     index_set cc;
//     choose_atoms_on_path(deler, failed, cc);
//     for (index_type j = 0; j < cc.size(); j++) {
//       index_set c;
// #ifdef VARIANT_B
//       meta_atom_to_atom_set(p, meta_atom_map, c);
//       meta_atom_to_atom_set(cc[j], meta_atom_map, c);
// #else
//       c.insert(p);
//       c.insert(cc[j]);
// #endif
//       acs.append_if_new(c);
//     }
//   }
//   else {
//     index_type ncd;
//     estimate_conflict_weight(index_pair(deler, failed), ncd);
//     assert((ncd != no_such_index) && (ncd < (rp_actions.size() + 1)));
//     // std::cerr << "deler = " << deler << ", failed = " << failed
//     //	<< ", ncd = " << ncd << std::endl;
//     index_set cc1;
//     choose_atoms_on_path(deler, ncd, cc1);
//     // std::cerr << "atom set 1 = " << cc1 << std::endl;
//     index_set cc2;
//     choose_atoms_on_path(failed, ncd, cc2);
//     // std::cerr << "atom set 2 = " << cc2 << std::endl;
//     cc2.insert(p);
//     assert(!cc1.have_common_element(cc2));
//     for (index_type j1 = 0; j1 < cc1.size(); j1++)
//       for (index_type j2 = 0; j2 < cc2.size(); j2++) {
// 	index_set c;
// #ifdef VARIANT_B
// 	meta_atom_to_atom_set(cc1[j1], meta_atom_map, c);
// 	meta_atom_to_atom_set(cc2[j2], meta_atom_map, c);
// #else
// 	c.insert(cc1[j1]);
// 	c.insert(cc2[j2]);
// #endif
// 	acs.append_if_new(c);
//       }
//   }
// }

// void ILB::conflicts
// (index_type failed,
//  lvector<bool_vec>& state,
//  lvector<bool_vec>& rstate,
//  const index_set& goal,
//  index_type s_fail,
//  index_set_vec& seq,
//  pair_set& cs)
// {
//   const index_set& c_fail =
//     (failed < rp_actions.size() ? ins.actions[rp_actions[failed]].pre : goal);
//   if (verbose_level > 1) {
//     std::cerr << "extracting conflicts: " << failed << ", c = ";
//     ins.write_atom_set(std::cerr, c_fail);
//     std::cerr << ", failed at step " << s_fail << std::endl;
//   }
//   assert(s_fail > 0);
//   // for every atom (p) in the failed condition that is not true in
//   // the current state (at step s_fail)...
//   for (index_type i = 0; i < c_fail.size(); i++) {
//     index_type p = c_fail[i];
//     if (!state[s_fail][p]) {
//       if (verbose_level > 2) {
// 	std::cerr << "p = " << ins.atoms[p].name << std::endl;
//       }
//       // find the last step where the atom was true (s_true)
//       index_type s_true = s_fail;
//       while ((s_true > 0) && !state[s_true - 1][p])
// 	s_true -= 1;
//       assert(s_true > 0);
//       s_true -= 1;
//       assert(s_true < s_fail);
//       assert(state[s_true][p]);
//       if (verbose_level > 2) {
// 	std::cerr << "s_true = " << s_true << std::endl;
//       }
//       // rem stores "future" relaxed plan actions
//       bool_vec rem(true, rp_actions.size());
//       for (index_type s = 0; s < s_true; s++)
// 	for (index_type k = 0; k < seq[s].size(); k++)
// 	  rem[seq[s][k]] = false;
//       // for every action taking place between s_true and current step...
//       for (index_type s = s_true; s < s_fail; s++) {
// 	if (verbose_level > 2) {
// 	  std::cerr << "step " << s << ":";
// 	  for (index_type k = 0; k < seq[s].size(); k++)
// 	    std::cerr << " " << ins.actions[rp_actions[seq[s][k]]].name;
// 	  std::cerr << std::endl;
// 	}
// 	for (index_type k = 0; k < seq[s].size(); k++)
// 	  rem[seq[s][k]] = false;
// 	for (index_type k = 0; k < seq[s].size(); k++) {
// 	  index_type a = rp_actions[seq[s][k]];
// 	  // if this action deletes p...
// 	  if (ins.actions[a].del.contains(p)) {
// 	    if (verbose_level > 2) {
// 	      std::cerr << seq[s][k] << " = " << ins.actions[a].name
// 			<< " deletes p" << std::endl;
// 	    }
// 	    if (rp_dg.reachable(seq[s][k], failed)) {
// 	      if (verbose_level > 2) {
// 		std::cerr << "case 1: " << seq[s][k]
// 			  << " is a predecessor of " << failed << std::endl;
// 	      }
// 	      index_set cc;
// 	      index_type l =
// 	      	rp_dg.union_of_edges_on_path(seq[s][k], failed, cc);
// 	      assert(l != no_such_index);
// 	      if (verbose_level > 2) {
// 	      	std::cerr << "chain conditions: ";
// 	      	ins.write_atom_set(std::cerr, cc);
// 	      	std::cerr << std::endl;
// 	      }
// 	      for (index_type j = 0; j < cc.size(); j++) {
// 	      	index_pair c(p, cc[j]);
// 	      	c.sort_ascending();
// 	      	cs.insert(c);
// 	      }
// 	      return;
// 	    }
// 	    else {
// 	      bool_vec ncd_set;
// 	      rp_dg.nearest_common_descendants(seq[s][k], failed, ncd_set);
// 	      index_type ncd = ncd_set.first(true);
// 	      assert(ncd < (rp_actions.size() + 1));
// 	      if (verbose_level > 2) {
// 		std::cerr << "case 2: " << seq[s][k]
// 			  << " and " << failed
// 			  << " have ncds " << index_set(ncd_set)
// 			  << std::endl;
// 		std::cerr << "chose ncd " << ncd << " = ";
// 		if (ncd < rp_actions.size())
// 		  std::cerr << ins.actions[rp_actions[ncd]].name << std::endl;
// 		else
// 		  std::cerr << "goal" << std::endl;
// 	      }
// 	      index_set cc1;
// 	      rp_dg.union_of_edges_on_path(seq[s][k], ncd, cc1);
// 	      index_set cc2;
// 	      rp_dg.union_of_edges_on_path(failed, ncd, cc2);
// 	      cc2.insert(p);
// 	      if (verbose_level > 2) {
// 		std::cerr << "conditions on chain 1: ";
// 		ins.write_atom_set(std::cerr, cc1);
// 		std::cerr << std::endl << "conditions on chain 2: ";
// 		ins.write_atom_set(std::cerr, cc2);
// 		std::cerr << std::endl;
// 	      }
// 	      assert(!cc1.have_common_element(cc2));
// 	      for (index_type j1 = 0; j1 < cc1.size(); j1++)
// 		for (index_type j2 = 0; j2 < cc2.size(); j2++) {
// 		  index_pair c(cc1[j1], cc2[j2]);
// 		  c.sort_ascending();
// 		  cs.insert(c);
// 		}
// 	      return;
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
//   assert(0);
// }
// 
// void ILB::conflicts_simple_min
// (index_type failed,
//  lvector<bool_vec>& state,
//  lvector<bool_vec>& rstate,
//  const index_set& goal,
//  index_type s_fail,
//  index_set_vec& seq,
//  pair_set& cs)
// {
//   const index_set& c_fail =
//     (failed < rp_actions.size() ? ins.actions[rp_actions[failed]].pre : goal);
//   if (verbose_level > 1) {
//     std::cerr << "extracting conflicts: " << failed << ", c = ";
//     ins.write_atom_set(std::cerr, c_fail);
//     std::cerr << ", failed at step " << s_fail << std::endl;
//   }
//   assert(s_fail > 0);
//   bool conflict_found = false;
//   pair_set best_cs;
//   index_type best_cs_size = 0;
//   // for every atom (p) in the failed condition that is not true in
//   // the current state (at step s_fail)...
//   for (index_type i = 0; i < c_fail.size(); i++) {
//     index_type p = c_fail[i];
//     if (!state[s_fail][p]) {
//       if (verbose_level > 2) {
// 	std::cerr << "p = " << ins.atoms[p].name << std::endl;
//       }
//       // find the last step where the atom was true (s_true)
//       index_type s_true = s_fail;
//       while ((s_true > 0) && !state[s_true - 1][p])
// 	s_true -= 1;
//       assert(s_true > 0);
//       s_true -= 1;
//       assert(s_true < s_fail);
//       assert(state[s_true][p]);
//       if (verbose_level > 2) {
// 	std::cerr << "s_true = " << s_true << std::endl;
//       }
//       // rem stores "future" relaxed plan actions
//       bool_vec rem(true, rp_actions.size());
//       for (index_type s = 0; s < s_true; s++)
// 	for (index_type k = 0; k < seq[s].size(); k++)
// 	  rem[seq[s][k]] = false;
//       // for every action taking place between s_true and current step...
//       for (index_type s = s_true; s < s_fail; s++) {
// 	if (verbose_level > 2) {
// 	  std::cerr << "step " << s << ":";
// 	  for (index_type k = 0; k < seq[s].size(); k++)
// 	    std::cerr << " " << ins.actions[rp_actions[seq[s][k]]].name;
// 	  std::cerr << std::endl;
// 	}
// 	for (index_type k = 0; k < seq[s].size(); k++)
// 	  rem[seq[s][k]] = false;
// 	for (index_type k = 0; k < seq[s].size(); k++) {
// 	  index_type a = rp_actions[seq[s][k]];
// 	  // if this action deletes p...
// 	  if (ins.actions[a].del.contains(p)) {
// 	    if (verbose_level > 2) {
// 	      std::cerr << seq[s][k] << " = " << ins.actions[a].name
// 			<< " deletes p" << std::endl;
// 	    }
// 	    if (rp_dg.reachable(seq[s][k], failed)) {
// 	      if (verbose_level > 2) {
// 		std::cerr << "case 1: " << seq[s][k]
// 			  << " is a predecessor of " << failed << std::endl;
// 	      }
// 	      index_set cc;
// 	      index_type l =
// 	      	rp_dg.union_of_edges_on_path(seq[s][k], failed, cc);
// 	      assert(l != no_such_index);
// 	      if (verbose_level > 2) {
// 	      	std::cerr << "chain conditions: ";
// 	      	ins.write_atom_set(std::cerr, cc);
// 	      	std::cerr << std::endl;
// 	      }
// 	      pair_set new_cs;
// 	      for (index_type j = 0; j < cc.size(); j++) {
// 	      	index_pair c(p, cc[j]);
// 	      	c.sort_ascending();
// 	      	new_cs.insert(c);
// 	      }
// 	      if (!conflict_found) {
// 		best_cs = new_cs;
// 		best_cs_size = new_cs.size() - new_cs.count_common(cs);
// 		conflict_found = true;
// 	      }
// 	      else if ((new_cs.size() - new_cs.count_common(cs))
// 		       < best_cs_size) {
// 		best_cs = new_cs;
// 		best_cs_size = new_cs.size() - new_cs.count_common(cs);
// 	      }
// 	    }
// 	    else {
// 	      bool_vec ncd_set;
// 	      rp_dg.nearest_common_descendants(seq[s][k], failed, ncd_set);
// 	      index_type ncd = ncd_set.first(true);
// 	      assert(ncd < (rp_actions.size() + 1));
// 	      if (verbose_level > 2) {
// 		std::cerr << "case 2: " << seq[s][k]
// 			  << " and " << failed
// 			  << " have ncds " << index_set(ncd_set)
// 			  << std::endl;
// 		std::cerr << "chose ncd " << ncd << " = ";
// 		if (ncd < rp_actions.size())
// 		  std::cerr << ins.actions[rp_actions[ncd]].name << std::endl;
// 		else
// 		  std::cerr << "goal" << std::endl;
// 	      }
// 	      index_set cc1;
// 	      rp_dg.union_of_edges_on_path(seq[s][k], ncd, cc1);
// 	      index_set cc2;
// 	      rp_dg.union_of_edges_on_path(failed, ncd, cc2);
// 	      cc2.insert(p);
// 	      if (verbose_level > 2) {
// 		std::cerr << "conditions on chain 1: ";
// 		ins.write_atom_set(std::cerr, cc1);
// 		std::cerr << std::endl << "conditions on chain 2: ";
// 		ins.write_atom_set(std::cerr, cc2);
// 		std::cerr << std::endl;
// 	      }
// 	      assert(!cc1.have_common_element(cc2));
// 	      pair_set new_cs;
// 	      for (index_type j1 = 0; j1 < cc1.size(); j1++)
// 		for (index_type j2 = 0; j2 < cc2.size(); j2++) {
// 		  index_pair c(cc1[j1], cc2[j2]);
// 		  c.sort_ascending();
// 		  new_cs.insert(c);
// 		}
// 	      if (!conflict_found) {
// 		best_cs = new_cs;
// 		best_cs_size = new_cs.size() - new_cs.count_common(cs);
// 		conflict_found = true;
// 	      }
// 	      else if ((new_cs.size() - new_cs.count_common(cs))
// 		       < best_cs_size) {
// 		best_cs = new_cs;
// 		best_cs_size = new_cs.size() - new_cs.count_common(cs);
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//   }
//   assert(conflict_found);
//   cs.insert(best_cs);
// }
// 
// bool ILB::exec_rp
// (bool_vec& rem_rp,
//  lvector<bool_vec>& state,
//  lvector<bool_vec>& rstate,
//  const index_set& goal,
//  index_type step,
//  index_set_vec& seq,
//  pair_set& cs)
// {
//   //std::cerr << " -> " << seq << std::endl;
//   if (rem_rp.count(true) == 0) {
//     assert(rstate[step].contains(goal));
//     if (state[step].contains(goal)) {
//       return true;
//     }
//     else {
//       conflicts_simple_min
// 	(rp_actions.size(), state, rstate, goal, step, seq, cs);
//       return false;
//     }
//   }
//   index_set_vec pnext;
//   partition_applicable(rem_rp, rstate[step], pnext);
//   assert(pnext.length() > 0);
//   for (index_type k = 0; k < pnext.length(); k++) {
//     state[step + 1] = state[step];
//     rstate[step + 1] = rstate[step];
//     bool ok = true;
//     for (index_type i = 0; i < pnext[k].size(); i++) {
//       if (!state[step].contains(ins.actions[rp_actions[pnext[k][i]]].pre)) {
// 	conflicts_simple_min(pnext[k][i], state, rstate, goal, step, seq, cs);
// 	ok = false;
//       }
//       else {
// 	state[step + 1].subtract(ins.actions[rp_actions[pnext[k][i]]].del);
// 	state[step + 1].insert(ins.actions[rp_actions[pnext[k][i]]].add);
// 	rstate[step + 1].insert(ins.actions[rp_actions[pnext[k][i]]].add);
// 	rem_rp[pnext[k][i]] = false;
//       }
//     }
//     if (ok) {
//       seq[step] = pnext[k];
//       bool solved =
// 	exec_rp(rem_rp, state, rstate, goal, step + 1, seq, cs);
//       if (solved) return true;
//       if (stats.break_signal_raised()) return false;
//     }
//     for (index_type i = 0; i < pnext[k].size(); i++)
//       rem_rp[pnext[k][i]] = true;
//   }
//   return false;
// }
// 
// bool ILB::check_rp
// (const index_set& init,
//  const index_set& goal,
//  index_set_vec& seq,
//  pair_set& cs)
// {
//   compute_non_redundant_rp(init, goal);
//   compute_dependency_graph(init, goal, rp_dg);
//   rp_dg.transitive_reduction();
//   if (verbose_level > 2) {
//     std::cerr << "dependency graph: " << rp_dg << std::endl;
//     rp_dg.write_digraph(std::cerr, "DG");
//   }
//   bool_vec rem_rp(true, rp_actions.size());
//   lvector<bool_vec>
//     state(bool_vec(init, ins.n_atoms()), rp_actions.size() + 1);
//   lvector<bool_vec>
//     rstate(bool_vec(init, ins.n_atoms()), rp_actions.size() + 1);
//   seq.assign_value(EMPTYSET, rp_actions.size());
//   bool solved = exec_rp(rem_rp, state, rstate, goal, 0, seq, cs);
//   if (!solved) assert(!cs.empty());
//   return solved;
// }

// bool ILB::propagate(cl_vec& cl, graph& prec)
// {
//   bool done = false;
//   while (!done) {
//     done = true;
//     for (index_type k = 0; k < cl.size(); k++)
//       for (index_type i = 0; i < cl[k].size(); i++) {
// 	index_type j = 0;
// 	while (j < cl[k][i].size()) {
// 	  if (prec.adjacent(k, cl[k][i][j])) {
// 	    cl[k][i].remove(j);
// 	  }
// 	  else {
// 	    j += 1;
// 	  }
// 	}
// 	if (cl[k][i].empty()) return false;
// 	if (cl[k][i].size() == 1)
// 	  if (!prec.adjacent(cl[k][i][0], k)) {
// 	    pair_set e;
// 	    prec.add_edge_to_transitive_closure(cl[k][i][0], k, e);
// 	    done = false;
// 	  }
//       }
//   }
//   return true;
// }
// 
// void ILB::print_rp
// (std::ostream& s, const index_set& goal, const cl_vec& cl, const graph& prec)
// {
//   s << "steps:" << std::endl;
//   for (index_type k = 0; k < cl.size(); k++) {
//     s << k << ". ";
//     if (k < rp_actions.size())
//       s << ins.actions[rp_actions[k]].name;
//     else
//       s << "goal";
//     s << std::endl;
//     for (index_type i = 0; i < cl[k].size(); i++) {
//       index_type p = (k < rp_actions.size() ?
// 		      ins.actions[rp_actions[k]].pre[i] : goal[i]);
//       s << " " << ins.atoms[p].name << ": " << cl[k][i] << std::endl;
//     }
//   }
//   s << cl.size() << ". init" << std::endl;
//   s << "precedence graph:" << std::endl;
//   prec.write_edge_set(s);
//   s << std::endl;
// }
// 
// index_type ILB::threat::count_options(const cl_vec& cl, const graph& prec)
// {
//   index_type opts = 0;
//   // can the deleter be ordered after the consumer?
//   if (!prec.adjacent(deler, c_node)) opts += 1;
//   // for each possible establisher, can the threat be ordered before it?
//   for (index_type i = 0; i < cl[c_node][c_pre].size(); i++)
//     if (!prec.adjacent(cl[c_node][c_pre][i], deler)) opts += 1;
//   return opts;
// }
// 
// bool ILB::threat::safe(const cl_vec& cl, const graph& prec)
// {
//   // a threat is safe (i.e., in fact not a threat) if the deleter
//   // is ordered after the consumer, or before at least one establisher
//   if (prec.adjacent(c_node, deler)) return true;
//   for (index_type i = 0; i < cl[c_node][c_pre].size(); i++)
//     if (prec.adjacent(deler, cl[c_node][c_pre][i])) return true;
//   return false;
// }
// 
// bool resolve_threats(cl_vec& cl, graph& prec, threat_vec& threats)
// {
//   // count resolve options and check if we have any unresolvable threat
//   index_vec opts(0, threats.size());
//   for (index_type k = 0; k < threats.size(); k++) {
//     opts[k] = threats[k].count_options(cl, prec);
//     if (opts[k] == 0) {
//       std::cerr << "unresolvable threat: " << threats[k] << std::endl;
//       return false;
//     }
//   }
//   // if not, take care of the ones that have only one resolver
//   for (index_type k = 0; k < threats.size(); k++)
//     if (opts[k] == 1) {
//       // ...
//     }
// }
// 
// bool ILB::check_rp
// (const index_set& init, const index_set& goal)
// {
//   index_type rp_size = rp_actions.size() + 2;
//   index_type goal_node = rp_actions.size();
//   index_type init_node = rp_actions.size() + 1;
// 
//   // precedence graph
//   graph prec(rp_size);
//   for (index_type k = 0; k < rp_actions.size(); k++) {
//     prec.add_edge(init_node, k);
//     prec.add_edge(k, goal_node);
//   }
// 
//   // potential causal links
//   cl_vec cl(index_set_vec(), rp_size - 1);
//   for (index_type k = 0; k < rp_actions.size(); k++) {
//     index_type a = rp_actions[k];
//     cl[k].assign_value(EMPTYSET, ins.actions[a].pre.size());
//     for (index_type i = 0; i < ins.actions[a].pre.size(); i++) {
//       index_type p = ins.actions[a].pre[i];
//       for (index_type j = 0; j < rp_actions.size(); j++)
// 	if ((j != k) && ins.actions[rp_actions[j]].add.contains(p))
// 	  cl[k][i].insert(j);
//       if (init.contains(p))
// 	cl[k][i].insert(init_node);
//       assert(!cl[k][i].empty());
//     }
//   }
//   cl[goal_node].assign_value(EMPTYSET, goal.size());
//   for (index_type i = 0; i < goal.size(); i++) {
//     index_type p = goal[i];
//     for (index_type j = 0; j < rp_actions.size(); j++)
//       if (ins.actions[rp_actions[j]].add.contains(p))
// 	cl[goal_node][i].insert(j);
//     if (init.contains(p))
//       cl[goal_node][i].insert(init_node);
//     assert(!cl[goal_node][i].empty());
//   }
// 
//   // do initial propagation
//   bool ok = propagate(cl, prec);
//   assert(ok);
//   print_rp(std::cerr, goal, cl, prec);
// 
//   // find the threats
//   threat_vec threats;
//   for (index_type k = 0; k < rp_actions.size(); k++)
//     for (index_type j = 0; j < cl.size(); j++)
//       for (index_type i = 0; i < cl[j].size(); i++) {
// 	index_type p = (j < rp_actions.size() ?
// 			ins.action[rp_actions[j]].pre[i] : goal[i]);
// 	if (ins.actions[k].del.contains(p)) {
// 	  threat t(k, j, i);
// 	  if (!t.safe(cl, prec))
// 	    threats.append(t);
// 	}
//       }
// 
//   // and check if they can be resolved
//   ok = resolve_threats(cl, prec, threats);
// }

//NTYPE ILB::compute_relaxed_plan_by_search
//(const index_set& init, const index_set& goal)
//{
//  h1->compute_H1(UnitACF());
//  if (INFINITE(h1->eval(goal))) {
//    return POS_INF;
//  }
//  rp_actions.clear();
//  Statistics sstats(&stats);
//  //RP_BFS_SearchResult res;
//  //RegressionLMCut h_lmc(ins, cost, sstats);
//  //RPState root(ins, cost, h_lmc, *h1, goal);
//  //BFS search(sstats, res);
//  //search.set_trace_level(verbose_level);
//  //search.start(root);
//
//  index_vec action_choice_order;
//  index_set zero_cost_actions;
//
//  precondition_cost_order o(ins, *h1);
//  for (index_type k = 0; k < ins.n_actions(); k++) {
//    if (IS_ZERO(cost(k)))
//      zero_cost_actions.insert(k);
//    else
//      action_choice_order.insert_ordered(k, o);
//  }
//
//  for (index_type k = 0; k < action_choice_order.size(); k++) {
//    std::cerr << action_choice_order[k] << ". "
//	      << ins.actions[action_choice_order[k]].name << ", "
//	      << PRINT_NTYPE(h1->eval_precondition(ins.actions[action_choice_order[k]]))
//	      << std::endl;
//  }
//
//  ForwardLMCut h_lmc(ins, goal, cost, sstats);
//  RPState2 root(ins, cost, action_choice_order, zero_cost_actions,
//		h_lmc, init, goal);
//  ActionSequenceSet rp_set;
//  Result res(&rp_set);
//  res.set_stop_condition(Result::stop_at_first);
//  BFS search(sstats, res, 1000007);
//  search.set_trace_level(verbose_level);
//  NTYPE v = search.start(root);
//
//  assert(search.solved());
//  assert(rp_set.size() == 1);
//  rp_actions.insert(rp_set[0]);
//  rp_cost = cost.sum(rp_actions);
//  assert(rp_cost == v);
//}
//
//ILB::RPState::RPState
//(Instance& i,
// const ACF& c,
// RegressionLMCut& _hlmc,
// CostTable& _h1,
// const index_set& goal)
//  : ins(i), cost(c), hlmc(_hlmc), h1(_h1),
//    d(0), e_new(no_such_index), delta(0), prec(i.n_actions() + 1)
//{
//  for (index_type k = 0; k < goal.size(); k++)
//    if (!ins.atoms[goal[k]].init)
//      oc.insert(index_pair(goal[k], ins.n_actions()));
//}
//
//ILB::RPState::RPState(const RPState& s)
//  : ins(s.ins), cost(s.cost), hlmc(s.hlmc), h1(s.h1), d(s.d), e_new(s.e_new),
//    delta(s.delta), est(s.est), plan(s.plan), prec(s.prec), oc(s.oc)
//{
//  // done
//}
//
//ILB::RPState::~RPState()
//{
//  // done
//}
//
//Transition* ILB::RPState::transition()
//{
//  return new RPTransition(e_new);
//}
//
//NTYPE ILB::RPState::delta_cost()
//{
//  return delta;
//}
//
//NTYPE ILB::RPState::acc_cost()
//{
//  return cost.sum(plan);
//}
//
//index_type ILB::RPState::depth() const
//{
//  return d;
//}
//
//NTYPE ILB::RPState::est_cost()
//{
//  return est;
//}
//
//bool  ILB::RPState::is_final()
//{
//  return oc.empty();
//}
//
//bool  ILB::RPState::is_max()
//{
//  return false;
//}
//
//void ILB::RPState::reevaluate()
//{
//  index_set g;
//  for (index_type i = 0; i < oc.size(); i++)
//    g.insert(oc[i].first);
//  DiscountACF hcost(cost, plan, ins.n_actions());
//  est = hlmc.compute(g, hcost);
//}
//
//NTYPE ILB::RPState::expand(Search& s, NTYPE bound)
//{
//  assert(!oc.empty());
//  // choose an open condition with min h^1 estimate
//  index_type i_min = 0;
//  NTYPE c_min = h1.eval(oc[0].first);
//  for (index_type i = 1; i < oc.size(); i++)
//    if (h1.eval(oc[i].first) < c_min) {
//      i_min = i;
//      c_min = h1.eval(oc[i].first);
//    }
//  index_type goal = oc[i_min].first;
//  index_type node = oc[i_min].second;
//  oc.remove(i_min);
//  index_pair s_edge = e_new;
//  NTYPE s_delta = delta;
//  NTYPE s_est = est;
//  c_min = POS_INF;
//  // first, check if we can establish it by adding only a causal link
//  for (index_type k = 0; k < plan.size(); k++)
//    if (ins.actions[plan[k]].add.contains(goal) &&
//	!prec.adjacent(node, plan[k])) {
//      pair_set es;
//      prec.add_edge_to_transitive_closure(plan[k], node, es);
//      e_new = index_pair(plan[k], node);
//      delta = 0;
//      reevaluate();
//      NTYPE c_new = s.new_state(*this, bound);
//      c_min = MIN(c_min, c_new);
//      prec.remove_edges(es);
//      if (s.done()) {
//	oc.insert(index_pair(goal, node));
//	e_new = s_edge;
//	delta = s_delta;
//	est = s_est;
//	return c_new;
//      }
//    }
//  // otherwise, we have to insert a new action
//  for (index_type k = 0; k < ins.atoms[goal].add_by.size(); k++) {
//    index_type a = ins.atoms[goal].add_by[k];
//    if (!plan.contains(a) && (cost(a) <= bound)) {
//      plan.insert(a);
//      pair_set es;
//      prec.add_edge_to_transitive_closure(a, node, es);
//      pair_set s_oc(oc);
//      index_type i = 0;
//      while (i < oc.size()) {
//	if (prec.adjacent(a, oc[i].second) &&
//	    ins.actions[a].add.contains(oc[i].first)) {
//	  oc.remove(i);
//	}
//	else {
//	  i += 1;
//	}
//      }
//      pair_set n_oc;
//      for (i = 0; i < ins.actions[a].pre.size(); i++)
//	if (!ins.atoms[ins.actions[a].pre[i]].init)
//	  oc.insert(index_pair(ins.actions[a].pre[i], a));
//      e_new = index_pair(a, node);
//      delta = cost(a);
//      reevaluate();
//      if ((delta + est) <= bound) {
//	NTYPE c_new = s.new_state(*this, bound - delta) + delta;
//	c_min = MIN(c_min, c_new);
//      }
//      else {
//	c_min = MIN(c_min, delta + est);
//      }
//      prec.remove_edges(es);
//      plan.subtract(a);
//      oc.assign_copy(s_oc);
//      if (s.done()) {
//	oc.insert(index_pair(goal, node));
//	delta = s_delta;
//	est = s_est;
//	return c_min;
//      }
//    }
//  }
//  oc.insert(index_pair(goal, node));
//  delta = s_delta;
//  est = s_est;
//  return c_min;
//}
//
//int ILB::RPState::compare(const State& s)
//{
//  const RPState& rps = (const RPState&)s;
//  if (oc < rps.oc)
//    return -1;
//  else if (oc > rps.oc)
//    return 1;
//  else if (plan < rps.plan)
//    return -1;
//  else if (plan > rps.plan)
//    return 1;
//  else
//    return prec.compare(rps.prec);
//}
//
//index_type ILB::RPState::hash()
//{
//  index_type h = 0;
//  for (index_type i = 0; i < oc.size(); i++) {
//    h = ILB::set_hash::f(oc[i].first, h);
//    h = ILB::set_hash::f(oc[i].second, h);
//  }
//  return h + ILB::set_hash::f(plan) + prec.hash(ILB::set_hash::f);
//}
//
//State* ILB::RPState::copy()
//{
//  return new RPState(*this);
//}
//
//void ILB::RPState::write(::std::ostream& s)
//{
//  s << "oc = " << oc << std::endl;
//  s << "plan = " << plan << ", ";
//  prec.write_edge_set(s);
//  s << std::endl;
//}
//
//ILB::RPTransition::RPTransition(const index_pair& e)
//  : edge(e)
//{
//  // done
//}
//
//ILB::RPTransition::RPTransition(const Transition& t)
//  : Transition(t)
//{
//  const RPTransition& rpt = (const RPTransition&)t;
//  edge = rpt.edge;
//}
//
//int ILB::RPTransition::compare(const Transition& t)
//{
//  const RPTransition& rpt = (const RPTransition&)t;
//  if (edge < rpt.edge)
//    return -1;
//  else if (edge > rpt.edge)
//    return 1;
//  else
//    return 0;
//}
//
//void ILB::RPTransition::insert(Plan& p)
//{
//  // does nada
//}
//
//void ILB::RPTransition::insert_path(Plan& p)
//{
//  // does nada
//}
//
//void ILB::RPTransition::write(::std::ostream& s) const
//{
//  s << edge.first << "->" << edge.second;
//}
//
//ILB::RP_BFS_SearchResult::RP_BFS_SearchResult()
//  : solved(false)
//{
//  // done
//}
//
//ILB::RP_BFS_SearchResult::~RP_BFS_SearchResult()
//{
//  // done
//}
//
//void ILB::RP_BFS_SearchResult::solution(State& s, NTYPE cost)
//{
//  const RPState& rps = (const RPState&)s;
//  plan = rps.plan;
//  assert(cost == rps.cost.sum(plan));
//  solved = true;
//}
//
//void ILB::RP_BFS_SearchResult::solution(State& s, Transition* p, NTYPE cost)
//{
//  const RPState& rps = (const RPState&)s;
//  plan = rps.plan;
//  assert(cost == rps.cost.sum(plan));
//  solved = true;
//}
//
//bool ILB::RP_BFS_SearchResult::more()
//{
//  return !solved;
//}
//
//
//ILB::RPState2::RPState2
//(Instance& i,
// const ACF& c,
// const index_vec& o,
// const index_vec& z,
// ForwardLMCut& h,
// const index_set& s0,
// const index_set& g)
//  : ins(i), cost(c), action_choice_order(o), zero_cost_actions(z),
//    hlmc(h), goal(g)
//{
//  holds.assign_value(false, ins.n_atoms());
//  holds.insert(s0);
//  last_choice = no_such_index;
//  apply_fixpoint();
//  reevaluate();
//}
//
//ILB::RPState2::RPState2(const RPState2& s)
//  : ins(s.ins), cost(s.cost), action_choice_order(s.action_choice_order),
//    zero_cost_actions(s.zero_cost_actions), hlmc(s.hlmc), goal(s.goal),
//    holds(s.holds), last_choice(s.last_choice), est(s.est)
//{
//  // done
//}
//
//ILB::RPState2::~RPState2()
//{
//  // done
//}
//
//void ILB::RPState2::apply_fixpoint()
//{
//  bool done = false;
//  while (!done) {
//    done = true;
//    for (index_type k = 0; k < zero_cost_actions.size(); k++)
//      if (holds.contains(ins.actions[zero_cost_actions[k]].pre))
//	if (!holds.contains(ins.actions[zero_cost_actions[k]].add)) {
//	  holds.insert(ins.actions[zero_cost_actions[k]].add);
//	  done = false;
//	}
//  }
//}
//
//void ILB::RPState2::reevaluate()
//{
//  bool_vec allowed(true, ins.n_actions());
//  if (last_choice != no_such_index) {
//    for (index_type k = 0; k <= last_choice; k++)
//      allowed[action_choice_order[k]] = false;
//  }
//  est = hlmc.compute(holds, goal, cost, &allowed);
//}
//
//ILB::RPState2* ILB::RPState2::apply_choice(index_type k)
//{
//  assert(k < action_choice_order.length());
//  index_type a = action_choice_order[k];
//  RPState2* new_s = new RPState2(*this);
//  new_s->holds.insert(ins.actions[a].add);
//  new_s->apply_fixpoint();
//  new_s->last_choice = k;
//  reevaluate();
//  new_s->State::set_predecessor(this);
//  return new_s;
//}
//
//Transition* ILB::RPState2::transition()
//{
//  if (last_choice != no_such_index) {
//    assert(last_choice < action_choice_order.size());
//    return new SeqProgTrans(action_choice_order[last_choice], 1);
//  }
//  else {
//    return 0;
//  }
//}
//
//NTYPE ILB::RPState2::delta_cost()
//{
//  if (last_choice == no_such_index) return 0;
//  return cost(action_choice_order[last_choice]);
//}
//
//NTYPE ILB::RPState2::est_cost()
//{
//  return est;
//}
//
//bool  ILB::RPState2::is_final()
//{
//  return holds.contains(goal);
//}
//
//bool  ILB::RPState2::is_max()
//{
//  return false;
//}
//
//NTYPE ILB::RPState2::expand(Search& s, NTYPE bound)
//{
//  index_type next_choice =
//    (last_choice == no_such_index ? 0 : last_choice + 1);
//  NTYPE c_min = POS_INF;
//  for (index_type k = next_choice; k < action_choice_order.size(); k++) {
//    index_type a = action_choice_order[k];
//    if (holds.contains(ins.actions[a].pre)) {
//      RPState2* new_s = apply_choice(k);
//      if (FINITE(new_s->est_cost())) {
//	if ((new_s->delta_cost() + new_s->est_cost()) <= bound) {
//	  NTYPE c_new = s.new_state(*new_s, bound - new_s->delta_cost())
//	    + new_s->delta_cost();
//	  c_min = MIN(c_min, c_new);
//	}
//	else {
//	  c_min = MIN(c_min, new_s->delta_cost() + new_s->est_cost());
//	}
//      }
//      delete new_s;
//      if (s.done()) {
//	return c_min;
//      }
//    }
//  }
//  return c_min;
//}
//
//int ILB::RPState2::compare(const State& s)
//{
//  const RPState2& rps = (const RPState2&)s;
//  if (last_choice < rps.last_choice)
//    return -1;
//  else if (last_choice > rps.last_choice)
//    return 1;
//  else
//    return holds.compare(rps.holds);
//}
//
//index_type ILB::RPState2::hash()
//{
//  return ILB::set_hash::f(holds) + last_choice;
//}
//
//State* ILB::RPState2::copy()
//{
//  return new RPState2(*this);
//}
//
//void ILB::RPState2::write(::std::ostream& s)
//{
//  s << index_set(holds) << " / " << last_choice;
//}

NTYPE ForwardHPlus::eval(const index_set& s)
{
  return ilb.hplus(s, goal, POS_INF, 0);
}

NTYPE ForwardHPlus::eval(const bool_vec& s)
{
  index_set s1(s);
  return ilb.hplus(s1, goal, POS_INF, 0);
}

void ForwardHPlus::print_stats(std::ostream& s) const
{
  ilb.print_hplus_stats(s);
}

END_HSPS_NAMESPACE
