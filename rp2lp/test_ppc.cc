
#include "parser.h"
#include "preprocess.h"
#include "exec.h"
#include "hypergraph.h"
#include "cost_table.h"
#include "enumerators.h"
#include "pcc.h"

#include "forward.h"
#include "seq_reg.h"
#include "bfs.h"
#include "plans.h"

#include <sstream>
#include <fstream>

#include <queue>
#include <hash_set>

// #define PCC_TPTP_SYNTAX
#define DO_TEST_PCC
// #define DEBUG_VERBOSE
#define NEW_VERSION_OF_OLD_PCC

#define NEW_PCC

#define MVS_USE_CPLEX
#define MVS_USE_VARIABLE_ORDER
#define MVS_USE_GREEDY_INITIAL
#define MVS_USE_IMPROVED_BOUND

#ifndef HAVE_CPLEX
#undef MVS_USE_CPLEX
#endif

#ifdef MVS_USE_CPLEX
#include <ilcplex/ilocplex.h>
#endif

BEGIN_HSPS_NAMESPACE

/// global stats variables

Statistics stats;
count_type n_hm_built = 0;
count_type n_hm_tests = 0;
Statistics hm_test_stats(&stats);
count_type n_pcc_tests = 0;
count_type n_pcc_proven = 0;
Statistics pcc_stats(&stats);
Stopwatch pcc_prep_stats(&stats);
count_type n_cega_tests = 0;
count_type n_cega_decided = 0;
Statistics cega_stats(&stats);
Statistics plan_stats(&stats);
count_type n_planner_calls = 0;
count_type n_plans_found = 0;
Statistics hg_stats(&stats);
count_type n_max_open_sets = 0;
count_type n_improving_sets = 0;

/// ppc representation and basic utils

// struct ppc {
//   const Name* name;
//   NTYPE     weight;
//   plan_constraint_type pct;
//   index_set s_c;
//   index_set s_t;
//   PDDL_Base::Preference* src;
// 
//   ppc()
//     : name(0), weight(0), pct(pc_at_end), src(0) { };
//   ppc(Name* n, NTYPE w, plan_constraint_type t, const index_set& a_c, const index_set& a_t, PDDL_Base::Preference* s)
//     : name(n), weight(w), pct(t), s_c(a_c), s_t(a_t), src(s) { };
// 
//   ppc& operator=(const ppc& p) {
//     name = p.name;
//     weight = p.weight;
//     pct = p.pct;
//     s_c.assign_copy(p.s_c);
//     s_t.assign_copy(p.s_t);
//     src = p.src;
//   };
// };
// 
// void write_ppc(std::ostream& s, const Instance& ins, const ppc& p)
// {
//   if (p.name)
//     s << p.name << ": ";
//   switch (p.pct) {
//   case pc_at_end:
//     s << "at-end /\\";
//     ins.write_atom_set(s, p.s_c);
//     break;
//   case pc_always:
//     s << "always /\\";
//     ins.write_atom_set(s, p.s_c);
//     break;
//   case pc_always_disjunction:
//     s << "always \\/";
//     ins.write_atom_set(s, p.s_c);
//     break;
//   case pc_sometime:
//     s << "sometime /\\";
//     ins.write_atom_set(s, p.s_c);
//     break;
//   case pc_at_most_once:
//     s << "at-most-once /\\";
//     ins.write_atom_set(s, p.s_c);
//     break;
//   case pc_sometime_before:
//     s << "sometime-before /\\";
//     ins.write_atom_set(s, p.s_t);
//     s << " /\\";
//     ins.write_atom_set(s, p.s_c);
//     break;
//   case pc_sometime_after:
//     s << "sometime-after /\\";
//     ins.write_atom_set(s, p.s_t);
//     s << " /\\";
//     ins.write_atom_set(s, p.s_c);
//     break;
//   default:
//     s << "ERROR!";
//     exit(255);
//   }
//   s << " (" << p.weight << ")";
// }
// 
// void write_PDDL_constraint(std::ostream& s, const Instance& ins, const ppc& p)
// {
//   switch (p.pct) {
//   case pc_at_end:
//     s << "(at end (and";
//     for (index_type i = 0; i < p.s_c.size(); i++)
//       s << " " << ins.atoms[p.s_c[i]].name;
//     s << "))";
//     break;
//   case pc_always:
//     s << "(always (and";
//     for (index_type i = 0; i < p.s_c.size(); i++)
//       s << " " << ins.atoms[p.s_c[i]].name;
//     s << "))";
//     break;
//   case pc_always_disjunction:
//     s << "(always (or";
//     for (index_type i = 0; i < p.s_c.size(); i++)
//       s << " " << ins.atoms[p.s_c[i]].name;
//     s << "))";
//     break;
//   case pc_sometime:
//     s << "(sometime (and";
//     for (index_type i = 0; i < p.s_c.size(); i++)
//       s << " " << ins.atoms[p.s_c[i]].name;
//     s << "))";
//     break;
//   case pc_sometime_disjunction:
//     s << "(sometime (or";
//     for (index_type i = 0; i < p.s_c.size(); i++)
//       s << " " << ins.atoms[p.s_c[i]].name;
//     s << "))";
//     break;
//   case pc_at_most_once:
//     s << "(at-most-once (and";
//     for (index_type i = 0; i < p.s_c.size(); i++)
//       s << " " << ins.atoms[p.s_c[i]].name;
//     s << "))";
//     break;
//   case pc_sometime_before:
//     s << "(sometime-before (and";
//     for (index_type i = 0; i < p.s_t.size(); i++)
//       s << " " << ins.atoms[p.s_t[i]].name;
//     s << ") (and";
//     for (index_type i = 0; i < p.s_c.size(); i++)
//       s << " " << ins.atoms[p.s_c[i]].name;
//     s << "))";
//     break;
//   case pc_sometime_after:
//     s << "(sometime-after (and";
//     for (index_type i = 0; i < p.s_t.size(); i++)
//       s << " " << ins.atoms[p.s_t[i]].name;
//     s << ") (and";
//     for (index_type i = 0; i < p.s_c.size(); i++)
//       s << " " << ins.atoms[p.s_c[i]].name;
//     s << "))";
//     break;
//   default:
//     s << "ERROR!";
//     exit(255);
//   }
// }
// 
// void write_PDDL_preference(std::ostream& s, const Instance& ins, const ppc& p)
// {
//   s << "(preference ";
//   if (p.name) s << p.name << " ";
//   write_PDDL_constraint(s, ins, p);
//   s << ")";
// }
// 
// void write_ppc_as_pcc(std::ostream& s, const Instance& ins, const ppc& p)
// {
//   switch (p.pct) {
// #ifdef PCC_TPTP_SYNTAX
//   case pc_at_end:
//     assert(p.s_c.length() == 1);
//     s << "fof(goal" << p.s_c[0] << ",axiom,at_end(";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << "))." << std::endl;
//     break;
//   case pc_always:
//     assert(p.s_c.length() == 1);
//     s << "fof(ac" << p.s_c[0] << ",axiom,always(";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << "))." << std::endl;
//     break;
//   case pc_sometime:
//     assert(p.s_c.length() == 1);
//     s << "fof(ec" << p.s_c[0] << ",axiom,sometime(";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << "))." << std::endl;
//     break;
//   case pc_at_most_once:
//     assert(p.s_c.length() == 1);
//     s << "fof(amoc" << p.s_c[0] << ",axiom,at_most_once(";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << "))." << std::endl;
//     break;
//   case pc_sometime_before:
//     assert(p.s_t.length() == 1);
//     assert(p.s_c.length() == 1);
//     s << "fof(sbc" << p.s_t[0] << "_" << p.s_c[0]
//       << ",axiom,sometime_before(";
//     ins.atoms[p.s_t[0]].name->write(s, Name::NC_INSTANCE);
//     s << ",";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << "))." << std::endl;
//     break;
//   case pc_sometime_after:
//     assert(p.s_t.length() == 1);
//     assert(p.s_c.length() == 1);
//     s << "fof(sac" << p.s_t[0] << "_" << p.s_c[0]
//       << ",axiom,sometime_after(";
//     ins.atoms[p.s_t[0]].name->write(s, Name::NC_INSTANCE);
//     s << ",";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << "))." << std::endl;
//     break;
// #else
//   case pc_at_end:
//     assert(p.s_c.length() == 1);
//     s << "AtEnd(";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << ")." << std::endl;
//     break;
//   case pc_always:
//     assert(p.s_c.length() == 1);
//     s << "Always(";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << ")." << std::endl;
//     break;
//   case pc_sometime:
//     assert(p.s_c.length() == 1);
//     s << "Sometime(";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << ")." << std::endl;
//     break;
//   case pc_at_most_once:
//     assert(p.s_c.length() == 1);
//     s << "AtMostOnce(";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << ")." << std::endl;
//     break;
//   case pc_sometime_before:
//     assert(p.s_t.length() == 1);
//     assert(p.s_c.length() == 1);
//     s << "SometimeBefore(";
//     ins.atoms[p.s_t[0]].name->write(s, Name::NC_INSTANCE);
//     s << ",";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << ")." << std::endl;
//     break;
//   case pc_sometime_after:
//     assert(p.s_t.length() == 1);
//     assert(p.s_c.length() == 1);
//     s << "SometimeAfter(";
//     ins.atoms[p.s_t[0]].name->write(s, Name::NC_INSTANCE);
//     s << ",";
//     ins.atoms[p.s_c[0]].name->write(s, Name::NC_INSTANCE);
//     s << ")." << std::endl;
//     break;
// #endif
//   default:
//     s << "ERROR!";
//     exit(255);
//   }
// }
// 
// void enforce_ppc(const ppc& p, Instance& ins, index_vec& map)
// {
//   switch (p.pct) {
//   case pc_always:
//     ins.enforce_pc_always_conjunction(p.s_c, p.name, map);
//     break;
//   case pc_always_disjunction:
//     ins.enforce_pc_always_disjunction(p.s_c, p.name, map);
//     break;
//   case pc_sometime:
//     ins.enforce_pc_sometime_conjunction(p.s_c, p.name, map);
//     break;
//   case pc_at_most_once:
//     ins.enforce_pc_at_most_once_conjunction(p.s_c, p.name, map);
//     break;
//   case pc_sometime_before:
//     ins.enforce_pc_sometime_before_cc(p.s_t, p.s_c, p.name, map);
//     break;
//   default:
//     std::cerr << "ERROR!";
//     exit(255);
//   }
// }
// 
// index_type compile_ppc(const ppc& p, Instance& ins)
// {
//   switch (p.pct) {
//   case pc_always:
//     return ins.compile_pc_always_conjunction(p.s_c, p.name);
//   case pc_always_disjunction:
//     return ins.compile_pc_always_disjunction(p.s_c, p.name);
//   case pc_sometime:
//     return ins.compile_pc_sometime_conjunction(p.s_c, p.name);
//   case pc_at_most_once:
//     return ins.compile_pc_at_most_once_conjunction(p.s_c, p.name);
//   case pc_sometime_before:
//     return ins.compile_pc_sometime_before_cc(p.s_t, p.s_c, p.name);
//   default:
//     std::cerr << "ERROR!";
//     exit(255);
//   }
// }
// 
// bool test_ppc(const ppc& p, ExecTrace* t)
// {
//   switch (p.pct) {
//   case pc_always:
//     return t->test_always_conjunction(p.s_c);
//   case pc_always_disjunction:
//     return t->test_always_disjunction(p.s_c);
//   case pc_sometime:
//     return t->test_sometime_conjunction(p.s_c);
//   case pc_at_most_once:
//     return t->test_at_most_once_conjunction(p.s_c);
//   case pc_sometime_before:
//     if (p.s_c.length() != 1) return false;
//     if (p.s_t.length() != 1) return false;
//     return t->test_sometime_before(p.s_t[0], p.s_c[0]);
//   case pc_sometime_after:
//     if (p.s_c.length() != 1) return false;
//     if (p.s_t.length() != 1) return false;
//     return t->test_sometime_after(p.s_t[0], p.s_c[0]);
//   default:
//     std::cerr << "ERROR!";
//     exit(255);
//   }
// }
// 
// typedef lvector<ppc> ppc_vec;
 
NTYPE ppc_set_value(const ppc_vec& pv, const index_set& set)
{
  NTYPE v = 0;
  for (index_type k = 0; k < set.length(); k++)
    v += pv[set[k]].weight;
  return v;
}

void write_ppc_set(std::ostream& s, const ppc_vec& pv, const index_set& set)
{
  s << "{";
  for (index_type k = 0; k < set.length(); k++) {
    if (k > 0) s << ", ";
    s << pv[set[k]].name;
  }
  s << "}";
}

// void write_ppc_set_as_pcc
// (std::ostream& s, Instance& ins, const ppc_vec& pv, const index_set& set,
//  bool has_val, NTYPE set_val, NTYPE set_pen)
// {
//   s << "%% set {";
//   for (index_type i = 0; i < set.length(); i++) {
//     if (i > 0) s << ", ";
//     s << pv[set[i]].name;
//   }
//   if (has_val)
//     s << "}: value = " << PRINT_NTYPE(set_val)
//       << ", penalty = " << PRINT_NTYPE(set_pen)
//       << std::endl;
//   else
//     s << "}" << std::endl;
//   for (index_type i = 0; i < set.length(); i++)
//     write_ppc_as_pcc(s, ins, pv[set[i]]);
//   for (index_type i = 0; i < ins.goal_atoms.length(); i++) {
// #ifdef PCC_TPTP_SYNTAX
//     s << "fof(goal" << ins.goal_atoms[i] << ",axiom,at_end(";
//     ins.atoms[ins.goal_atoms[i]].name->write(s, Name::NC_INSTANCE);
//     s << "))." << std::endl;
// #else
//     s << "AtEnd(";
//     ins.atoms[ins.goal_atoms[i]].name->write(s, Name::NC_INSTANCE);
//     s << ")." << std::endl;
// #endif
//   }
// }

// bool append_ppc
// (PDDL_Base* b,
//  PDDL_Base::Preference* p,
//  PDDL_Base::Goal* g1,
//  PDDL_Base::Goal* g2,
//  Instance& ins,
//  plan_constraint_type t,
//  index_type k,
//  ppc_vec& ppcv)
// {
//   PDDL_Base::signed_atom_vec a;
//   bool g1_is_dis;
//   index_set s1;
//   if (b->goal_to_atom_vec(g1, a, g1_is_dis)) if (!g1_is_dis) {
//     Name* n = (p->name ?
// 	       (Name*)new StringName(p->name->print_name) :
// 	       (Name*)new EnumName("ppc", k));
//     NTYPE w = (b->metric ? p->value(b->metric_type, b->metric) : 1);
//     b->instantiate_atom_set(ins, a, s1);
//     for (index_type i = 0; i < s1.length(); i++)
//       ins.atoms[s1[i]].goal = true;
//     if (g2) {
//       PDDL_Base::signed_atom_vec a2;
//       bool g2_is_dis;
//       index_set s2;
//       if (b->goal_to_atom_vec(g2, a2, g2_is_dis)) if (!g2_is_dis) {
// 	b->instantiate_atom_set(ins, a2, s2);
// 	for (index_type i = 0; i < s2.length(); i++)
// 	  ins.atoms[s2[i]].goal = true;
// 	ppcv.append(ppc(n, w, t, s1, s2, p));
//       }
//       else {
// 	std::cerr << "error: can't handle goal ";
// 	p->print(std::cerr);
// 	std::cerr << std::endl;
// 	return false;
//       }
//     }
//     else {
//       ppcv.append(ppc(n, w, t, s1, EMPTYSET, p));
//     }
//     return true;
//   }
//   else {
//     std::cerr << "error: can't handle goal ";
//     p->print(std::cerr);
//     std::cerr << std::endl;
//     return false;
//   }
// }
// 
// bool extract_ppcs(PDDL_Base* b, Instance& ins, ppc_vec& ppcv)
// {
//   bool ok = true;
//   for (index_type k = 0; k < b->dom_preferences.length(); k++) {
//     if (b->dom_preferences[k]->goal->g_class == PDDL_Base::goal_always) {
//       if (!append_ppc(b, b->dom_preferences[k],
// 		      ((PDDL_Base::SimpleSequenceGoal*)
// 		       b->dom_preferences[k]->goal)->constraint,
// 		      0, ins, pc_always, k, ppcv))
// 	ok = false;
//     }
//     else if (b->dom_preferences[k]->goal->g_class == PDDL_Base::goal_sometime) {
//       if (!append_ppc(b, b->dom_preferences[k],
// 		      ((PDDL_Base::SimpleSequenceGoal*)
// 		       b->dom_preferences[k]->goal)->constraint,
// 		      0, ins, pc_sometime, k, ppcv))
// 	ok = false;
//     }
//     else if (b->dom_preferences[k]->goal->g_class
// 	     == PDDL_Base::goal_at_most_once) {
//       if (!append_ppc(b, b->dom_preferences[k],
// 		      ((PDDL_Base::SimpleSequenceGoal*)
// 		       b->dom_preferences[k]->goal)->constraint,
// 		      0, ins, pc_at_most_once, k, ppcv))
// 	ok = false;
//     }
//     else if (b->dom_preferences[k]->goal->g_class
// 	     == PDDL_Base::goal_sometime_before) {
//       if (!append_ppc(b, b->dom_preferences[k],
// 		      ((PDDL_Base::TriggeredSequenceGoal*)
// 		       b->dom_preferences[k]->goal)->constraint,
// 		      ((PDDL_Base::TriggeredSequenceGoal*)
// 		       b->dom_preferences[k]->goal)->trigger,
// 		      ins, pc_sometime_before, k, ppcv))
// 	ok = false;
//     }
//     else {
//       if (!append_ppc(b, b->dom_preferences[k], b->dom_preferences[k]->goal,
// 		      0, ins, pc_at_end, k, ppcv))
// 	ok = false;
//     }
//   }
//   return ok;
// }

/// PCC implementations

class BasicPCC {
  // subconstructors
  void init_premises(ppc_vec& v, const index_set& s);

protected:
  Instance* ins;
  Statistics& stats;
  StaticMutex* mx;
  set_edge_vec* trlm;
  bool_vec rem_trlm;
  bool print_proof;

  // input constraints
  index_set c_always;
  index_set c_sometime;
  index_set c_at_most_once;
  pair_set  c_sometime_before;

  // derived constraints
  bool_vec  allowed_actions;
  bool      aa_changed;
  bool_vec  possible_atoms;
  bool      pa_changed;
  bool_vec  active_sometime;
  bool      as_changed;

  // note: sometime_before holds the inverse of the sometime-before relation,
  // i.e. p -> q iff sometime-before(q, p); p is the one that comes "before",
  // and q is the trigger.
  graph     sometime_before;
  bool      sb_changed;
  const graph* never_after;

  void init_derived(const graph* lmg, const graph* nag);
  void update();

  void assert_always(index_type i);
  void assert_never(index_type i);
  void assert_sometime(index_type i);
  // note: assert follows the PDDL3 convention, i.e., sb(i, j)
  void assert_sometime_before(index_type i, index_type j);
  void infer_never_from_sb_cycle();
  void infer_never_from_na_and_sb();
  //void infer_na_from_sb();
  //bool check_asb_na();
  bool check_sometime_never();
  bool check_as_na_cycle();
  bool check_amo_count(const index_set& amo_atoms);

public:
  BasicPCC(Instance& i, Statistics& st, ppc_vec& v, const index_set& s,
	   bool opt_print_proof);
  BasicPCC(Instance& i, Statistics& st, ppc_vec& v, const index_set& s,
	   StaticMutex* m, const graph* lmg, set_edge_vec* t, const graph* nag,
	   bool opt_print_proof);
  ~BasicPCC();
  void print_premises();
  void print_graph(const graph* lmg, std::ostream& to);
  bool run();
};

void BasicPCC::init_premises(ppc_vec& v, const index_set& s)
{
  for (index_type k = 0; k < s.length(); k++) {
    if (v[s[k]].pct == pc_always) {
      assert(v[s[k]].s_c.length() == 1);
      c_always.insert(v[s[k]].s_c);
    }
    else if (v[s[k]].pct == pc_sometime) {
      assert(v[s[k]].s_c.length() == 1);
      c_sometime.insert(v[s[k]].s_c);
    }
    else if (v[s[k]].pct == pc_at_most_once) {
      assert(v[s[k]].s_c.length() == 1);
      c_at_most_once.insert(v[s[k]].s_c);
    }
    else if (v[s[k]].pct == pc_sometime_before) {
      assert(v[s[k]].s_t.length() == 1);
      assert(v[s[k]].s_c.length() == 1);
      c_sometime_before.insert(index_pair(v[s[k]].s_c[0],v[s[k]].s_t[0]));
    }
  }
}

void BasicPCC::init_derived(const graph* lmg, const graph* nag)
{
  allowed_actions.assign_value(true, ins->n_actions());
  possible_atoms.assign_value(true, ins->n_atoms());
  active_sometime.assign_value(false, ins->n_atoms());
  if (lmg)
    sometime_before.copy(*lmg);
  else
    sometime_before.init(ins->n_atoms());
  never_after = nag;
  aa_changed = false;
  pa_changed = false;
  as_changed = false;
  sb_changed = false;
}

void BasicPCC::print_premises()
{
  std::cout << "PREMISES:" << std::endl;
  for (index_type k = 0; k < c_always.size(); k++)
    std::cout << "Always(" << ins->atoms[c_always[k]].name << ")"
	      << std::endl;
  for (index_type k = 0; k < c_sometime.size(); k++)
    std::cout << "Sometime(" << ins->atoms[c_sometime[k]].name << ")"
	      << std::endl;
  for (index_type k = 0; k < c_at_most_once.size(); k++)
    std::cout << "AtMostOnce(" << ins->atoms[c_at_most_once[k]].name << ")"
	      << std::endl;
  for (index_type k = 0; k < c_sometime_before.size(); k++)
    std::cout << "SometimeBefore("
	      << ins->atoms[c_sometime_before[k].second].name << ","
	      << ins->atoms[c_sometime_before[k].first].name << ")"
	      << std::endl;
}

void BasicPCC::print_graph(const graph* lmg, std::ostream& to)
{
  to << "digraph PROOF {" << std::endl;
  for (index_type i = 0; i < ins->n_atoms(); i++) {
    const char* shape = "circle";
    const char* color = "white";
    int periph = 1;
    if (c_always.contains(i)) {
      if (active_sometime[i])
	shape = "octagon";
      else
	shape = "box";
    }
    else if (active_sometime[i])
      shape = "diamond";
    if (c_sometime.contains(i) || ins->goal_atoms.contains(i))
      periph = 2;
    if (!possible_atoms[i])
      color = "red";
    to << i << "[shape=" << shape << ",peripheries=" << periph
       << ",style=filled,fillcolor=" << color << "];" << std::endl;
  }
  for (index_type i = 0; i < ins->n_atoms(); i++)
    for (index_type j = 0; j < ins->n_atoms(); j++) {
      if (sometime_before.adjacent(i, j)) {
	to << i << " -> " << j << "[";
	bool first = true;
	if (c_sometime_before.contains(index_pair(i, j))) {
	  to << "color=green";
	  first = false;
	}
	else if (lmg) {
	  if (lmg->adjacent(i, j)) {
	    to << "color=blue";
	    first = false;
	  }
	}
	if (never_after) {
	  if (never_after->adjacent(i, j)) {
	    if (!first) to << ",";
	    to << "style=dashed";
	  }
	}
	to << "];" << std::endl;
      }
      else if (never_after) {
	if (never_after->adjacent(i, j)) {
	  to << i << " -> " << j << "[style=dashed];" << std::endl;
	}
      }
    }
  to << "}" << std::endl;
}

BasicPCC::BasicPCC
(Instance& i, Statistics& st, ppc_vec& v, const index_set& s,
 bool opt_print_proof)
  : ins(&i), stats(st), mx(0), trlm(0), print_proof(opt_print_proof)
{
  init_premises(v, s);
  init_derived(0, 0);
}

BasicPCC::BasicPCC
(Instance& i, Statistics& st, ppc_vec& v, const index_set& s,
 StaticMutex* m, const graph* lmg, set_edge_vec* t, const graph* nag,
 bool opt_print_proof)
  : ins(&i), stats(st), mx(m), trlm(t), print_proof(opt_print_proof)
{
  init_premises(v, s);
  init_derived(lmg, nag);
  if (trlm) {
    rem_trlm.assign_value(true, trlm->size());
  }
}

BasicPCC::~BasicPCC()
{
  // nothing to do
}

void BasicPCC::update()
{
  if (!trlm) return;
  while (aa_changed) {
    aa_changed = false;
    for (index_type k = 0; k < trlm->size(); k++)
      if (rem_trlm[k]) {
	if (!(*trlm)[k].first.have_common_element(allowed_actions)) {
	  if (print_proof) {
	    std::cout << "Triggered landmark: (";
	    for (index_type i = 0; i < (*trlm)[k].first.size(); i++) {
	      if (i > 0) std::cout << " & ";
	      std::cout << "Disallowed("
			<< ins->actions[(*trlm)[k].first[i]].name << ")";
	    }
	    std::cout << ") => SometimeBefore("
		      << ins->atoms[(*trlm)[k].second.second].name
		      << ", " << ins->atoms[(*trlm)[k].second.first].name
		      << ")" << std::endl;
	  }
	  assert_sometime_before((*trlm)[k].second.second,
				 (*trlm)[k].second.first);
	  rem_trlm[k] = false;
	}
      }
    for (index_type k = 0; k < ins->n_atoms(); k++)
      if (!ins->atoms[k].init && possible_atoms[k])
	if (!(ins->atoms[k].add_by.have_common_element(allowed_actions))) {
	  if (print_proof) {
	    std::cout << "Triggered unreachable: (";
	    for (index_type i = 0; i < ins->atoms[k].add_by.size(); i++) {
	      if (i > 0) std::cout << " & ";
	      std::cout << "Disallowed("
			<< ins->actions[ins->atoms[k].add_by[i]].name << ")";
	    }
	    std::cout << ") => Never(" << ins->atoms[k].name << ")"
		      << std::endl;
	  }
	  assert_never(k);
	}
  }
}

void BasicPCC::assert_always(index_type i)
{
  for (index_type k = 0; k < ins->atoms[i].del_by.length(); k++)
    if (allowed_actions[ins->atoms[i].del_by[k]]) {
      if (print_proof) {
	std::cout << "Always(" << ins->atoms[i].name << ") => Disallowed("
		  << ins->actions[ins->atoms[i].del_by[k]].name << ")"
		  << std::endl;
      }
      allowed_actions[ins->atoms[i].del_by[k]] = false;
      aa_changed = true;
    }
#ifdef NEW_VERSION_OF_OLD_PCC
  if (mx) {
    for (index_type j = 0; j < ins->n_atoms(); j++)
      if (mx->mutex(i, j)) {
	if (print_proof) {
	  std::cout << "Always(" << ins->atoms[i].name << ") => Never("
		    << ins->atoms[j].name << ")" << std::endl;
	}
	assert_never(j);
      }
  }
#endif
}

void BasicPCC::assert_never(index_type i)
{
  if (possible_atoms[i]) {
    possible_atoms[i] = false;
    pa_changed = true;
    for (index_type k = 0; k < ins->atoms[i].add_by.length(); k++)
      if (allowed_actions[ins->atoms[i].add_by[k]]) {
	if (print_proof) {
	  std::cout << "Never(" << ins->atoms[i].name << ") => Disallowed("
		    << ins->actions[ins->atoms[i].add_by[k]].name << ")"
		    << std::endl;
	}
	allowed_actions[ins->atoms[i].add_by[k]] = false;
	aa_changed = true;
      }
    for (index_type k = 0; k < ins->atoms[i].req_by.length(); k++)
      if (allowed_actions[ins->atoms[i].req_by[k]]) {
	if (print_proof) {
	  std::cout << "Never(" << ins->atoms[i].name << ") => Disallowed("
		    << ins->actions[ins->atoms[i].req_by[k]].name << ")"
		    << std::endl;
	}
	allowed_actions[ins->atoms[i].req_by[k]] = false;
	aa_changed = true;
      }
    const index_set& a = sometime_before.successors(i);
    for (index_type k = 0; k < a.length(); k++)
      if (possible_atoms[a[k]]) {
	if (print_proof) {
	  std::cout << "Never(" << ins->atoms[i].name
		    << ") & SometimeBefore(" << ins->atoms[a[k]].name
		    << "," << ins->atoms[i].name << ") => Never("
		    << ins->atoms[a[k]].name << ")" << std::endl;
	}
	assert_never(a[k]);
      }
  }
}

void BasicPCC::assert_sometime(index_type i)
{
  if (!active_sometime[i]) {
    active_sometime[i] = true;
    as_changed = true;
    const index_set& b = sometime_before.predecessors(i);
    for (index_type k = 0; k < b.length(); k++)
      if (!active_sometime[b[k]]) {
	if (print_proof) {
	  std::cout << "Sometime(" << ins->atoms[i].name
		    << ") & SometimeBefore(" << ins->atoms[i].name
		    << "," << ins->atoms[b[k]].name << ") => Sometime("
		    << ins->atoms[b[k]].name << ")" << std::endl;
	}
	assert_sometime(b[k]);
      }
  }
}

void BasicPCC::assert_sometime_before(index_type i, index_type j)
{
  if (!sometime_before.adjacent(j, i)) {
    pair_set ie;
    sometime_before.add_edge_to_transitive_closure(j, i, ie);
    for (index_type k = 0; k < ie.size(); k++) {
      if (active_sometime[ie[k].second] && !active_sometime[ie[k].first]) {
	if (print_proof) {
	  std::cout << "Sometime(" << ins->atoms[ie[k].second].name
		    << ") & SometimeBefore(" << ins->atoms[ie[k].second].name
		    << "," << ins->atoms[ie[k].first].name << ") => Sometime("
		    << ins->atoms[ie[k].first].name << ")" << std::endl;
	}
	assert_sometime(ie[k].first);
      }
      if (!possible_atoms[ie[k].first] && possible_atoms[ie[k].second]) {
	if (print_proof) {
	  std::cout << "Never(" << ins->atoms[ie[k].first].name
		    << ") & SometimeBefore(" << ins->atoms[ie[k].second].name
		    << "," << ins->atoms[ie[k].first].name << ") => Never("
		    << ins->atoms[ie[k].second].name << ")" << std::endl;
	}
	assert_never(ie[k].second);
      }
    }
    sb_changed = true;
  }
}

void BasicPCC::infer_never_from_sb_cycle()
{
  for (index_type i = 0; i < ins->n_atoms(); i++)
    for (index_type j = i + 1; j < ins->n_atoms(); j++)
      if (sometime_before.adjacent(i, j) && sometime_before.adjacent(j, i)) {
	if (possible_atoms[i] || possible_atoms[j]) {
	  if (print_proof) {
	    std::cout << "SometimeBefore(" << ins->atoms[i].name << ","
		      << ins->atoms[j].name << ") & SometimeBefore("
		      << ins->atoms[j].name << "," << ins->atoms[i].name
		      << ") => Never(" << ins->atoms[i].name << ") & Never("
		      << ins->atoms[j].name << ")" << std::endl;
	  }
	  assert_never(i);
	  assert_never(j);
	}
      }
}

void BasicPCC::infer_never_from_na_and_sb()
{
  for (index_type i = 0; i < ins->n_atoms(); i++)
    for (index_type j = i + 1; j < ins->n_atoms(); j++)
      if (sometime_before.adjacent(i, j) && never_after->adjacent(i, j)) {
	if (possible_atoms[i] || possible_atoms[j]) {
	  if (print_proof) {
	    std::cout << "SometimeBefore(" << ins->atoms[j].name << ","
		      << ins->atoms[i].name << ") & NeverAfter("
		      << ins->atoms[j].name << "," << ins->atoms[i].name
		      << ") => Never(" << ins->atoms[i].name << ")"
		      << std::endl;
	  }
	  assert_never(i);
	}
      }
}

// void BasicPCC::infer_na_from_sb()
// {
//   for (index_type i = 0; i < ins->n_atoms(); i++)
//     for (index_type j = 0; j < ins->n_atoms(); j++)
//       if (never_after.adjacent(i, j)) {
// 	const index_set& a = sometime_before.successors(i);
// 	for (index_type k = 0; k < a.size(); k++)
// 	  if (!never_after.adjacent(a[k], j)) {
// 	    if (print_proof) {
// 	      std::cout << "NeverAfter(" << ins->atoms[i].name << ","
// 			<< ins->atoms[j].name << ") & SometimeBefore("
// 			<< ins->atoms[i].name << "," << ins->atoms[a[k]].name
// 			<< ") => NeverAfter(" << ins->atoms[a[k]].name
// 			<< "," << ins->atoms[j].name << ")" << std::endl;
// 	    }
// 	    never_after.add_edge(a[k], j);
// 	  }
//       }
// }

bool BasicPCC::check_sometime_never()
{
  for (index_type k = 0; k < ins->n_atoms(); k++)
    if (active_sometime[k] && !possible_atoms[k]) {
      if (print_proof) {
	std::cout << "Sometime(" << ins->atoms[k].name << ") & Never("
		  << ins->atoms[k].name << ") => CONTRADICTION"
		  << std::endl;
      }
      return true;
    }
  return false;
}

// bool BasicPCC::check_asb_na()
// {
//   for (index_type i = 0; i < ins->n_atoms(); i++)
//     if (active_sometime[i]) {
//       for (index_type j = 0; j < ins->n_atoms(); j++)
// 	if (never_after.adjacent(j, i) && sometime_before.adjacent(j, i)) {
// 	  if (print_proof) {
// 	    std::cout << "Sometime(" << ins->atoms[i].name
// 		      << ") & SometimeBefore(" << ins->atoms[i].name
// 		      << ", " << ins->atoms[j].name
// 		      << ") & NeverAfter(" << ins->atoms[j].name
// 		      << ", " << ins->atoms[i].name
// 		      << ") => CONTRADICTION"
// 		      << std::endl;
// 	  }
// 	  return true;
// 	}
//     }
//   return false;
// }

bool BasicPCC::check_as_na_cycle()
{
  for (index_type i = 0; i < ins->n_atoms(); i++)
    if (active_sometime[i])
      for (index_type j = i + 1; j < ins->n_atoms(); j++)
	if (active_sometime[j] &&
	    never_after->adjacent(i, j) && never_after->adjacent(j, i)) {
	  index_type a = 0;
	  index_type b = 0;
	  bool have_common_add = false;
	  while ((a < ins->atoms[i].add_by.size()) &&
		 (b < ins->atoms[j].add_by.size()) &&
		 !have_common_add) {
	    if (ins->atoms[i].add_by[a] < ins->atoms[j].add_by[b])
	      a += 1;
	    else if (ins->atoms[i].add_by[a] > ins->atoms[j].add_by[b])
	      b += 1;
	    else if (allowed_actions[ins->atoms[i].add_by[a]])
	      have_common_add = true;
	  }
	  if (!have_common_add) {
	    if (print_proof) {
	      std::cout << "Sometime(" << ins->atoms[i].name
			<< ") & Sometime(" << ins->atoms[j].name
			<< ") & NeverAfter(" << ins->atoms[i].name
			<< ", " << ins->atoms[j].name
			<< ") & NeverAfter(" << ins->atoms[j].name
			<< ", " << ins->atoms[i].name
			<< ") & NoCommonAdd(" << ins->atoms[i].name
			<< ", " << ins->atoms[j].name
			<< ") => CONTRADICTION"
			<< std::endl;
	    }
	    return true;
	  }
	}
  return false;
}

bool BasicPCC::check_amo_count(const index_set& amo_atoms)
{
  // find set of actions whose number of executions are limited
  // (to one) by each at-most-once constraint; also count how
  // many of the atoms contribute at least one action to this set
  index_set_vec amo_acts(EMPTYSET, 2 * amo_atoms.size());
  for (index_type i = 0; i < amo_atoms.size(); i++) {
    // actions that consume the atom (i.e., atom in pre and in del)
    for (index_type k = 0; k < ins->atoms[amo_atoms[i]].del_by.size(); k++)
      if (allowed_actions[ins->atoms[amo_atoms[i]].del_by[k]] &&
	  ins->actions[ins->atoms[amo_atoms[i]].del_by[k]].pre
	  .contains2(amo_atoms[i]))
	amo_acts[2*i].insert(ins->atoms[amo_atoms[i]].del_by[k]);
    if (print_proof) {
      std::cout << "AtMostOnce(" << ins->atoms[amo_atoms[i]].name
		<< ") => AtMostOnce(";
      ins->write_action_set(std::cout, amo_acts[2*i]);
      std::cout << ")" << std::endl;
    }
    // actions that add the atom and whose prec. is mutex with it
    if (mx) {
      for (index_type k = 0; k < ins->atoms[amo_atoms[i]].add_by.size(); k++)
	if (allowed_actions[ins->atoms[amo_atoms[i]].add_by[k]] &&
	    mx->mutex(ins->actions[ins->atoms[amo_atoms[i]].add_by[k]].pre,
		      amo_atoms[i]))
	  amo_acts[(2*i)+1].insert(ins->atoms[amo_atoms[i]].add_by[k]);
      if (print_proof) {
	std::cout << "AtMostOnce(" << ins->atoms[amo_atoms[i]].name
		  << ") => AtMostOnce(";
	ins->write_action_set(std::cout, amo_acts[(2*i)+1]);
	std::cout << ")" << std::endl;
      }
    }
  }
  // find set of not initially true active sometime atoms whose (still
  // allowed) adders all belong to sets of actions thus limited:
  // amo_sets and req_sets are parallel vectors (pairs):
  // req_sets[i] contains a set of atoms such that all allowed adders
  // of each one of them is contained in the union of the amo_acts
  // sets indexed by amo_sets[i], and, for each amo atom p, amo_sets[i]
  // does not contain both the set of consumers (even) and producers
  // (odd) of p.
  index_set_vec amo_sets;
  index_set_vec req_sets;
  for (index_type i = 0; i < ins->n_atoms(); i++)
    if (!ins->atoms[i].init && active_sometime[i]) {
      index_set s;
      bool ok = true;
      for (index_type k = 0; k < ins->atoms[i].add_by.size(); k++)
	if (allowed_actions[ins->atoms[i].add_by[k]]) {
	  index_type j = amo_acts.first_contains(ins->atoms[i].add_by[k]);
	  if (j == no_such_index) {
	    ok = false;
	  }
	  else if ((j % 2) == 1) {
	    if (s.contains(j - 1))
	      ok = false;
	    else
	      s.insert(j);
	  }
	  else {
	    if (s.contains(j + 1))
	      ok = false;
	    else
	      s.insert(j);
	  }
	}
      if (ok) {
	assert(req_sets.size() == amo_sets.size());
	bool app = true;
	for (index_type j = 0; j < amo_sets.size(); j++)
	  if (amo_sets[j].contains(s)) {
	    req_sets[j].insert(i);
	    if (amo_sets[j] == s)
	      app = false;
	  }
	if (app) {
	  amo_sets.append(s);
	  req_sets.append().assign_singleton(i);
	}
	// if the number of active sometime atoms that require one of the
	// actions in a union of amo_act sets to make it true exceeds the
	// maximum number of actions from this set that can be executed
	// (which is the number of amo_act sets in the union), we may have
	// a contradiction, but to be sure, we have to check for actions
	// that add more than one of the atoms
	for (index_type j = 0; j < amo_sets.size(); j++)
	  if (req_sets[j].size() > amo_sets[j].size()) {
	    // construct a graph over atoms in req_set, with an edge between
	    // any pair of atoms that are both added by the same action
	    graph req_ca(req_sets[j].size());
	    for (index_type i1 = 0; i1 < req_sets[j].size(); i1++)
	      for (index_type i2 = i1 + 1; i2 < req_sets[j].size(); i2++)
		if (ins->atoms[req_sets[j][i1]].add_by
		    .first_common(ins->atoms[req_sets[j][i2]].add_by)
		    != index_pair(no_such_index))
		  req_ca.add_undirected_edge(i1, i2);
	    // find an independent set in this graph: if the size of this
	    // independent set is greater than amo_sets[j], then we have
	    // a contradiction
	    index_set ind_req_set;
	    req_ca.apx_independent_set(ind_req_set);
	    if (ind_req_set.size() > amo_sets[j].size()) {
	      if (print_proof) {
		std::cout << "|";
		index_set s2(req_sets[j], ind_req_set);
		ins->write_atom_set(std::cout, s2);
		std::cout << "| > " << amo_sets[j].size()
			  << " => CONTRADICTION"
			  << std::endl;
	      }
	      return true;
	    }
	  }
      }
    }
  return false;
}

bool BasicPCC::run()
{
  stats.start();
  if (print_proof) {
    print_premises();
    std::cout << "BEGIN PROOF" << std::endl;
  }
  //std::cerr << "asserting given always constraints..." << std::endl;
  for (index_type k = 0; k < c_always.length(); k++)
    assert_always(c_always[k]);
  //std::cerr << "asserting given sometime-before constraints..." << std::endl;
  for (index_type k = 0; k < c_sometime_before.length(); k++)
    sometime_before.add_edge(c_sometime_before[k].first,
			     c_sometime_before[k].second);
  //std::cerr << "infering never from cyclic sb's..." << std::endl;
  infer_never_from_sb_cycle();
#ifdef NEW_VERSION_OF_OLD_PCC
  if (never_after)
    infer_never_from_na_and_sb();
#endif
  while (aa_changed) {
    sb_changed = false;
    //std::cerr << "updating landmarks..." << std::endl;
    update();
    if (sb_changed) {
      //std::cerr << "infering never from cyclic sb's..." << std::endl;
      infer_never_from_sb_cycle();
#ifdef NEW_VERSION_OF_OLD_PCC
      if (never_after)
	infer_never_from_na_and_sb();
#endif
    }
  }
  //std::cerr << rem_trlm.count(false) << " triggered landmarks" << std::endl;
  //std::cerr << "asserting goals..." << std::endl;
  for (index_type k = 0; k < ins->goal_atoms.length(); k++)
    assert_sometime(ins->goal_atoms[k]);
  //std::cerr << "asserting given sometime constraints..." << std::endl;
  for (index_type k = 0; k < c_sometime.length(); k++)
    assert_sometime(c_sometime[k]);
  bool con = check_sometime_never();
  if (con) {
    stats.stop();
    return true;
  }
  //std::cerr << "inferring more na's from sb's..." << std::endl;
#ifndef NEW_VERSION_OF_OLD_PCC
  infer_na_from_sb();
  con = check_asb_na();
  if (con) {
    stats.stop();
    return true;
  }
#endif
  if (never_after) {
    con = check_as_na_cycle();
    if (con) {
      stats.stop();
      return true;
    }
  }
  con = check_amo_count(c_at_most_once);
  if (con) {
    stats.stop();
    return true;
  }
  if (print_proof) {
    std::cout << "END PROOF" << std::endl;
  }
  stats.stop();
  //std::cerr << "failed to prove contradiction" << std::endl;
  return false;
}


class StrongPCC : public BasicPCC {
  bool own_ins;
  bool opt_H2;

  // landmark graph
  graph lm;

  void update_landmark_graph();
  void update();
  void update_never_after();

public:
  StrongPCC(Instance& i, Statistics& st, ppc_vec& v, const index_set& s,
	    bool optH2, bool opt_compile_amo, bool opt_print_proof);
  ~StrongPCC();
  bool run();
};

StrongPCC::StrongPCC
(Instance& i, Statistics& st, ppc_vec& v, const index_set& s,
 bool optH2, bool opt_compile_amo, bool opt_print_proof)
  : BasicPCC(i, st, v, s, opt_print_proof), opt_H2(optH2)
{
  if (opt_compile_amo && !c_at_most_once.empty()) {
    ins = i.copy();
    for (index_type k = 0; k < s.length(); k++)
      if (v[s[k]].pct == pc_at_most_once) {
	mapping map(ins->n_actions());
	ins->enforce_pc_at_most_once_conjunction
	  (v[s[k]].s_c, v[s[k]].name, map);
      }
    ins->cross_reference();
    own_ins = true;
  }
  else {
    ins = &i;
    own_ins = false;
  }
  mx = new StaticMutex(*ins, false);
  init_derived(0, 0);
  never_after = new graph(ins->n_atoms());
}

StrongPCC::~StrongPCC()
{
  delete mx;
  if (own_ins)
    delete ins;
  delete never_after;
}

void StrongPCC::update_landmark_graph()
{
  lm.init(ins->n_atoms());
  for (index_type i = 0; i < ins->n_atoms(); i++)
    if (possible_atoms[i]) {
      Instance* test_ins = new Instance(*ins);
      test_ins->atoms[i].init = false;
      for (index_type k = 0; k < test_ins->n_actions(); k++)
	test_ins->actions[k].add.subtract(i);
      test_ins->cross_reference();
      StaticMutex* test_mutex = new StaticMutex(*test_ins, false);
      test_mutex->recompute(allowed_actions);
      for (index_type j = 0; j < ins->n_atoms(); j++)
	if ((i != j) && possible_atoms[j]) {
	  if (test_mutex->unreachable(j))
	    lm.add_edge(i, j);
	}
      delete test_mutex;
      delete test_ins;
    }
}

void StrongPCC::update()
{
  mx->recompute(allowed_actions);
  for (index_type k = 0; k < ins->n_atoms(); k++)
    if (possible_atoms[k] && mx->unreachable(k)) {
      if (print_proof)
	std::cout << "-Reachable(" << ins->atoms[k].name << ")" << std::endl;
      assert_never(k);
    }
  for (index_type k = 0; k < ins->n_actions(); k++)
    if (allowed_actions[k] && mx->mutex(ins->actions[k].pre)) {
      if (print_proof)
	std::cout << "-Reachable(" << ins->actions[k].name << ")" << std::endl;
      allowed_actions[k] = false;
    }
  update_landmark_graph();
  for (index_type i = 0; i < ins->n_atoms(); i++)
    for (index_type j = 0; j < ins->n_atoms(); j++) {
      if (lm.adjacent(i, j) && !sometime_before.adjacent(i, j)) {
	if (print_proof) {
	  std::cout << "Landmark(" << ins->atoms[i].name
		    << "," << ins->atoms[j].name << ") => SometimeBefore("
		    << ins->atoms[i].name << "," << ins->atoms[j].name
		    << ")" << std::endl;
	}
	assert_sometime_before(j, i);
      }
    }
  aa_changed = false;
}

void StrongPCC::update_never_after()
{
  graph* nag = (graph*)never_after;
  nag->init(ins->n_atoms());
  for (index_type i = 0; i < ins->n_atoms(); i++)
    if (possible_atoms[i]) {
      index_set init;
      for (index_type j = 0; j < ins->n_atoms(); j++)
	if (possible_atoms[j] && !mx->mutex(i, j))
	  init.insert(j);
      StaticMutex m2(*ins, init);
      for (index_type j = 0; j < ins->n_atoms(); j++)
	if (possible_atoms[j] && m2.unreachable(j)) {
	  if (print_proof) {
	    std::cout << "NeverAfter(" << ins->atoms[i].name << ","
		      << ins->atoms[j].name << ")" << std::endl;
	  }
	  nag->add_edge(i, j);
	}
    }
  never_after = nag;
}

bool StrongPCC::run()
{
  stats.start();
  if (print_proof) {
    print_premises();
    std::cout << "BEGIN PROOF" << std::endl;
  }
  std::cerr << "asserting given always constraints..." << std::endl;
  for (index_type k = 0; k < c_always.length(); k++)
    assert_always(c_always[k]);
  std::cerr << "asserting given sometime-before constraints..." << std::endl;
  for (index_type k = 0; k < c_sometime_before.length(); k++)
    sometime_before.add_edge(c_sometime_before[k].first,
			     c_sometime_before[k].second);
  
  std::cerr << "initialising reachability/landmarks..." << std::endl;
  update();

  std::cerr << "infering never from cyclic sb's..." << std::endl;
  infer_never_from_sb_cycle();
  infer_never_from_na_and_sb();
  while (aa_changed) {
    std::cerr << "updating reachability/landmarks..." << std::endl;
    update();
    infer_never_from_sb_cycle();
    infer_never_from_na_and_sb();
  }

  std::cerr << "asserting goals..." << std::endl;
  for (index_type k = 0; k < ins->goal_atoms.length(); k++)
    assert_sometime(ins->goal_atoms[k]);
  std::cerr << "asserting given sometime constraints..." << std::endl;
  for (index_type k = 0; k < c_sometime.length(); k++)
    assert_sometime(c_sometime[k]);

  bool con = check_sometime_never();
  if (con) {
    stats.stop();
    return true;
  }

  std::cerr << "computing never-after..." << std::endl;
  update_never_after();

  con = check_as_na_cycle();
  if (con) {
    stats.stop();
    return true;
  }

  con = check_amo_count(c_at_most_once);
  if (con) {
    stats.stop();
    return true;
  }

  if (print_proof) {
    std::cout << "END PROOF" << std::endl;
  }
  stats.stop();
  std::cerr << "failed to prove contradiction" << std::endl;
  return false;
}

/// CEGA

bool cega
(Instance& ins,
 Statistics& stats,
 long test_time_limit,
 const ppc_vec& ppcv,
 const index_set& s,
 const StaticMutex& mx,
 const graph& lmg,
 bool opt_verbose,
 Schedule*& soln)
{
  bool dec = false;
  Statistics my_stats(&stats);
  if (opt_verbose) {
    std::cerr << "CEGA: testing set ";
    write_ppc_set(std::cerr, ppcv, s);
    std::cerr << std::endl;
  }
  if (test_time_limit > 0) {
    if (opt_verbose) {
      std::cerr << "CEGA: time limit is " << test_time_limit << " sec."
		<< std::endl;
    }
    my_stats.enable_time_out(test_time_limit);
  }
  my_stats.start();
  Instance* c_ins = ins.copy();
  mapping map(c_ins->n_actions());
  for (index_type k = 0; k < s.size(); k++) {
    ppcv[s[k]].enforce(*c_ins, map);
  }
  c_ins->cross_reference();
  // construct initial abstraction
  // step 1: collect atoms mentioned in goal and atoms created
  // by compiling in selected ppcs.
  index_set abs0(c_ins->goal_atoms);
  // include predecessors of goals in lmg - is this a good idea?
  //for (index_type i = 0; i < ins.goal_atoms.size(); i++)
  //  abs0.insert(lmg.predecessors(ins.goal_atoms[i]));
  for (index_type i = ins.n_atoms(); i < c_ins->n_atoms(); i++)
    abs0.insert(i);
  if (opt_verbose) {
    std::cerr << "CEGA: atoms of interest = ";
    c_ins->write_atom_set(std::cerr, abs0);
    std::cerr << std::endl;
  }
  // step 2: connect atoms of interest
  graph cg;
  c_ins->causal_graph(cg);
  cg.make_undirected();
  weighted_graph adist(abs0.size());
  for (index_type i = 0; i < abs0.size(); i++)
    for (index_type j = i + 1; j < abs0.size(); j++) {
      index_type dmin = cg.distance(abs0[i], abs0[j]);
      adist.add_undirected_edge
	(i, j, (dmin == no_such_index ? POS_INF : I_TO_N(dmin)));
    }
  graph amst;
  adist.minimum_spanning_tree(amst);
  index_set abs(abs0);
  for (index_type i = 0; i < abs0.size(); i++)
    for (index_type j = i + 1; j < abs0.size(); j++)
      if (amst.adjacent(i, j)) {
	index_vec p;
	index_type l = cg.shortest_path(abs0[i], abs0[j], p);
	abs.insert(p);
      }
  if (opt_verbose) {
    std::cerr << "CEGA: initial abstraction = ";
    c_ins->write_atom_set(std::cerr, abs);
    std::cerr << std::endl;
  }
  index_type iteration = 0;
  while (!dec && !my_stats.break_signal_raised()) {
    iteration += 1;
    index_vec atm_map;
    index_vec act_map;
    Instance* abs_ins = new Instance(c_ins->name);
    abs_ins->abstracted_copy(*c_ins, abs, atm_map, act_map);
    abs_ins->cross_reference();
    std::cerr << "CEGA: abstraction " << iteration << " has "
	      << abs_ins->n_atoms() << "/" << c_ins->n_atoms() << " atoms and "
	      << abs_ins->n_actions() << "/" << c_ins->n_actions() << " actions"
	      << std::endl;
    UnitACF acf;
    //CostTable h(*abs_ins, my_stats);
    //h.compute_H2(acf);
    //SeqRegState s0(*abs_ins, h, acf, abs_ins->goal_atoms);
    ForwardFF h(*abs_ins, acf, abs_ins->goal_atoms, my_stats);
    SeqProgState s0(*abs_ins, h, acf, abs_ins->init_atoms);
    ActionSequenceSet abs_plan_set;
    Result search_res(&abs_plan_set);
    search_res.set_stop_condition(Result::stop_at_first);
    BFS search(my_stats, search_res, 1000007);
    search.start(s0);
    std::cerr << "CEGA: search " << iteration << " done, "
	      << my_stats.total_time() << " sec. and "
	      << my_stats.total_nodes() << " nodes total"
	      << std::endl;
    if (search.solved()) {
      std::cerr << "CEGA: " << abs_plan_set.size()
		<< " abstract plan(s) found" << std::endl;
      PrintActions plan_printer(*abs_ins, std::cerr);
      abs_plan_set.output(plan_printer);
      std::cerr << std::endl;
      index_set new_abs;
      for (index_type k = 0; (k < abs_plan_set.size()) && !dec; k++) {
	ExecState as(*abs_ins, abs_ins->init_atoms);
	ExecErrorSet ases;
	bool ok = abs_plan_set[k].simulate(*abs_ins, as, 0, &ases, false);
	if (!ok) {
	  std::cerr << "error: abstract plan failed in abstraction!"
		    << std::endl << "errors: ";
	  ases.write(std::cerr);
	  std::cerr << std::endl;
	  exit(255);
	}
	index_set add_to_abs_k;
	ExecState vs(*c_ins, c_ins->init_atoms);
	//ExecTrace tr(*c_ins);
	//tr.append((ExecState*)vs.copy());
	Schedule vp(ins);
	for (index_type p = 0; p < abs_plan_set[k].size(); p++) {
	  index_set ca;
	  index_type ca_min = no_such_index;
	  index_set up_min;
	  mapping::inverse_map_image(act_map, abs_plan_set[k][p], ca);
	  //std::cerr << "abstract action " << p << ": " << std::endl;
	  //abs_ins->print_action(std::cerr,
	  //			abs_ins->actions[abs_plan_set[k][p]]);
	  //std::cerr << "candidate actions: ";
	  //c_ins->write_action_set(std::cerr, ca);
	  //std::cerr << std::endl << "current state: " << vs;
	  for (index_type i = 0; i < ca.size(); i++) {
	    ExecErrorSet es;
	    vs.applicable(ca[i], &es, p);
	    index_set up;
	    for (index_type j = 0; j < es.size(); j++)
	      if (es[j]->type_of_error() ==
		  ExecError::error_unsatisfied_precondition) {
		UnsatisfiedPreconditionError* u =
		  (UnsatisfiedPreconditionError*)es[j];
		assert(!abs.contains(u->precondition().index));
		up.insert(u->precondition().index);
	      }
	    if ((ca_min == no_such_index) || (up.size() < up_min.size())) {
	      ca_min = ca[i];
	      up_min.assign_copy(up);
	    }
	  }
	  assert(ca_min != no_such_index);
	  assert((map[ca_min] != no_such_index) &&
		 (map[ca_min] < ins.n_actions()));
	  //std::cerr << "selected action: " << std::endl;
	  //c_ins->print_action(std::cerr, c_ins->actions[ca_min]);
	  //std::cerr << "unsatisfied preconditions: ";
	  //c_ins->write_atom_set(std::cerr, up_min);
	  //std::cerr << std::endl;
	  vp.insert(map[ca_min]);
	  vs.apply(ca_min, 0, p);
	  NTYPE d = vs.min_delta();
	  //std::cerr << "state after apply: " << vs;
	  //std::cerr << "min delta = " << d << std::endl;
	  vp.advance(d);
	  vs.advance(d, 0);
	  //std::cerr << "state after advance: " << vs;
	  //tr.append((ExecState*)vs.copy());
	  add_to_abs_k.insert(up_min);
	}
	vp.end();
	if (add_to_abs_k.empty()) {
	  std::cerr << "CEGA: abstract plan appears to be valid!" << std::endl;
	  soln = new Schedule(vp);
	  //std::ofstream dump("cega-dump.txt");
	  //Instance::always_write_parameters = true;
	  //c_ins->write_domain(dump);
	  //c_ins->write_problem(dump);
	  //soln->write(dump, Name::NC_INSTANCE);
	  //tr.write(dump);
	  dec = true;
	}
	else {
	  new_abs.insert(add_to_abs_k);
	}
      }
      if (!dec) {
	assert(!new_abs.empty());
	std::cerr << "CEGA: extending abstraction with ";
	c_ins->write_atom_set(std::cerr, new_abs);
	std::cerr << std::endl;
	assert(!abs.contains(new_abs));
	abs.insert(new_abs);
      }
    }
    else if (!my_stats.break_signal_raised()) {
      std::cerr << "CEGA: abstract problem is unsolvable" << std::endl;
      soln = 0;
      dec = true;
    }
    else if (opt_verbose) {
      std::cerr << "CEGA: search interrupted" << std::endl;
    }
    delete abs_ins;
  }
  delete c_ins;
  my_stats.stop();
  return dec;
}

class CEGATest : public ConsistencyTest {
  Instance& instance;
  Statistics& cega_stats;
  count_type& n_cega_tests;
  const StaticMutex& mutex;
  const graph& landmark_graph;

 public:
  bool verbose;

  CEGATest(const ppc_vec& constraints,
	   Instance& ins,
	   Statistics& test_stats,
	   count_type& n_tests,
	   const StaticMutex& mx,
	   const graph& lmg);
  ~CEGATest();

  virtual bool test(const index_set& s, bool find_conflict, bool min_conflict);
};

CEGATest::CEGATest
(const ppc_vec& constraints,
 Instance& ins,
 Statistics& test_stats,
 count_type& n_tests,
 const StaticMutex& mx,
 const graph& lmg)
  : ConsistencyTest(constraints),
    instance(ins), cega_stats(test_stats), n_cega_tests(n_tests),
    mutex(mx), landmark_graph(lmg)
{
  // done
}

CEGATest::~CEGATest()
{
  // done
}

bool CEGATest::test
(const index_set& s, bool find_conflict, bool min_conflict)
{
  Schedule* soln = NULL;
  bool dec = cega(instance, cega_stats, 0, ppcs, s, mutex, landmark_graph,
		  verbose, soln);
  n_cega_tests += 1;
  if (dec && (soln == NULL)) {
    conflict.assign_copy(s);
    return true;
  }
  else {
    if (soln != NULL)
      delete soln;
    return false;
  }
}

/// basic planner (FF)

bool plan_compiled
(Instance& cins,
 Statistics& stats,
 long test_time_limit,
 bool opt_verbose,
 Schedule*& soln)
{
  Statistics my_stats(&stats);
  if (test_time_limit > 0) {
    my_stats.enable_time_out(test_time_limit);
  }
  UnitACF acf;
  ForwardFF h(cins, acf, cins.goal_atoms, my_stats);
  // SeqProgState s0(cins, h, acf, cins.goal_atoms);
  RedSeqProgState s0(cins, h, acf, cins.init_atoms, h);
  ActionSequenceSet plan_set;
  Result search_res(&plan_set);
  search_res.set_stop_condition(Result::stop_at_first);
  BFS search(my_stats, search_res);
  search.set_trace_level(0);
  search.greedy = true;
  search.start(s0);
  if (search.solved()) {
    std::cerr << "planner: plan found!" << std::endl;
    assert(plan_set.size() > 0);
    PrintActions plan_printer(cins, std::cerr);
    plan_set.output(plan_printer);
    std::cerr << std::endl;
  }
  else {
    std::cerr << "planner: search interrupted or exhausted" << std::endl;
  }
  my_stats.stop();
  return search.solved();
}

bool plan_compiled_lds
(Instance& cins,
 Statistics& stats,
 index_type d_limit,
 long test_time_limit,
 bool opt_verbose,
 Schedule*& soln)
{
  Statistics my_stats(&stats);
  if (test_time_limit > 0) {
    my_stats.enable_time_out(test_time_limit);
  }
  UnitACF acf;
  ForwardFF h(cins, acf, cins.goal_atoms, my_stats);
  ActionSequenceSet plan_set;
  Result search_res(&plan_set);
  search_res.set_stop_condition(Result::stop_at_first);
  BFS search(my_stats, search_res);
  search.greedy = true;
  index_type d = 0;
  while (!search.solved() && (d <= d_limit) && !my_stats.break_signal_raised()) {
    std::cerr << "trying d = " << d << "..." << std::endl;
    RedSeqProgState s0(cins, h, acf, cins.init_atoms, h, d);
    search.start(s0);
    if (search.solved()) {
      std::cerr << "planner: plan found!" << std::endl;
      assert(plan_set.size() > 0);
      PrintActions plan_printer(cins, std::cerr);
      plan_set.output(plan_printer);
      std::cerr << std::endl;
    }
    else {
      d += 1;
    }
  }
  my_stats.stop();
  if (!search.solved()) {
    std::cerr << "planner: search interrupted or exhausted" << std::endl;
  }
  return search.solved();
}

bool plan
(Instance& ins,
 const ppc_vec& ppcv,
 const index_set& sel,
 Statistics& stats,
 long test_time_limit,
 bool opt_verbose,
 Schedule*& soln)
{
  Instance* c_ins = ins.copy();
  mapping map(c_ins->n_actions());
  for (index_type k = 0; k < sel.size(); k++) {
    ppcv[sel[k]].enforce(*c_ins, map);
  }
  c_ins->cross_reference();
  bool solved =
    plan_compiled(*c_ins, stats, test_time_limit, opt_verbose, soln);
  // bool solved =
  //   plan_compiled_lds(*c_ins, stats, 3, test_time_limit, opt_verbose, soln);
  delete c_ins;
  return solved;
}

/// single h^m test

bool hm_test
(Instance& ins,
 const ppc_vec& ppcv,
 const index_set& sel,
 Statistics& stats,
 bool opt_H3,
 bool opt_H2,
 bool opt_verbose)
{
  bool inc = false;
  stats.start();
  Instance* test_ins = ins.copy();
  mapping map(test_ins->n_actions());
  for (index_type i = 0; i < sel.length(); i++) {
    if (opt_verbose) {
      std::cerr << "compiling ";
      ppcv[sel[i]].write_PDDL_constraint(std::cerr, *test_ins);
      std::cerr << "..." << std::endl;
    }
    ppcv[sel[i]].enforce(*test_ins, map);
    if (opt_verbose) {
      std::cerr << test_ins->n_actions() << " actions" << std::endl;
    }
  }
  test_ins->cross_reference();
  Heuristic* h0 = 0;
  if (opt_H3) {
    CostTable* hm = new CostTable(*test_ins, stats);
    hm->compute_H3(ZeroACF());
    h0 = hm;
  }
  else if (opt_H2) {
    StaticMutex* mx = new StaticMutex(*test_ins);
    h0 = mx;
  }
  else {
    h0 = new Reachability(*test_ins);
  }
  if (INFINITE(h0->eval(test_ins->goal_atoms))) {
    std::cerr << "set " << sel << " = ";
    write_ppc_set(std::cerr, ppcv, sel);
    std::cerr << " proven unachievable" << std::endl;
    inc = true;
  }
  delete h0;
  delete test_ins;
  stats.stop();
  return inc;
}

/// max value computation

#ifdef MVS_USE_CPLEX

NTYPE max_value
(const ppc_vec& ppcv, const hypergraph& h_inc, index_set& s_max)
{
  NTYPE sum = 0;
  IloEnv env;
  try {
    env.setOut(env.getNullStream());
    //env.setWarning(env.getNullStream());
    IloModel model(env);

    IloBoolVarArray vars(env, ppcv.size());
    IloNumArray value(env);
    for (index_type k = 0; k < ppcv.size(); k++) {
      value.add(N_TO_D(ppcv[k].weight));
    }
    model.add(vars);
    model.add(IloObjective(env, IloScalProd(value, vars))); 

    for (index_type k = 0; k < h_inc.n_edges(); k++) {
      IloExpr expr(env);
      for (index_type i = 0; i < h_inc.edge(k).size(); i++) {
	IloBoolVar &var = vars[h_inc.edge(k)[i]];
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

    s_max.fill(ppcv.size());
    for (index_type k = 0; k < ppcv.size(); k++) {
      double v = cplex.getValue(vars[k]);
      if (v > 0.5)
	s_max.subtract(k);
      else
	sum += ppcv[k].weight;
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
  std::cerr << "max val = " << sum << std::endl;
  // return optimal value
  return sum;
}

#else // !MVS_USE_CPLEX

NTYPE mvs_bound
(const ppc_vec& ppcv, index_vec& ord, const hypergraph& h_inc,
 index_type p, bool_vec& in, bool_vec& not_in, bool_vec& active_edge)
{
  NTYPE ub1 = 0;
  for (index_type k = p + 1; k < ppcv.size(); k++)
    if (!not_in[ord[k]])
      ub1 += ppcv[ord[k]].weight;
#ifdef MVS_USE_IMPROVED_BOUND
  index_set ae_set(active_edge);
  weighted_graph aeg(ae_set.size());
  for (index_type k = 0; k < ae_set.size(); k++) {
    NTYPE w_min = POS_INF;
    for (index_type i = 0; i < h_inc.edge(ae_set[k]).size(); i++) {
      //if(not_in[h_inc.edge(ae_set[k])[i]]) {
      //	std::cerr << "assert fail:" << std::endl;
      //	std::cerr << "active edges = " << active_edge << " = " << ae_set
      //		  << std::endl;
      //	std::cerr << "active edge #" << k << " = " << ae_set[k]
      //		  << " = " << h_inc.edge(ae_set[k])
      //		  << std::endl;
      //	std::cerr << "not_in = " << not_in << std::endl;
      //}
      assert(!not_in[h_inc.edge(ae_set[k])[i]]);
      if (!in[h_inc.edge(ae_set[k])[i]]) {
	w_min = MIN(w_min, ppcv[h_inc.edge(ae_set[k])[i]].weight);
	for (index_type j = k + 1; j < ae_set.size(); j++)
	  if (h_inc.edge(ae_set[j]).contains(h_inc.edge(ae_set[k])[i]))
	    aeg.add_undirected_edge(k, j);
      }
    }
    assert(FINITE(w_min));
    aeg.set_weight(k, w_min);
  }
  index_set s;
  NTYPE ws = aeg.apx_weighted_independent_set(s);
  ub1 -= ws;
#endif
  return ub1;
}

void max_value
(const ppc_vec& ppcv, index_vec& ord, const hypergraph& h_inc,
 index_type p, bool_vec& in, NTYPE val, bool_vec& not_in,
 bool_vec& active_edge, index_vec& n_in_edge,
 NTYPE& best_val, index_set& best_set)
{
  if (p < ppcv.size()) {
    bool any_edge_changed = false;
    // if possible, try include ppc #p
    if (!not_in[ord[p]]) {
      in[ord[p]] = true;
      NTYPE new_val = val + ppcv[ord[p]].weight;
      bool_vec ch_edge(false, h_inc.n_edges());
      bool_vec ch_active(false, h_inc.n_edges());
      bool_vec ch_not_in(false, ppcv.size());
      for (index_type i = 0; i < h_inc.n_edges(); i++)
	if (active_edge[i] && h_inc.edge(i).contains(ord[p])) {
	  assert(n_in_edge[i] > 1);
	  n_in_edge[i] -= 1;
	  ch_edge[i] = true;
	  any_edge_changed = true;
	  if (n_in_edge[i] == 1) {
	    index_type o = no_such_index;
	    for (index_type j = 0;
		 (j < h_inc.edge(i).size()) && (o == no_such_index); j++)
	      if (!in[h_inc.edge(i)[j]])
		o = h_inc.edge(i)[j];
	    assert(o != no_such_index);
	    assert(o < ppcv.size());
	    if (!not_in[o]) {
	      not_in[o] = true;
	      ch_not_in[o] = true;
	      for (index_type j = 0; j < h_inc.adjacent_edges(o).size(); j++)
		if (active_edge[h_inc.adjacent_edges(o)[j]]) {
		  active_edge[h_inc.adjacent_edges(o)[j]] = false;
		  ch_active[h_inc.adjacent_edges(o)[j]] = true;
		}
	    }
	  }
	}
      //for (index_type k = 0; k < h_inc.n_edges(); k++)
      //	if (active_edge[k]) {
      //	  if (h_inc.edge(k).have_common_element(not_in)) {
      //	    std::cerr << "assert fail #1:" << std::endl;
      //	    std::cerr << "active edges = " << active_edge << std::endl;
      //	    std::cerr << "edge #" << k << " = " << h_inc.edge(k)
      //		      << std::endl;
      //	    std::cerr << "not_in = " << not_in << std::endl;
      //	    exit(1);
      //	  }
      //	}
      //std::cerr << "before call 1: active_edge = " << active_edge << std::endl;
      //std::cerr << "before call 1: not_in = " << not_in << std::endl;
      NTYPE ub = mvs_bound(ppcv, ord, h_inc, p, in, not_in, active_edge);
      if ((new_val + ub) > best_val) {
	max_value(ppcv, ord, h_inc, p + 1, in, new_val, not_in,
		  active_edge, n_in_edge, best_val, best_set);
      }
      // reset
      for (index_type k = 0; k < ppcv.size(); k++)
	if (ch_not_in[k]) not_in[k] = false;
      for (index_type i = 0; i < h_inc.n_edges(); i++)
	if (ch_edge[i]) n_in_edge[i] += 1;
      for (index_type i = 0; i < h_inc.n_edges(); i++)
	if (ch_active[i]) active_edge[i] = true;
      in[ord[p]] = false;
    }
    // if not possible to include #p, or if including #p changed
    // any edge, try excluding ppc #p
    if (not_in[ord[p]] || any_edge_changed) {
      bool ch_not_in = !not_in[ord[p]];
      not_in[ord[p]] = true;
      bool_vec ch_active(false, h_inc.n_edges());
      for (index_type j = 0; j < h_inc.adjacent_edges(ord[p]).size(); j++)
	if (active_edge[h_inc.adjacent_edges(ord[p])[j]]) {
	  active_edge[h_inc.adjacent_edges(ord[p])[j]] = false;
	  ch_active[h_inc.adjacent_edges(ord[p])[j]] = true;
	}
      //for (index_type k = 0; k < h_inc.n_edges(); k++)
      //	if (active_edge[k]) {
      //	  if (h_inc.edge(k).have_common_element(not_in)) {
      //	    std::cerr << "assert fail #2:" << std::endl;
      //	    std::cerr << "active edges = " << active_edge << std::endl;
      //	    std::cerr << "edge #" << k << " = " << h_inc.edge(k)
      //		      << std::endl;
      //	    std::cerr << "not_in = " << not_in << std::endl;
      //	    exit(1);
      //	  }
      //	}
      //std::cerr << "before call 2: active_edge = " << active_edge << std::endl;
      //std::cerr << "before call 2: not_in = " << not_in << std::endl;
      NTYPE ub = mvs_bound(ppcv, ord, h_inc, p, in, not_in, active_edge);
      if ((val + ub) > best_val) {
	max_value(ppcv, ord, h_inc, p + 1, in, val, not_in,
		  active_edge, n_in_edge, best_val, best_set);
      }
      if (ch_not_in) not_in[ord[p]] = false;
      for (index_type i = 0; i < h_inc.n_edges(); i++)
	if (ch_active[i]) active_edge[i] = true;
    }
  }
  else {
    assert(val > best_val);
    best_val = val;
    best_set = index_set(in);
  }
}

class max_value_search_order : public lvector<index_type>::order {
  const ppc_vec& ppcv;
  const hypergraph& h_inc;
public:
  max_value_search_order(const ppc_vec& v, const hypergraph& h)
    : ppcv(v), h_inc(h) { };
//   bool operator()(const index_type& v0, const index_type& v1) const {
//     assert(v0 < ppcv.size());
//     assert(v1 < ppcv.size());
//     return (ppcv[v0].weight > ppcv[v1].weight);
//   };
  bool operator()(const index_type& v0, const index_type& v1) const {
    assert(v0 < h_inc.size());
    assert(v1 < h_inc.size());
    return (h_inc.adjacent_edges(v0) > h_inc.adjacent_edges(v1));
  };
};

index_type greedy_select_next
(const ppc_vec& ppcv, const bool_vec& in, const bool_vec ex)
{
  index_type p = no_such_index;
  NTYPE w_max = 0;
  for (index_type k = 0; k < ppcv.size(); k++)
    if (!in[k] && !ex[k])
      if (ppcv[k].weight > w_max) {
	w_max = ppcv[k].weight;
	p = k;
      }
  return p;
}

NTYPE greedy_max_value
(const ppc_vec& ppcv, const hypergraph& h_inc, index_set& s_max)
{
  bool_vec in(false, ppcv.size());
  bool_vec ex(false, ppcv.size());
  for (index_type k = 0; k < h_inc.n_edges(); k++)
    if (h_inc.edge(k).size() == 1)
      ex[h_inc.edge(k)[0]] = true;
  index_type p = greedy_select_next(ppcv, in, ex);
  while (p != no_such_index) {
    in[p] = true;
    for (index_type k = 0; k < h_inc.n_edges(); k++)
      if (h_inc.edge(k).count_common(in) == (h_inc.edge(k).size() - 1)) {
	index_set s(h_inc.edge(k));
	s.subtract(in);
	assert(s.size() == 1);
	assert(!in[s[0]]);
	ex[s[0]] = true;
      }
    p = greedy_select_next(ppcv, in, ex);
  }
  NTYPE val = 0;
  s_max.clear();
  for (index_type k = 0; k < ppcv.size(); k++)
    if (in[k]) {
      val += ppcv[k].weight;
      s_max.insert(k);
    }
  return val;
}

NTYPE max_value
(const ppc_vec& ppcv, const hypergraph& h_inc, index_set& s_max)
{
  bool_vec in(false, ppcv.size());
  bool_vec not_in(false, ppcv.size());
  bool_vec active_edge(true, h_inc.n_edges());
  index_vec n_in_edge(0, h_inc.n_edges());
  for (index_type i = 0; i < h_inc.n_edges(); i++) {
    n_in_edge[i] = h_inc.edge(i).size();
    if (h_inc.edge(i).size() == 1) {
      index_type o = h_inc.edge(i)[0];
      not_in[o] = true;
      for (index_type j = 0; j < h_inc.adjacent_edges(o).size(); j++)
	if (active_edge[h_inc.adjacent_edges(o)[j]])
	  active_edge[h_inc.adjacent_edges(o)[j]] = false;
    }
  }
#ifdef MVS_USE_GREEDY_INITIAL
  NTYPE v_max = greedy_max_value(ppcv, h_inc, s_max);
  //std::cerr << "greedy: v_max = " << v_max << ", s_max = ";
  //write_ppc_set(std::cerr, ppcv, s_max);
  //std::cerr << std::endl;
#else
  NTYPE v_max = 0;
  s_max.clear();
#endif
  index_vec ord;
#ifdef MVS_USE_VARIABLE_ORDER
  max_value_search_order mvso(ppcv, h_inc);
  for (index_type k = 0; k < ppcv.size(); k++)
    ord.insert_ordered(k, mvso);
#else
  index_vec_util::fill(ord, ppcv.size());
#endif
  //std::cerr << "ord = " << ord << std::endl;
  max_value(ppcv, ord, h_inc, 0, in, 0, not_in,
	    active_edge, n_in_edge, v_max, s_max);
  return v_max;
}

#endif // !MVS_USE_CPLEX

/// main bits

bool save_split = false;
count_type file_write_limit = 1000;
count_type files_written = 0;

typedef lvector<plan_vec> plan_set_vec;

void analyse_plan
(Instance& ins,
 ppc_vec& ppcv,
 Schedule* plan,
 bool untimed,
 index_set_vec& max_sat,
 plan_set_vec& max_plan,
 bool opt_verbose)
{
  ExecTrace* trace = new ExecTrace(ins);
  ExecErrorSet* errors = new ExecErrorSet();
  bool ok = plan->simulate(trace, errors, false);
  if (!ok) {
    std::cerr << "plan " << plan->plan_name() << " not executable: ";
    errors->write(std::cerr);
    std::cerr << std::endl;
  }
  else {
    if (untimed) {
      ExecTrace* t1 = trace;
      trace = t1->stable_trace();
      delete t1;
    }
    index_set sat;
    NTYPE val = 0;
    for (index_type i = 0; i < ppcv.length(); i++) {
      if (ppcv[i].test(trace)) {
	if (opt_verbose) {
	  std::cerr << "plan " << plan->plan_name()
		    << " satisfies #" << i << ": ";
	  ppcv[i].write(std::cerr, ins);
	  //std::cerr << std::endl;
	}
	sat.insert(i);
	val += ppcv[i].weight;
      }
      else if (opt_verbose) {
	  std::cerr << "plan " << plan->plan_name()
		    << " violates #" << i << ": ";
	  ppcv[i].write(std::cerr, ins);
	  //std::cerr << std::endl;
      }
    }
    std::cerr << "plan " << plan->plan_name()
	      << ": final value = " << PRINT_NTYPE(val)
	      << std::endl;
    // check if set of satisfied constraints is a strict subset
    // of an already known maximal set
    if (max_sat.first_strict_superset(sat) == no_such_index) {
      // if not, check if it equals a known maximal set
      index_type e = max_sat.first(sat);
      if (e != no_such_index) {
	// if so, add the plan to the planset for that set
	assert(max_plan.size() == max_sat.size());
	max_plan[e].append(plan);
      }
      else {
	// if not, this is a new maximal set: check if it dominates
	// some of the previous ones
	bool_vec d(false, max_sat.size());
	for (index_type k = 0; k < max_sat.size(); k++)
	  if (sat.contains(max_sat[k])) d[k] = true;
	max_sat.remove(d);
	max_plan.remove(d);
	max_sat.append(sat);
	max_plan.append().append(plan);
	assert(max_plan.size() == max_sat.size());
      }
    }
  }
  delete errors;
  delete trace;
}

void analyse_input_plans
(PDDL_Base* b,
 Instance& ins,
 ppc_vec& ppcv,
 const index_vec& action_map,
 bool untimed,
 index_set_vec& max_sat,
 plan_set_vec& max_plan,
 bool opt_verbose)
{
  NTYPE v_max_plan = 0;
  for (index_type k = 0; k < b->n_plans(); k++) {
    Schedule* plan = new Schedule(ins);
    b->export_plan(k, ins, action_map, *plan);
    analyse_plan(ins, ppcv, plan, untimed, max_sat, max_plan, opt_verbose);
  }
}

void save_problem
(Parser* rd, Instance& ins, const ppc_vec& ppcs, NTYPE v_sum,
 const index_set& sel, const char* desc, index_type num)
{
  if (files_written++ > file_write_limit) {
    std::cerr << "error: writing too many files!" << std::endl;
    exit(1);
  }
  std::ostringstream fname;
  char* b = rd->problem_file_basename();
  if (b) fname << b << "-";
  fname << desc << num << ".pddl";
  std::cerr << "writing file " << fname.str() << "..." << std::endl;
  std::ofstream s_out(fname.str().c_str());
  s_out << ";; " << ins.name << ", problem "
	<< desc << " #" << num << std::endl;
  s_out << ";; selected constraints: {";
  NTYPE val = 0;
  for (index_type i = 0; i < sel.length(); i++) {
    if (i > 0) s_out << ", ";
    s_out << ppcs[sel[i]].name;
    val += ppcs[sel[i]].weight;
  }
  s_out << "}" << std::endl;
  s_out << ";; value = " << PRINT_NTYPE(val)
	<< ", penalty = " << PRINT_NTYPE(v_sum - val)
	<< std::endl;
  rd->write_problem_begin(s_out);
  rd->write_objects(s_out, true);
  rd->write_init(s_out);
  s_out << " (:goal (and";
  for (index_type i = 0; i < rd->dom_goals.length(); i++) {
    if (rd->dom_goals[i]->is_state() &&
	rd->dom_goals[i]->is_propositional()) {
      s_out << " ";
      rd->dom_goals[i]->print(s_out);
    }
  }
  s_out << "))" << std::endl << " (:constraints (and";
  for (index_type i = 0; i < sel.length(); i++) {
    s_out << " ";
    ppcs[sel[i]].src->goal->print(s_out);
  }
  s_out << "))" << std::endl << ")" << std::endl;
  s_out.close();
}

void save_compiled_problem
(Parser* rd, Instance& ins, const ppc_vec& ppcs, NTYPE v_sum,
 Preprocessor& prep,
 const index_set& sel, Instance* c_ins, const char* desc, index_type num)
{
  bool own = false;
  if (c_ins == 0) {
    c_ins = ins.copy();
    mapping map(c_ins->n_actions());
    for (index_type i = 0; i < sel.length(); i++)
      ppcs[sel[i]].enforce(*c_ins, map);
    own = true;
  }
  if (files_written++ > file_write_limit) {
    std::cerr << "error: writing too many files!" << std::endl;
    exit(1);
  }
  std::ostringstream fname;
  char* b = rd->problem_file_basename();
  if (b) fname << b << "-";
  if (save_split)
    fname << desc << num << "-compiled.domain.pddl";
  else
    fname << desc << num << "-compiled.pddl";
  std::cerr << "writing file " << fname.str() << "..." << std::endl;
  std::ofstream s_out(fname.str().c_str());
  s_out << ";; " << ins.name << ", problem "
	<< desc << " #" << num << std::endl;
  s_out << ";; selected constraints: {";
  NTYPE val = 0;
  for (index_type i = 0; i < sel.length(); i++) {
    if (i > 0) s_out << ", ";
    s_out << ppcs[sel[i]].name;
    val += ppcs[sel[i]].weight;
  }
  s_out << "}" << std::endl;
  s_out << ";; value = " << PRINT_NTYPE(val)
	<< ", penalty = " << PRINT_NTYPE(v_sum - val)
	<< std::endl;
  c_ins->write_domain_header(s_out);
  c_ins->write_domain_declarations(s_out);
  c_ins->write_domain_actions(s_out);
  if (HSPS::Instance::write_DKEL)
    c_ins->write_domain_DKEL_items(s_out);
  if (HSPS::Instance::write_extra) {
    name_vec set_names(0, 0);
    index_set_vec sets;
    rd->export_action_partitions(set_names, sets);
    s_out << ";; action partitions" << std::endl;
    for (index_type k = 0; k < sets.length(); k++) {
      ins.remap_set(sets[k], prep.action_map);
      ins.write_domain_action_set(s_out, sets[k], set_names[k]);
    }
  }
  s_out << ")";
  if (save_split) {
    s_out.close();
    std::ostringstream pname;
    if (b) pname << b << "-";
    pname << desc << num << "-compiled.problem.pddl";
    std::ofstream p_out(pname.str().c_str());
    c_ins->write_problem(p_out);
    p_out.close();
  }
  else {
    c_ins->write_problem(s_out);
    s_out.close();
  }
  if (own) {
    delete c_ins;
  }
}

// find all minimal and maximal open sets:
// - h_inc is the hypergraph of known-to-be inconsistent sets
// - max_sat is the list of maximal known-to-be satisfiable sets
// - n is the number of elements
// a set is "open" if (a) it constains no h_inc edge (i.e., it is
//  an independent set wrt the hypergraph) and (b) it is not
//  contained in any max_sat set.

void compute_minimal_open
(const hypergraph& h_inc,
 const index_set_vec& max_sat,
 index_type n,
 index_set_vec& min_open)
{
  // make a hypergraph whose edges are complements of max_sat sets
  hypergraph h_max(n);
  for (index_type k = 0; k < max_sat.size(); k++) {
    index_set e;
    e.fill(n);
    e.subtract(max_sat[k]);
    h_max.add_edge(e);
  }
  // find minimal edge-covering sets in h_max
  index_set_vec cs;
  h_max.covering_sets_2(cs);

  // every h_max edge-cover cs[k] is a set that contains at least one
  // element not in each max_sat set; thus, it is not a subset of any
  // max_sat set. moreover, these are all *minimal* sets with this
  // property. thus, if cs[k] does not contain any edge of h_inc, it
  // is a minimal open set.
  for (index_type k = 0; k < cs.size(); k++)
    if (h_inc.edges().first_subset(cs[k]) == no_such_index)
      min_open.append(cs[k]);
}

void compute_maximal_open
(const hypergraph& h_inc,
 const index_set_vec& max_sat,
 index_type n,
 index_set_vec& max_open)
{
  // find the maximal independent sets of h_inc
  index_set_vec non_inc;
  h_inc.independent_sets(non_inc);

  // each independent set is a non-superset of every known-to-be
  // inconsistent set (i.e., h_inc edge), and it is maximal; thus, if
  // it is not contained in any known satisfiable set, it's a maximal
  // open set.
  for (index_type k = 0; k < non_inc.size(); k++)
    if (max_sat.first_superset(non_inc[k]) == no_such_index)
      max_open.append(non_inc[k]);
}

// // OLD IMPL.
// void compute_minimal_open
// (const hypergraph& h_inc,
//  const index_set_vec& max_sat,
//  index_type n,
//  index_set_vec& min_open)
// {
//   // find minimal edge-covering sets in h_inc
//   index_set_vec cs;
//   h_inc.covering_sets(cs);
//   //std::cerr << "#cs = " << cs.size() << std::endl;
//   for (index_type k = 0; k < cs.length(); k++) {
//     // the complement of each edge-covering set is an independent set
//     // for h_inc; if it is not contained in any max_sat set, it is an
//     // open set, and moreover it is a *maximal* open set (and if it is
//     // contained in one of the max_sat sets, it's not open, and there's
//     // no way to make it open).
//     index_set is;
//     is.fill(n);
//     is.subtract(cs[k]);
//     if (max_sat.first_superset(is) == no_such_index) {
//       //std::cerr << "cs #" << k << " is a candidate" << std::endl;
//       // but we're looking for *minimal* open sets, so we have to consider
//       // the subsets of the independent set that are not contained in any
//       // max_sat set. to do this, we build a "complementary" hypergraph,
//       // h_cand, where the i:th edge is (a) a subset of the independent set
//       // and (b) contains only elements not max_sat[i]; then we find the
//       // (minimal) edge covers of this hypergraph.
//       hypergraph h_cand(n);
//       for (index_type i = 0; i < max_sat.length(); i++) {
// 	index_set e(is);
// 	e.subtract(max_sat[i]);
// 	assert(!e.empty());
// 	h_cand.add_edge(e);
//       }
//       index_set_vec ccs;
//       h_cand.covering_sets(ccs);
//       //std::cerr << "#ccs = " << ccs.size() << std::endl;
//       for (index_type k = 0; k < ccs.length(); k++)
// 	min_open.insert_minimal(ccs[k]);
//     }
//   }
// }

index_type generate_constraints
(PDDL_Base* b,
 Instance& ins,
 Preprocessor& prep,
 ppc_vec& ppcv,
 const index_vec& action_map,
 bool untimed,
 index_type md)
{
  HSPS::ScheduleSet plans(ins);
  for (index_type k = 0; k < b->n_plans(); k++) {
    Schedule* s = new Schedule(ins);
    bool ok = b->export_plan(k, ins, action_map, *s);
    if (ok) {
      plans.add_schedule(s);
    }
    else {
      std::cerr << "error: failed to export plan #" << k;
      if (b->input_plans[k]->name)
	std::cerr << ": " << b->input_plans[k]->name;
      std::cerr << std::endl;
    }
  }
  TraceSet* trajs = plans.trace_set(true, true, untimed);
  std::cerr << trajs->length() << " input traces" << std::endl;
  if (trajs->length() == 0) return 0;
  StaticMutex* mx = prep.inconsistency();
  index_set_vec common_always;
  trajs->find_common_always(mx, md, common_always);
  std::cerr << "found " << common_always.length()
	    << " common A-constraints" << std::endl;
  for (index_type i = 0; i < common_always.size(); i++) {
    if (common_always[i].size() == 1)
      ppcv.append(ppc(new NameWithContext(new EnumName("A", i), Name::NC_INSTANCE, 0),
		      1, pc_always, common_always[i], false, NULL));
    else
      ppcv.append(ppc(new NameWithContext(new EnumName("A", i), Name::NC_INSTANCE, 0),
		      1, pc_always, common_always[i], true, NULL));
  }
  HSPS::index_set common_sometime;
  HSPS::pair_set common_sb;
  HSPS::pair_set common_sa;
  trajs->find_common_order_constraints(prep.necessary_sb_graph(),
				       common_sometime, common_sb, common_sa);
  std::cerr << "found "
	    << common_sometime.length() << " common E-constraints, "
	    << common_sb.length() << " common SB-constraints, and "
	    << common_sa.length() << " common SA-constraints"
	    << std::endl;
  // create common constraints
  for (index_type i = 0; i < common_sometime.size(); i++) {
    index_set s(common_sometime[i]);
    ppcv.append(ppc(new NameWithContext(new EnumName("E", i),
					Name::NC_INSTANCE, 0),
		    1, pc_sometime, s, false, NULL));
  }
  for (index_type i = 0; i < common_sb.size(); i++) {
    index_set sc(common_sb[i].second);
    index_set st(common_sb[i].first);
    ppcv.append(ppc(new NameWithContext(new EnumName("SB", i),
					Name::NC_INSTANCE, 0),
		    1, pc_sometime_before, sc, false, st, false, NULL));
  }
  for (index_type i = 0; i < common_sa.size(); i++) {
    index_set sc(common_sa[i].second);
    index_set st(common_sa[i].first);
    ppcv.append(ppc(new NameWithContext(new EnumName("SA", i),
					Name::NC_INSTANCE, 0),
		    1, pc_sometime_after, sc, false, st, false, NULL));
  }
  index_type n_common = ppcv.size();
  // create negations of the common constraints
  for (index_type i = 0; i < common_always.size(); i++) {
    ins.complete_atom_negations(common_always[i]);
    index_set s;
    ins.negation_atom_set(common_always[i], s);
    if (!mx->mutex(s)) {
      ppcv.append(ppc(new NameWithContext(new EnumName("notA", i),
					  Name::NC_INSTANCE, 0),
		      1, pc_sometime, s, false, NULL));
    }
  }
  for (index_type i = 0; i < common_sometime.size(); i++) {
    index_type a = ins.complete_atom_negation(common_sometime[i]);
    if (ins.atoms[a].init) {
      index_set s(a);
      ppcv.append(ppc(new NameWithContext(new EnumName("notE", i),
					  Name::NC_INSTANCE, 0),
		      1, pc_always, s, false, NULL));
    }
  }
  for (index_type i = 0; i < common_sb.size(); i++) {
    // SB(b,a) is unsatisfiable if b is initially true
    if (!ins.atoms[common_sb[i].second].init) {
      index_set sc(common_sb[i].first);
      index_set st(common_sb[i].second);
      ppcv.append(ppc(new NameWithContext(new EnumName("notSB", i),
					  Name::NC_INSTANCE, 0),
		      1, pc_sometime_before, sc, false, st, false, NULL));
    }
  }
  for (index_type i = 0; i < common_sa.size(); i++) {
    // the reverse of an SA-constraint might be satisfied
    bool ok = true;
    for (index_type k = 0; (k < trajs->length()) && ok; k++) {
      if ((*trajs)[k].first->test_sometime_after(common_sa[i].second,
						 common_sa[i].first))
	ok = false;
    }
    if (ok) {
      index_set sc(common_sa[i].first);
      index_set st(common_sa[i].second);
      ppcv.append(ppc(new NameWithContext(new EnumName("notSA", i),
					  Name::NC_INSTANCE, 0),
		      1, pc_sometime_after, sc, false, st, false, NULL));
    }
  }
  if (!ins.cross_referenced())
    ins.cross_reference();
  return n_common;
}

void save_preference_set
(Instance& ins, const ppc_vec& ppcv, const index_set& s, std::string name)
{
  std::string filename = name + ".pddl";
  std::ofstream s_out(filename.c_str());
  s_out << "(define (problem " << name << ")" << std::endl;
  s_out << "  (:constraints (and" << std::endl;
  for (index_type k = 0; k < s.size(); k++) {
    s_out << "    ";
    ppcv[s[k]].write(s_out, ins);
    s_out << std::endl;
  }
  s_out << "  ))" << std::endl;
  s_out << "  (:metric minimize (+ " << std::endl;
  for (index_type k = 0; k < s.size(); k++) {
    s_out << "    (* (is-violated " << ppcv[s[k]].name << ") "
	  << PRINT_NTYPE(ppcv[s[k]].weight) << ")" << std::endl;
  }
  s_out << "  ))";
  s_out << " )" << std::endl;
  s_out.close();
}

void save_constraint_set
(Instance& ins, const ppc_vec& ppcv, const index_set& s, std::string name)
{
  std::string filename = name + ".pddl";
  std::ofstream s_out(filename.c_str());
  s_out << "(define (problem " << name << ")" << std::endl;
  s_out << "  (:constraints (and" << std::endl;
  for (index_type k = 0; k < s.size(); k++) {
    s_out << "    ";
    ppcv[s[k]].write_PDDL_constraint(s_out, ins);
    s_out << std::endl;
  }
  s_out << "  ))" << std::endl;
  s_out << " )" << std::endl;
  s_out.close();
}

void variety
(Instance& ins,
 Preprocessor& prep,
 ppc_vec& ppcv,
 index_type n_common,
 hypergraph& h_inc,
 index_type variations_limit,
 index_type failure_limit,
 Statistics& stats,
 RNG& rng)
{
  assert(ppcv.size() > n_common);

  // prepare for PCC
  set_edge_vec trlm;
  StaticMutex* mx = prep.inconsistency();
  graph* lmg = prep.necessary_sb_graph();
  landmark_graph_triggered_edges(ins, *lmg, trlm);
#ifdef NEW_PCC
  PCC test(ins, stats, ppcv, mx, lmg, &trlm, prep.never_after_graph());
#endif

  index_set common_unary;
  for (index_type k = 0; k < n_common; k++)
    if (((ppcv[k].pct == pc_always) || (ppcv[k].pct == pc_sometime)) &&
	(ppcv[k].s_c.size() == 1))
      common_unary.insert(k);
  std::cerr << common_unary.size() << " common unary constraints" << std::endl;
  if (common_unary.size() < 2) return;
  index_type n_variations = 0;
  index_type n_failures = 0;
  while (n_variations < variations_limit) {
    index_set cset;
    index_set eatoms;
    index_set rs;
    rng.select_fixed_set(rs, common_unary.size() / 2, common_unary.size());
    for (index_type i = 0; i < rs.size(); i++) {
      index_type c = common_unary[rs[i]];
      cset.insert(c);
      if (ppcv[c].pct == pc_sometime)
	eatoms.insert(ppcv[c].s_c[0]);
    }
    index_set common_binary;
    for (index_type k = 0; k < n_common; k++)
      if (((ppcv[k].pct == pc_sometime_before) ||
	   (ppcv[k].pct == pc_sometime_after)) &&
	  (ppcv[k].s_t.size() == 1) && (ppcv[k].s_c.size() == 1) &&
	  eatoms.contains(ppcv[k].s_t[0]) &&
	  eatoms.contains(ppcv[k].s_c[0]))
	common_binary.insert(k);
    std::cerr << common_binary.size() << " common binary constraints"
	      << std::endl;
    if (common_binary.size() < cset.size()) {
      n_failures += 1;
      if (n_failures >= failure_limit) {
	std::cerr << "failure limit exceeded!" << std::endl;
	return;
      }
    }
    else {
      rng.select_fixed_set(rs, cset.size(), common_binary.size());
      for (index_type i = 0; i < rs.size(); i++) {
	index_type c = common_binary[rs[i]];
	cset.insert(c);
      }
#ifdef NEW_PCC
      test.set_tracker(new InferenceTracker(test));
      bool common_unsat = test.run(cset);
      if (common_unsat) {
	test.tracker()->print_proof(std::cout);
	save_preference_set(ins, ppcv, cset, "failure");
      }
      assert(!common_unsat);
      test.set_tracker(NULL);
#endif
      index_set out;
      for (index_type k = n_common; k < h_inc.n_edges(); k++) {
	if (h_inc.edge(k).size() == 2) {
	  if (cset.contains(h_inc.edge(k)[0]))
	    out.insert(h_inc.edge(k)[1]);
	  else if (cset.contains(h_inc.edge(k)[1]))
	    out.insert(h_inc.edge(k)[0]);
	}
	else if (h_inc.edge(k).size() == 3) {
	  if (cset.contains(h_inc.edge(k)[0]) &&
	      cset.contains(h_inc.edge(k)[1]))
	    out.insert(h_inc.edge(k)[2]);
	  else if (cset.contains(h_inc.edge(k)[0]) &&
		   cset.contains(h_inc.edge(k)[2]))
	    out.insert(h_inc.edge(k)[1]);
	  else if (cset.contains(h_inc.edge(k)[1]) &&
		   cset.contains(h_inc.edge(k)[2]))
	    out.insert(h_inc.edge(k)[0]);
	}
      }
      index_set possible;
      for (index_type k = n_common; k < ppcv.size(); k++)
	possible.insert(k);
      possible.subtract(out);
      bool ok = false;
      while (!ok) {
	index_set cset1(cset);
	rng.select_fixed_set(rs, ((possible.size() / 2) < (cset.size() / 10) ?
				  (possible.size() / 2) : (cset.size() / 10)),
			     possible.size());
	for (index_type i = 0; i < rs.size(); i++) {
	  index_type c = (n_common + rs[i]);
	  cset1.insert(c);
	}
	std::cerr << cset1.size() << " constraints selected: " << cset1
		  << std::endl;
#ifdef NEW_PCC
	bool cp = test.run(cset1);
#else
	BasicPCC test(ins, stats, ppcv, cset1,
		      prep.inconsistency(), prep.landmark_graph(),
		      &trlm, prep.never_after_graph(), false);
	bool cp = test.run();
#endif
	if (cp) {
	  std::cerr << "selected constraints are not satisfiable!" << std::endl;
	  n_failures += 1;
	  if (n_failures >= failure_limit) {
	    std::cerr << "failure limit exceeded!" << std::endl;
	    return;
	  }
	}
	else {
	  ok = true;
	  cset.assign_copy(cset1);
	}
      }
      assert(ok);
      std::cerr << "success!" << std::endl;
      n_variations += 1;
      std::cerr << "variation #" << n_variations << std::endl;
      for (index_type k = 0; k < cset.size(); k++) {
	std::cerr << " - ";
	ppcv[cset[k]].write(std::cerr, ins);
	std::cerr << std::endl;
      }
      std::ostringstream vname1;
      vname1 << "variation" << n_variations;
      std::cerr << "writing file " << vname1.str() << "..." << std::endl;
      save_constraint_set(ins, ppcv, cset, vname1.str());
      n_failures = 0;
    }
  }
}

class set_value_order {
  ppc_vec& ppcv;
public:
  set_value_order(ppc_vec& v) : ppcv(v) { };
  NTYPE value(const index_set& set) const;
  bool operator()(const index_set& a, const index_set& b) const;
};

NTYPE set_value_order::value(const index_set& set) const
{
  NTYPE sum = 0;
  for (index_type i = 0; i < set.size(); i++) {
    assert(set[i] < ppcv.size());
    sum += ppcv[set[i]].weight;
  }
  return sum;
}

bool set_value_order::operator()(const index_set& a, const index_set& b) const
{
  return (value(a) < value(b));
}

void minimise_conflict
(ppc_vec& ppcv,
 index_set& cset,
 PCC& tester)
{
  index_set s(cset);
  for (index_type k = 0; k < cset.size(); k++) {
    s.subtract(cset[k]);
    bool contradiction = tester.run(s);
    if (!contradiction)
      s.insert(cset[k]);
    n_pcc_tests += 1;
  }
  cset = s;
}

bool cdastar_expand
(ppc_vec& ppcv,
 const index_set& s,
 hypergraph& conflicts,
 index_type next,
 __gnu_cxx::hash_set<index_set, set_hash_function>& v,
 std::priority_queue<index_set, std::vector<index_set>, set_value_order>& q,
 set_value_order& svo,
 PCC& tester)
{
  if (next < conflicts.n_edges()) {
    if (s.contains(conflicts.edge(next))) {
      for (index_type i = 0; i < conflicts.edge(next).size(); i++) {
	index_set new_s(s);
	new_s.subtract(conflicts.edge(next)[i]);
	// std::cerr << "(old conflict) new candidate: ";
	// write_ppc_set(std::cerr, ppcv, new_s);
	// std::cerr << ", value = " << PRINT_NTYPE(svo.value(new_s))
	// 	  << std::endl;
	__gnu_cxx::hash_set<index_set, set_hash_function>::const_iterator t =
	  v.find(new_s);
	if (t == v.end()) {
	  // std::cerr << "new candidate added to queue" << std::endl;
	  v.insert(new_s);
	  q.push(new_s);
	}
      }
      return false;
    }
    else {
      return cdastar_expand(ppcv, s, conflicts, next + 1, v, q, svo, tester);
    }
  }
  // s resolves all known conflicts
  else {
    // std::cerr << "testing ";
    // write_ppc_set(std::cerr, ppcv, s);
    // std::cerr << "..." << std::endl;
    InferenceTracker* tracker = new InferenceTracker(tester);
    tester.set_tracker(tracker);
    bool fail = tester.run(s);
    if (fail) {
      // std::cerr << "contradiction found!" << std::endl;
      tracker->print_proof(std::cerr);
      index_set new_conflict;
      InferenceTracker::assertion a(pc_contradiction);
      bool ok = tracker->extract_proof_premises(a, new_conflict);
      assert(ok);
      assert(new_conflict.size() > 0);
      // std::cerr << "new_conflict = " << new_conflict
      // 		<< ", max = " << ppcv.size() << std::endl;
      new_conflict.remove_greater_than(ppcv.size() - 1);
      std::cerr << "new conflict: ";
      write_ppc_set(std::cerr, ppcv, new_conflict);
      std::cerr << std::endl;
      tester.set_tracker(NULL);
      minimise_conflict(ppcv, new_conflict, tester);
      std::cerr << "minimal conflict: ";
      write_ppc_set(std::cerr, ppcv, new_conflict);
      std::cerr << std::endl;
      // state s should exhibit the new conflict
      assert(s.contains(new_conflict));
      // new conflict should not be one of the known conflicts
      assert(conflicts.is_independent(new_conflict));
      for (index_type i = 0; i < new_conflict.size(); i++) {
	index_set new_s(s);
	new_s.subtract(new_conflict[i]);
	// std::cerr << "new candidate: ";
	// write_ppc_set(std::cerr, ppcv, new_s);
	// std::cerr << ", value = " << PRINT_NTYPE(svo.value(new_s))
	// 	  << std::endl;
	__gnu_cxx::hash_set<index_set, set_hash_function>::const_iterator t =
	  v.find(new_s);
	if (t == v.end()) {
	  // std::cerr << "new candidate added to queue" << std::endl;
	  v.insert(new_s);
	  q.push(new_s);
	}
      }
      conflicts.add_edge(new_conflict);
    }
    tester.set_tracker(NULL);
    n_pcc_tests += 1;
    return !fail;
  }
}

NTYPE cdastar
(Instance& ins,
 Preprocessor& prep,
 ppc_vec& ppcv,
 hypergraph& hinc,
 index_set& soln)
{
  // prepare for PCC
  set_edge_vec trlm;
  pcc_prep_stats.start();
  StaticMutex* mx = prep.inconsistency();
  graph* lmg = prep.necessary_sb_graph();
  landmark_graph_triggered_edges(ins, *lmg, trlm);
  pcc_prep_stats.stop();
  PCC tester(ins, pcc_stats, ppcv, mx, lmg, &trlm, prep.never_after_graph());
  set_hash_function shfn(ppcv.size());
  __gnu_cxx::hash_set<index_set, set_hash_function> v(31337, shfn);
  set_value_order svo(ppcv);
  std::priority_queue<index_set, std::vector<index_set>, set_value_order>
   q(svo);
  index_set s0(ppcv.size(), EMPTYSET); // init s = full set of ppcv.size()
  v.insert(s0);
  q.push(s0);
  while (!q.empty()) {
    index_set s(q.top());
    NTYPE v_s = svo.value(s);
    q.pop();
    std::cerr << "best candidate value = " << PRINT_NTYPE(v_s) << std::endl;
    if (stats.break_signal_raised())
      return v_s;
    bool solved = cdastar_expand(ppcv, s, hinc, 0, v, q, svo, tester);
    if (solved) {
      std::cerr << "solved!" << std::endl;
      soln = s;
      return svo.value(s);
    }
    std::cerr << "|q| = " << q.size() << std::endl;
  }
  return POS_INF;
}

NTYPE incremental_conflicts
(ppc_vec& ppcv,
 ConsistencyTest& tester,
 bool opt_min_conflict,
 hypergraph& h_inc,
 index_set& best)
{
  hg_stats.start();
  NTYPE val = max_value(ppcv, h_inc, best);
  hg_stats.stop();
  bool fail = true;
  while (fail && !stats.break_signal_raised()) {
    std::cerr << "best value = " << PRINT_NTYPE(val) << std::endl;
    std::cerr << "testing ";
    write_ppc_set(std::cerr, ppcv, best);
    std::cerr << "..." << std::endl;
    fail = tester.test(best, true, opt_min_conflict);
    if (fail) {
      const index_set& new_conflict = tester.last_conflict();
      assert(new_conflict.size() > 0);
      std::cerr << "new conflict: ";
      write_ppc_set(std::cerr, ppcv, new_conflict);
      std::cerr << std::endl;
      // state s should exhibit the new conflict
      assert(best.contains(new_conflict));
      // new conflict should not be one of the known conflicts
      assert(h_inc.is_independent(new_conflict));
      h_inc.add_edge(new_conflict);
    }
    hg_stats.start();
    val = max_value(ppcv, h_inc, best);
    hg_stats.stop();
  }
  return val;
}

int main(int argc, char *argv[]) {
  bool     opt_preprocess = true;
  bool     opt_complete_negation = false;
  bool     opt_generate = false;
  index_type opt_generate_md = 2;
  bool     opt_floor = false;
  long     opt_scale_up = 1;
  bool     opt_integrify = false;
  bool     opt_hm_test = false;
  bool     opt_pcc_test = false;
  bool     opt_cda = false;
  bool     opt_icg = false;
  bool     opt_min_conflict = true;
  bool     opt_max = true;
  bool     opt_max_sets = false;
  bool     opt_open = false;
  bool     opt_iterate = false;
  bool     opt_min_inc = false;
  bool     opt_max_sat = false;
  bool     opt_H2 = true;
  bool     opt_H3 = false;
  bool     opt_split = false;
  bool     opt_precheck = false;
  bool     opt_load = true;
  bool     opt_pddl = false;
  bool     opt_pddl_all = false;
  bool     opt_untimed = false;
  bool     opt_save_test = false;
  bool     opt_save_min = false;
  bool     opt_save_max = false;
  bool     opt_save_open = false;
  bool     opt_save_ext = false;
  bool     opt_save_compiled = false;
  bool     opt_compile_to_action_cost = false;
  bool     opt_compile_KG = false;
  NTYPE    opt_slack = 0;
  bool     opt_pcc = false;
  bool     opt_nsb_H2 = false;
  bool     opt_nag_H2 = false;
  bool     opt_cega = false;
  long     cega_time_limit = 0;
  bool     opt_plan = false;
  long     plan_time_limit = 0;
  bool     opt_strong_pcc = false;
  bool     opt_verbose = false;
  bool     opt_print_proof = false;
  bool     opt_print_graph = false;
  index_type opt_min_d = 1;
  index_type opt_max_d = 0;
  unsigned int random_seed = 0;
  long     main_time_limit = 0;
  long     memory_limit = 0;

  PDDL_Base::warning_level = 0;
  PDDL_Base::compile_away_plan_constraints = false;
  PDDL_Base::make_types_from_static_predicates = false;

  Instance::write_PDDL2 = false;
  Instance::write_PDDL3 = false;
  Instance::write_metric = false;
  Instance::write_time = false;
  Instance::neg_atom_name = "not";

  stats.enable_interrupt();

  StringTable symbols(50, lowercase_map);
  Parser* reader = new Parser(symbols);

  for (int k = 1; k < argc; k++) {
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      Instance::default_trace_level = atoi(argv[++k]);
      Preprocessor::default_trace_level = Instance::default_trace_level;
      opt_verbose = true;
    }
    else if (strcmp(argv[k],"-proof") == 0) {
      opt_print_proof = true;
    }
    else if (strcmp(argv[k],"-no-prep") == 0) {
      opt_preprocess = false;
    }
    else if (strcmp(argv[k],"-complete-negation") == 0) {
      opt_complete_negation = true;
    }
    else if (strcmp(argv[k],"-generate") == 0) {
      opt_generate = true;
      Preprocessor::prep_remove_useless_actions = false;
    }
    else if ((strcmp(argv[k],"-generate-md") == 0) && (k < argc - 1)) {
      opt_generate = true;
      Preprocessor::prep_remove_useless_actions = false;
      opt_generate_md = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-floor") == 0) {
      opt_floor = true;
    }
    else if ((strcmp(argv[k],"-scale-up") == 0) && (k < argc - 1)) {
      opt_scale_up = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-integrify") == 0) {
      opt_integrify = true;
    }
    else if ((strcmp(argv[k],"-catc") == 0) && (k < argc - 1)) {
      HSPS::PDDL_Name::catc = *argv[++k];
    }
    else if (strcmp(argv[k],"-3") == 0) {
      opt_H3 = true;
    }
    else if (strcmp(argv[k],"-1") == 0) {
      opt_H2 = false;
    }
    // input options
    else if (strcmp(argv[k],"-load") == 0) {
      opt_load = true;
    }
    else if (strcmp(argv[k],"-no-load") == 0) {
      opt_load = false;
    }
    else if (strcmp(argv[k],"-untimed") == 0) {
      opt_untimed = true;
    }
    // hm-test and related options
    else if (strcmp(argv[k],"-hm-test") == 0) {
      opt_hm_test = true;
      opt_pddl = true; // make sure we output the result!
    }
    else if ((strcmp(argv[k],"-d-min") == 0) && (k < argc - 1)) {
      opt_min_d = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-d-max") == 0) && (k < argc - 1)) {
      opt_max_d = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-split") == 0) {
      opt_split = true;
    }
    else if (strcmp(argv[k],"-precheck") == 0) {
      opt_precheck = true;
    }
    else if (strcmp(argv[k],"-pcc-test") == 0) {
      opt_pcc_test = true;
      opt_pddl = true;
    }
    else if (strcmp(argv[k],"-cda") == 0) {
      opt_cda = true;
      opt_pddl = true;
    }
    else if (strcmp(argv[k],"-icg") == 0) {
      opt_icg = true;
      opt_pddl = true;
    }
    else if (strcmp(argv[k],"-non-min-conflict") == 0) {
      opt_min_conflict = false;
    }
    else if (strcmp(argv[k],"-min-inc-hm") == 0) {
      opt_hm_test = true;
      opt_min_inc = true;
    }
    else if (strcmp(argv[k],"-max-sat-hm") == 0) {
      opt_hm_test = true;
      opt_max_sat = true;
    }
    else if (strcmp(argv[k],"-sb2") == 0) {
      opt_nsb_H2 = true;
    }
    else if (strcmp(argv[k],"-na2") == 0) {
      opt_nag_H2 = true;
    }
    // maximal consistent sets and related options
    else if (strcmp(argv[k],"-no-max") == 0) {
      opt_max = false;
    }
    else if (strcmp(argv[k],"-max-sets") == 0) {
      opt_max_sets = true;
    }
    else if ((strcmp(argv[k],"-slack") == 0) && (k < argc - 1)) {
      opt_slack = A_TO_N(argv[++k]);
    }
    else if (strcmp(argv[k],"-max-slack") == 0) {
      opt_slack = POS_INF;
    }
    // open sets
    else if (strcmp(argv[k],"-open") == 0) {
      opt_open = true;
    }
    else if (strcmp(argv[k],"-iterate") == 0) {
      opt_iterate = true;
    }
    else if (strcmp(argv[k],"-plan") == 0) {
      opt_plan = true;
    }
    else if ((strcmp(argv[k],"-plan-time") == 0) && (k < argc - 1)) {
      plan_time_limit = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-max-sat-plan") == 0) {
      opt_plan = true;
      opt_max_sat = true;
    }
    else if (strcmp(argv[k],"-min-inc-plan") == 0) {
      opt_plan = true;
      opt_min_inc = true;
    }
    // output options
    else if (strcmp(argv[k],"-pddl") == 0) {
      opt_pddl = true;
    }
    else if (strcmp(argv[k],"-no-pddl") == 0) {
      opt_pddl = false;
    }
    else if (strcmp(argv[k],"-all") == 0) {
      opt_pddl_all = true;
    }
    else if (strcmp(argv[k],"-pddl-all") == 0) {
      opt_pddl = true;
      opt_pddl_all = true;
    }
    // save/format options
    else if ((strcmp(argv[k],"-save-test") == 0) ||
	     (strcmp(argv[k],"-save-tests") == 0)) {
      opt_save_test = true;
    }
    else if ((strcmp(argv[k],"-save-ext") == 0) ||
	     (strcmp(argv[k],"-save-exts") == 0)) {
      opt_save_ext = true;
    }
    else if (strcmp(argv[k],"-save-max") == 0) {
      opt_max_sets = true;
      opt_save_max = true;
    }
    else if (strcmp(argv[k],"-save-open") == 0) {
      opt_save_open = true;
    }
    else if (strcmp(argv[k],"-save-min") == 0) {
      opt_save_min = true;
    }
    else if (strcmp(argv[k],"-save-compiled") == 0) {
      opt_save_compiled = true;
    }
    else if (strcmp(argv[k],"-save-split") == 0) {
      save_split = true;
    }
    else if (strcmp(argv[k],"-no-dkel") == 0) {
      HSPS::Instance::write_DKEL = false;
    }
    else if (strcmp(argv[k],"-no-extra") == 0) {
      HSPS::Instance::write_extra = false;
    }
    else if (strcmp(argv[k],"-no-negation") == 0) {
      HSPS::Instance::write_negation = false;
    }
    else if ((strcmp(argv[k],"-file-limit") == 0) && (k < argc - 1)) {
      file_write_limit = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-cac") == 0) {
      opt_compile_to_action_cost = true;
    }
    else if (strcmp(argv[k],"-kg") == 0) {
      opt_compile_KG = true;
    }
    // pcc test
    else if (strcmp(argv[k],"-pcc") == 0) {
      opt_pcc = true;
      opt_pddl = true;
    }
    else if (strcmp(argv[k],"-strong-pcc") == 0) {
      opt_strong_pcc = true;
      opt_pcc = true;
      opt_pddl = true;
    }
    else if (strcmp(argv[k],"-min-inc-pcc") == 0) {
      opt_pcc = true;
      opt_min_inc = true;
    }
    else if (strcmp(argv[k],"-max-sat-pcc") == 0) {
      opt_pcc = true;
      opt_max_sat = true;
    }
    else if (strcmp(argv[k],"-print-graph") == 0) {
      opt_print_graph = true;
    }
    // cega option
    else if (strcmp(argv[k],"-cega") == 0) {
      opt_cega = true;
      opt_pddl = true;
    }
    else if ((strcmp(argv[k],"-cega-time") == 0) && (k < argc - 1)) {
      cega_time_limit = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-min-inc-cega") == 0) {
      opt_cega = true;
      opt_min_inc = true;
    }
    else if (strcmp(argv[k],"-max-sat-cega") == 0) {
      opt_cega = true;
      opt_max_sat = true;
    }
    // misc. options
    else if ((strcmp(argv[k],"-t") == 0) && (k < argc - 1)) {
      main_time_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-y") == 0) && (k < argc - 1)) {
      memory_limit = atoi(argv[++k]);
    }
    else if (((strcmp(argv[k],"-rnd") == 0) ||
	      (strcmp(argv[k],"-r") == 0)) &&
	     (k < argc - 1)) {
      random_seed = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-dba-semantics") == 0) {
      HSPS::PDDL_Base::del_before_add_semantics = true;
    }
    else if (*argv[k] != '-') {
      reader->read(argv[k], false);
    }
  }

  if (main_time_limit > 0) {
    stats.enable_time_out(main_time_limit);
  }
  if (memory_limit > 0) {
    stats.enable_memory_limit(memory_limit);
  }
  stats.start();

  Instance instance;
  ppc_vec  ppcs;

  reader->post_process();

  reader->instantiate(instance);
  index_set hard_goals;
  for (index_type k = 0; k < instance.n_atoms(); k++)
    if (instance.atoms[k].goal)
      hard_goals.insert(k);

  if (!opt_generate) {
    bool ok = extract_ppcs(reader, instance, ppcs);
    if (!ok) {
      std::cerr << "error extracting ppcs!" << std::endl;
      exit(255);
    }
  }

  if (opt_scale_up > 1) {
    for (index_type k = 0; k < ppcs.length(); k++)
      ppcs[k].weight = (ppcs[k].weight * opt_scale_up);
  }
  if (opt_floor) {
    for (index_type k = 0; k < ppcs.length(); k++)
      ppcs[k].weight = FLOOR(ppcs[k].weight);
  }

  Preprocessor prep(instance, stats);

  if (opt_preprocess) {
    prep.preprocess();
    for (index_type k = 0; k < ppcs.length(); k++) {
      instance.remap_set(ppcs[k].s_c, prep.atom_map);
      instance.remap_set(ppcs[k].s_t, prep.atom_map);
    }
    instance.remap_set(hard_goals, prep.atom_map);
  }
  else {
    instance.cross_reference();
  }
  if (opt_complete_negation) {
    std::cerr << "adding atom negations..." << std::endl;
    instance.complete_atom_negations();
    if (!instance.cross_referenced()) {
      std::cerr << "re-cross referencing..." << std::endl;
      instance.cross_reference();
    }
  }
  instance.set_goal(hard_goals);

  std::cerr << "instantiation and preprocessing finished in "
	    << stats.time() << " seconds" << std::endl;

  // graph* g0 = prep.landmark_graph();
  // graph g1;
  // prep.nc_landmark_graph(g1);
  // for (index_type i = 0; i < instance.init_atoms.size(); i++)
  //   for (index_type j = 0; j < instance.n_atoms(); j++)
  //     if (!instance.atoms[j].init) {
  // 	assert(instance.init_atoms[i] != j);
  // 	if (!g1.adjacent(instance.init_atoms[i], j))
  // 	  g1.add_edge(instance.init_atoms[i], j);
  //     }
  // for (index_type i = 0; i < instance.n_atoms(); i++)
  //   for (index_type j = 0; j < instance.n_atoms(); j++)
  //     if (g0->adjacent(i, j) && !g1.adjacent(i, j)) {
  // 	std::cerr << "A-not-B: " << instance.atoms[i].name << " -> "
  // 		  << instance.atoms[j].name << std::endl;
  //     }
  //     else if (!g0->adjacent(i, j) && g1.adjacent(i, j)) {
  // 	std::cerr << "B-not-A: " << instance.atoms[i].name << " -> "
  // 		  << instance.atoms[j].name << std::endl;
  //     }

  index_type n_common = 0;
  if (opt_generate) {
    if (reader->n_plans() == 0) {
      std::cerr << "error: can't generate constraints with no input plans!"
		<< std::endl;
      exit(255);
    }
    n_common = generate_constraints(reader, instance, prep, ppcs,
				    prep.action_map, opt_untimed,
				    opt_generate_md);
    std::cerr << "generated " << n_common << " satisfied and "
	      << (ppcs.size() - n_common) << " violated constraints"
	      << std::endl;
  }
  else {
    // this is only relevant if we're loading a set of previously
    // generated constraints.
    index_type k = 0;
    bool done = (k >= ppcs.size());
    while (!done) {
      std::string s = ppcs[k].name->to_string();
      std::string s1(s, 0, 3);
      if (s1 == std::string("not")) {
	n_common = k;
	done = true;
      }
      else {
	k += 1;
	done = (k >= ppcs.size());
      }
    }
  }

  NTYPE v_sum = 0;
  for (index_type k = 0; k < ppcs.length(); k++)
    v_sum += ppcs[k].weight;
  NTYPE v_max = v_sum;
  NTYPE v_max_sat = 0;

  std::cerr << ppcs.length() << " preferences over plan constraints"
	    << std::endl;
  std::cerr << "sum penalty = " << v_sum << std::endl;
  std::cerr << instance.n_atoms() << " atoms and "
	    << instance.n_actions() << " actions"
	    << std::endl;

  if (opt_verbose)
    for (index_type k = 0; k < ppcs.length(); k++) {
      std::cerr << k << ". ";
      ppcs[k].write(std::cerr, instance);
      //std::cerr << std::endl;
    }

#ifdef BUILD_WITH_SOFT
  if (opt_compile_to_action_cost) {
    // make a soft instance by compiling ppc's to soft goals
    SoftInstance* s_ins = new SoftInstance(instance.name);
    s_ins->copy(instance);
    s_ins->hard = hard_goals;
    s_ins->cross_reference();
    mapping map(s_ins->n_actions());
    for (index_type k = 0; k < ppcs.length(); k++) {
      std::cerr << "compiling trajectory constraint ";
      ppcs[k].write(std::cerr, instance);
      std::cerr << "..." << std::endl;
      index_type pk = ppcs[k].compile(*s_ins, map);
      SoftInstance::SoftGoal& sgk = s_ins->new_soft_goal();
      sgk.name = ppcs[k].name;
      sgk.atoms.assign_singleton(pk);
      sgk.weight = ppcs[k].weight;
      Stopwatch::seconds();
      std::cerr << "now " << s_ins->n_actions() << " actions, "
		<< Stopwatch::peak_memory() << "k" << std::endl;
    }
    s_ins->null_value = 0;
    if (!s_ins->cross_referenced())
      s_ins->cross_reference();

    long divr = (opt_integrify ? s_ins->integrify_weights() : 1);

    //  Instance::write_PDDL2 = true;
    //  Instance::write_PDDL3 = true;
    //  Instance::write_metric = true;
    //  std::cerr << s_ins->n_hard() << " hard and " << s_ins->n_soft()
    //     << " soft goals, max action cost = " << s_ins->max_cost
    //     << std::endl;
    //  s_ins->write_problem_header(std::cerr);
    //  s_ins->write_problem_init(std::cerr);
    //  s_ins->write_problem_goal(std::cerr);
    //  s_ins->write_problem_metric(std::cerr);
    //  std::cerr << ")" << std::endl;
    //  Instance::write_PDDL3 = false;

    // make a regular instance by compiling away the soft goals
    Instance* c_ins = new Instance(instance.name);
    if (opt_compile_KG) {
      s_ins->compile_KG(*c_ins);
    }
    else {
      s_ins->compile_direct_cost(*c_ins);
    }
    c_ins->cross_reference();
    c_ins->assign_unique_action_names(false);

    // c_ins->print(std::cerr);

    // set options for "nice" output
    Instance::write_negation = false;
    Instance::write_PDDL2 = true;
    Instance::write_metric = true;
    Instance::write_time = false;
    Instance::write_DKEL = false;
    Instance::write_extra = false;
    Instance::always_write_parameters = true;
    Instance::always_write_requirements = true;
    Instance::always_write_precondition = true;
    Instance::always_write_effect = false;
    Instance::always_write_conjunction = true;

    // save domain
    char* b = reader->problem_file_basename();
    std::ostringstream fname1;
    if (b) fname1 << b << "-";
    fname1 << "ac.domain.pddl";
    std::cerr << "writing file " << fname1.str() << "..." << std::endl;
    std::ofstream s1_out(fname1.str().c_str());
    if (opt_integrify) {
      s1_out << ";; divisor = " << divr << std::endl;
    }
    c_ins->write_domain_header(s1_out);
    c_ins->write_domain_declarations(s1_out);
    c_ins->write_domain_actions(s1_out);
    c_ins->write_domain_DKEL_items(s1_out);
    s1_out << ")";
    s1_out.close();

    // save problem
    std::ostringstream fname2;
    if (b) fname2 << b << "-";
    fname2 << "ac.problem.pddl";
    std::cerr << "writing file " << fname2.str() << "..." << std::endl;
    std::ofstream s2_out(fname2.str().c_str());
    c_ins->write_problem(s2_out);
    s2_out.close(); 

    delete c_ins;
    delete s_ins;
    exit(0);
  }
#endif

  // hypergraph over ppcs
  hypergraph h_inc(ppcs.length());
  index_set_vec min_inc;

  if (opt_load) {
    std::cerr << "loading inconsistent sets..." << std::endl;
    PDDL_Base::warning_level = (opt_verbose ? 1 : 0);
    name_vec ppc_name(0, ppcs.length());
    for (index_type k = 0; k < ppcs.length(); k++)
      ppc_name[k] = ppcs[k].name;
    index_set_vec h_load;
    reader->export_sets(ppc_name, h_load);
    std::cerr << h_load.length() << " sets loaded" << std::endl;
    for (index_type k = 0; k < h_load.length(); k++)
      min_inc.insert_minimal(h_load[k]);
    std::cerr << min_inc.length()
	      << " minimal inconsistent sets read from input"
	      << std::endl;
  }

  for (index_type k = 0; k < min_inc.length(); k++)
    if (h_inc.is_independent(min_inc[k]))
      h_inc.add_edge(min_inc[k]);

  // maximal known-to-be-satisfiable sets of ppcs, and corresponding plans
  index_set_vec max_sat;
  plan_set_vec max_plan;

  if (reader->n_plans() > 0) {
    std::cerr << "analysing " << reader->n_plans() << " plans..." << std::endl;
    analyse_input_plans(reader, instance, ppcs, prep.action_map,
			opt_untimed, max_sat, max_plan, opt_verbose);
    if (opt_verbose)
      std::cerr << "maximal satisfied sets of constraints:" << std::endl;
    for (index_type k = 0; k < max_sat.length(); k++) {
      NTYPE val = ppc_set_value(ppcs, max_sat[k]);
      if (val > v_max_sat) v_max_sat = val;
      if (opt_verbose) {
	std::cerr << k + 1 << ": ";
	write_ppc_set(std::cerr, ppcs, max_sat[k]);
	std::cerr << " (value: " << PRINT_NTYPE(val) << ", plans:";
	for (index_type i = 0; i < max_plan[k].size(); i++)
	  std::cerr << " " << max_plan[k][i]->plan_name();
	std::cerr << ")" << std::endl;
      }
    }
    std::cerr << "max plan value = " << PRINT_NTYPE(v_max_sat)
	      << std::endl
	      << "min plan penalty = " << PRINT_NTYPE(v_sum - v_max_sat)
	      << std::endl;
  }

  set_edge_vec trlm;

  // test PCC (int. or ext.) with known cases
  if (opt_pcc && (opt_min_inc || opt_max_sat)) {
    if (opt_pcc && !opt_strong_pcc) {
      std::cerr << "preparing for PCC..." << std::endl;
      pcc_prep_stats.start();
      StaticMutex* mx = prep.inconsistency();
      graph* lmg = prep.necessary_sb_graph(opt_nsb_H2);
      landmark_graph_triggered_edges(instance, *lmg, trlm);
      if (opt_verbose) {
	std::cerr << "necessary SB's:" << std::endl;
	for (index_type i = 0; i < instance.n_atoms(); i++)
	  for (index_type j = 0; j < instance.n_atoms(); j++)
	    if (lmg->adjacent(i, j)) {
	      std::cerr << i << "." << instance.atoms[i].name << " SB "
			<< j << "." << instance.atoms[j].name << std::endl;
	    }
	std::cerr << trlm.size() << " triggered landmarks:" << std::endl;
	for (index_type k = 0; k < trlm.size(); k++) {
	  instance.write_action_set(std::cout, trlm[k].first);
	  std::cerr << " => (" << instance.atoms[trlm[k].second.first].name
		    << " -> " << instance.atoms[trlm[k].second.second].name
		    << ")" << std::endl;
	}
      }
      graph* nag = prep.never_after_graph(opt_nag_H2);
      if (opt_verbose) {
	std::cerr << "never after:" << std::endl;
	for (index_type i = 0; i < instance.n_atoms(); i++)
	  for (index_type j = 0; j < instance.n_atoms(); j++)
	    if (nag->adjacent(i, j)) {
	      std::cerr << j << "." << instance.atoms[j].name << " NA "
			<< i << "." << instance.atoms[i].name << std::endl;
	    }
      }
      pcc_prep_stats.stop();
      std::cerr << "finished in " << pcc_prep_stats.total_time()
		<< " seconds" << std::endl;
      std::cerr << trlm.size() << " potential triggered landmarks"
		<< std::endl;
    }
#ifdef NEW_PCC
    PCC test(instance, pcc_stats, ppcs,
	     prep.inconsistency(), prep.necessary_sb_graph(opt_nsb_H2),
	     &trlm, prep.never_after_graph(opt_nag_H2));
#endif
    for (index_type k = 0; (k < min_inc.length()) && opt_min_inc; k++) {
#ifdef DO_TEST_PCC
      std::cout << "testing inconsistent set ";
      write_ppc_set(std::cout, ppcs, min_inc[k]);
      std::cout << std::endl;
      bool cp = false;
#ifdef NEW_PCC
      if (opt_print_proof) {
	InferenceTracker* proof = new InferenceTracker(test);
	proof->set_verbose(opt_verbose);
	proof->set_debug_mode(opt_verbose);
	test.set_tracker(proof);
      }
      cp = test.run(min_inc[k]);
#else
      if (opt_strong_pcc) {
	StrongPCC test(instance, pcc_stats, ppcs, min_inc[k],
		       opt_H2, true, opt_verbose);
	cp = test.run();
      }
      else {
	BasicPCC test(instance, pcc_stats, ppcs, min_inc[k],
		      prep.inconsistency(), prep.landmark_graph(),
		      &trlm, prep.never_after_graph(), opt_verbose);
	cp = test.run();
	if (cp && opt_print_graph) {
	  test.print_graph(prep.landmark_graph(), std::cout);
	}
      }
#endif // NEW_PCC
      if (cp) {
	std::cout << "contradiction proven" << std::endl;
	n_pcc_proven += 1;
#ifdef NEW_PCC
	if (test.tracker()) {
	  test.tracker()->print_proof(std::cout);
	}
#endif
      }
      else {
	std::cout << "contradiction not proven" << std::endl;
      }
#ifdef NEW_PCC
      test.set_tracker(NULL);
#endif      
      n_pcc_tests += 1;
#else
      write_ppc_set_as_pcc(std::cout, instance, ppcs, min_inc[k], false, 0, 0);
#endif
    }
    for (index_type k = 0; (k < max_sat.length()) && opt_max_sat; k++) {
#ifdef DO_TEST_PCC
      std::cout << "testing consistent set ";
      write_ppc_set(std::cout, ppcs, max_sat[k]);
      std::cout << std::endl;
      bool cp = false;
#ifdef NEW_PCC
      if (opt_print_proof) {
	InferenceTracker* proof = new InferenceTracker(test);
	proof->set_verbose(opt_verbose);
	proof->set_debug_mode(opt_verbose);
	test.set_tracker(proof);
      }
      cp = test.run(max_sat[k]);
#else
      if (opt_strong_pcc) {
	StrongPCC test(instance, pcc_stats, ppcs, max_sat[k],
		       opt_H2, true, opt_verbose);
	cp = test.run();
      }
      else {
	BasicPCC test(instance, pcc_stats, ppcs, max_sat[k],
		      prep.inconsistency(), prep.landmark_graph(),
		      &trlm, prep.never_after_graph(), opt_verbose);
	cp = test.run();
      }
#endif
      if (cp) {
	std::cout << "error: contradiction found!" << std::endl;
#ifdef NEW_PCC
	if (test.tracker()) {
	  test.tracker()->print_proof(std::cout);
	}
#endif
	exit(255);
      }
      else {
	std::cout << "ok, no contradiction" << std::endl;
      }
#ifdef NEW_PCC
      test.set_tracker(NULL);
#endif      
      n_pcc_tests += 1;
#else
      write_ppc_set_as_pcc(std::cout, instance, ppcs, max_sat[k], false, 0, 0);
#endif
    }
#ifdef DO_TEST_PCC
    std::cout << n_pcc_tests << " PCC tests in "
	      << pcc_prep_stats.total_time() << " (prep) + "
	      << pcc_stats.total_time() << " seconds (avg. "
	      << pcc_stats.total_time() / n_pcc_tests
	      << " sec/test)" << std::endl;
    std::cout << n_pcc_proven << " of " << min_inc.length()
	      << " contradictions proven" << std::endl;
#endif
    exit(0);
  }

  // test CEGA with known cases
  if (opt_cega && (opt_min_inc || opt_max_sat)) {
    count_type n_pos = 0;
    count_type n_neg = 0;
    if (opt_min_inc)
      for (index_type k = 0; (k < min_inc.length()) &&
	     !stats.break_signal_raised(); k++) {
	n_cega_tests += 1;
	Schedule* soln = 0;
	bool dec = cega(instance, cega_stats, cega_time_limit,
			ppcs, min_inc[k], *(prep.inconsistency()),
			*(prep.landmark_graph()), true, soln);
	if (dec) {
	  n_cega_decided += 1;
	  n_neg += 1;
	}
	if (dec && (soln != 0)) {
	  std::cerr << "error: plan found for inconsistent set!" << std::endl;
	  exit(255);
	}
      }
    if (opt_max_sat)
      for (index_type k = 0; (k < max_sat.length()) &&
	     !stats.break_signal_raised(); k++) {
	n_cega_tests += 1;
	Schedule* soln = 0;
	bool dec = cega(instance, cega_stats, cega_time_limit, ppcs,
			max_sat[k], *(prep.inconsistency()),
			*(prep.landmark_graph()), true, soln);
	if (dec) {
	  n_cega_decided += 1;
	  n_pos += 1;
	}
	if (dec && (soln == 0)) {
	  std::cerr << "error: consistent set unsolvable!" << std::endl;
	  exit(255);
	}
      }
    std::cout << n_cega_tests << " CEGA tests in "
	      << cega_stats.total_time() << " seconds (avg. "
	      << cega_stats.total_time() / n_cega_tests
	      << " sec/test), "
	      << n_cega_decided << " of " << n_cega_tests
	      << " decided (" << n_neg << " negative, "
	      << n_pos << " positive)" << std::endl;
    exit(0);
  }

  // test h^m on known cases
  if (opt_hm_test && (opt_max_sat || opt_min_inc)) {
    if (opt_max_sat) {
      for (index_type k = 0; k < max_sat.length(); k++) {
	std::cout << "testing consistent set ";
	write_ppc_set(std::cout, ppcs, max_sat[k]);
	std::cout << std::endl;
	n_hm_tests += 1;
	bool inc = hm_test(instance, ppcs, max_sat[k], hm_test_stats,
			   opt_H3, opt_H2, opt_verbose);
	if (inc) {
	  std::cerr << "error: set proven unsatisfiable!"
		    << std::endl;
	  exit(255);
	}
      }
    }
    count_type n_hm_proven = 0;
    if (opt_min_inc) {
      for (index_type k = 0; k < min_inc.length(); k++) {
	n_hm_tests += 1;
	bool inc = hm_test(instance, ppcs, min_inc[k], hm_test_stats,
			   opt_H3, opt_H2, opt_verbose);
	if (inc) n_hm_proven += 1;
      }
    }
    std::cout << n_hm_tests << " h^m tests in "
	      << hm_test_stats.total_time() << " seconds (avg. "
	      << hm_test_stats.total_time() / n_hm_tests
	      << " sec/test), " << n_hm_proven << " proven inconsistent"
	      << std::endl;
    exit(0);
  }

  // test planner with known cases
  if (opt_plan && (opt_max_sat || opt_min_inc)) {
    if (opt_max_sat) {
      for (index_type k = 0; k < max_sat.length(); k++) {
	n_planner_calls += 1;
	Schedule* soln = 0;
	bool solved = plan(instance, ppcs, max_sat[k], plan_stats,
			   plan_time_limit, false, soln);
	if (solved) n_plans_found += 1;
      }
    }
    if (opt_min_inc) {
      for (index_type k = 0; k < min_inc.length(); k++) {
	n_planner_calls += 1;
	Schedule* soln = 0;
	bool solved = plan(instance, ppcs, min_inc[k], plan_stats,
			   plan_time_limit, false, soln);
	if (solved) n_plans_found += 1;
      }
    }
    std::cout << n_planner_calls << " planner calls in "
	      << plan_stats.total_time() << " seconds (avg. "
	      << plan_stats.total_time() / n_planner_calls
	      << " sec/test), "
	      << n_plans_found << " of " << n_planner_calls
	      << " solved" << std::endl;
    exit(0);
  }

  // h^m & pcc test-related stuff

  // first, put all constraints in the ca set (because pcc uses only that)
  index_set ce;
  index_set ca;
  ca.fill(ppcs.length());

  if (opt_cda) {
    std::cerr << "running CDA*..." << std::endl;
    index_set soln;
    NTYPE best = cdastar(instance, prep, ppcs, h_inc, soln);
    std::cerr << "best candidate value = " << best << std::endl;
    std::cerr << h_inc.n_edges() << " inconsistent sets found ("
	      << n_pcc_tests << " tests, " << pcc_stats << ")"
	      << std::endl;
  }

  else if (opt_icg) {
    std::cerr << "running ICG..." << std::endl;
    index_set s_max;
    ConsistencyTest* tester;
    if (opt_cega) {
      CEGATest* cega_tester =
	new CEGATest(ppcs, instance, cega_stats, n_cega_tests,
		     *(prep.inconsistency()), *(prep.landmark_graph()));
      cega_tester->verbose = true;
      tester = cega_tester;
    }
    else if (opt_hm_test) {
      hmTest* hm_tester =
	new hmTest(ppcs, instance, hm_test_stats, n_hm_built, n_hm_tests);
      hm_tester->verbose = true;
      hm_tester->opt_H2 = opt_H2;
      if (opt_pcc_test) {
	PCCTest* pcc_tester =
	  new PCCTest(ppcs, instance, prep, pcc_prep_stats, pcc_stats,
		      n_pcc_tests);
	pcc_tester->verbose1 = true;
	pcc_tester->verbose2 = opt_verbose;
	tester = new SerialTest(ppcs, *pcc_tester, *hm_tester);
      }
      else {
	tester = hm_tester;
      }
    }
    else {
      PCCTest* pcc_tester =
	new PCCTest(ppcs, instance, prep, pcc_prep_stats, pcc_stats,
		    n_pcc_tests);
      pcc_tester->verbose1 = true;
      pcc_tester->verbose2 = opt_verbose;
      tester = pcc_tester;
    }
    NTYPE v_max = incremental_conflicts(ppcs, *tester, opt_min_conflict,
					h_inc, s_max);
    std::cerr << "best candidate: ";
    write_ppc_set(std::cerr, ppcs, s_max);
    std::cerr << std::endl;
    std::cerr << "best candidate value = " << v_max << std::endl;
    std::cerr << h_inc.n_edges() << " inconsistent sets found ("
	      << n_pcc_tests << " tests, " << pcc_stats << ")"
	      << std::endl;
  }

  if ((opt_max_d >= opt_min_d) && opt_pcc_test) {
    std::cerr << "preparing for PCC..." << std::endl;
    pcc_prep_stats.start();
    StaticMutex* mx = prep.inconsistency();
    graph* lmg = prep.necessary_sb_graph(opt_nsb_H2);
    landmark_graph_triggered_edges(instance, *lmg, trlm);
    graph* nag = prep.never_after_graph(opt_nag_H2);
    pcc_prep_stats.stop();
    std::cerr << "finished in " << pcc_prep_stats.total_time()
	      << " seconds" << std::endl;
    std::cerr << trlm.size() << " potential triggered landmarks"
	      << std::endl;

    pcc_stats.start();
#ifdef NEW_PCC
    PCC test(instance, pcc_stats, ppcs,
	     prep.inconsistency(), prep.necessary_sb_graph(opt_nsb_H2),
	     &trlm, prep.never_after_graph(opt_nag_H2));
#endif
    index_type d = opt_min_d;
    while ((d <= opt_max_d) && (d <= ca.length())) {
      std::cerr << "checking " << d << "-consistency..." << std::endl;
      mSubsetEnumerator e(ca.length(), d);
      bool more = e.first();
      while (more) {
	assert(e.current_set_size() == d);
	index_set sel;
	e.current_set(ca, sel);
	if (h_inc.is_independent(sel) &&
	    (max_sat.first_superset(sel) == no_such_index)) {
	  if (opt_verbose) {
	    std::cerr << "testing ";
	    write_ppc_set(std::cerr, ppcs, sel);
	    std::cerr << " (" << h_inc.n_edges() << " inc. sets, "
		      << n_pcc_tests << " tests, "
		      << pcc_stats.time() << " sec., "
		      << pcc_stats.peak_memory()
		      << "k)..." << std::endl;
	  }
#ifdef NEW_PCC
	  if (opt_print_proof) {
	    InferenceTracker* proof = new InferenceTracker(test);
	    proof->set_verbose(opt_verbose);
	    test.set_tracker(proof);
	  }
	  bool cp = test.run(sel);
#else
	  BasicPCC test(instance, pcc_stats, ppcs, sel,
			prep.inconsistency(), prep.landmark_graph(),
			&trlm, prep.never_after_graph(), opt_verbose);
	  bool cp = test.run();
#endif
	  if (cp) {
	    std::cerr << "set " << sel << " = ";
	    write_ppc_set(std::cerr, ppcs, sel);
	    std::cerr << " is unachievable" << std::endl;
	    h_inc.add_edge(sel);
#ifdef NEW_PCC
	    if (test.tracker()) {
	      test.tracker()->print_proof(std::cout);
	    }
#endif
	  }
#ifdef NEW_PCC
	  test.set_tracker(NULL);
#endif
	  n_pcc_tests += 1;
	}
	more = e.next();
      }
      std::cerr << h_inc.n_edges() << " inconsistent sets found ("
		<< n_pcc_tests << " tests, " << pcc_stats << ")"
		<< std::endl;
      d += 1;
    }
    pcc_stats.stop();
  }

  // if split option enabled, now we can move E-constraints to ce.
  if (opt_split) {
    for (index_type k = 0; k < ppcs.length(); k++) {
      if (ppcs[k].pct == pc_sometime)
	ce.insert(k);
    }
    ca.subtract(ce);
    std::cerr << ce.length() << " type-E constraints: " << ce << std::endl;
    std::cerr << ca.length() << " type-A constraints: " << ca << std::endl;
  }

  if ((opt_max_d >= opt_min_d) && opt_hm_test) {
    hm_test_stats.start();
    if (opt_precheck && (opt_H2 || opt_H3)) {
      std::cerr << "pre-checking..." << std::endl;
      Heuristic* h0 = 0;
      if (opt_H3) {
	CostTable* hm = new CostTable(instance, hm_test_stats);
	hm->compute_H3(ZeroACF());
	h0 = hm;
      }
      else {
	StaticMutex* mx = new StaticMutex(instance);
	h0 = mx;
      }
      for (index_type i = 0; i < ppcs.length(); i++)
	if (ppcs[i].pct == pc_always)
	  for (index_type j = 0; j < ppcs.length(); j++)
	    if (ppcs[j].pct == pc_sometime) {
	      index_set s(ppcs[i].s_c);
	      s.insert(ppcs[j].s_c);
	      n_hm_tests += 1;
	      if (INFINITE(h0->eval(s))) {
		index_set sel;
		sel.insert(i);
		sel.insert(j);
		std::cerr << "set " << sel << " = ";
		write_ppc_set(std::cerr, ppcs, sel);
		std::cerr << " is unachievable" << std::endl;
		h_inc.add_edge(sel);
	      }
	    }
      delete h0;
      n_hm_built += 1;
    }

    if (ce.length() > 0) {
      std::cerr << "checking type-E constraints..." << std::endl;
      Instance* test_ins = instance.copy();
      index_vec ace(no_such_index, ce.length());
      for (index_type k = 0; k < ce.length(); k++) {
	ace[k] =
	  test_ins->compile_pc_sometime_conjunction(ppcs[ce[k]].s_c, ppcs[ce[k]].name);
	test_ins->atoms[ace[k]].goal = false;
	std::cerr << k + 1 << " of " << ce.length() << " compiled, "
		  << test_ins->n_actions() << " actions" << std::endl;
      }
      test_ins->cross_reference();
      Heuristic* h0 = 0;
      if (opt_H3) {
	CostTable* hm = new CostTable(*test_ins, hm_test_stats);
	hm->compute_H3(ZeroACF());
	h0 = hm;
      }
      else if (opt_H2) {
	StaticMutex* mx = new StaticMutex(*test_ins);
	h0 = mx;
      }
      else {
	//CostTable* hm = new CostTable(*test_ins, hm_test_stats);
	//hm->compute_H1(ZeroACF());
	//h0 = hm;
	h0 = new Reachability(*test_ins);
      }
      bool_vec test_pos(false, ce.length());
      for (index_type k = 0; k < ce.length(); k++) {
	if (opt_verbose) {
	  std::cerr << "checking " << ce[k] << "..." << std::endl;
	}
	n_hm_tests += 1;
	if (INFINITE(h0->incremental_eval(test_ins->goal_atoms, ace[k]))) {
	  std::cerr << "constraint " << ppcs[ce[k]].name
		    << " is unachievable" << std::endl;
	  test_pos[k] = true;
	  h_inc.add_singleton_edge(ce[k]);
	}
      }
      if ((opt_H2 || opt_H3) && (opt_max_d > 1)) {
	for (index_type k = 0; k < ce.length(); k++)
	  if (!test_pos[k]) {
	    index_set g1(test_ins->goal_atoms);
	    g1.insert(ace[k]);
	    for (index_type i = k + 1; i < ce.length(); i++) {
	      if (opt_verbose) {
		std::cerr << "checking " << ce[k] << " & " << ce[i] << "..."
			  << std::endl;
	      }
	      n_hm_tests += 1;
	      if (INFINITE(h0->incremental_eval(g1, ace[i]))) {
		std::cerr << "constraints " << ppcs[ce[k]].name
			  << " and " << ppcs[ce[i]].name
			  << " are jointly unachievable" << std::endl;
		h_inc.add_binary_edge(ce[k], ce[i]);
	      }
	      else if (opt_H3 && (opt_max_d >= 2)) {
		index_set g2(g1);
		g2.insert(ace[i]);
		for (index_type j = i + 1; j < ce.length(); j++) {
		  if (opt_verbose) {
		    std::cerr << "checking " << ce[k] << " & " << ce[i]
			      << " & " << ce[j] << "..." << std::endl;
		  }
		  n_hm_tests += 1;
		  if (INFINITE(h0->incremental_eval(g2, ace[j]))) {
		    std::cerr << "constraints " << ppcs[ce[k]].name
			      << ", " << ppcs[ce[i]].name
			      << " and " << ppcs[ce[j]].name
			      << " are jointly unachievable" << std::endl;
		    index_set e;
		    e.insert(ce[k]);
		    e.insert(ce[i]);
		    e.insert(ce[j]);
		    h_inc.add_edge(e);
		  }
		}
	      }
	    }
	  }
      }
      delete h0;
      delete test_ins;
      n_hm_built += 1;
    }

    index_type d = opt_min_d;
    while ((d <= opt_max_d) && (d <= ca.length()) && !stats.break_signal_raised()) {
      std::cerr << "checking " << d << "-consistency..." << std::endl;
      mSubsetEnumerator e(ca.length(), d);
      bool more = e.first();
      while (more && !stats.break_signal_raised()) {
	assert(e.current_set_size() == d);
	index_set sel;
	e.current_set(ca, sel);
	if (h_inc.is_independent(sel) &&
	    (max_sat.first_superset(sel) == no_such_index)) {
	  if (opt_verbose)
	    std::cerr << "testing " << sel << " ("
		      << h_inc.n_edges() << " inc. sets, "
		      << n_hm_built << "/" << n_hm_tests << " tests, "
		      << hm_test_stats.time() << " sec., "
		      << hm_test_stats.peak_memory()
		      << "k)..." << std::endl;
	  Instance* test_ins = instance.copy();
	  mapping map(test_ins->n_actions());
	  for (index_type i = 0; i < sel.length(); i++)
	    ppcs[sel[i]].enforce(*test_ins, map);
	  index_vec ace;
	  if (ce.length() > 0) {
	    ace.assign_value(no_such_index, ce.length());
	    for (index_type k = 0; k < ce.length(); k++) {
	      ace[k] = test_ins->compile_pc_sometime_conjunction(ppcs[ce[k]].s_c, ppcs[ce[k]].name);
	      test_ins->atoms[ace[k]].goal = false;
	    }
	  }
	  test_ins->cross_reference();
	  Heuristic* h0 = 0;
	  if (opt_H3) {
	    CostTable* hm = new CostTable(*test_ins, hm_test_stats);
	    hm->compute_H3(ZeroACF());
	    h0 = hm;
	  }
	  else if (opt_H2) {
	    StaticMutex* mx = new StaticMutex(*test_ins);
	    h0 = mx;
	  }
	  else {
	    //CostTable* hm = new CostTable(*test_ins, hm_test_stats);
	    //hm->compute_H1(ZeroACF());
	    //h0 = hm;
	    h0 = new Reachability(*test_ins);
	  }
	  n_hm_tests += 1;
	  if (INFINITE(h0->eval(test_ins->goal_atoms))) {
	    std::cerr << "set " << sel << " = ";
	    write_ppc_set(std::cerr, ppcs, sel);
	    std::cerr << " is unachievable" << std::endl;
	    h_inc.add_edge(sel);
	  }
	  else {
	    if (opt_save_test) {
	      save_compiled_problem(reader, instance, ppcs, v_sum, prep,
				    sel, test_ins, "test", n_hm_tests + 1);
	    }
	    if (ce.length() > 0) {
	      bool_vec test_pos(false, ce.length());
	      for (index_type k = 0; k < ce.length(); k++) {
		if (opt_verbose) {
		  std::cerr << "checking " << sel << " & " << ce[k] << "..."
			    << std::endl;
		}
		n_hm_tests += 1;
		if (INFINITE(h0->incremental_eval(test_ins->goal_atoms, ace[k]))) {
		  test_pos[k] = true;
		  index_set sel1(sel);
		  sel1.insert(ce[k]);
		  if (h_inc.is_independent(sel1)) {
		    std::cerr << "set " << sel1 << " = ";
		    write_ppc_set(std::cerr, ppcs, sel1);
		    std::cerr << " is unachievable" << std::endl;
		    h_inc.add_edge(sel1);
		  }
		}
	      }
	      if ((opt_H2 || opt_H3) && (opt_max_d > 1)) {
		for (index_type k = 0; k < ce.length(); k++)
		  if (!test_pos[k]) {
		    index_set g1(test_ins->goal_atoms);
		    g1.insert(ace[k]);
		    for (index_type i = k + 1; i < ce.length(); i++) {
		      if (opt_verbose) {
			std::cerr << "checking " << sel << " & " << ce[k]
				  << " & " << ce[i] << "..." << std::endl;
		      }
		      n_hm_tests += 1;
		      if (INFINITE(h0->incremental_eval(g1, ace[i]))) {
			index_set sel2(sel);
			sel2.insert(ce[k]);
			sel2.insert(ce[i]);
			if (h_inc.is_independent(sel2)) {
			  std::cerr << "set " << sel2 << " = ";
			  write_ppc_set(std::cerr, ppcs, sel2);
			  std::cerr << " is unachievable" << std::endl;
			  h_inc.add_edge(sel2);
			}
		      }
		      else if (opt_H3 && (opt_max_d >= 2)) {
			index_set g2(g1);
			g2.insert(ace[i]);
			for (index_type j = i + 1; j < ce.length(); j++) {
			  if (opt_verbose) {
			    std::cerr << "checking " << sel << " & "
				      << ce[k] << " & " << ce[i] << " & "
				      << ce[j] << "..." << std::endl;
			  }
			  n_hm_tests += 1;
			  if (INFINITE(h0->incremental_eval(g2, ace[j]))) {
			    index_set sel3(sel);
			    sel3.insert(ce[k]);
			    sel3.insert(ce[i]);
			    sel3.insert(ce[j]);
			    if (h_inc.is_independent(sel3)) {
			      std::cerr << "set " << sel3 << " = ";
			      write_ppc_set(std::cerr, ppcs, sel3);
			      std::cerr << " is unachievable" << std::endl;
			      h_inc.add_edge(sel3);
			    }
			  }
			}
		      }
		    }
		  }
	      }
	    }
	  }
	  delete h0;
	  delete test_ins;
	  n_hm_built += 1;
	}
	more = e.next();
      }
      std::cerr << h_inc.n_edges() << " inconsistent sets found ("
		<< n_hm_built << "/" << n_hm_tests << " tests, "
		<< hm_test_stats << ")"	<< std::endl;
      d += 1;
    }
    hm_test_stats.stop();
    //h_inc.reduce_to_minimal();
    //std::cerr << h_inc.n_edges() << " minimal inconsistent sets" << std::endl;
  } // if (opt_hm_test) ...

  if (opt_save_min) {
    for (index_type k = 0; k < h_inc.n_edges(); k++)
      if (h_inc.is_minimal_edge(k)) {
	if (opt_save_compiled) {
	  save_compiled_problem(reader, instance, ppcs, v_sum, prep,
				h_inc.edge(k), 0, "min", k + 1);
	}
	else {
	  save_problem(reader, instance, ppcs, v_sum,
		       h_inc.edge(k), "min", k + 1);
	}
      }
  }

  if ((opt_open || opt_max_sets) && opt_pcc && !opt_strong_pcc) {
    std::cerr << "preparing for PCC..." << std::endl;
    pcc_prep_stats.start();
    StaticMutex* mx = prep.inconsistency();
    graph* lmg = prep.necessary_sb_graph(opt_nsb_H2);
    landmark_graph_triggered_edges(instance, *lmg, trlm);
    graph* nag = prep.never_after_graph(opt_nag_H2);
    pcc_prep_stats.stop();
    std::cerr << "finished in " << pcc_prep_stats.total_time()
	      << " seconds" << std::endl;
    std::cerr << trlm.size() << " potential triggered landmarks"
	      << std::endl;
  }

  if (opt_open && !stats.break_signal_raised()) {
    std::cerr << "computing minimal open sets..." << std::endl;
    index_set_vec min_open;
    hg_stats.start();
    compute_minimal_open(h_inc, max_sat, ppcs.size(), min_open);
    hg_stats.stop();
    std::cerr << min_open.length() << " minimal open sets found" << std::endl;
    index_type first_open = 0;
    index_type n_inc = h_inc.n_edges();
    count_type n_hm_proven = 0;
#ifdef NEW_PCC
    PCC test(instance, pcc_stats, ppcs,
	     prep.inconsistency(), prep.necessary_sb_graph(opt_nsb_H2),
	     &trlm, prep.never_after_graph(opt_nag_H2));
#endif
    while ((first_open < min_open.length()) && !stats.break_signal_raised()) {
      for (index_type k = first_open; (k < min_open.length()) && !stats.break_signal_raised(); k++) {
	bool unsolvable = false;
	bool solved = false;

	if (opt_pcc) {
	  std::cerr << "testing open set ";
	  write_ppc_set(std::cerr, ppcs, min_open[k]);
#ifdef NEW_PCC
	  if (opt_print_proof) {
	    InferenceTracker* proof = new InferenceTracker(test);
	    proof->set_verbose(opt_verbose);
	    test.set_tracker(proof);
	  }
	  unsolvable = test.run(min_open[k]);
#else
	  if (opt_strong_pcc) {
	    std::cerr << " with strong PCC" << std::endl;
	    StrongPCC test(instance, pcc_stats, ppcs, min_open[k],
			   opt_H2, true, opt_verbose);
	    unsolvable = test.run();
	  }
	  else {
	    std::cerr << " with weak PCC" << std::endl;
	    BasicPCC test(instance, pcc_stats, ppcs, min_open[k],
			  prep.inconsistency(), prep.landmark_graph(),
			  &trlm, prep.never_after_graph(), opt_verbose);
	    unsolvable = test.run();
	  }
#endif
	  if (unsolvable) {
	    std::cerr << "contradiction proven" << std::endl;
#ifdef NEW_PCC
	    if (test.tracker()) {
	      test.tracker()->print_proof(std::cout);
	    }
#endif
	    h_inc.add_edge(min_open[k]);
	    n_pcc_proven += 1;
	  }
	  else {
	    std::cerr << "contradiction not proven" << std::endl;
	  }
#ifdef NEW_PCC
	  test.set_tracker(NULL);
#endif
	  n_pcc_tests += 1;
	}

	if ((opt_hm_test || opt_plan) && !unsolvable) {
	  std::cerr << "compiling open set ";
	  write_ppc_set(std::cerr, ppcs, min_open[k]);
	  std::cerr << "..." << std::endl;
	  Instance* test_ins = instance.copy();
	  mapping map(test_ins->n_actions());
	  for (index_type i = 0; i < min_open[k].length(); i++)
	    ppcs[min_open[k][i]].enforce(*test_ins, map);
	  test_ins->cross_reference();
	  if (opt_hm_test) {
	    std::cerr << "testing open set ";
	    write_ppc_set(std::cerr, ppcs, min_open[k]);
	    std::cerr << " with h^" << (opt_H3 ? 3 : (opt_H2 ? 2 : 1))
		      << std::endl;
	    Heuristic* h0 = 0;
	    if (opt_H3) {
	      CostTable* hm = new CostTable(*test_ins, hm_test_stats);
	      hm->compute_H3(ZeroACF());
	      h0 = hm;
	    }
	    else if (opt_H2) {
	      hm_test_stats.start();
	      StaticMutex* mx = new StaticMutex(*test_ins);
	      hm_test_stats.stop();
	      h0 = mx;
	    }
	    else {
	      hm_test_stats.start();
	      h0 = new Reachability(*test_ins);
	      hm_test_stats.stop();
	    }
	    if (INFINITE(h0->eval(test_ins->goal_atoms))) {
	      std::cerr << "set proven unachievable" << std::endl;
	      h_inc.add_edge(min_open[k]);
	      n_hm_proven += 1;
	      unsolvable = true;
	    }
	    else {
	      std::cerr << "unachievability not proven" << std::endl;
	    }
	    delete h0;
	    n_hm_tests += 1;
	  }
	  if (opt_plan && !unsolvable) {
	    std::cerr << "searching for a plan..." << std::endl;
	    n_planner_calls += 1;
	    Schedule* soln = 0;
	    solved = plan_compiled(*test_ins, plan_stats, plan_time_limit,
				   false, soln);
	    if (solved) {
	      std::cerr << "set proven solvable!" << std::endl;
	      n_plans_found += 1;
	      max_sat.append(min_open[k]);
	      plan_vec dummy;
	      max_plan.append(dummy);
	    }
	  }
	  delete test_ins;
	}

	if (opt_cega && !unsolvable && !solved) {
	  n_cega_tests += 1;
	  Schedule* soln = 0;
	  bool dec = cega(instance, cega_stats, cega_time_limit, ppcs,
			  min_open[k], *(prep.inconsistency()),
			  *(prep.landmark_graph()), true, soln);
	  if (dec) {
	    n_cega_decided += 1;
	    if (soln) {
	      max_sat.append(min_open[k]);
	      plan_vec nsv;
	      nsv.append(soln);
	      max_plan.append(nsv);
	    }
	    else {
	      h_inc.add_edge(min_open[k]);
	    }
	  }
	}

	if (!opt_pcc && !opt_hm_test && !opt_cega) {
	  std::cerr << k + 1 << ". ";
	  write_ppc_set(std::cerr, ppcs, min_open[k]);
	  std::cerr << std::endl;
	}
	if (opt_save_open && !solved && !unsolvable) {
	  if (h_inc.is_independent(min_open[k])) {
	    if (opt_save_compiled) {
	      save_compiled_problem(reader, instance, ppcs, v_sum, prep,
				    min_open[k], 0, "open", k + 1);
	    }
	    else {
	      save_problem(reader, instance, ppcs, v_sum,
			   min_open[k], "open", k + 1);
	    }
	  }
	}
      }

      first_open = min_open.length();
      if (opt_iterate && !stats.break_signal_raised()) {
	if (h_inc.n_edges() > n_inc) {
	  n_inc = h_inc.n_edges();
	  index_set_vec new_min_open;
	  std::cerr << "computing minimal open sets..." << std::endl;
	  hg_stats.start();
	  compute_minimal_open(h_inc, max_sat, ppcs.size(), new_min_open);
	  hg_stats.stop();
	  for (index_type i = 0; i < new_min_open.length(); i++)
	    if (!min_open.contains(new_min_open[i]))
	      min_open.append(new_min_open[i]);
	  std::cerr << min_open.length() - first_open
		    << " new minimal open sets found" << std::endl;
	}
      }
    }

    // some summary stats at end of open set testing
    std::cerr << n_pcc_tests << " PCC tests in "
	      << pcc_prep_stats.total_time() << " (prep) + "
	      << pcc_stats.total_time() << " seconds (avg. "
	      << pcc_stats.total_time() / n_pcc_tests
	      << " sec/test), " << n_pcc_proven << " proven inconsistent"
	      << std::endl;
    std::cerr << n_hm_built << "/" << n_hm_tests << " h^m tests in "
	      << hm_test_stats.total_time() << " seconds (avg. "
	      << hm_test_stats.total_time() / n_hm_tests
	      << " sec/test), " << n_hm_proven << " proven inconsistent"
	      << std::endl;
    std::cerr << n_cega_tests << " CEGA tests in "
	      << cega_stats.total_time() << " seconds (avg. "
	      << cega_stats.total_time() / n_cega_tests
	      << " sec/test), "
	      << n_cega_decided << " of " << n_cega_tests
	      << " decided" << std::endl;
    std::cerr << n_planner_calls << " planner calls in "
	      << plan_stats.total_time() << " seconds (avg. "
	      << plan_stats.total_time() / n_planner_calls
	      << " sec/test), "
	      << n_plans_found << " of " << n_planner_calls
	      << " solved" << std::endl;
  }

  // LC_RNG rng;
  // rng.seed_with_pid();
  // variety(instance, prep, ppcs, n_common, h_inc, 1000, 1000, stats, rng);

  if (opt_max_sets && !stats.break_signal_raised()) {
    std::cerr << "computing maximal consistent sets..." << std::endl;
    index_set_vec mcs;
    hg_stats.start();
    h_inc.independent_sets(mcs);
    std::cerr << mcs.length() << " maximal consistent sets ("
	      << hg_stats.time() << " sec)" << std::endl;
    n_max_open_sets = mcs.length();
    hg_stats.stop();
    cost_vec set_value(0, mcs.length());
    cost_vec set_penalty(0, mcs.length());
    for (index_type k = 0; k < mcs.length(); k++) {
      if (opt_verbose) std::cerr << "{";
      for (index_type i = 0; i < mcs[k].length(); i++) {
	if (opt_verbose) {
	  if (i > 0) std::cerr << ", ";
	  std::cerr << ppcs[mcs[k][i]].name;
	}
	set_value[k] += ppcs[mcs[k][i]].weight;
      }
      set_penalty[k] = v_sum - set_value[k];
      if (opt_verbose) 
	std::cerr << "}: value = " << PRINT_NTYPE(set_value[k])
		  << ", penalty = " << PRINT_NTYPE(set_penalty[k])
		  << std::endl;
      bool unsolvable = false;
      if (opt_pcc) {
	if (opt_strong_pcc) {
	  std::cerr << "testing with strong PCC" << std::endl;
	  StrongPCC test(instance, pcc_stats, ppcs, mcs[k],
			 opt_H2, true, opt_verbose);
	  unsolvable = test.run();
	}
	else {
	  std::cerr << "testing with weak PCC" << std::endl;
	  BasicPCC test(instance, pcc_stats, ppcs, mcs[k],
			prep.inconsistency(), prep.landmark_graph(),
			&trlm, prep.never_after_graph(), opt_verbose);
	  unsolvable = test.run();
	}
	if (unsolvable) {
	  std::cerr << "contradiction proven" << std::endl;
	  h_inc.add_edge(mcs[k]);
	  n_pcc_proven += 1;
	}
	n_pcc_tests += 1;
      }
      if (opt_plan && !unsolvable) {
	std::cerr << "searching for a plan..." << std::endl;
	n_planner_calls += 1;
	Schedule* soln = 0;
	bool solved = plan(instance, ppcs, mcs[k], plan_stats,
			   plan_time_limit, false, soln);
	if (solved) {
	  std::cerr << "set proven solvable!" << std::endl;
	  n_plans_found += 1;
	  max_sat.append(mcs[k]);
	  plan_vec dummy;
	  max_plan.append(dummy);
	  if (set_value[k] > v_max_sat)
	    v_max_sat = set_value[k];
	}
      }
      if (set_value[k] > v_max_sat)
	n_improving_sets += 1;
    }

    v_max = cost_vec_util::max(set_value);

    std::cerr << "max value = " << PRINT_NTYPE(v_max) << std::endl;
    std::cerr << "min penalty = " << PRINT_NTYPE(v_sum - v_max) << std::endl;
    std::cerr << n_improving_sets << " strictly improving open sets"
	      << std::endl;

    if (opt_save_max) {
      for (index_type k = 0; k < mcs.length(); k++)
	if (set_value[k] >= (v_max - opt_slack)) {
	  if (opt_save_compiled) {
	    save_compiled_problem(reader, instance, ppcs, v_sum, prep,
				  mcs[k], 0, "max", k + 1);
	  }
	  else {
	    save_problem(reader, instance, ppcs, v_sum, mcs[k], "max", k + 1);
	  }
	}
    }
  }

  else if (opt_max && !stats.break_signal_raised()) {
    std::cerr << "computing maximum value..." << std::endl;
    index_set s_max;
    hg_stats.start();
    v_max = max_value(ppcs, h_inc, s_max);
    hg_stats.stop();
    std::cerr << "max value = " << PRINT_NTYPE(v_max) << std::endl;
    std::cerr << "min penalty = " << PRINT_NTYPE(v_sum - v_max) << std::endl;
    if (opt_verbose) {
      std::cerr << "max value set: ";
      write_ppc_set(std::cerr, ppcs, s_max);
      std::cerr << std::endl;
    }
  }

  if (opt_save_ext && (max_sat.length() > 0) && !stats.break_signal_raised()) {
    index_set_vec exts;
    for (index_type k = 0; k < max_sat.length(); k++) {
      bool_vec in(max_sat[k], ppcs.length());
      for (index_type i = 0; i < ppcs.length(); i++) if (!in[i]) {
	in[i] = true;
	if (h_inc.is_independent(in)) {
	  exts.insert_maximal(index_set(in));
	}
	in[i] = false;
      }
    }
    std::cerr << exts.length() << " consistent 1-step extensions..."
	      << std::endl;
    for (index_type k = 0; k < exts.length(); k++) {
      if (opt_save_compiled) {
	save_compiled_problem(reader, instance, ppcs, v_sum, prep,
			      exts[k], 0, "ext", k + 1);
      }
      else {
	save_problem(reader, instance, ppcs, v_sum, exts[k], "ext", k + 1);
      }
    }
  }

  // print final report (to stdout)

  std::cout << ";; " << n_hm_built << "/" << n_hm_tests << " h^m tests, "
	    << hm_test_stats.total_time()
	    << " seconds" << std::endl;
  std::cout << ";; " << n_pcc_tests << " PCC tests, "
	    << pcc_prep_stats.total_time() << " + " << pcc_stats.total_time()
	    << " seconds" << std::endl;
  std::cout << ";; " << n_cega_tests << " CEGA tests in "
	    << cega_stats.total_time() << " seconds" << std::endl;
  std::cout << ";; " << n_planner_calls << " planner calls in "
	    << plan_stats.total_time() << " seconds" << std::endl;

  std::cout << ";; " << hg_stats.total_time()
	    << " seconds computing consistent/open sets"
	    << std::endl;
  std::cout << ";; " << stats.total_time() << " seconds total"
	    << std::endl;
  std::cout << ";; " << max_sat.size() << " known-to-be satisfiable and "
	    << h_inc.n_edges() << " known-to-be inconsistent sets"
	    << std::endl;
  std::cout << ";; achieved: max value = " << PRINT_NTYPE(v_max_sat)
	    << " (" << v_max_sat
	    << "), min penalty = " << PRINT_NTYPE(v_sum - v_max_sat)
	    << " (" << (v_sum - v_max_sat)
	    << ")" << std::endl;
  std::cout << ";; possible: max value = " << PRINT_NTYPE(v_max)
	    << " (" << v_max
	    << "), min penalty = " << PRINT_NTYPE(v_sum - v_max)
	    << " (" << (v_sum - v_max)
	    << ")" << std::endl;
  if (opt_max_sets) {
    std::cout << ";; " << n_max_open_sets << " maximal open sets"
	      << std::endl;
    std::cout << ";; " << n_improving_sets << " strictly improving open sets"
	      << std::endl;
  }

  if (opt_pddl) {
    std::cout << "(define (problem " << reader->problem_name << ")"
	      << std::endl;
    if (opt_generate) {
      std::cout << " ;; generated constraints:" << std::endl;
      std::cout << " (:constraints (and" << std::endl;
      for (index_type k = 0; k < ppcs.size(); k++) {
	std::cout << "   ";
	ppcs[k].write(std::cout, instance);
	//std::cout << std::endl;
      }
      std::cout << "))" << std::endl;
    }
    std::cout << " ;; inconsistent preference sets:" << std::endl;
    for (index_type k = 0; k < h_inc.n_edges(); k++) {
      if (!opt_load || !min_inc.contains(h_inc.edge(k)) || opt_pddl_all) {
	std::cout << " (:set";
	for (index_type i = 0; i < h_inc.edge(k).length(); i++)
	  std::cout << " " << ppcs[h_inc.edge(k)[i]].name;
	std::cout << ")" << std::endl;
      }
    }
    std::cout << " ;; maximal preference sets/plans:" << std::endl;
    Schedule::write_traits = false;
    for (index_type k = 0; k < max_sat.size(); k++) {
      std::cout << " ;; set ";
      write_ppc_set(std::cout, ppcs, max_sat[k]);
      NTYPE val = ppc_set_value(ppcs, max_sat[k]);
      std::cout << ", value = " << PRINT_NTYPE(val)
		<< ", penalty = " << PRINT_NTYPE(v_sum - val)
		<< std::endl;
      for (index_type i = 0; i < max_plan[k].size(); i++)
	if ((max_plan[k][i]->find_trait("SourceFileName") == 0) ||
	    opt_pddl_all) {
	  max_plan[k][i]->write(std::cout);
	}
    }
    std::cout << ")" << std::endl;
  }

  return 0;
}

END_HSPS_NAMESPACE

#ifdef USE_HSPS_NAMESPACE

int main(int argc, char *argv[])
{
  return HSPS::main(argc, argv);
}

#endif
