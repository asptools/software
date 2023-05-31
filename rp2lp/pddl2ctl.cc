
#include "parser.h"
#include "preprocess.h"
#include "exec.h"
#include "hypergraph.h"
#include "cost_table.h"
#include "enumerators.h"

#include <sstream>
#include <fstream>

int main(int argc, char *argv[]) {
  bool     opt_preprocess = true;
  bool     opt_verbose = false;
  bool     opt_single_line = false;
  bool     opt_bi_ssa = false;
  unsigned int random_seed = 0;

  HSPS::PDDL_Base::warning_level = 0;
  HSPS::PDDL_Base::compile_away_plan_constraints = false;

  HSPS::Instance::write_PDDL2 = false;
  HSPS::Instance::write_PDDL3 = false;
  HSPS::Instance::write_metric = false;
  HSPS::Instance::write_time = false;

  HSPS::Statistics stats;
  stats.enable_interrupt();
  stats.start();

  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  HSPS::Parser* reader = new HSPS::Parser(symbols);

  for (int k = 1; k < argc; k++) {
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      HSPS::Instance::default_trace_level = atoi(argv[++k]);
      opt_verbose = true;
    }
    else if (strcmp(argv[k],"-no-prep") == 0) {
      opt_preprocess = false;
    }
    else if (strcmp(argv[k],"-s") == 0) {
      opt_single_line = true;
    }
    else if (strcmp(argv[k],"-b") == 0) {
      opt_bi_ssa = true;
    }
    else if ((strcmp(argv[k],"-catc") == 0) && (k < argc - 1)) {
      HSPS::PDDL_Name::catc = *argv[++k];
    }
    else if (((strcmp(argv[k],"-rnd") == 0) ||
	      (strcmp(argv[k],"-r") == 0)) &&
	     (k < argc - 1)) {
      random_seed = atoi(argv[++k]);
    }
    else if (*argv[k] != '-') {
      reader->read(argv[k], false);
    }
  }

  HSPS::Instance instance;
//   ppc_vec  ppcs;

  reader->post_process();

  stats.start();
  reader->instantiate(instance);
  HSPS::index_set hard_goals;
  for (HSPS::index_type k = 0; k < instance.n_atoms(); k++)
    if (instance.atoms[k].goal)
      hard_goals.insert(k);

//   bool ok = true;
//   for (index_type k = 0; k < reader->dom_preferences.length(); k++) {
//     if (reader->dom_preferences[k]->goal->g_class
// 	== PDDL_Base::goal_always) {
//       if (!append_ppc(reader,
// 		      reader->dom_preferences[k],
// 		      ((PDDL_Base::SimpleSequenceGoal*)
// 		       reader->dom_preferences[k]->goal)->constraint,
// 		      0,
// 		      instance,
// 		      pc_always,
// 		      k,
// 		      ppcs))
// 	ok = false;
//     }
//     else if (reader->dom_preferences[k]->goal->g_class
// 	     == PDDL_Base::goal_sometime) {
//       if (!append_ppc(reader,
// 		      reader->dom_preferences[k],
// 		      ((PDDL_Base::SimpleSequenceGoal*)
// 		       reader->dom_preferences[k]->goal)->constraint,
// 		      0,
// 		      instance,
// 		      pc_sometime,
// 		      k,
// 		      ppcs))
// 	ok = false;
//     }
//     else if (reader->dom_preferences[k]->goal->g_class
// 	     == PDDL_Base::goal_at_most_once) {
//       if (!append_ppc(reader,
// 		      reader->dom_preferences[k],
// 		      ((PDDL_Base::SimpleSequenceGoal*)
// 		       reader->dom_preferences[k]->goal)->constraint,
// 		      0,
// 		      instance,
// 		      pc_at_most_once,
// 		      k,
// 		      ppcs))
// 	ok = false;
//     }
//     else if (reader->dom_preferences[k]->goal->g_class
// 	     == PDDL_Base::goal_sometime_before) {
//       if (!append_ppc(reader,
// 		      reader->dom_preferences[k],
// 		      ((PDDL_Base::TriggeredSequenceGoal*)
// 		       reader->dom_preferences[k]->goal)->constraint,
// 		      ((PDDL_Base::TriggeredSequenceGoal*)
// 		       reader->dom_preferences[k]->goal)->trigger,
// 		      instance,
// 		      pc_sometime_before,
// 		      k,
// 		      ppcs))
// 	ok = false;
//     }
//     else {
//       if (!append_ppc(reader,
// 		      reader->dom_preferences[k],
// 		      reader->dom_preferences[k]->goal,
// 		      0,
// 		      instance,
// 		      pc_always,
// 		      k,
// 		      ppcs))
// 	ok = false;
//     }
//   }
//
//   if (!ok) exit(255);

  HSPS::Preprocessor prep(instance, stats);

  if (opt_preprocess) {
    prep.preprocess();
//     for (index_type k = 0; k < ppcs.length(); k++) {
//       instance.remap_set(ppcs[k].s_c, prep.atom_map);
//       instance.remap_set(ppcs[k].s_t, prep.atom_map);
//     }
    instance.remap_set(hard_goals, prep.atom_map);
  }
  else {
    instance.cross_reference();
  }
  instance.set_goal(hard_goals);

  stats.stop();
  std::cerr << "instantiation and preprocessing finished in "
	    << stats.time() << " seconds" << std::endl;

  const char* separator = (opt_single_line ? " & " : "\n");

  // action preconditions: AG (act -> pre)
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
    if (!opt_single_line)
      std::cout << std::endl << "%% pre of " << instance.actions[k].name;
    if (k > 0) std::cout << separator;
    std::cout << "AG (";
    instance.actions[k].name->write
      (std::cout, HSPS::Name::NC_INSTANCE);
    std::cout << " => (";
    for (HSPS::index_type i = 0; i < instance.actions[k].pre.size(); i++) {
      if (i > 0) std::cout << " & ";
      instance.atoms[instance.actions[k].pre[i]].name->write
	(std::cout, HSPS::Name::NC_INSTANCE);
    }
    std::cout << "))";
  }

  // action mutex constraints
  if (!opt_single_line)
    std::cout << std::endl << "%% action exclusion constraints";
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
    for (HSPS::index_type j = k + 1; j < instance.n_actions(); j++)
      if (!instance.non_interfering(k, j)) {
	std::cout << separator << "AG (~";
	instance.actions[k].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << " | ~";
	instance.actions[j].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << ")";
      }

  // successor state axioms: AG ((AX p) <-> (add(p) | (p & ~del(p))))
  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++) {
    if (!opt_single_line)
      std::cout << std::endl << "%% ssa of " << instance.atoms[i].name;

    // bi-imp form
    if (opt_bi_ssa) {
      std::cout << separator << "AG ((AX ";
      instance.atoms[i].name->write
	(std::cout, HSPS::Name::NC_INSTANCE);
      std::cout << ") <=> (";
      if (!instance.atoms[i].add_by.empty()) {
	for (HSPS::index_type k = 0; k < instance.atoms[i].add_by.size(); k++) {
	  instance.actions[instance.atoms[i].add_by[k]].name->write
	    (std::cout, HSPS::Name::NC_INSTANCE);
	  std::cout << " | ";
	}
      }
      if (!instance.atoms[i].del_by.empty()) {
	std::cout << "(";
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	for (HSPS::index_type k = 0; k < instance.atoms[i].del_by.size(); k++) {
	  std::cout << " & ~";
	  instance.actions[instance.atoms[i].del_by[k]].name->write
	    (std::cout, HSPS::Name::NC_INSTANCE);
	}
	std::cout << ")";
      }
      else {
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
      }
      std::cout << "))";
    }

    // 2 single-imp form
    else {
      // positive
      std::cout << separator << "AG ((";
      if (!instance.atoms[i].add_by.empty()) {
	for (HSPS::index_type k = 0; k < instance.atoms[i].add_by.size(); k++) {
	  instance.actions[instance.atoms[i].add_by[k]].name->write
	    (std::cout, HSPS::Name::NC_INSTANCE);
	  std::cout << " | ";
	}
      }
      if (!instance.atoms[i].del_by.empty()) {
	std::cout << "(";
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	for (HSPS::index_type k = 0; k < instance.atoms[i].del_by.size(); k++) {
	  std::cout << " & ~";
	  instance.actions[instance.atoms[i].del_by[k]].name->write
	    (std::cout, HSPS::Name::NC_INSTANCE);
	}
	std::cout << ")";
      }
      else {
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
      }
      std::cout << ") => (AX ";
      instance.atoms[i].name->write
	(std::cout, HSPS::Name::NC_INSTANCE);
      std::cout << "))";

      // negative
      std::cout << separator << "AG ((";
      if (!instance.atoms[i].del_by.empty()) {
	for (HSPS::index_type k = 0; k < instance.atoms[i].del_by.size(); k++) {
	  instance.actions[instance.atoms[i].del_by[k]].name->write
	    (std::cout, HSPS::Name::NC_INSTANCE);
	  std::cout << " | ";
	}
      }
      if (!instance.atoms[i].add_by.empty()) {
	std::cout << "(~";
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	for (HSPS::index_type k = 0; k < instance.atoms[i].add_by.size(); k++) {
	  std::cout << " & ~";
	  instance.actions[instance.atoms[i].add_by[k]].name->write
	    (std::cout, HSPS::Name::NC_INSTANCE);
	}
	std::cout << ")";
      }
      else {
	std::cout << "~";
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
      }
      std::cout << ") => (AX ~";
      instance.atoms[i].name->write
	(std::cout, HSPS::Name::NC_INSTANCE);
      std::cout << "))";
    }
  }

  // initial state
  if (!opt_single_line)
    std::cout << std::endl << "%% initial state";
  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++) {
    std::cout << separator;
    if (!instance.atoms[i].init) std::cout << "~";
    instance.atoms[i].name->write(std::cout, HSPS::Name::NC_INSTANCE);
  }

  // goal
  if (!opt_single_line)
    std::cout << std::endl << "%% goal";
  std::cout << separator << "EF (";
  for (HSPS::index_type i = 0; i < hard_goals.size(); i++) {
    if (i > 0) std::cout << " & ";
    instance.atoms[hard_goals[i]].name->write
      (std::cout, HSPS::Name::NC_INSTANCE);
  }
  std::cout << ")" << std::endl;
}
