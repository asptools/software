
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

  // action descriptions
  std::cout << "AXIOMS" << std::endl;
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
    for (HSPS::index_type i = 0; i < instance.actions[k].pre.size(); i++) {
      std::cout << "(<";
      instance.actions[k].name->write
	(std::cout, HSPS::Name::NC_INSTANCE);
      std::cout << "> True) => ";
      instance.atoms[instance.actions[k].pre[i]].name->write
	(std::cout, HSPS::Name::NC_INSTANCE);
      std::cout << std::endl;
    }
    for (HSPS::index_type i = 0; i < instance.n_atoms(); i++) {
      if (instance.actions[k].add.contains(i)) {
	std::cout << "[";
	instance.actions[k].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << "] ";
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << std::endl;
      }
      else if (instance.actions[k].del.contains(i)) {
	std::cout << "[";
	instance.actions[k].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << "] ~";
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << std::endl;
      }
      else {
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << " => ([";
	instance.actions[k].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << "] ";
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << ")" << std::endl;
	// negation
	std::cout << "~";
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << " => ([";
	instance.actions[k].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << "] ~";
	instance.atoms[i].name->write
	  (std::cout, HSPS::Name::NC_INSTANCE);
	std::cout << ")" << std::endl;
      }
    }
  }

  std::cout << "QUERY" << std::endl;
  // initial state
  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++) {
    if (instance.atoms[i].init) {
      instance.atoms[i].name->write
	(std::cout, HSPS::Name::NC_INSTANCE);
      std::cout << std::endl;
    }
    else {
      std::cout << "~";
      instance.atoms[i].name->write
	(std::cout, HSPS::Name::NC_INSTANCE);
      std::cout << std::endl;
    }
  }

  // goal
  std::cout << "<*(";
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
    if (k > 0) std::cout << " + ";
    instance.actions[k].name->write
      (std::cout, HSPS::Name::NC_INSTANCE);
  }
  std::cout << ")> (";
  for (HSPS::index_type i = 0; i < hard_goals.size(); i++) {
    if (i > 0) std::cout << " & ";
    instance.atoms[hard_goals[i]].name->write
      (std::cout, HSPS::Name::NC_INSTANCE);
  }
  std::cout << ")" << std::endl;
}
