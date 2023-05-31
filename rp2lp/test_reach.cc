
#include "problem.h"
#include "preprocess.h"
#include "parser.h"
#include "cost_table.h"
#include "heuristic.h"

extern "C" {
  typedef unsigned int reach_fluent_id;
  typedef unsigned int reach_action_id;
  struct reach_struct;

  typedef enum {
    no_segment,
    reading,
    generate,
    solve,
    landmark,
    hittingset,
    time_segments
  } segment;

  typedef struct {
    int startsec;
    int startusec;
    int last_call;
    int latest[time_segments];
    int cumulative[time_segments];
  } reach_time_record;

  reach_struct* reach_getProblem();
  void reach_setVerbosity(int v, reach_struct *r);
  reach_fluent_id reach_addFluent(char *s, reach_struct *r);
  reach_action_id reach_addAction(reach_fluent_id *prec,
				  reach_fluent_id *postc,
				  reach_struct *r);
  void reach_set_initial(reach_fluent_id initial[],
			 reach_struct *r);
  void reach_set_goal(reach_fluent_id goal[],
		      reach_struct *r);
  int reach_checkReachable(reach_struct *r);
  int hplus(reach_struct *r);
  reach_action_id *reach_hittingSet(reach_struct *r, reach_time_record *t);
}

int main(int argc, char *argv[]) {
  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  int         verbose_level = 1;

  HSPS::Statistics  stats;
  //stats.enable_interrupt();
  stats.start();
  HSPS::Statistics  prep_stats(&stats);
  HSPS::Statistics  h_stats(&stats);
  HSPS::Statistics  reach_stats(&stats);

  HSPS::Parser* reader = new HSPS::Parser(symbols);

  for (int k = 1; k < argc; k++) {
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      verbose_level = atoi(argv[++k]);
    }
    else if (*argv[k] != '-') {
      reader->read(argv[k], false);
    }
  }

  HSPS::Heuristic::default_trace_level = verbose_level;
  HSPS::Instance::default_trace_level = verbose_level;
  HSPS::Preprocessor::default_trace_level = verbose_level;
  if (verbose_level <= 0) HSPS::PDDL_Base::warning_level = 0;
  if (verbose_level > 1) HSPS::PDDL_Base::write_info = true;

  prep_stats.start();
  std::cerr << "instantiating..." << std::endl;
  HSPS::PDDL_Base::name_instance_by_problem_file = true;
  HSPS::Instance instance;
  reader->instantiate(instance);

  HSPS::Preprocessor prep(instance, prep_stats);
  std::cerr << "preprocessing..." << std::endl;
  prep.preprocess(false);
  prep.compute_irrelevant_atoms();
  prep.remove_irrelevant_atoms();
  if (!instance.cross_referenced()) {
    std::cerr << "re-cross referencing..." << std::endl;
    instance.cross_reference();
  }

  prep_stats.stop();
  std::cerr << "instance " << instance.name << " built in "
	    << prep_stats.total_time() << " seconds" << std::endl;
  std::cout << instance.name << ": "
	    << instance.n_atoms() << " atoms, "
	    << instance.n_actions() << " actions, "
	    << std::endl;

  HSPS::CostTable* h = new HSPS::CostTable(instance, h_stats);
  h->compute_H1(HSPS::UnitACF());
  NTYPE c1 = h->eval(instance.goal_atoms);
  std::cout << "h^1 = " << PRINT_NTYPE(c1) << std::endl;

  std::cerr << "initialising reach struct..." << std::endl;
  reach_struct* p = reach_getProblem();
  reach_setVerbosity(verbose_level, p);

  reach_fluent_id atom_rfi[instance.n_atoms()];
  for (HSPS::index_type k = 0; k < instance.n_atoms(); k++) {
    atom_rfi[k] = reach_addFluent(instance.atoms[k].name->to_cstring(), p);
  }

  HSPS::lvector<reach_action_id> action_id(0, instance.n_actions());
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
    reach_fluent_id prec[instance.actions[k].pre.size() + 2];
    for (HSPS::index_type i = 0; i < instance.actions[k].pre.size(); i++)
      prec[i] = atom_rfi[instance.actions[k].pre[i]];
    prec[instance.actions[k].pre.size()] = 0;
    reach_fluent_id adds[instance.actions[k].add.size() + 2];
    for (HSPS::index_type i = 0; i < instance.actions[k].add.size(); i++)
      adds[i] = atom_rfi[instance.actions[k].add[i]];
    adds[instance.actions[k].add.size()] = 0;
    action_id[k] = reach_addAction(prec, adds, p);
  }

  reach_fluent_id r_init[instance.init_atoms.size() + 2];
  for (HSPS::index_type i = 0; i < instance.init_atoms.size(); i++)
    r_init[i] = atom_rfi[instance.init_atoms[i]];
  r_init[instance.init_atoms.size()] = 0;
  reach_set_initial(r_init, p);

  reach_fluent_id r_goal[instance.goal_atoms.size() + 2];
  for (HSPS::index_type i = 0; i < instance.goal_atoms.size(); i++)
    r_goal[i] = atom_rfi[instance.goal_atoms[i]];
  r_goal[instance.goal_atoms.size()] = 0;
  reach_set_goal(r_goal, p);

  std::cerr << "checking reachability..." << std::endl;
  int reach_ok = reach_checkReachable(p);
  if (!reach_ok) {
    std::cerr << "goal not reachable!" << std::endl;
    return 1;
  }

  std::cerr << "computing relaxed plan..." << std::endl;
  reach_time_record t;
  reach_stats.start();
  reach_action_id* rplan = reach_hittingSet(p, &t);
  reach_stats.stop();
  std::cerr << "done" << std::endl;

  HSPS::index_type k = 0;
  while (rplan[k] != 0) {
    HSPS::index_type a = action_id.first(rplan[k]);
    if (a == HSPS::no_such_index) {
      std::cerr << "error: invalid index " << rplan[k] << " in plan"
		<< std::endl;
      exit(1);
    }
    std::cout << instance.actions[a].name << std::endl;
    k += 1;
  }
  std::cout << "hplus = " << k << std::endl;
  stats.stop();
  std::cout << "time: " << prep_stats.total_time()
	    << " sec. preprocessing," << std::endl
	    << "      " << reach_stats.total_time()
	    << " sec. computing relaxed plan, " << std::endl
	    << "      " << stats.total_time()
	    << " sec. total" << std::endl;
  std::cout << "reach internal timing:" << std::endl
	    << "solving: " << (t.cumulative[solve] / 1000.0) << std::endl
	    << "landmark: " << (t.cumulative[landmark] / 1000.0) << std::endl
	    << "hitting set: " << (t.cumulative[hittingset] / 1000.0)
	    << std::endl;

  return 0;
}
