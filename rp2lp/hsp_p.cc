
#include "problem.h"
#include "preprocess.h"
#include "enumerators.h"
#include "parser.h"
#include "cost_table.h"
#include "additive.h"
#include "seq_reg.h"
#include "plans.h"
#include "ida.h"
#include "bfs.h"
#include "bfhs.h"
#include "soft.h"

BEGIN_HSPS_NAMESPACE

int main(int argc, char *argv[]) {
  StringTable symbols(50, lowercase_map);
  bool        opt_H1 = false;
  bool        opt_AH = false;
  bool        opt_load_partition = false;
  bool        opt_pia1 = false;
  bool        opt_pia2a = false;
  bool        opt_round = false;
  bool        opt_find_invariants = false;
  bool        opt_apply_invariants = false;
  bool        opt_resource = false;
  bool        opt_compose_resources = false;
  index_type  composite_resource_size = 2;
  bool        opt_R2 = false;
  bool        opt_apply_cuts = true;
  bool        opt_preprocess = true;
  bool        opt_rm_irrelevant = false;
  bool        opt_print_plan = true;
  bool        opt_pddl = false;
  bool        opt_ipc = false;
  bool        opt_cc = false;
  index_type  tt_size = 100007;
  long        time_limit = 0;
  long        memory_limit = 0;
  count_type  node_limit = 0;
  NTYPE       cost_limit = 0;
  int         verbose_level = 1;

  Statistics stats;
  Statistics parse_stats(&stats);
  Statistics search_stats(&stats);
  Parser* reader = new Parser(symbols);

  for (int k = 1; k < argc; k++) {
    // verbose level
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      verbose_level = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-print-options-max") == 0) && (k < argc - 1)) {
      MaxValueSearch::print_options_max = atoi(argv[++k]);
    }

    // problem input/tranformation options
    else if (strcmp(argv[k],"-use-strict-borrow") == 0) {
      PDDL_Base::use_strict_borrow_definition = true;
    }
    else if (strcmp(argv[k],"-use-extended-borrow") == 0) {
      PDDL_Base::use_extended_borrow_definition = true;
    }
    else if (strcmp(argv[k],"-dba-semantics") == 0) {
      PDDL_Base::del_before_add_semantics = true;
    }

    // search space (problem type) selection
    else if (strcmp(argv[k],"-res") == 0) {
      opt_resource = true;
    }

    // additional search space options
    else if (strcmp(argv[k],"-cut") == 0) {
      opt_apply_cuts = true;
    }
    else if (strcmp(argv[k],"-no-cut") == 0) {
      opt_apply_cuts = false;
    }
    else if (strcmp(argv[k],"-inv") == 0) {
      opt_apply_invariants = true;
    }

    // preprocessing options
    else if (strcmp(argv[k],"-prep") == 0) {
      opt_preprocess = true;
    }
    else if (strcmp(argv[k],"-no-prep") == 0) {
      opt_preprocess = false;
    }
    else if (strcmp(argv[k],"-rm") == 0) {
      opt_rm_irrelevant = true;
    }
    else if (strcmp(argv[k],"-find") == 0) {
      opt_find_invariants = true;
    }
    else if (strcmp(argv[k],"-no-find") == 0) {
      opt_find_invariants = false;
    }
    else if ((strcmp(argv[k],"-compose") == 0) && (k < argc - 1)) {
      opt_compose_resources = true;
      composite_resource_size = atoi(argv[++k]);
    }

    // heuristic options
    else if (strcmp(argv[k],"-1") == 0) {
      opt_H1 = true;
    }
    else if (strcmp(argv[k],"-AH") == 0) {
      opt_AH = true;
    }
    else if (strcmp(argv[k],"-use-lse") == 0) {
      AH::use_linear_scan_eval = true;
    }
    else if (strcmp(argv[k],"-pia1") == 0) {
      opt_pia1 = true;
      opt_find_invariants = true;
    }
    else if (strcmp(argv[k],"-pia2a") == 0) {
      opt_pia2a = true;
    }
    else if (strcmp(argv[k],"-load-p") == 0) {
      opt_load_partition = true;
    }
    else if (strcmp(argv[k],"-R2") == 0) {
      opt_R2 = true;
    }

    // search algorithm selection & options
    else if (strcmp(argv[k],"-cc") == 0) {
      opt_cc = true;
    }
    else if ((strcmp(argv[k],"-tt-size") == 0) && (k < argc - 1)) {
      tt_size = atoi(argv[++k]);
    }

    // limit-setting options
    else if ((strcmp(argv[k],"-t") == 0) && (k < argc - 1)) {
      time_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-y") == 0) && (k < argc - 1)) {
      memory_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-c") == 0) && (k < argc - 1)) {
      cost_limit = A_TO_N(argv[++k]);
    }
    else if ((strcmp(argv[k],"-n") == 0) && (k < argc - 1)) {
      node_limit = atoi(argv[++k]);
    }

    // output/formatting options
    else if (strcmp(argv[k],"-no-plan") == 0) {
      opt_print_plan = false;
    }
    else if (strcmp(argv[k],"-pddl") == 0) {
      opt_print_plan = true;
      opt_pddl = true;
    }
    else if (strcmp(argv[k],"-ipc") == 0) {
      opt_print_plan = true;
      opt_ipc = true;
    }
    else if (strcmp(argv[k],"-nsn") == 0) {
      Instance::write_atom_set_with_symbolic_names = false;
      Instance::write_action_set_with_symbolic_names = false;
    }
    else if (strcmp(argv[k],"-print-rational") == 0) {
      Print::decimal_time = false;
    }

    // misc. options
    else if ((strcmp(argv[k],"-exclude") == 0) && (k < argc - 1)) {
      char* tag = argv[++k];
      if (strcmp(tag, "all") == 0) {
	PDDL_Base::exclude_all_dkel_items = true;
      }
      else {
	const StringTable::Cell* c = symbols.find(tag);
	if (c) PDDL_Base::excluded_dkel_tags.insert(c->text);
      }
    }
    else if ((strcmp(argv[k],"-require") == 0) && (k < argc - 1)) {
      const StringTable::Cell* c = symbols.find(argv[++k]);
      if (c) PDDL_Base::required_dkel_tags.insert(c->text);
    }

    // input file
    else if (*argv[k] != '-') {
      parse_stats.start();
      reader->read(argv[k], false);
      parse_stats.stop();
    }
  }

  SearchAlgorithm::default_trace_level = verbose_level;
  Heuristic::default_trace_level = verbose_level - 1;
  Instance::default_trace_level = verbose_level - 1;
  Preprocessor::default_trace_level = verbose_level - 1;
  if (verbose_level < 1) opt_print_plan = false;
  if (verbose_level <= 0) PDDL_Base::warning_level = 0;
  if (verbose_level > 1) PDDL_Base::write_info = true;

  SoftInstance instance;
  Store       result(instance);

  stats.enable_interrupt();
  if (time_limit > 0) stats.enable_time_out(time_limit);
  if (memory_limit > 0) stats.enable_memory_limit(memory_limit);

  stats.start();
  std::cerr << "instantiating..." << std::endl;
  reader->instantiate(instance);
  reader->instantiate_soft(instance);

  if (opt_resource && opt_compose_resources && (instance.n_resources() > 1)) {
    mSubsetEnumerator crs(instance.n_resources(), composite_resource_size);
    bool more = crs.first();
    while (more) {
      index_set s;
      crs.current_set(s);
      instance.create_composite_resource(s);
      more = crs.next();
    }
  }

  Preprocessor prep(instance, stats);
  if (opt_preprocess) {
    std::cerr << "preprocessing..." << std::endl;
    prep.preprocess();
    if (opt_rm_irrelevant) {
      prep.compute_irrelevant_atoms();
      prep.remove_irrelevant_atoms();
      if (!instance.cross_referenced()) {
	std::cerr << "re-cross referencing..." << std::endl;
	instance.cross_reference();
      }
    }
    instance.remap_hard_goals(prep.atom_map);
    instance.remap_soft_goals(prep.atom_map);
  }
  else {
    std::cerr << "cross referencing..." << std::endl;
    instance.cross_reference();
  }
  if (opt_find_invariants) {
    prep.bfs_find_invariants();
  }
  instance.verify_invariants();
  stats.stop();

  std::cerr << "instance " << instance.name << " built in "
	    << stats.time() << " seconds" << std::endl;
  std::cerr << instance.n_atoms() << " atoms, "
	    << instance.n_hard() << " hard goals, "
	    << instance.n_soft() << " soft goals, "
	    << instance.n_resources() << " resources ("
	    << instance.n_reusable_resources() << " reusable, "
	    << instance.n_consumable_resources() << " consumable), "
	    << instance.n_actions() << " actions, "
	    << instance.n_invariants() << " invariants"
	    << std::endl;

  Heuristic* h = 0;
  CostACF f(instance);

  stats.start();
  if (opt_AH) {
    std::cerr << "computing AH heuristic..." << std::endl;
    AH* ah = new AH(instance, stats);
    if (opt_load_partition) {
      stats.start();
      name_vec pnames(0, 0);
      index_set_vec partition;
      reader->export_action_partitions(pnames, partition);
      instance.remap_sets(partition, prep.action_map);
      ah->compute_additive(CostACF(instance), partition, !opt_H1);
      stats.stop();
    }
    else if (opt_pia1) {
      ah->compute_with_iterative_assignment_1
	(f, instance.goal_atoms, !opt_H1, false, false, instance.goal_atoms);
    }
    else if (opt_pia2a) {
      index_set a;
      a.fill(instance.n_atoms());
      ah->compute_with_iterative_assignment_2
	(f, a, !opt_H1, false, instance.goal_atoms);
    }
    else {
      ah->compute_with_iterative_assignment_2
	(f, instance.goal_atoms, !opt_H1, false, instance.goal_atoms);
    }
    // if (opt_round) {
    //   NTYPE c = f.cost_gcd(instance.n_actions());
    //   std::cerr << "rounding up to 1/" << c.divisor() << std::endl;
    //   RoundUp* rh = new RoundUp(instance, *ah, c.divisor());
    //   h = rh;
    // }
    // else {
    //   h = ah;
    // }
    h = ah;
  }
  else if (opt_H1) {
    CostTable* h1 = new CostTable(instance, stats);
    h1->compute_H1(f);
    h = h1;
  }
  else {
    CostTable* h2 = new CostTable(instance, stats);
    h2->compute_H2(f);
    h = h2;
  }
  stats.stop();

  index_type init_n = 0;
  NTYPE      init_best = 0;
  index_type final_n = 0;
  NTYPE      final_best = 0;
  index_type final_best_size = 0;
  bool       solved = false;

  if (!stats.break_signal_raised()) {
    std::cerr << "heuristic computed in " << stats.time() << " seconds"
	      << std::endl;

    RegressionResourceState* rs = 0;
    if (opt_resource && (instance.n_resources() > 0)) {
      estimator_vec* rest = new estimator_vec(0, instance.n_resources());
      for (index_type k = 0; k < instance.n_resources(); k++) {
	CostTable* c = new CostTable(instance, stats);
	if (opt_R2)
	  c->compute_H2(ResourceConsACF(instance, k));
	else
	  c->compute_H1(ResourceConsACF(instance, k));
	(*rest)[k] = c;
      }
      rs = new RegressionResourceState(instance, *rest);
    }

    MaxNetBenefit mnb(instance, f, search_stats, result, *h, rs, tt_size, opt_cc);
    mnb.init();
    init_n = mnb.n_options();
    init_best = mnb.best_option_estimated_value();
    std::cerr << "initial best value estimate: " << init_best << std::endl;
    std::cerr << "searching..." << std::endl;
    final_best = mnb.main();
    final_best_size = mnb.best_option_size();
    solved = mnb.solved();
  }

  if (opt_ipc || opt_pddl || (verbose_level == 0)) {
    if (opt_ipc) {
      std::cout << "; Time " << Stopwatch::seconds() << std::endl;
      std::cout << "; ParsingTime " << parse_stats.total_time() << std::endl;
      if (solved) {
	std::cout << "; MetricValue " << PRINT_NTYPE(final_best)
		  << std::endl;
	PrintIPC print_plan(instance, std::cout);
	result.output(print_plan);
      }
      else {
	std::cout << "; Not Solved" << std::endl;
      }
    }

    else if (opt_pddl) {
      if (solved) {
	PrintPDDL print_plan(instance, std::cout);
	result.output(print_plan);
      }
    }

    if (opt_ipc || opt_pddl)
      std::cout << ";; stats: ";
    std::cout << instance.name
	      << ' ' << (solved ? 1 : 0)
	      << ' ' << PRINT_NTYPE(init_best)
	      << ' ' << PRINT_NTYPE(final_best)
	      << ' ' << stats.total_nodes()
	      << ' ' << stats.total_time()
	      << ' ' << stats.peak_memory()
#ifdef RSS_FROM_PSINFO
	      << ' ' << stats.peak_total_size()
#endif
	      << ' ' << stats.time()
	      << ' ' << init_n
	      << ' ' << final_n
	      << ' ' << final_best_size
	      << ' ' << stats.flags()
	      << std::endl;
  }

  else {
    if (solved) {
      std::cout << "solution value: " << final_best << std::endl;
      if (opt_print_plan) {
	Print print_plan(instance, std::cout);
	result.output(print_plan);
      }
    }
    else {
      std::cout << "no solution found" << std::endl;
    }
    stats.print_total(std::cout);
    std::cout << "search: " << search_stats << std::endl;
    double total_t = Stopwatch::seconds();
    std::cout << total_t << " seconds total (" << total_t - stats.total_time()
	      << " sec. not accounted for)" << std::endl;
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
