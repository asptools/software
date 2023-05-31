
#include "problem.h"
#include "preprocess.h"
#include "parser.h"
#include "cost_table.h"
#include "seq_reg.h"
#include "para_reg.h"
#include "temporal.h"
#include "plans.h"
#include "idao.h"
#include "ida.h"
#include <fstream>

BEGIN_HSPS_NAMESPACE

int main(int argc, char *argv[])
{
  bool        opt_round = false;
  bool        opt_round_up = false;
  bool        opt_round_down = false;
  bool        opt_H1 = false;
  bool        opt_sequential = false;
  bool        opt_temporal = false;
  bool        opt_resource = false;
  bool        opt_cost = false;
  bool        opt_preprocess = true;
  bool        opt_rm_irrelevant = false;
  bool        opt_cc = false;
  bool        opt_tt = false;
  index_type  opt_tt_size = 31337;
  index_type  opt_st_size = 10007;
  bool        opt_bst = true;
  index_type  opt_bst_size = 10007;
  bool        opt_boost = true;
  bool        opt_init_boost = true;
  bool        opt_dynamic_wps = false;
  bool        opt_save_h = false;
  index_type  max_apx_limit = no_such_index;
  count_type  boost_wps_limit = count_type_max;
  long        time_limit = 0;
  int         verbose_level = 1;
  bool        opt_print_plan = true;
  bool        opt_pddl = false;

  Statistics  stats;
  stats.start();
  Statistics  prep_stats(&stats);
  Statistics  hb_stats(&stats);
  Statistics  a_stats(&stats);
  Statistics  s_stats(&stats);

  StringTable symbols(50, lowercase_map);
  Parser* reader = new Parser(symbols);

  for (int k = 1; k < argc; k++) {
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      verbose_level = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-name-by-file") == 0) {
      HSPS::PDDL_Base::name_instance_by_problem_file = true;
    }
    else if (strcmp(argv[k],"-dba-semantics") == 0) {
      HSPS::PDDL_Base::del_before_add_semantics = true;
    }
    else if (strcmp(argv[k],"-use-strict-borrow") == 0) {
      PDDL_Base::use_strict_borrow_definition = true;
    }
    else if (strcmp(argv[k],"-use-extended-borrow") == 0) {
      PDDL_Base::use_extended_borrow_definition = true;
    }
    else if (strcmp(argv[k],"-no-compile") == 0) {
      PDDL_Base::compile_away_disjunctive_preconditions = false;
      PDDL_Base::compile_away_conditional_effects = false;
    }
    else if (strcmp(argv[k],"-no-compact") == 0) {
      PDDL_Base::compact_resource_effects = false;
    }
    else if (strcmp(argv[k],"-round") == 0) {
      opt_round = true;
      opt_round_up = false;
      opt_round_down = false;
    }
    else if (strcmp(argv[k],"-round-up") == 0) {
      opt_round_up = true;
      opt_round = false;
      opt_round_down = false;
    }
    else if (strcmp(argv[k],"-round-down") == 0) {
      opt_round_down = true;
      opt_round_up = false;
      opt_round = false;
    }
    else if (strcmp(argv[k],"-seq") == 0) {
      opt_sequential = true;
    }
    else if (strcmp(argv[k],"-time") == 0) {
      opt_temporal = true;
    }
    else if (strcmp(argv[k],"-res") == 0) {
      opt_resource = true;
    }
    else if (strcmp(argv[k],"-cost") == 0) {
      opt_cost = true;
    }
    else if (strcmp(argv[k],"-no-prep") == 0) {
      opt_preprocess = false;
    }
    else if (strcmp(argv[k],"-rm") == 0) {
      opt_rm_irrelevant = true;
    }
    else if (strcmp(argv[k],"-no-boost") == 0) {
      opt_init_boost = false;
      opt_boost = false;
    }
    else if (strcmp(argv[k],"-no-init-boost") == 0) {
      opt_init_boost = false;
    }
    else if (strcmp(argv[k],"-dwps") == 0) {
      opt_dynamic_wps = true;
    }
    else if (strcmp(argv[k],"-1") == 0) {
      opt_H1 = true;
    }
    else if (strcmp(argv[k],"-tt") == 0) {
      opt_tt = true;
    }
    else if ((strcmp(argv[k],"-tt-size") == 0) && (k < argc - 1)) {
      opt_tt_size = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-st-size") == 0) && (k < argc - 1)) {
      opt_st_size = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-bst-size") == 0) && (k < argc - 1)) {
      opt_bst_size = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-no-bst") == 0) {
      opt_bst = false;
    }
    else if (strcmp(argv[k],"-cc") == 0) {
      opt_cc = true;
    }
    else if ((strcmp(argv[k],"-t") == 0) && (k < argc - 1)) {
      time_limit = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-no-plan") == 0) {
      opt_print_plan = false;
    }
    else if (strcmp(argv[k],"-pddl") == 0) {
      opt_pddl = true;
    }
    else if ((strcmp(argv[k],"-m") == 0) && (k < argc - 1)) {
      max_apx_limit = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-save") == 0) {
      opt_save_h = true;
    }
    else if (*argv[k] != '-') reader->read(argv[k], false);
  }

  SearchAlgorithm::default_trace_level = verbose_level;
  Heuristic::default_trace_level = verbose_level - 1;
  Instance::default_trace_level = verbose_level - 1;
  Preprocessor::default_trace_level = verbose_level - 1;
  if (verbose_level < 1) opt_print_plan = false;

  Instance    instance;

  NTYPE      root_est_cost = 0;
  bool       solved = false;
  NTYPE      current_est_cost = 0;
  index_type current_apx_limit;

  stats.enable_interrupt();
  if (time_limit > 0) stats.enable_time_out(time_limit);

  prep_stats.start();
  std::cerr << "instantiating..." << std::endl;
  reader->instantiate(instance);
  if (verbose_level > 0) {
    std::cout << "instance: " << instance.name << std::endl;
  }
  if (opt_round_up) instance.round_durations_up();
  else if (opt_round_down) instance.round_durations_down();
  else if (opt_round) instance.round_durations();

  Preprocessor prep(instance, prep_stats);
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
  }
  else {
    std::cerr << "cross referencing..." << std::endl;
    instance.cross_reference();
  }
  prep_stats.stop();

  std::cerr << "instance " << instance.name << " built in "
	    << prep_stats.total_time() << " seconds" << std::endl;
  std::cerr << instance.n_atoms() << " atoms, "
	    << instance.n_resources() << " resources ("
	    << instance.n_reusable_resources() << " reusable, "
	    << instance.n_consumable_resources() << " consumable), "
	    << instance.n_actions() << " actions" << std::endl;

  if (!stats.break_signal_raised()) {
    CostTableBoost cost_tab(instance, hb_stats);
    estimator_vec resource_cons_est(0, instance.n_resources());

    std::cerr << "computing heuristic..." << std::endl;
    if (opt_H1) {
      if (opt_temporal)
	cost_tab.compute_H1(MakespanACF(instance));
      else if (opt_sequential && opt_cost)
	cost_tab.compute_H1(CostACF(instance));
      else
	cost_tab.compute_H1(UnitACF());
    }
    else if (opt_temporal) {
      cost_tab.compute_H2C(MakespanACF(instance), opt_resource);
    }
    else if (opt_sequential) {
      if (opt_cost)
	cost_tab.compute_H2(CostACF(instance));
      else
	cost_tab.compute_H2(UnitACF());
    }
    else {
      cost_tab.compute_H2C(UnitACF(), false);
    }
    if (!stats.break_signal_raised()) {
      std::cerr << "heuristic computed in " << hb_stats.total_time()
		<< " seconds" << std::endl;
    }

#ifdef APPLY_NCW_NOOP_TRICK
    if (opt_temporal && !stats.break_signal_raised()) {
      std::cerr << "computing n.c.w. sets..." << std::endl;
      prep.compute_ncw_sets(cost_tab);
    }
#endif

    if (opt_resource && !stats.break_signal_raised()) {
      std::cerr << "computing resource estimators..." << std::endl;
      hb_stats.start();
      for (index_type k = 0; k < instance.n_resources(); k++) {
	CostTable* rce = new CostTable(instance, hb_stats);
	rce->compute_H1(ResourceConsACF(instance, k));
	resource_cons_est[k] = rce;
      }
      hb_stats.stop();
    }

    if (!stats.break_signal_raised()) { // after H^2 and n.c.w.
      root_est_cost = cost_tab.eval(instance.goal_atoms);
      std::cerr << "estimated goal cost: " << root_est_cost << std::endl;

      StateFactory* b_root = 0;
      RegressionResourceState* rcs =
	(opt_resource ? new RegressionResourceState(instance, resource_cons_est) : 0);

      if (opt_temporal) {
	b_root = new TemporalRSRegState(instance, cost_tab, rcs);
      }
      else if (opt_sequential) {
	ACF* acf =
	  (opt_cost ? (ACF*)new CostACF(instance) : (ACF*)new UnitACF());
	b_root = new SeqCRegState(instance, cost_tab, *acf, rcs);
      }
      else {
	b_root = new ParaRegState(instance, cost_tab);
      }

      HashTable* bstab = (opt_bst ? new HashTable(opt_bst_size) : 0);

      current_est_cost = root_est_cost;
      current_apx_limit = (opt_H1 ? 2 : 3);
      bool apx_converged = false;
      bool done = false;

      if (opt_init_boost) {
	hb_stats.start();
	CostTable::Entry* list = cost_tab.boostable_entries();
	std::cerr << "boosting heuristic ("
		  << (list ? list->list_length() : 0)
		  << " entries)..." << std::endl;
	list = cost_tab.boost(list, *b_root, bstab, current_est_cost,
			      instance.goal_atoms, boost_wps_limit, false);
	if (list) list->delete_list();
	if (!stats.break_signal_raised()) {
	  std::cerr << "boost completed in " << hb_stats.time()
		    << " seconds (" << stats << ", TUF: "
		    << (bstab ? bstab->TUF() : 0) << ", HCF: "
		    << (bstab ? bstab->HCF() : 0) << ")" << std::endl;
	}
	else {
	  done = true;
	}
	hb_stats.stop();
	current_est_cost = cost_tab.eval(instance.goal_atoms);
      }

      State* search_root = 0;
      rcs = (opt_resource? new RegressionResourceState(instance, resource_cons_est) : 0);
      if (opt_temporal) {
	search_root = new TemporalRSRegState(instance, cost_tab,
					     instance.goal_atoms, rcs);
      }
      else if (opt_sequential) {
	ACF* acf =
	  (opt_cost ? (ACF*)new CostACF(instance) : (ACF*)new UnitACF());
	search_root = new SeqCRegState(instance, cost_tab, *acf,
				       instance.goal_atoms, rcs);
      }
      else {
	search_root = new ParaRSRegState(instance, cost_tab,
					 instance.goal_atoms);
      }

      Result result;
      if (opt_print_plan) {
	result.set_plan_set(new Print(instance, std::cout));
      }
      IDA* search = 0;
      if (opt_tt) {
	std::cerr << "using IDA* with transposition table..." << std::endl;
	HashTable* tt = new HashTable(opt_tt_size);
	search = new IDA(s_stats, result, tt);
      }
      else {
	std::cerr << "using IDA*..." << std::endl;
	search = new IDA(s_stats, result);
      }
      search->set_cycle_check(opt_cc);
      s_stats.enable_iteration_limit(1);

      ApxResult apx_res;
      HashTable* sol_tab = new HashTable(opt_st_size);
      IDAO* apx_search = new IDAO(a_stats, apx_res, sol_tab);
      apx_search->set_store_cost(true);
      apx_search->set_cycle_check(opt_cc);
      bool first_try_at_current_apx_limit = true;

      std::cerr << "current est. cost: " << current_est_cost << std::endl;
      std::cerr << "searching..." << std::endl;
      current_est_cost = search->start(*search_root);
      solved = search->solved();
      count_type search_work = s_stats.evaluations();
      if (search->solved() || stats.break_signal_raised()) done = true;

      while (!done) {
	if (!done && !apx_converged) {
	  std::cerr << "current est. cost: " << current_est_cost << std::endl;
	  std::cerr << "computing recursive " << current_apx_limit
		    << "-approximation..." << std::endl;
	  AtomSet::max_set_size_encountered = 0;

	  State* apx_root = 0;
	  RegressionResourceState* rcs =
	    (opt_resource? new RegressionResourceState(instance, resource_cons_est) : 0);
	  if (opt_temporal) {
	    apx_root = new ApxTemporalRegState(instance, cost_tab,
					       instance.goal_atoms, rcs,
					       current_apx_limit);
	  }
	  else if (opt_sequential) {
	    ACF* acf =
	      (opt_cost ? (ACF*)new CostACF(instance) : (ACF*)new UnitACF());
	    apx_root = new ApxSeqRegState(instance, cost_tab, *acf,
					  instance.goal_atoms, rcs,
					  current_apx_limit);
	  }
	  else {
	    apx_root = new ApxParaRegState(instance, cost_tab,
					   instance.goal_atoms,
					   current_apx_limit);
	  }

	  NTYPE new_cost = 0;
	  apx_search->set_cost_limit(current_est_cost);
	  if (first_try_at_current_apx_limit) {
	    new_cost = apx_search->start(*apx_root);
	    search_work += a_stats.evaluations();
	    first_try_at_current_apx_limit = false;
	  }
	  else {
	    count_type work_before = a_stats.evaluations();
	    new_cost = apx_search->resume(*apx_root, apx_root->est_cost());
	    search_work += (a_stats.evaluations() - work_before);
	  }

	  if (!a_stats.break_signal_raised()) {
	    if (apx_search->solved()) {
	      std::cerr << current_apx_limit << "-apxroximate solution: "
			<< new_cost << " (max depth = "
			<< apx_res.solution_depth()
			<< ", TUF: " << sol_tab->TUF()
			<< ", HCF: " << sol_tab->HCF()
			<< ", " << stats << ")" << std::endl;
	      if (apx_res.min_solution()) apx_converged = true;
	      if (current_apx_limit == max_apx_limit) {
		apx_converged = true;
	      }
	      else {
		current_apx_limit += 1;
		first_try_at_current_apx_limit = true;
	      }
	    }
	    else {
	      std::cerr << "no " << current_apx_limit
			<< "-apxroximate solution within current bound (TUF: "
			<< sol_tab->TUF() << ", HCF: " << sol_tab->HCF()
			<< ", " << stats << ")" << std::endl;
	    }
	  }

	  current_est_cost = MAX(current_est_cost, new_cost);
	  if (stats.break_signal_raised()) done = true;
	  delete apx_root;
	}

	if (!done && opt_boost) {
	  hb_stats.start();
	  CostTable::Entry* list = cost_tab.boostable_entries();
	  std::cerr << "current est. cost: " << current_est_cost << std::endl;
	  if (list) {
	    std::cerr << "boosting heuristic (" << list->list_length()
		      << " entries)..." << std::endl;
	    if (opt_dynamic_wps) {
	      boost_wps_limit = (10*search_work)/list->list_length();
	    }
	    list = cost_tab.boost(list, *b_root, bstab, current_est_cost,
				  instance.goal_atoms, boost_wps_limit,
				  !opt_dynamic_wps);
	    if (list) list->delete_list();
	    if (!stats.break_signal_raised()) {
	      std::cerr << "boost completed in " << hb_stats.time()
			<< " seconds (" << stats << ", TUF: "
			<< (bstab ? bstab->TUF() : 0) << ", HCF: "
			<< (bstab ? bstab->HCF() : 0) << ")" << std::endl;
	    }
	    else {
	      done = true;
	    }
	  }
	  current_est_cost =
	    MAX(cost_tab.eval(instance.goal_atoms), current_est_cost);
	  hb_stats.stop();
	}

	if (!done) {
	  std::cerr << "current est. cost: " << current_est_cost << std::endl;
	  std::cerr << "searching..." << std::endl;
	  s_stats.enable_iteration_limit(s_stats.total_complete_iterations() + 1);
	  current_est_cost = search->resume(*search_root, current_est_cost);
	  solved = search->solved();
	  search_work = s_stats.evaluations();
	  if (search->solved() || stats.break_signal_raised()) done = true;
	}
      }

      if (opt_save_h) {
	std::ofstream h_file("h.pddl");
	cost_tab.write_pddl(h_file, instance);
	h_file.close();
      }
    }
  }

  if (verbose_level > 0) {
    if (solved) {
      std::cout << "solution cost: " << current_est_cost << std::endl;
    }
    else {
      std::cout << "no solution found" << std::endl;
    }
    stats.print_total(std::cout);
    std::cout << "preprocessing: " << prep_stats << std::endl;
    std::cout << "heuristic & boosting: " << hb_stats << std::endl;
    std::cout << "approximate search: " << a_stats << std::endl;
    std::cout << "main search: " << s_stats << std::endl;
    double total_t = Stopwatch::seconds();
    std::cout << total_t << " seconds total (" << total_t - stats.total_time()
	      << " sec. not accounted for)" << std::endl;
  }

  else {
    std::cout << instance.name
	      << ' ' << (solved ? 1 : 0)
	      << ' ' << PRINT_NTYPE(root_est_cost)
	      << ' ' << PRINT_NTYPE(current_est_cost)
	      << ' ' << stats.total_nodes()
	      << ' ' << stats.total_time()
	      << ' ' << stats.peak_memory()
	      << ' ' << s_stats.total_nodes()
	      << ' ' << s_stats.total_time()
	      << ' ' << current_apx_limit
	      << ' ' << a_stats.total_min_nodes()
	      << ' ' << a_stats.total_max_nodes()
	      << ' ' << a_stats.total_time()
	      << ' ' << hb_stats.total_nodes()
	      << ' ' << hb_stats.total_time()
	      << ' ' << CostTableBoost::n_entries_boosted
	      << ' ' << CostTableBoost::n_entries_solved
	      << ' ' << CostTableBoost::n_entries_discarded
	      << ' ' << CostTableBoost::n_boost_searches
#ifdef EVAL_EXTRA_STATS
	      << ' ' << CostNode::eval_count
	      << ' ' << CostNode::eval_rec_count
#endif
#ifdef SEARCH_EXTRA_STATS
	      << ' ' << (rminx_size/(double)rminx_count)
	      << ' ' << (rminc_succ_size_ratio/(double)rminc_count)
	      << ' ' << (rminx_succ/(double)rminx_count)
	      << ' ' << (rmaxx_size/(double)rmaxx_count)
	      << ' ' << (rmaxx_succ/(double)rmaxx_count)
	      << ' ' << (trie_count/(double)rminx_count)
	      << ' ' << (tries_applicable/(double)trie_count)
	      << ' ' << (tries_within_bound/(double)trie_count)
#endif
	      << std::endl;
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
