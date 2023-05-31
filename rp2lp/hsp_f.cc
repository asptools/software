
#include "problem.h"
#include "preprocess.h"
#include "parser.h"
#include "additive.h"
#include "lmcut.h"
#include "pdb_construction.h"
#include "forward.h"
#include "ext_state.h"
#include "temporal.h"
#include "plans.h"
#include "ida.h"
#include "bfs.h"
#include "bfhs.h"
#include "bb.h"
#include "bss.h"
#include "pop.h"
#include "enumerators.h"
#include "explore.h"
#include "ilb.h"
#include <fstream>

int main(int argc, char *argv[]) {
  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  bool        opt_forward_rc = false;
  bool        opt_forward_H1 = false;
  bool        opt_forward_H2 = false;
  // bool        opt_forward_sumx = false;
  bool        opt_lm_cut = false;
  bool        opt_lm_cut_2 = false;
  bool        opt_hplus = false;
  bool        opt_reverse_H1 = false;
  bool        opt_reverse_H2 = false;
  bool        opt_reverse_AH2 = false;
  bool        opt_fpia = false;
  bool        opt_bu_partition = false;
  bool        opt_bu_glue = false;
  bool        opt_bu2_partition = false;
  NTYPE       bu2_ratio = R_TO_N(1, 4);
  bool        opt_kcut_partition = false;
  HSPS::index_type  opt_kcuts = 2;
  bool        opt_compare_h = false;
  bool        opt_H0 = false;
  bool        opt_HAlmost0 = false;
  bool        opt_pdb = true;
  bool        opt_pdb_load = false;
  bool        opt_ipdb = false;
  bool        opt_fast_pdb = false;
  HSPS::index_type  opt_repeat = 1;
  bool        opt_pdb_incremental = false;
  bool        opt_pdb_bin = false;
  bool        opt_pdb_random_bin = false;
  bool        opt_pdb_random_independent_bin = false;
  bool        opt_pdb_spanning = true;
  HSPS::index_type  opt_pdb_span_search_limit = HSPS::no_such_index;
  HSPS::index_type  opt_pdb_random_bin_swaps = 10000;
  bool        opt_pdb_ip = false;
  bool        opt_pdb_weighted_independent_bin = false;
  bool        opt_pdb_collapse = true;
  HSPS::index_type  opt_pdb_size = 1000000;
  HSPS::index_type  opt_total_size = 10000000;
  HSPS::index_type  opt_pdb_set_size = HSPS::index_type_max;
  bool        opt_pdb_add = true;
  bool        opt_pdb_ext_add = false;
  bool        opt_pdb_add_all = false;
  HSPS::rational    ext_add_threshold = 0;
  bool        opt_pdb_inc = true;
  bool        opt_maximal_add = true;
  bool        opt_sas_min = false;
  bool        opt_sas_safe = false;
  bool        opt_sas_select = false;
  bool        opt_apx_clique = false;
  bool        opt_extended = false;
  bool        opt_relaxed = false;
  bool        opt_resource = false;
  bool        opt_res_ipdb = false;
  bool        opt_compose_resources = false;
  HSPS::index_type  composite_resource_size = 2;
  bool        opt_R2 = false;
  bool        opt_sum_time = false;
  bool        opt_cost = false;
  bool        opt_zero = false;
  NTYPE       opt_scale_zero_cost = 1;
  bool        opt_apply_cuts = true;
  bool        opt_preprocess = true;
  bool        opt_preprocess_2 = true;
  bool        opt_rm_irrelevant = false;
  bool        opt_find_invariants = true;
  bool        opt_quick_find_invariants = false;
  bool        opt_verify_invariants = true;
  bool        opt_extend_goal = true;
  bool        opt_find_all = false;
  bool        opt_all_different = false;
  bool        opt_exhaustive = false;
  bool        opt_print_plan = true;
  bool        opt_print_plan_on_store = false;
  bool        opt_save_plan_to_file = false;
  bool        opt_schedule = false;
  bool        opt_post_op = false;
  bool        opt_validate = false;
  bool        opt_pddl = false;
  bool        opt_ipc = false;
  bool        opt_strict_ipc = false;
#ifdef NTYPE_RATIONAL
  NTYPE       epsilon = HSPS::rational(1,100);
#else
  NTYPE       epsilon = 0.001;
#endif

  HSPS::index_type  ipdb_param_d_skip = 1;
  HSPS::index_type  ipdb_param_s_max = HSPS::index_type_max;
  double      ipdb_param_i_min = 0;
  HSPS::index_type  ipdb_param_n_trials = 10;
  HSPS::index_type  ipdb_param_n_samples = 100;

  bool        opt_save = false;
  bool        opt_cc = false;
  bool        opt_tt = false;
  HSPS::index_type  opt_tt_size = 31337;
  bool        opt_bfs = false;
  bool        opt_bfs_px = false;
  NTYPE       px_threshold = 0;
  bool        opt_bss = false;
  bool        opt_bb = false;
  bool        opt_dfs = false;
  bool        opt_bfhs = false;
  bool        opt_bfida = false;
  bool        opt_bfxd = false;
  bool        opt_test = false;
  long        time_limit = 0;
  long        memory_limit = 0;
  HSPS::index_type  iteration_limit = 0;
  HSPS::count_type  node_limit = 0;
  NTYPE       cost_limit = -1;
  int         verbose_level = 1;
  bool        opt_bfs_write_graph = false;
  bool        opt_bfs_write_stats = false;
  bool        opt_explore_to_depth = false;
  HSPS::index_type  depth_limit = 0;
  bool        opt_explore_to_bound = false;
  bool        opt_ff = false;
  bool        opt_ff_P2 = false;
  bool        opt_cce = false;

  HSPS::Statistics  stats;
  stats.start();
  HSPS::Statistics parse_stats(&stats);
  HSPS::Statistics prep_stats(&stats);
  HSPS::Statistics h_stats(&stats);
  HSPS::Statistics search_stats(&stats);

  HSPS::Parser* reader = new HSPS::Parser(symbols);

  HSPS::LC_RNG rng;

  for (int k = 1; k < argc; k++) {
    // verbose level
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      verbose_level = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-no-warnings") == 0) {
      HSPS::PDDL_Base::warning_level = 0;
    }
    else if (strcmp(argv[k],"-no-info") == 0) {
      HSPS::PDDL_Base::write_info = false;
    }

    // problem input/tranformation options
    else if (strcmp(argv[k],"-dba-semantics") == 0) {
      HSPS::PDDL_Base::del_before_add_semantics = true;
    }
    else if (strcmp(argv[k],"-use-strict-borrow") == 0) {
      HSPS::PDDL_Base::use_strict_borrow_definition = true;
    }
    else if (strcmp(argv[k],"-use-extended-borrow") == 0) {
      HSPS::PDDL_Base::use_extended_borrow_definition = true;
    }
    else if (strcmp(argv[k],"-no-compile") == 0) {
      HSPS::PDDL_Base::create_all_atoms = true;
      HSPS::PDDL_Base::compile_away_disjunctive_preconditions = false;
      HSPS::PDDL_Base::compile_away_conditional_effects = false;
    }
    else if (strcmp(argv[k],"-cce") == 0) {
      HSPS::PDDL_Base::compile_away_conditional_effects = false;
      opt_cce = true;
    }
    else if (strcmp(argv[k],"-no-compact") == 0) {
      HSPS::PDDL_Base::compact_resource_effects = false;
    }
    else if (strcmp(argv[k],"-test") == 0) {
      opt_test = true;
    }

    // search space (problem type) selection
    else if (strcmp(argv[k],"-x") == 0) {
      opt_extended = true;
    }
    else if (strcmp(argv[k],"-relax") == 0) {
      opt_relaxed = true;
    }
    else if (strcmp(argv[k],"-res") == 0) {
      opt_resource = true;
    }
    else if (strcmp(argv[k],"-res-ipdb") == 0) {
      opt_res_ipdb = true;
    }
    else if ((strcmp(argv[k],"-compose") == 0) && (k < argc - 1)) {
      opt_compose_resources = true;
      composite_resource_size = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-sum-time") == 0) {
      opt_sum_time = true;
    }
    else if (strcmp(argv[k],"-cost") == 0) {
      opt_cost = true;
    }
    else if (strcmp(argv[k],"-zero") == 0) {
      opt_zero = true;
    }
    else if ((strcmp(argv[k],"-scale-zero-cost") == 0) && (k < argc - 1)) {
      opt_scale_zero_cost = A_TO_N(argv[++k]);
    }
    else if (strcmp(argv[k],"-sua") == 0) {
      HSPS::SeqProgState::separate_update_actions = true;
    }

    // additional search space options
    else if (strcmp(argv[k],"-cut") == 0) {
      opt_apply_cuts = true;
    }
    else if (strcmp(argv[k],"-no-cut") == 0) {
      opt_apply_cuts = false;
    }

    // preprocessing options
    else if (strcmp(argv[k],"-prep") == 0) {
      opt_preprocess = true;
    }
    else if (strcmp(argv[k],"-prep-1") == 0) {
      opt_preprocess = true;
      opt_preprocess_2 = false;
    }
    else if (strcmp(argv[k],"-no-prep") == 0) {
      opt_preprocess = false;
    }
    else if (strcmp(argv[k],"-extend") == 0) {
      opt_extend_goal = true;
    }
    else if (strcmp(argv[k],"-no-extend") == 0) {
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-rm") == 0) {
      opt_rm_irrelevant = true;
    }
    else if (strcmp(argv[k],"-find") == 0) {
      opt_find_invariants = true;
    }
    else if (strcmp(argv[k],"-quick-find") == 0) {
      opt_quick_find_invariants = true;
    }
    else if (strcmp(argv[k],"-no-find") == 0) {
      opt_find_invariants = false;
      opt_quick_find_invariants = false;
    }
    else if (strcmp(argv[k],"-verify") == 0) {
      opt_verify_invariants = true;
    }
    else if (strcmp(argv[k],"-no-verify") == 0) {
      opt_verify_invariants = false;
    }

    // heuristic options
    else if (strcmp(argv[k],"-r1") == 0) {
      opt_reverse_H1 = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-r") == 0) {
      opt_reverse_H2 = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-rAH") == 0) {
      opt_reverse_AH2 = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-fpia") == 0) {
      opt_fpia = true;
    }
    else if (strcmp(argv[k],"-use-lse") == 0) {
      HSPS::AH::use_linear_scan_eval = true;
    }
    else if (strcmp(argv[k],"-pbu") == 0) {
      opt_bu_partition = true;
    }
    else if (strcmp(argv[k],"-pbug") == 0) {
      opt_bu_partition = true;
      opt_bu_glue = true;
    }
    else if ((strcmp(argv[k],"-pbu2") == 0) && (k < argc - 1)) {
      opt_bu2_partition = true;
      bu2_ratio = A_TO_N(argv[++k]);
    }
    else if ((strcmp(argv[k],"-pkc") == 0) && (k < argc - 1)) {
      opt_kcut_partition = true;
      opt_kcuts = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-f") == 0) {
      opt_forward_rc = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-f1") == 0) {
      opt_forward_H1 = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-lm-cut") == 0) {
      opt_lm_cut = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-lm-cut-2") == 0) {
      opt_lm_cut_2 = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-hplus") == 0) {
      opt_hplus = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-f2") == 0) {
      opt_forward_H2 = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    // else if (strcmp(argv[k],"-fx") == 0) {
    //   opt_forward_sumx = true;
    //   opt_pdb = false;
    //   opt_find_invariants = false;
    //   opt_verify_invariants = false;
    // }
    else if (strcmp(argv[k],"-0") == 0) {
      opt_H0 = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-almost-blind") == 0) {
      opt_HAlmost0 = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
      opt_extend_goal = false;
    }
    else if (strcmp(argv[k],"-compare") == 0) {
      opt_compare_h = true;
    }
    else if (strcmp(argv[k],"-pdb-load") == 0) {
      opt_pdb_load = true;
      opt_pdb = false; // AARGH! (only to test use of PDBCollection instead
                       // of old PDBHeuristic...)
    }
    else if (strcmp(argv[k],"-pdb-fast") == 0) {
      opt_fast_pdb = true;
    }
    else if (strcmp(argv[k],"-ipdb") == 0) {
      opt_ipdb = true;
      opt_pdb = false; // looks stupid, but it is correct...
      opt_extend_goal = false;
      ipdb_param_d_skip = 1;
      ipdb_param_i_min = 0.01;
      ipdb_param_n_trials = 10;
      ipdb_param_n_samples = 100;
    }
    else if ((strcmp(argv[k],"-ipdb-param") == 0) && (k < argc - 5)) {
      // opt_ipdb = true;
      // opt_pdb = false; // looks stupid, but it is correct...
      // opt_extend_goal = false;
      ipdb_param_d_skip = atoi(argv[++k]);
      ipdb_param_s_max = atoi(argv[++k]);
      ipdb_param_i_min = atof(argv[++k]);
      ipdb_param_n_trials = atoi(argv[++k]);
      ipdb_param_n_samples = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-repeat") == 0) && (k < argc - 1)) {
      opt_repeat = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-pdb-i") == 0) {
      opt_pdb_incremental = true;
      opt_extend_goal = false;
      ipdb_param_d_skip = 1;
      ipdb_param_n_samples = 25;
    }
    else if ((strcmp(argv[k],"-pdb-i-param") == 0) && (k < argc - 2)) {
      opt_pdb_incremental = true;
      opt_extend_goal = false;
      ipdb_param_d_skip = atoi(argv[++k]);
      ipdb_param_n_samples = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-pdb-bin") == 0) {
      opt_pdb_bin = true;
    }
    else if (strcmp(argv[k],"-pdb-rbin") == 0) {
      opt_pdb_random_bin = true;
    }
    else if (strcmp(argv[k],"-pdb-rib") == 0) {
      opt_pdb_random_independent_bin = true;
    }
    else if (strcmp(argv[k],"-pdb-strict-rib") == 0) {
      opt_pdb_random_independent_bin = true;
      opt_pdb_collapse = false;
    }
    else if ((strcmp(argv[k],"-rbin-swaps") == 0) && (k < argc - 1)) {
      opt_pdb_random_bin_swaps = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-pdb-no-span") == 0) {
      opt_pdb_spanning = false;
    }
    else if ((strcmp(argv[k],"-pdb-span-search-limit") == 0) && (k < argc - 1)) {
      opt_pdb_span_search_limit = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-pdb-ip") == 0) {
      opt_pdb_ip = true;
    }
    else if (strcmp(argv[k],"-pdb-windbin") == 0) {
      opt_pdb_weighted_independent_bin = true;
    }
    else if (strcmp(argv[k],"-pdb-strict-windbin") == 0) {
      opt_pdb_weighted_independent_bin = true;
      opt_pdb_collapse = false;
    }
    else if ((strcmp(argv[k],"-pdb-size") == 0) && (k < argc - 1)) {
      opt_pdb_size = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-pdb-total-size") == 0) && (k < argc - 1)) {
      opt_total_size = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-pdb-set-size") == 0) && (k < argc - 1)) {
      opt_pdb_set_size = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-add") == 0) {
      opt_pdb_add = true;
    }
    else if (strcmp(argv[k],"-no-add") == 0) {
      opt_pdb_add = false;
    }
    else if (strcmp(argv[k],"-apx-add") == 0) {
      opt_pdb_add = true;
      opt_maximal_add = false;
    }
    else if (strcmp(argv[k],"-add-all") == 0) {
      opt_pdb_add = true;
      opt_pdb_add_all = true;
    }
    else if ((strcmp(argv[k],"-xadd") == 0) && (k < argc - 1)) {
      opt_pdb_add = true;
      opt_pdb_ext_add = true;
      ext_add_threshold = HSPS::rational::ator(argv[++k]);
    }
    else if (strcmp(argv[k],"-inc") == 0) {
      opt_pdb_inc = true;
    }
    else if (strcmp(argv[k],"-no-inc") == 0) {
      opt_pdb_inc = false;
    }
    else if (strcmp(argv[k],"-R2") == 0) {
      opt_R2 = true;
    }

    // search algorithm selection & options
    else if (strcmp(argv[k],"-bfs") == 0) {
      opt_bfs = true;
      opt_apply_cuts = false;
    }
    else if (strcmp(argv[k],"-bfida") == 0) {
      opt_bfida = true;
      opt_apply_cuts = false;
    }
    else if (strcmp(argv[k],"-bfhs") == 0) {
      opt_bfhs = true;
      opt_apply_cuts = false;
    }
    else if (strcmp(argv[k],"-bfxd") == 0) {
      opt_bfxd = true;
      opt_apply_cuts = false;
    }
    else if ((strcmp(argv[k],"-px") == 0) && (k < argc - 1)) {
      opt_bfs_px = true;
      px_threshold = A_TO_N(argv[++k]);
      opt_apply_cuts = false;
    }
    else if (strcmp(argv[k],"-bss") == 0) {
      opt_bss = true;
      opt_apply_cuts = false;
    }
    else if (strcmp(argv[k],"-bb") == 0) {
      opt_bb = true;
    }
    else if (strcmp(argv[k],"-dfs") == 0) {
      opt_dfs = true;
    }
    else if (strcmp(argv[k],"-cc") == 0) {
      opt_cc = true;
    }
    else if (strcmp(argv[k],"-tt") == 0) {
      opt_tt = true;
    }
    else if ((strcmp(argv[k],"-tt-size") == 0) && (k < argc - 1)) {
      opt_tt_size = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-all") == 0) {
      opt_find_all = true;
    }
    else if (strcmp(argv[k],"-all-different") == 0) {
      opt_find_all = true;
      opt_all_different = true;
    }
    else if (strcmp(argv[k],"-ex") == 0) {
      opt_find_all = false;
      opt_exhaustive = true;
    }
    else if (strcmp(argv[k],"-ff") == 0) {
      opt_ff = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
    }
    else if (strcmp(argv[k],"-ff-P2") == 0) {
      opt_ff = true;
      opt_ff_P2 = true;
      opt_pdb = false;
      opt_find_invariants = false;
      opt_verify_invariants = false;
    }

    // limit-setting options
    else if ((strcmp(argv[k],"-t") == 0) && (k < argc - 1)) {
      time_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-y") == 0) && (k < argc - 1)) {
      memory_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-i") == 0) && (k < argc - 1)) {
      iteration_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-c") == 0) && (k < argc - 1)) {
      cost_limit = A_TO_N(argv[++k]);
    }
    else if ((strcmp(argv[k],"-n") == 0) && (k < argc - 1)) {
      node_limit = atoi(argv[++k]);
    }

    // output/formatting options
    else if (strcmp(argv[k],"-val") == 0) {
      opt_validate = true;
    }
    else if (strcmp(argv[k],"-schedule") == 0) {
      opt_schedule = true;
    }
    else if (strcmp(argv[k],"-post") == 0) {
      opt_schedule = true;
      opt_post_op = true;
    }
    else if (strcmp(argv[k],"-no-plan") == 0) {
      opt_print_plan = false;
    }
    else if (strcmp(argv[k],"-print-on-store") == 0) {
      opt_print_plan = false;
      opt_print_plan_on_store = true;
    }
    else if (strcmp(argv[k],"-save-plan-to-file") == 0) {
      opt_print_plan = false;
      opt_save_plan_to_file = true;
    }
    else if (strcmp(argv[k],"-pddl") == 0) {
      opt_print_plan = true;
      opt_pddl = true;
      opt_ipc = false;
    }
    else if (strcmp(argv[k],"-ipc") == 0) {
      opt_print_plan = true;
      opt_ipc = true;
      opt_pddl = false;
    }
    else if (strcmp(argv[k],"-strict-ipc") == 0) {
      opt_print_plan = true;
      opt_ipc = true;
      opt_strict_ipc = true;
      opt_pddl = false;
    }
    else if ((strcmp(argv[k],"-epsilon") == 0) && (k < argc - 1)) {
      epsilon = A_TO_N(argv[++k]);
    }
    else if (strcmp(argv[k],"-nsn") == 0) {
      HSPS::Instance::write_atom_set_with_symbolic_names = false;
      HSPS::Instance::write_action_set_with_symbolic_names = false;
    }

    // misc. options
    else if (strcmp(argv[k],"-save") == 0) {
      opt_save = true;
    }
    else if (strcmp(argv[k],"-sas-select") == 0) {
      opt_sas_select = true;
    }
    else if (strcmp(argv[k],"-sas-min") == 0) {
      opt_sas_min = true;
    }
    else if (strcmp(argv[k],"-sas-safe") == 0) {
      opt_sas_safe = true;
    }
    else if (strcmp(argv[k],"-ac") == 0) {
      opt_apx_clique = true;
      opt_maximal_add = false;
    }
    else if (((strcmp(argv[k],"-rnd") == 0) ||
	      (strcmp(argv[k],"-r") == 0)) &&
	     (k < argc - 1)) {
      rng.seed(atoi(argv[++k]));
    }
    else if ((strcmp(argv[k],"-rnd-pid") == 0) ||
	     (strcmp(argv[k],"-rp") == 0)) {
      rng.seed_with_pid();
    }
    else if ((strcmp(argv[k],"-rnd-time") == 0) ||
	     (strcmp(argv[k],"-rt") == 0)) {
      rng.seed_with_time();
    }
    else if ((strcmp(argv[k],"-exclude") == 0) && (k < argc - 1)) {
      char* tag = argv[++k];
      if (strcmp(tag, "all") == 0) {
	HSPS::PDDL_Base::exclude_all_dkel_items = true;
      }
      else {
	const HSPS::StringTable::Cell* c = symbols.find(tag);
	if (c) HSPS::PDDL_Base::excluded_dkel_tags.insert(c->text);
      }
    }
    else if ((strcmp(argv[k],"-require") == 0) && (k < argc - 1)) {
      const HSPS::StringTable::Cell* c = symbols.find(argv[++k]);
      if (c) HSPS::PDDL_Base::required_dkel_tags.insert(c->text);
    }
    else if (strcmp(argv[k],"-bfs-write-stats") == 0) {
      opt_bfs_write_stats = true;
    }
    else if (strcmp(argv[k],"-bfs-write-graph") == 0) {
      opt_bfs_write_graph = true;
    }
    else if (((strcmp(argv[k],"-e") == 0) ||
	      (strcmp(argv[k],"-ed") == 0)) &&
	     (k < argc - 1)) {
      opt_explore_to_depth = true;
      depth_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-eb") == 0) && (k < argc - 1)) {
      opt_explore_to_bound = true;
      cost_limit = A_TO_N(argv[++k]);
    }

    // input file
    else if (*argv[k] != '-') {
      parse_stats.start();
      reader->read(argv[k], false);
      parse_stats.stop();
    }
  }

  HSPS::SearchAlgorithm::default_trace_level = verbose_level;
  HSPS::Heuristic::default_trace_level = verbose_level - 1;
  HSPS::Instance::default_trace_level = verbose_level - 1;
  HSPS::Preprocessor::default_trace_level = verbose_level - 1;
  if (verbose_level < 1) opt_print_plan = false;
  if (verbose_level <= 0) HSPS::PDDL_Base::warning_level = 0;
  if (verbose_level > 1) HSPS::PDDL_Base::write_info = true;

  HSPS::Instance    instance;
  HSPS::cost_vec    saved_dur;
  HSPS::cost_vec    saved_res;
  HSPS::count_type  pre_search_nodes = 0;
  HSPS::Preprocessor prep(instance, prep_stats);
  NTYPE       root_est_cost = 0;
  bool        solved = false;
  bool        optimally = false;
  NTYPE       solution_cost = 0;
  HSPS::Store       store(instance);
  HSPS::HashTable* tt = 0; // so we can access table use stats at end

  stats.enable_interrupt();
  if (time_limit > 0) stats.enable_time_out(time_limit);
  if (memory_limit > 0) stats.enable_memory_limit(memory_limit);
  if (node_limit > 0) stats.enable_node_limit(node_limit);
  if (iteration_limit > 0) stats.enable_iteration_limit(iteration_limit);

  prep_stats.start();
  std::cerr << "instantiating..." << std::endl;
  reader->instantiate(instance);
  if (opt_resource && opt_compose_resources && (instance.n_resources() > 1)) {
    HSPS::mSubsetEnumerator crs(instance.n_resources(),
				composite_resource_size);
    bool more = crs.first();
    while (more) {
      HSPS::index_set s;
      crs.current_set(s);
      instance.create_composite_resource(s);
      more = crs.next();
    }
  }
  if (opt_cce) {
    if (opt_preprocess) {
      prep.preprocess_with_ce(false, true);
    }
    instance.compile_conditional_effects(prep_stats, true);
  }
  else if (opt_preprocess) {
    std::cerr << "preprocessing..." << std::endl;
    prep.preprocess(opt_preprocess_2 && !opt_rm_irrelevant);
    if (opt_rm_irrelevant) {
      prep.compute_irrelevant_atoms();
      prep.remove_irrelevant_atoms();
      if (opt_preprocess_2)
	prep.preprocess(true);
    }
    if (!instance.cross_referenced()) {
      std::cerr << "re-cross referencing..." << std::endl;
      instance.cross_reference();
    }
  }
  else {
    std::cerr << "cross referencing..." << std::endl;
    instance.cross_reference();
  }
#ifdef COMMENTED_OUT
  if (opt_schedule) {
    HSPS::bool_vec n_act(false, instance.n_actions());
    for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
      for (HSPS::index_type r = 0; r < instance.n_reusable(); r++)
	if (instance.actions[k].use[r] > instance.reusables[r].init)
	  n_act[k] = true;
    if (n_act.count(true) > 0) {
      HSPS::index_vec n_map;
      std::cerr << "removing " << n_act.count(true)
		<< " non-executable actions..."
		<< std::endl;
      instance.remove_actions(n_act, n_map);
      instance.clear_cross_reference();
      instance.cross_reference();
    }
  }
#endif
  if (opt_quick_find_invariants) {
    HSPS::graph* g_inc = prep.inconsistency_graph();
    prep.find_inconsistent_set_invariants(*g_inc);
    instance.add_missing_negation_invariants();
  }
  else if (opt_find_invariants) {
    prep.bfs_find_invariants();
  }
  if (opt_verify_invariants) {
    prep.verify_invariants(*(prep.inconsistency()));
    prep.remove_unverified_invariants();
  }
  HSPS::index_set original_goal_atoms(instance.goal_atoms);
  if (opt_extend_goal) {
    HSPS::index_set new_goals;
    prep.implied_atom_set(instance.goal_atoms,new_goals,*prep.inconsistency());
    std::cerr << new_goals.length() << " implied goals found" << std::endl;
    if (new_goals.length() > 0) {
      HSPS::index_set new_goal(instance.goal_atoms);
      new_goal.insert(new_goals);
      instance.set_goal(new_goal);
    }
  }
  // scaling of action costs
  if (opt_scale_zero_cost > 1) {
    for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
      instance.actions[k].cost =
	(instance.actions[k].cost * opt_scale_zero_cost);
      if (instance.actions[k].cost == 0)
	instance.actions[k].cost = 1;
    }
  }
  prep_stats.stop();

  std::cerr << "instance " << instance.name << " built in "
	    << prep_stats.total_time() << " seconds" << std::endl;
  std::cerr << instance.n_atoms() << " atoms ("
	    << instance.goal_atoms.length() << " goals), "
	    << instance.n_resources() << " resources ("
	    << instance.n_reusable_resources() << " reusable, "
	    << instance.n_consumable_resources() << " consumable), "
	    << instance.n_actions() << " actions, "
	    << instance.n_invariants() << " invariants"
	    << std::endl;

  HSPS::Heuristic* cost_est = 0;
  HSPS::estimator_vec resource_est(0, 0);

  if (!stats.break_signal_raised()) { // instance construction finished
    h_stats.start();

    HSPS::ACF* cost = (opt_sum_time ?
		       (HSPS::ACF*)new HSPS::MakespanACF(instance) :
		       (opt_cost ?
			(HSPS::ACF*)new HSPS::CostACF(instance) :
			(opt_zero ? (HSPS::ACF*)new HSPS::ZeroACF() :
			 (HSPS::ACF*)new HSPS::UnitACF())));

    if (opt_H0) {
      cost_est = new HSPS::ZeroHeuristic(instance);
    }

    else if (opt_HAlmost0) {
      cost_est = new HSPS::FwdUnitHeuristic(instance);
    }

    else if (opt_forward_rc) {
      HSPS::ForwardReachabilityCheck* fh =
	new HSPS::ForwardReachabilityCheck(instance, instance.goal_atoms);
      cost_est = fh;
    }

    else if (opt_forward_H1) {
      HSPS::ForwardH1* fh =
	new HSPS::ForwardH1(instance, instance.goal_atoms, *cost, h_stats);
      cost_est = fh;
    }

    else if (opt_forward_H2) {
      HSPS::ForwardH2* fh =
	new HSPS::ForwardH2(instance, instance.goal_atoms, *cost, h_stats);
      cost_est = fh;
    }

    else if (opt_lm_cut) {
      HSPS::ForwardLMCut* h =
	new HSPS::ForwardLMCut(instance, instance.goal_atoms, *cost, h_stats);
      cost_est = h;
    }

    else if (opt_lm_cut_2) {
      HSPS::ForwardLMCut2* h =
	new HSPS::ForwardLMCut2(instance, instance.goal_atoms, *cost, h_stats);
      cost_est = h;
    }

    else if (opt_hplus) {
      HSPS::ForwardHPlus* h =
	new HSPS::ForwardHPlus(instance, *cost, instance.goal_atoms, h_stats);
      cost_est = h;
    }

    else if (opt_ff) {
      if (opt_ff_P2) {
	HSPS::Instance* i2 = new HSPS::Instance(instance.name);
	HSPS::s2index   s2;
	HSPS::index_vec am;
	i2->makeP2(instance, prep.inconsistency(), s2, am);
	i2->cross_reference();
	// note: using the same ACF for P^2 problem; this will probably
	// crash for any ACF other than unit cost
	HSPS::ForwardFF* h =
	  new HSPS::ForwardFF(*i2, *cost, i2->goal_atoms, h_stats);
	HSPS::ToP2Adapter* a = new HSPS::ToP2Adapter(instance, s2, *h);
	cost_est = a;
      }
      else {
	HSPS::ForwardFF* h =
	  new HSPS::ForwardFF(instance, *cost, instance.goal_atoms, h_stats);
	cost_est = h;
      }
    }

    // else if (opt_forward_sumx) {
    //   SumX* h = new SumX(instance, *cost, stats);
    //   cost_est = h;
    // }

    else if (opt_reverse_H1 || opt_reverse_H2 || opt_reverse_AH2) {
      std::cerr << "constructing reversed instance..." << std::endl;
      HSPS::Instance* c_instance = instance.copy();
      c_instance->complete_atom_negations();
      c_instance->cross_reference();
      HSPS::pair_vec pn;
      for (HSPS::index_type k = 0; k < c_instance->n_atoms(); k++)
	if (k < c_instance->atoms[k].neg)
	  pn.append(HSPS::index_pair(k, c_instance->atoms[k].neg));
      assert((pn.length() * 2) == c_instance->n_atoms());
      HSPS::mapping bm(c_instance->n_atoms());
      for (HSPS::index_type k = instance.n_atoms(); k < c_instance->n_atoms(); k++)
	bm[k] = HSPS::no_such_index;
      HSPS::AtomMapAdapter c_inc(*c_instance, bm, *(prep.inconsistency()));
      for (HSPS::index_type i = 0; i < pn.length(); i++) {
	if (INFINITE(c_inc.incremental_eval(c_instance->goal_atoms, pn[i].first))) {
	  c_instance->atoms[pn[i].second].goal = true;
	  c_instance->goal_atoms.insert(pn[i].second);
	}
	else if (INFINITE(c_inc.incremental_eval(c_instance->goal_atoms, pn[i].second))) {
	  c_instance->atoms[pn[i].first].goal = true;
	  c_instance->goal_atoms.insert(pn[i].first);
	}
      }
      for (HSPS::index_type k = 0; k < c_instance->n_actions(); k++) {
	for (HSPS::index_type i = 0; i < c_instance->actions[k].add.length(); i++)
	  if (INFINITE(c_inc.incremental_eval(c_instance->actions[k].pre,
					      c_instance->actions[k].add[i])))
	    c_instance->actions[k].pre.insert(c_instance->atoms[c_instance->actions[k].add[i]].neg);
	for (HSPS::index_type i = 0; i < c_instance->actions[k].del.length(); i++)
	  if (INFINITE(c_inc.incremental_eval(c_instance->actions[k].pre,
					      c_instance->atoms[c_instance->actions[k].del[i]].neg)))
	    c_instance->actions[k].pre.insert(c_instance->actions[k].del[i]);
      }

      HSPS::Instance* r_instance = new HSPS::Instance(instance.name);
      r_instance->reverse_copy(*c_instance);
      r_instance->cross_reference();
      std::cerr << "reversed instance: "
		<< r_instance->n_atoms() << " atoms ("
		<< r_instance->goal_atoms.length() << " goals), "
		<< r_instance->n_actions() << " actions, "
		<< std::endl;
      delete c_instance;
      std::cerr << "computing heuristic..." << std::endl;
      if (opt_reverse_AH2) {
	HSPS::AH* h = new HSPS::AH(*r_instance, h_stats);
	if (opt_kcut_partition) {
	  h->compute_with_k_cuts(*cost, r_instance->goal_atoms, opt_kcuts, true);
	  h->compute_max(*cost, true);
	}
	else if (opt_bu2_partition) {
	  h->compute_bottom_up_2
	    (*cost, r_instance->goal_atoms,
	     (HSPS::index_type)CEIL_TO_INT(bu2_ratio * r_instance->n_atoms()), true);
	  // implicit: h->compute_max(*cost, true);
	}
	else if (opt_bu_partition) {
	  h->compute_bottom_up(*cost, r_instance->goal_atoms, opt_bu_glue, true);
	  // implicit: h->compute_max(*cost, true);
	}
	else {
	  h->compute_with_iterative_assignment_2
	    (*cost, r_instance->goal_atoms, true, opt_fpia, r_instance->goal_atoms);
	}
	HSPS::CompleteNegationAdapter* rh =
	  new HSPS::CompleteNegationAdapter(instance, pn, *h);
	cost_est = rh;
      }
      else {
	HSPS::CostTable* h = new HSPS::CostTable(*r_instance, h_stats);
	if (opt_reverse_H1)
	  h->compute_H1(*cost);
	else
	  h->compute_H2(*cost);
	HSPS::CompleteNegationAdapter* rh =
	  new HSPS::CompleteNegationAdapter(instance, pn, *h);
	cost_est = rh;
      }
    }

    else if (opt_pdb_load) {
      std::cerr << "constructing SAS+ instance..." << std::endl;
      HSPS::SASInstance* sas_instance =
	new HSPS::SASInstance(instance, opt_sas_select, opt_sas_min, opt_sas_safe);
      HSPS::index_set goal_variables;
      sas_instance->goal_state.defined_set(goal_variables);
      std::cerr << goal_variables.length() << " of "
		<< sas_instance->n_variables()
		<< " variables with goal value"	<< std::endl;
      // std::cerr << "SAS instance:" << std::endl;
      // sas_instance->write_domain(std::cerr);
      HSPS::Heuristic* inc = (opt_pdb_inc ? prep.inconsistency() : 0);
      HSPS::MDDNode* sinc =
	(opt_pdb_inc ? HSPS::makeMDD(prep.inconsistency(),
				     sas_instance->atom_map_defined(),
				     sas_instance->atom_map_n()) : 0);
      HSPS::name_vec vnames(0, 0);
      sas_instance->variable_names(vnames);
      HSPS::index_set_vec sets;
      reader->export_sets(vnames, sets);
      HSPS::PDBCollection* h_col =
	new HSPS::PDBCollection(*sas_instance, *cost, sinc, inc, h_stats);
      for (HSPS::index_type k = 0; k < sets.length(); k++)
	h_col->addProgressionPDB(sets[k], opt_fast_pdb);
      std::cerr << h_col->n_patterns() << " PDBs, total size = "
		<< h_col->total_size() << std::endl;
      if (opt_save) {
	h_col->write_collection(std::cout);
	h_col->write_PDB(std::cout);
	exit(0);
      }
      cost_est = new HSPS::FwdSASHAdapter(instance, *sas_instance, *h_col);
    }

    else if (opt_ipdb) {
      std::cerr << "constructing SAS+ instance..." << std::endl;
      HSPS::SASInstance* sas_instance =
	new HSPS::SASInstance(instance, opt_sas_select, opt_sas_min, opt_sas_safe);
      HSPS::index_set goal_variables;
      sas_instance->goal_state.defined_set(goal_variables);
      std::cerr << h_stats.time() << " seconds, "
		<< goal_variables.length() << " of "
		<< sas_instance->n_variables()
		<< " variables with goal value"	<< std::endl;
      // std::cerr << "SAS instance:" << std::endl;
      // sas_instance->write_domain(std::cerr);
      HSPS::Heuristic* inc = (opt_pdb_inc ? prep.inconsistency() : 0);
      HSPS::MDDNode* sinc = (opt_pdb_inc ? HSPS::makeMDD(prep.inconsistency(),
					     sas_instance->atom_map_defined(),
					     sas_instance->atom_map_n()) : 0);
      HSPS::index_set vars;
      if (opt_pdb_spanning) {
	// std::cerr << "searching for independent sets..." << std::endl;
	// HSPS::IndependentVariableSets iv(*sas_instance);
	// if (opt_apx_clique) {
	//   iv.compute_approximate_independent_sets();
	// }
	// else {
	//   iv.compute_maximal_independent_sets();
	// }
	// std::cerr << stats.time() << "seconds, "
	// 	  << iv.length() << " independent sets found" << std::endl;
	// std::cerr << "searching for spanning sets..."
	// 	  << std::endl;
	// iv.compute_spanning_sets(true, opt_pdb_span_search_limit);
	// iv.union_set(vars);
	// std::cerr << stats.time() << "seconds, "
	// 	  << iv.length() << " spanning sets, "
	// 	  << vars.length() << " of " << sas_instance->n_variables()
	// 	  << " variables span the state space"
	// 	  << std::endl;
	std::cerr << "searching for spanning set..." << std::endl;
	HSPS::SpanningVariableSet svs(*sas_instance);
	std::cerr << svs.length() << " of " << sas_instance->n_variables()
		  << " variables in spanning set ("
		  << h_stats.time() << " seconds)"
		  << std::endl;
	vars.assign_copy(svs);
      }
      else {
	vars.fill(sas_instance->n_variables());
      }
      HSPS::UnitACF dummy;
      std::cerr << "variables: ";
      sas_instance->write_variable_set(std::cerr, vars);
      std::cerr << std::endl << "goal: ";
      HSPS::partial_state sg(sas_instance->goal_state, vars);
      sas_instance->write_partial_state(std::cerr, sg);
      std::cerr << std::endl;
      std::cerr << "searching for patterns..." << std::endl;
      double i_out = N_TO_D(NEG_INF);
      HSPS::index_type s_out = HSPS::no_such_index;
      if (opt_repeat > 1) {
	HSPS::MaxH* h_max = new HSPS::MaxH();
	for (HSPS::index_type k = 0; k < opt_repeat; k++) {
	  std::cerr << "collection #" << k + 1 << "..." << std::endl;
	  HSPS::PDBCollection* h_col =
	    HSPS::build_collection(instance, dummy, *sas_instance, *cost,
				   sinc, inc, vars,
				   opt_pdb_size, opt_total_size,
				   ipdb_param_d_skip, ipdb_param_s_max,
				   ipdb_param_i_min,
				   ipdb_param_n_trials, ipdb_param_n_samples,
				   opt_fast_pdb, i_out, s_out, rng, h_stats);
	  h_max->new_component(h_col);
	}
	cost_est = new HSPS::FwdSASHAdapter(instance, *sas_instance, *h_max);
      }
      else {
	HSPS::PDBCollection* h_col =
	  HSPS::build_collection(instance, dummy, *sas_instance, *cost,
				 sinc, inc, vars, opt_pdb_size, opt_total_size,
				 ipdb_param_d_skip, ipdb_param_s_max,
				 ipdb_param_i_min,
				 ipdb_param_n_trials, ipdb_param_n_samples,
				 opt_fast_pdb, i_out, s_out, rng, h_stats);
	std::cerr << "built PDB collection: " << h_col->n_patterns()
		  << " patterns, total size = " << h_col->total_size()
		  << ", best skipped extension score = " << i_out
		  << ", size of smallest too large PDB = " << s_out
		  << std::endl;
	cost_est = new HSPS::FwdSASHAdapter(instance, *sas_instance, *h_col);
      }
    }

    if (opt_pdb || opt_compare_h) {
      std::cerr << "constructing SAS+ instance..." << std::endl;
      HSPS::SASInstance* sas_instance =
	new HSPS::SASInstance(instance, opt_sas_select, opt_sas_min, opt_sas_safe);
      HSPS::index_set goal_variables;
      sas_instance->goal_state.defined_set(goal_variables);
      std::cerr << goal_variables.length() << " of "
		<< sas_instance->n_variables()
		<< " variables with goal value"	<< std::endl;

      if (verbose_level > 2) {
	for (HSPS::index_type k = 0; k < goal_variables.size(); k++) {
	  std::cerr << goal_variables[k] << ". ";
	  sas_instance->write_variable(std::cerr, sas_instance->
				       variables[goal_variables[k]]);
	  std::cerr << std::endl;
	}
      }

      HSPS::ProgressionPDBSize set_size_fcn(sas_instance->signature);
      HSPS::PDBHeuristic* h_pdb = 0;
      HSPS::Heuristic* inc = (opt_pdb_inc ? prep.inconsistency() : 0);
      HSPS::MDDNode* sinc = (opt_pdb_inc ? HSPS::makeMDD(prep.inconsistency(),
					     sas_instance->atom_map_defined(),
					     sas_instance->atom_map_n()) : 0);
      HSPS::SASHeuristic* h_sas = 0;

      if (opt_pdb_incremental) {
	HSPS::index_set vars;
	if (opt_pdb_spanning) {
	  // std::cerr << "searching for independent sets..." << std::endl;
	  // HSPS::IndependentVariableSets iv(*sas_instance);
	  // if (opt_apx_clique) {
	  //   iv.compute_approximate_independent_sets();
	  // }
	  // else {
	  //   iv.compute_maximal_independent_sets();
	  // }
	  // std::cerr << iv.length() << " independent sets found"
	  //  << std::endl;
	  // std::cerr << "searching for spanning sets..."
	  // 	    << std::endl;
	  // iv.compute_spanning_sets(true, opt_pdb_span_search_limit);
	  // iv.union_set(vars);
	  // std::cerr << iv.length() << " spanning sets: "
	  // 	    << vars.length() << " of " << sas_instance->n_variables()
	  // 	    << " variables span the state space"
	  // 	    << std::endl;
	  std::cerr << "searching for spanning set..." << std::endl;
	  HSPS::SpanningVariableSet svs(*sas_instance);
	  std::cerr << svs.length() << " of " << sas_instance->n_variables()
		    << " variables in spanning set ("
		    << h_stats.time() << " seconds)"
		    << std::endl;
	  vars.assign_copy(svs);
	}
	else {
	  vars.fill(sas_instance->n_variables());
	}
	std::cerr << "computing PDB:s..." << std::endl;
	HSPS::IncrementalProgressionPDB* h_inc =
	  new HSPS::IncrementalProgressionPDB(*sas_instance, h_stats);
	h_inc->compute_sets(vars, opt_pdb_size, opt_pdb_set_size, *cost,
			    sinc, inc, ipdb_param_d_skip,
			    ipdb_param_n_samples, rng);
// 	std::cerr << "time: " << stats.time() << std::endl;
// 	std::cerr << StateAbstraction::n_walks << " random walks, "
// 		  << StateAbstraction::n_cut_walks << " cut short"
// 		  << std::endl;
// 	exit(0);
	if (!stats.break_signal_raised())
	  h_inc->compute_progression_PDB(*cost, sinc, inc);
	if ((opt_pdb_ext_add || opt_pdb_add) &&
	    !stats.break_signal_raised())
	  h_inc->compute_additive_groups(opt_maximal_add);
	h_pdb = h_inc;
	if (opt_save && !stats.break_signal_raised()) {
	  h_pdb->write(std::cout);
	  exit(0);
	}
      }

      else {
	if (opt_pdb_load) {
	  HSPS::name_vec vnames(0, 0);
	  sas_instance->variable_names(vnames);
	  HSPS::index_set_vec sets;
	  reader->export_sets(vnames, sets);
	  h_pdb = new HSPS::PDBHeuristic(*sas_instance, h_stats);
	  h_pdb->assign_sets(sets);
	  if (verbose_level > 1) {
	    std::cerr << "PDB sets: ";
	    h_pdb->write_variable_sets(std::cerr);
	    std::cerr << std::endl;
	  }
	  if (opt_pdb_ext_add || opt_pdb_add)
	    h_pdb->compute_additive_groups(opt_maximal_add);
	}

	else if (opt_pdb_random_bin) {
	  HSPS::RandomBinPDB* h_bin =
	    new HSPS::RandomBinPDB(*sas_instance, h_stats);
	  h_bin->compute_sets(goal_variables, opt_pdb_size, set_size_fcn,
			      rng, opt_pdb_random_bin_swaps);
	  if (verbose_level > 0) {
	    std::cerr << "random bin sets: ";
	    h_bin->write_variable_sets(std::cerr);
	    std::cerr << std::endl;
	  }
	  if (opt_pdb_ext_add || opt_pdb_add)
	    h_bin->compute_additive_groups(opt_maximal_add);
	  h_pdb = h_bin;
	}

	else if (opt_pdb_random_independent_bin) {
	  HSPS::RandomIndependentBinPDB* h_bin =
	    new HSPS::RandomIndependentBinPDB(*sas_instance, h_stats);
	  h_bin->compute_sets(opt_pdb_size, set_size_fcn,
			      opt_pdb_spanning, opt_pdb_collapse,
			      rng, opt_pdb_random_bin_swaps);
	  if (verbose_level > 0) {
	    std::cerr << "random bin sets: ";
	    h_bin->write_variable_sets(std::cerr);
	    std::cerr << std::endl;
	  }
	  if (opt_pdb_ext_add || opt_pdb_add)
	    h_bin->compute_additive_groups(opt_maximal_add);
	  h_pdb = h_bin;
	}

	else if (opt_pdb_weighted_independent_bin) {
	  std::cerr << "searching for independent sets..." << std::endl;
	  HSPS::IndependentVariableSets iv(*sas_instance, goal_variables);
	  if (opt_apx_clique) {
	    iv.compute_approximate_independent_sets();
	  }
	  else {
	    iv.compute_maximal_independent_sets();
	  }
	  std::cerr << iv.length() << " independent sets found" << std::endl;
	  if (opt_pdb_spanning) {
	    std::cerr << "searching for spanning sets..." << std::endl;
	    iv.compute_spanning_sets(false, opt_pdb_span_search_limit);
	    HSPS::index_set uss;
	    iv.union_set(uss);
	    std::cerr << iv.length() << " spanning sets: "
		      << uss.length() << " of " << goal_variables.length()
		      << " (goal) variables included"
		      << std::endl;
	  }
	  if (opt_pdb_collapse) {
	    HSPS::index_set small;
	    HSPS::bool_vec  tbr(false, iv.length());
	    for (HSPS::index_type k = 0; k < iv.length(); k++)
	      if (set_size_fcn(iv[k]) <= opt_pdb_size) {
		small.insert(iv[k]);
		tbr[k] = true;
	      }
	    if (tbr.count(true) > 1) {
	      std::cerr << "collapsing " << tbr.count(true)
			<< " small sets..." << std::endl;
	      iv.remove(tbr);
	      iv.append(small);
	    }
	  }
	  HSPS::GoalStateInterference* sv_gsi =
	    new HSPS::GoalStateInterference(*sas_instance);
	  HSPS::InverseCGFraction* sv_icg =
	    new HSPS::InverseCGFraction(*sas_instance);
	  HSPS::SetValueSum* sv_sum =
	    new HSPS::SetValueSum(*sv_gsi, *sv_icg);
	  HSPS::WeightedBinPDB* h_wib =
	    new HSPS::WeightedBinPDB(*sas_instance, h_stats, *sv_sum);
	  h_wib->compute_sets(iv, opt_pdb_size, set_size_fcn);
	  if (verbose_level > 0) {
	    std::cerr << "weighted independent bin sets: ";
	    h_wib->write_variable_sets(std::cerr);
	    std::cerr << std::endl;
	  }
	  if (opt_pdb_ext_add || opt_pdb_add)
	    h_wib->compute_additive_groups(opt_maximal_add);
	  h_pdb = h_wib;
	}

	else {
	  HSPS::Max1PDB* h_1 = new HSPS::Max1PDB(*sas_instance, h_stats);
	  h_1->compute_sets(goal_variables);
	  if (opt_pdb_ext_add || opt_pdb_add)
	    h_1->compute_additive_groups(opt_maximal_add);
	  h_pdb = h_1;
	}

	std::cerr << "computing PDB:s..." << std::endl;
	h_pdb->compute_progression_PDB(*cost, sinc, inc);
	if (opt_save && !stats.break_signal_raised()) {
	  // h_pdb->write(std::cout);
	  for (HSPS::index_type k = 0; k < h_pdb->n_components(); k++) {
	    std::cout << "PDB #" << k << " ";
	    sas_instance->
	      write_variable_set(std::cout, h_pdb->variable_set(k));
	    std::cout << std::endl;
	    ((HSPS::AbstractionHeuristic*)h_pdb->component(k))->
	      write(std::cout);
	    std::cout << std::endl;
	    ((HSPS::AbstractionHeuristic*)h_pdb->component(k))->
	      write_graph(std::cout, sas_instance->goal_state, HSPS::EMPTYSET,
			  false, true, true);
	    std::cout << std::endl;
	  }
	  exit(0);
	}
      }

      if (opt_pdb_add && !stats.break_signal_raised()) {
	std::cerr << "creating additive groups..." << std::endl;
	if (verbose_level > 0) {
	  std::cerr << "additive groups: ";
	  h_pdb->write_additive_groups(std::cerr);
	  std::cerr << std::endl;
	}
	HSPS::MaxAddH* h_ma = h_pdb->make_max_add();
	h_sas = h_ma;
      }
      else {
	h_sas = h_pdb;
      }

      if (opt_compare_h) {
	HSPS::Heuristic* h_alt =
	  new HSPS::FwdSASHAdapter(instance, *sas_instance, *h_sas);
	cost_est = new HSPS::CompareEval(instance, *cost_est, *h_alt);
      }
      else {
	cost_est = new HSPS::FwdSASHAdapter(instance, *sas_instance, *h_sas);
      }
    }

    if (opt_resource && opt_res_ipdb && !stats.break_signal_raised()) {
      std::cerr << "computing resource estimators..." << std::endl;
      resource_est.assign_value(0, instance.n_resources());

      std::cerr << "constructing SAS+ instance..." << std::endl;
      HSPS::SASInstance* sas_instance =
	new HSPS::SASInstance(instance, opt_sas_select, opt_sas_min, opt_sas_safe);
      HSPS::index_set goal_variables;
      sas_instance->goal_state.defined_set(goal_variables);
      std::cerr << h_stats.time() << " seconds, "
		<< goal_variables.length() << " of "
		<< sas_instance->n_variables()
		<< " variables with goal value"	<< std::endl;
      HSPS::Heuristic* inc = (opt_pdb_inc ? prep.inconsistency() : 0);
      HSPS::MDDNode* sinc =
	(opt_pdb_inc ? HSPS::makeMDD(prep.inconsistency(),
				     sas_instance->atom_map_defined(),
				     sas_instance->atom_map_n()) : 0);
      HSPS::index_set vars;
      if (opt_pdb_spanning) {
	// std::cerr << "searching for independent sets..." << std::endl;
	// HSPS::IndependentVariableSets iv(*sas_instance);
	// if (opt_apx_clique) {
	//   iv.compute_approximate_independent_sets();
	// }
	// else {
	//   iv.compute_maximal_independent_sets();
	// }
	// std::cerr << stats.time() << "seconds, "
	// 	  << iv.length() << " independent sets found" << std::endl;
	// std::cerr << "searching for spanning sets..."
	// 	  << std::endl;
	// iv.compute_spanning_sets(true, opt_pdb_span_search_limit);
	// iv.union_set(vars);
	// std::cerr << stats.time() << "seconds, "
	// 	  << iv.length() << " spanning sets, "
	// 	  << vars.length() << " of " << sas_instance->n_variables()
	// 	  << " variables span the state space"
	// 	  << std::endl;
	std::cerr << "searching for spanning set..." << std::endl;
	HSPS::SpanningVariableSet svs(*sas_instance);
	std::cerr << svs.length() << " of " << sas_instance->n_variables()
		  << " variables in spanning set ("
		  << h_stats.time() << " seconds)"
		  << std::endl;
	vars.assign_copy(svs);
      }
      else {
	vars.fill(sas_instance->n_variables());
      }
      std::cerr << "variables: ";
      sas_instance->write_variable_set(std::cerr, vars);
      std::cerr << std::endl << "goal: ";
      HSPS::partial_state sg(sas_instance->goal_state, vars);
      sas_instance->write_partial_state(std::cerr, sg);
      std::cerr << std::endl;

      for (HSPS::index_type k = 0; k < instance.n_resources(); k++) {
	std::cerr << "computing PDB collection for resource " << k
		  << ": " << instance.resources[k].name << "..."
		  << std::endl;
	HSPS::ResourceConsACF* rce = new HSPS::ResourceConsACF(instance, k);
	std::cerr << "searching for patterns..." << std::endl;
	double i_out = N_TO_D(NEG_INF);
	HSPS::index_type s_out = HSPS::no_such_index;
	HSPS::PDBCollection* h_col =
	  HSPS::build_collection(instance, *rce, *sas_instance, *rce,
				 sinc, inc, vars, opt_pdb_size, opt_total_size,
				 ipdb_param_d_skip, ipdb_param_s_max,
				 ipdb_param_i_min,
				 ipdb_param_n_trials, ipdb_param_n_samples,
				 opt_fast_pdb, i_out, s_out, rng, h_stats);
	std::cerr << "built PDB collection: " << h_col->n_patterns()
		  << " patterns, total size = " << h_col->total_size()
		  << ", best skipped extension score = " << i_out
		  << ", size of smallest too large PDB = " << s_out
		  << std::endl;
	resource_est[k] =
	  new HSPS::FwdSASHAdapter(instance, *sas_instance, *h_col);
      }
    }

    h_stats.stop();
    if (!stats.break_signal_raised()) {
      std::cerr << "heuristic computed in " << stats.time()
		<< " seconds (" << stats << ")" << std::endl;
    }

    pre_search_nodes = stats.total_nodes();

    if (!stats.break_signal_raised()) { // heuristic finished
      HSPS::State* search_root = 0;

      if (opt_extended) {
	search_root =
	  new HSPS::ExtendedSeqProgState(instance, *cost_est, *cost,
					 instance.init_atoms);
      }
      else if (opt_ff) {
	// HSPS::ForwardFF* h_ff = (HSPS::ForwardFF*)cost_est;
	// search_root =
	//   new HSPS::RedSeqProgState(instance, *cost_est, *cost,
	// 			    instance.init_atoms, *h_ff);
	search_root = new HSPS::SeqProgState(instance, *cost_est, *cost,
					     instance.init_atoms);
      }
      else if (opt_relaxed) {
	HSPS::index_set x_atoms;
	for (HSPS::index_type k = 0; k < instance.n_invariants(); k++)
	  if ((instance.invariants[k].set.length() > 1) &&
	      (instance.invariants[k].lim == 1) &&
	      instance.invariants[k].exact)
	    x_atoms.insert(instance.invariants[k].set);
	search_root = new HSPS::RelaxedSeqProgState(instance, x_atoms,
						    *cost_est, *cost,
						    instance.init_atoms);
      }
      else {
	HSPS::BasicResourceState* rcs =
	  (opt_resource ? new HSPS::BasicResourceState(instance, resource_est) : 0);
	if (opt_apply_cuts)
	  search_root = new HSPS::SeqCProgState(instance, *cost_est, *cost,
						instance.init_atoms, rcs);
	else
	  search_root = new HSPS::SeqProgState(instance, *cost_est, *cost,
					       instance.init_atoms, rcs);
      }

      std::cerr << "search root: " << *search_root << std::endl;
      std::cerr << "estimated goal cost: "
		<< cost_est->eval(instance.init_atoms)
		<< std::endl;

      if (opt_explore_to_depth || opt_explore_to_bound) {
	HSPS::TreeStatistics tree_stats;
	if (node_limit > 0) tree_stats.node_limit = node_limit;
	HSPS::Tree* search_tree =
	  new HSPS::Tree(*search_root, stats, tree_stats);
	if (opt_explore_to_depth) {
	  for (HSPS::index_type k = 1; k <= depth_limit; k++) {
	    std::cerr << "exploring to depth " << k << "..." << std::endl;
	    search_tree->build(k);
	    std::cerr << stats << std::endl;
	  }
	}
	else if (opt_explore_to_bound) {
	  std::cerr << "exploring to bound " << cost_limit << "..."
		    << std::endl;
	  search_tree->build(cost_limit);
	  std::cerr << stats << std::endl;
	}
	std::cerr << stats << std::endl;
	tree_stats.write(std::cout);
	std::ofstream st_out("search_tree.dot");
	search_tree->write_dot(st_out, false, false);
	st_out.close();
	return 0;
      }

      if (opt_save_plan_to_file) {
	store.copy_to = new HSPS::SaveAllToFile(instance, &stats, "plan");
      }
      else if (opt_print_plan_on_store) {
	if (opt_ipc) {
	  store.copy_to = new HSPS::PrintIPC(instance, std::cout);
	}
	else if (opt_pddl) {
	  store.copy_to = new HSPS::PrintPDDL(instance, std::cout);
	}
	else {
	  store.copy_to = new HSPS::Print(instance, std::cout);
	}
      }

      HSPS::SearchAlgorithm* search = 0;
      HSPS::SearchResult* result = &store;
      store.set_stop_condition(HSPS::Result::stop_at_first);
      if (opt_bb || opt_bss) {
	store.set_stop_condition(HSPS::Result::stop_at_optimal);
	result = new HSPS::StoreMinCost(store);
      }
      else if (opt_all_different) {
	result = new HSPS::StoreDistinct(instance, store);
      }

      if (opt_find_all) {
	store.set_stop_condition(HSPS::Result::stop_at_all_optimal);
      }
      else if (opt_exhaustive) {
	store.set_stop_condition(HSPS::Result::stop_at_all);
      }

      if (opt_bfs) {
	if (opt_bfs_px) {
	  std::cerr << "using partial expansion A*..." << std::endl;
	  search = new HSPS::BFS_PX(search_stats, *result, px_threshold);
	}
	else {
	  std::cerr << "using A*..." << std::endl;
	  search = new HSPS::BFS(search_stats, *result);
	}
      }
      else if (opt_ff) {
	HSPS::BFS* bfs = new HSPS::BFS(search_stats, *result);
	bfs->greedy = true;
	search = bfs;
      }
      else if (opt_bfhs) {
	std::cerr << "using BFHS..." << std::endl;
	//search = new HSPS::BFHS(search_stats, *result, 10007);
	search = new HSPS::BFHS(search_stats, *result, 2000003);
      }
      else if (opt_bfida) {
	std::cerr << "using BFIDA..." << std::endl;
	search = new HSPS::BFIDA(search_stats, *result, 10007);
      }
      else if (opt_bfxd) {
	std::cerr << "using BFHS with Exponential Deepening..." << std::endl;
	search = new HSPS::BFHS_XD(search_stats, *result, 10007);
      }
      else if (opt_dfs) {
	std::cerr << "exhaustive DFS to " << cost_limit << "..." << std::endl;
	HSPS::DFS* df_search = new HSPS::DFS(search_stats, *result);
	df_search->set_cycle_check(opt_cc);
	df_search->set_upper_bound(cost_limit);
	search = df_search;
      }
      else if (opt_bb) {
	std::cerr << "using DFS branch-and-bound..." << std::endl;
	HSPS::DFS_BB* bb_search = 0;
	if (opt_tt) {
	  tt = new HSPS::HashTable(opt_tt_size);
	  bb_search = new HSPS::DFS_BB(search_stats, *result, tt);
	}
	else {
	  bb_search = new HSPS::DFS_BB(search_stats, *result);
	}
	bb_search->set_cycle_check(opt_cc);
	if (cost_limit < 0) cost_limit = POS_INF;
	bb_search->set_upper_bound(cost_limit);
	search = bb_search;
      }
      else if (opt_bss) {
	std::cerr << "using beam-stack search..." << std::endl;
	HSPS::BeamStackSearch* bss =
	  new HSPS::BeamStackSearch(search_stats, *result, 500);
	if (cost_limit < 0) {
	  std::cerr << "error: -bss without -c <limit>" << std::endl;
	  exit(1);
	}
	//bss->set_cost_limit(cost_limit);
	search = bss;
      }
      else if (opt_tt) {
	std::cerr << "using IDA* with transposition table..." << std::endl;
	tt = new HSPS::HashTable(opt_tt_size);
	HSPS::IDA* i_search = new HSPS::IDA(search_stats, *result, tt);
	i_search->set_cycle_check(opt_cc);
	search = i_search;
      }
      else {
	std::cerr << "using IDA*..." << std::endl;
	HSPS::IDA* i_search = new HSPS::IDA(search_stats, *result);
	i_search->set_cycle_check(opt_cc);
	search = i_search;
      }

      if (cost_limit >= 0) search->set_cost_limit(cost_limit);

      root_est_cost = search_root->est_cost();
      std::cerr << "searching..." << std::endl;
      solution_cost = search->start(*search_root);

      if (!stats.break_signal_raised()) {
	std::cerr << "search complete (" << stats << ")" << std::endl;
      }
      solved = search->solved();
      optimally = search->optimal();
      if (!search->optimal())
	solution_cost = search->cost();

      if (opt_bfs && (opt_bfs_write_graph || opt_bfs_write_stats)) {
	HSPS::NodeSet& g = ((HSPS::BFS*)search)->state_space();
	std::cerr << g.root_nodes().size() << " root nodes" << std::endl;
	// g.backup_costs();
	g.mark_solved();
	if (opt_bfs_write_graph) {
	  g.write_graph(std::cout);
	  g.write_graph_compact(std::cout);
	}
	if (opt_bfs_write_stats)
	  g.write_short(std::cout, instance.name);
      }
    }
  }

  NTYPE min_makespan = POS_INF;

  if (solved && opt_schedule) {
    for (HSPS::index_type p = 0; p < store.n_solutions(); p++) {
      HSPS::Schedule* plan = store.plan(p);
      std::cerr << "processing plan " << p << "..." << std::endl;
      plan->write(std::cerr);
      HSPS::graph s_prec;
      plan->deorder(s_prec);
      HSPS::index_vec s_acts(plan->step_actions());
      std::cerr << "actions:";
      for (HSPS::index_type k = 0; k < s_acts.length(); k++)
	std::cerr << " " << k << ":" << instance.actions[s_acts[k]].name;
      std::cerr << std::endl;
      std::cerr << "precedence graph: " << s_prec << std::endl;
      plan->schedule(s_acts, s_prec);
      plan->write(std::cerr);
      min_makespan = MIN(min_makespan, plan->makespan());
      optimally = false;
    }
  }

  if (solved && opt_post_op) {
    std::cerr << "post-optimizing:" << std::endl;

    // compute makespan heuristic
    std::cerr << "computing makespan heuristic..." << std::endl;
    HSPS::CostTable* cost_tab = new HSPS::CostTable(instance, h_stats);
    cost_tab->compute_H2C(HSPS::MakespanACF(instance), opt_resource);
    std::cerr << "lower bound = " << cost_tab->eval(instance.goal_atoms)
	      << ", upper bound = " << min_makespan
	      << std::endl;

    if (opt_resource && !stats.break_signal_raised()) {
      std::cerr << "computing resource estimators..." << std::endl;
      h_stats.start();
      for (HSPS::index_type k = 0; k < instance.n_resources(); k++) {
	HSPS::CostTable* rce = new HSPS::CostTable(instance, stats);
	if (opt_R2) {
	  rce->compute_H2(HSPS::ResourceConsACF(instance, k));
	}
	else {
	  rce->compute_H1(HSPS::ResourceConsACF(instance, k));
	}
	if (verbose_level > 2) {
	  std::cout << "resource " << instance.resources[k].name
		    << " estimator:" << std::endl;
	  rce->write(std::cout);
	}
	resource_est[k] = rce;
      }
      h_stats.stop();
    }

    // construct new search root and run DFS-BB
    if (!stats.break_signal_raised()) {
      HSPS::RegressionResourceState* rcs =
	(opt_resource? new HSPS::RegressionResourceState(instance, resource_est) : 0);
      HSPS::TemporalRSRegState* search_root =
	new HSPS::TemporalRSRegState(instance, *cost_tab, instance.goal_atoms, rcs);
      std::cerr << "searching..." << std::endl;
      store.clear();
      store.set_stop_condition(HSPS::Result::stop_at_all_optimal);
      HSPS::StoreMinCost result(store);
      HSPS::HashTable* tt = new HSPS::HashTable(opt_tt_size);
      HSPS::DFS_BB op_search(search_stats, result, tt);
      op_search.set_cycle_check(opt_cc);
      op_search.start(*search_root, min_makespan);
      std::cerr << "final cost: " << op_search.cost()
		<< " (improvement: "
		<< solution_cost - op_search.cost()
		<< " (abs.) "
		<< (solution_cost - op_search.cost())/op_search.cost()
		<< " (rel.), " << stats << ")" << std::endl;
      solution_cost = op_search.cost();
      optimally = op_search.optimal();
    }
  }

  if (solved && opt_validate) {
    for (HSPS::index_type p = 0; p < store.n_solutions(); p++) {
      HSPS::Schedule* plan = store.plan(p);
      HSPS::ExecTrace trace(instance);
      std::cerr << "simulating plan " << p << "..." << std::endl;
      bool ok = plan->simulate(&trace);
      if (ok) {
	std::cout << "plan " << p << " ok" << std::endl;
	if (instance.n_resources() > 0) {
	  HSPS::amt_vec peak;
	  trace.peak_resource_use(peak);
	  for (HSPS::index_type k = 0; k < instance.n_resources(); k++) {
	    NTYPE ratio = (peak[k] / instance.resources[k].init);
	    std::cout << instance.resources[k].name
		      << '\t' << peak[k]
		      << '\t' << instance.resources[k].init
		      << '\t' << PRINT_NTYPE(ratio*100)
		      << std::endl;
	  }
	}
      }
      else {
	std::cerr << "plan " << p << " failed" << std::endl;
      }
    }
  }

  if (opt_hplus) {
    HSPS::ForwardHPlus* h = (HSPS::ForwardHPlus*)cost_est;
    //h->print_stats(std::cerr);
  }

  if (opt_ipc || (verbose_level == 0)) {
    if (opt_ipc) {
      std::cout << "; Time " << HSPS::Stopwatch::seconds() << std::endl;
      std::cout << "; ParsingTime " << parse_stats.total_time() << std::endl;
      if (solved) {
	if (opt_cost) {
	  std::cout << "; NrActions " << std::endl;
	  std::cout << "; MakeSpan " << std::endl;
	  std::cout << "; MetricValue " << PRINT_NTYPE(solution_cost)
		    << std::endl;
	}
	else {
	  std::cout << "; NrActions " << solution_cost << std::endl;
	  std::cout << "; MakeSpan " << std::endl;
	  std::cout << "; MetricValue " << std::endl;
	}
	HSPS::PrintIPC print_plan(instance, std::cout);
	store.output(print_plan);
      }
      else {
	std::cout << "; Not Solved (" << stats.flags() << ")" << std::endl;
      }
    }

    else if (opt_pddl) {
      if (solved && opt_print_plan && !opt_print_plan_on_store) {
	HSPS::PrintPDDL print_plan(instance, std::cout);
	store.output(print_plan);
      }
    }

    if (opt_ipc || opt_pddl)
      std::cout << ";; stats: ";
    std::cout << instance.name
	      << ' ' << (solved ? 1 : 0)
	      << ' ' << (optimally ? 1 : 0)
	      << ' ' << PRINT_NTYPE(root_est_cost)
	      << ' ' << PRINT_NTYPE(solution_cost)
	      << ' ' << stats.total_nodes()
	      << ' ' << stats.total_time()
	      << ' ' << stats.peak_memory()
	      << ' ' << stats.peak_stack_size()
#ifdef RSS_FROM_PSINFO
	      << ' ' << stats.peak_total_size()
#endif
	      << ' ' << stats.complete_iterations()
	      << ' ' << stats.nodes_at_max_lower_bound()
	      << ' ' << stats.time()
#ifdef EVAL_EXTRA_STATS
	      << ' ' << HSPS::CostNode::eval_count
	      << ' ' << HSPS::CostNode::eval_rec_count
#endif
#ifdef AH_EXTRA_STATS
	      << ' ' << HSPS::AH::Hsum_wins/((double)(HSPS::AH::Hsum_wins + HSPS::AH::Hmax_wins + HSPS::AH::draws))
	      << ' ' << HSPS::AH::draws/((double)(HSPS::AH::Hsum_wins + HSPS::AH::Hmax_wins + HSPS::AH::draws))
#endif
#ifdef SEARCH_EXTRA_STATS
	      << ' ' << (HSPS::rminx_size/(double)HSPS::rminx_count)
	      << ' ' << (HSPS::rminc_succ_size_ratio/(double)HSPS::rminc_count)
	      << ' ' << (HSPS::rminx_succ/(double)HSPS::rminx_count)
	      << ' ' << (HSPS::trie_count/(double)HSPS::rminx_count)
	      << ' ' << (HSPS::trie_count == 0 ? 0 :
			 (HSPS::tries_applicable/(double)HSPS::trie_count))
	      << ' ' << (HSPS::trie_count == 0 ? 0 :
			 (HSPS::tries_within_bound/(double)HSPS::trie_count))
#endif
	      << std::endl;
    if (opt_ipc || opt_pddl) {
      std::cout << ";; ";
      std::cout << "(:heuristic (";
      for (HSPS::index_type k = 0; k < instance.goal_atoms.length(); k++)
	std::cout << instance.atoms[instance.goal_atoms[k]].name;
      std::cout << ") " << stats.max_lower_bound() << ")" << std::endl;
    }
  }

  else { // verbose_level > 0 && !opt_ipc
    if (solved) {
      std::cout << "solution cost: " << solution_cost;
      if (optimally)
	std::cout << " (optimal)";
      else
	std::cout << " (upper bound)";
      std::cout << std::endl;
      std::cout << store.n_solutions() << " solutions";
      if (opt_all_different) {
	std::cout << " (" << HSPS::StoreDistinct::n_discarded
		  << " equivalent solutions discarded)";
      }
      std::cout << std::endl;
    }
    else {
      std::cout << "no solution found" << std::endl;
    }
    stats.print_total(std::cout);
    std::cout << "parsing: " << parse_stats << std::endl;
    std::cout << "preprocessing: " << prep_stats << std::endl;
    std::cout << "heuristic: " << h_stats << std::endl;
    std::cout << "search: " << search_stats << std::endl;
    std::cerr << "highest lower bound " << stats.max_lower_bound()
	      << " proven at " << stats.nodes_at_max_lower_bound()
	      << " nodes" << std::endl;
    double total_t = HSPS::Stopwatch::seconds();
    std::cout << total_t << " seconds total ("
	      << total_t - stats.total_time()
	      << " sec. not accounted for)" << std::endl;
#ifdef AH_EXTRA_STATS
    if (opt_reverse_AH2) {
      std::cout << "AH: " << HSPS::AH::Hmax_wins << " Hmax, "
		<< HSPS::AH::Hsum_wins << " Hsum, "
		<< HSPS::AH::draws << " equal ("
		<< (HSPS::AH::Hsum_wins/((double)(HSPS::AH::Hsum_wins + HSPS::AH::Hmax_wins + HSPS::AH::draws)))*100
		<< "% improved, "
		<< (HSPS::AH::draws/((double)(HSPS::AH::Hsum_wins + HSPS::AH::Hmax_wins + HSPS::AH::draws)))*100
		<< "% equal)" << std::endl;
    }
#endif
    if (opt_tt) {
      assert(tt);
      std::cout << "transposition table stats: TUF = " << tt->TUF()
		<< ", HCF: " << tt->HCF() << std::endl;
    }
#ifdef SEARCH_EXTRA_STATS
    std::cout << "search space extra stats:"
	      << ' ' << (HSPS::rminx_size/(double)HSPS::rminx_count)
	      << ' ' << (HSPS::rminc_succ_size_ratio/(double)HSPS::rminc_count)
	      << ' ' << (HSPS::rminx_succ/(double)HSPS::rminx_count)
	      << ' ' << (HSPS::trie_count/(double)HSPS::rminx_count)
	      << ' ' << (HSPS::trie_count == 0 ? 0 :
			 (HSPS::tries_applicable/(double)HSPS::trie_count))
	      << ' ' << (HSPS::trie_count == 0 ? 0 :
			 (HSPS::tries_within_bound/(double)HSPS::trie_count))
	      << std::endl;
#endif
    if (solved && opt_print_plan && !opt_print_plan_on_store) {
      if (opt_pddl) {
	HSPS::PrintPDDL print_plan(instance, std::cout);
	store.output(print_plan);
      }
      else {
	for (HSPS::index_type p = 0; p < store.n_solutions(); p++) {
	  HSPS::Schedule* plan = store.plan(p);
	  HSPS::State* path_end_state = store.solution(p);
	  std::cout << "plan #" << p << ":" << std::endl;
	  HSPS::Print print_plan(instance, std::cout);
	  plan->output(print_plan);
	}
      }
    }
  }

  return 0;
}
