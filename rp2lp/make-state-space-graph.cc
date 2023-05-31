
#include "problem.h"
#include "preprocess.h"
#include "parser.h"
#include "forward.h"
#include "ext_state.h"
#include "cost_table.h"
#include "plans.h"
#include "bfs.h"
#include "enumerators.h"
#include "ilb.h"
#include <fstream>

NTYPE eval_h_add
(const HSPS::index_set& g, const HSPS::cost_vec& h)
{
  NTYPE s = 0;
  for (HSPS::index_type i = 0; i < g.size(); i++)
    s += h[g[i]];
  return s;
}

void compute_h_add
(HSPS::Instance& ins, const HSPS::bool_vec& s, HSPS::cost_vec& h)
{
  h.assign_value(POS_INF, ins.n_atoms());
  for (HSPS::index_type i = 0; i < ins.n_atoms(); i++)
    if (s[i]) h[i] = 0;
  bool done = false;
  while (!done) {
    done = true;
    for (HSPS::index_type k = 0; k < ins.n_actions(); k++) {
      NTYPE newc = eval_h_add(ins.actions[k].pre, h) + 1;
      if (FINITE(newc)) {
	for (HSPS::index_type i = 0; i < ins.actions[k].add.size(); i++)
	  if (newc < h[ins.actions[k].add[i]]) {
	    h[ins.actions[k].add[i]] = newc;
	    done = false;
	  }
      }
    }
  }
}

int main(int argc, char *argv[]) {
  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  bool opt_apply_cuts = true;
  bool opt_preprocess = true;
  bool opt_preprocess_2 = true;
  bool opt_rm_irrelevant = false;
  int  heuristic = -1;
  bool opt_fvalue = false;
  int  verbose_level = 0;
  NTYPE w = 1;

  HSPS::Statistics  stats;
  stats.start();

  HSPS::Parser* reader = new HSPS::Parser(symbols);

  for (int k = 1; k < argc; k++) {
    // verbose level
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      verbose_level = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-h") == 0) && (k < argc - 1)) {
      heuristic = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-w") == 0) && (k < argc - 1)) {
      w = A_TO_N(argv[++k]);
    }
    else if (strcmp(argv[k],"-f") == 0) {
      opt_fvalue = true;
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
    else if (strcmp(argv[k],"-no-compact") == 0) {
      HSPS::PDDL_Base::compact_resource_effects = false;
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
    else if (strcmp(argv[k],"-rm") == 0) {
      opt_rm_irrelevant = true;
    }

    // input file
    else if (*argv[k] != '-') {
      reader->read(argv[k], false);
    }
  }

  HSPS::SearchAlgorithm::default_trace_level = verbose_level;
  HSPS::Heuristic::default_trace_level = verbose_level - 1;
  HSPS::Instance::default_trace_level = verbose_level - 1;
  HSPS::Preprocessor::default_trace_level = verbose_level - 1;
  if (verbose_level <= 0) HSPS::PDDL_Base::warning_level = 0;
  if (verbose_level > 1) HSPS::PDDL_Base::write_info = true;

  HSPS::Instance    instance;
  HSPS::Preprocessor prep(instance, stats);
  HSPS::Store       store(instance);

  stats.enable_interrupt();

  std::cerr << "instantiating..." << std::endl;
  reader->instantiate(instance);
  if (opt_preprocess) {
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

  std::cerr << "instance " << instance.name << " built in "
	    << stats.total_time() << " seconds" << std::endl;
  std::cerr << instance.n_atoms() << " atoms ("
	    << instance.goal_atoms.length() << " goals), "
	    << instance.n_resources() << " resources ("
	    << instance.n_reusable_resources() << " reusable, "
	    << instance.n_consumable_resources() << " consumable), "
	    << instance.n_actions() << " actions, "
	    << instance.n_invariants() << " invariants"
	    << std::endl;

  if (stats.break_signal_raised()) {
    exit(0);
  }

  HSPS::Heuristic* h0 = new HSPS::ZeroHeuristic(instance);
  HSPS::ACF* cost = new HSPS::UnitACF();

  HSPS::State* root =
    new HSPS::SeqProgState(instance, *h0, *cost, instance.init_atoms);

  store.set_stop_condition(HSPS::Result::stop_at_all);
  HSPS::SearchAlgorithm* search = new HSPS::BFS(stats, store);

  std::cerr << "searching..." << std::endl;
  NTYPE solution_cost = search->start(*root);

  if (stats.break_signal_raised()) {
    exit(0);
  }

  if (!stats.break_signal_raised()) {
    std::cerr << "search complete (" << stats << ")" << std::endl;
  }
  bool solved = search->solved();
  bool optimally = search->optimal();

  HSPS::NodeSet& g = ((HSPS::BFS*)search)->state_space();
  std::cerr << g.root_nodes().size() << " root nodes" << std::endl;
  g.mark_solved();

  std::cerr << "evaluating heuristic..." << std::endl;
  HSPS::node_vec v(0, 0);
  g.collect_nodes(v);
  if (heuristic == 5) {
    HSPS::cost_vec h_add;
    for (HSPS::index_type k = 0; k < v.size(); k++) {
      HSPS::SeqProgState* s = (HSPS::SeqProgState*)v[k]->state;
      assert(s != NULL);
      compute_h_add(instance, s->atom_set(), h_add);
      if (opt_fvalue)
	v[k]->est = (v[k]->acc + (w * eval_h_add(instance.goal_atoms, h_add)));
      else
	v[k]->est = eval_h_add(instance.goal_atoms, h_add);
    }
  }
  else if (heuristic > 0) {
    HSPS::Heuristic* h = 0;
    switch (heuristic) {
    case 1:
      h = new HSPS::ForwardH1(instance, instance.goal_atoms, *cost, stats);
      break;
    case 2:
      h = new HSPS::ForwardH2(instance, instance.goal_atoms, *cost, stats);
      break;
    case 3:
      h = new HSPS::ForwardFF(instance, *cost, instance.goal_atoms, stats);
      break;
    case 4:
      h = new HSPS::ForwardLMCut(instance, instance.goal_atoms, *cost, stats);
      break;
    default:
      std::cerr << "warning: invalid heuristic choice -- using h^1"
		<< std::endl;
      h = new HSPS::ForwardH1(instance, instance.goal_atoms, *cost, stats);
    }
    for (HSPS::index_type k = 0; k < v.size(); k++) {
      HSPS::SeqProgState* s = (HSPS::SeqProgState*)v[k]->state;
      assert(s != NULL);
      if (opt_fvalue)
	v[k]->est = (v[k]->acc + (w * h->eval(s->atom_set())));
      else
	v[k]->est = h->eval(s->atom_set());
    }
  }
  else if (heuristic == 0) {
    HSPS::ILB* hplus = new HSPS::ILB(instance, *cost, 0, stats);
    for (HSPS::index_type k = 0; k < v.size(); k++) {
      HSPS::SeqProgState* s = (HSPS::SeqProgState*)v[k]->state;
      assert(s != NULL);
      HSPS::index_set sk(s->atom_set());
      if (opt_fvalue)
	v[k]->est =
	  (v[k]->acc + (w * hplus->hplus(sk, instance.goal_atoms, POS_INF, 0)));
      else
	v[k]->est = hplus->hplus(sk, instance.goal_atoms, POS_INF, 0);
    }
  }
  else {
    std::cerr << "warning: no heuristic chosen -- using h*" << std::endl;
    for (HSPS::index_type k = 0; k < v.size(); k++) {
      if (opt_fvalue)
	v[k]->est = (v[k]->acc + v[k]->opt);
      else
	v[k]->est = v[k]->opt;
    }
  }
  //g.write_graph_rainbow(std::cout);
  g.write_graph(std::cout);

  return 0;
}
