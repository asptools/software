
#include "preprocess.h"
#include "parser.h"
#include "lmcut.h"
#include "forward.h"
#include "plans.h"
#include "bfs.h"

struct StepSearchResult : HSPS::SearchResult {
  HSPS::Instance& instance;
  HSPS::ActionSequence& global_plan;
  HSPS::index_set step_goals;
  HSPS::index_set lah_goals;
  HSPS::bool_vec end_state;
  HSPS::NodeSet* bfs_space;

  StepSearchResult(HSPS::Instance& ins, HSPS::ActionSequence& gp);
  virtual ~StepSearchResult() { };

  void setup(const HSPS::index_set& sg,
	     const HSPS::index_set& g,
	     HSPS::NodeSet* ss);

  virtual void solution(HSPS::State& s, NTYPE cost);
  virtual void solution(HSPS::State& s, HSPS::Transition* p, NTYPE cost);

  // return true iff more solutions are desired
  virtual bool  more() { return false; };
};


StepSearchResult::StepSearchResult
(HSPS::Instance& ins, HSPS::ActionSequence& gp)
  : instance(ins),
    global_plan(gp),
    end_state(false, ins.n_atoms()),
    bfs_space(NULL)
{
  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++)
    end_state[i] = instance.atoms[i].init;
}

void StepSearchResult::setup
(const HSPS::index_set& sg, const HSPS::index_set& g, HSPS::NodeSet* ss)
{
  step_goals.assign_copy(sg);
  lah_goals.assign_copy(g);
  bfs_space = ss;
}

void StepSearchResult::solution(HSPS::State& s, NTYPE cost)
{
  std::cerr << "step solution found!" << std::endl;
  assert(bfs_space);
  HSPS::Node* p_end = bfs_space->find_node(s);
  assert(p_end);
  HSPS::node_vec path;
  bfs_space->back_path_to_sequence(p_end, path);
  std::cerr << "step solution length is " << path.length() << std::endl;
  assert(path.length() > 0);
  HSPS::bool_matrix rel(false, path.length(), instance.n_atoms());
  HSPS::index_type p = path.length() - 1;
  for (HSPS::index_type i = 0; i < lah_goals.size(); i++)
    rel[p][lah_goals[i]] = true;
  while (p > 0) {
    HSPS::SeqProgTrans* t = (HSPS::SeqProgTrans*)path[p]->bp_trans;
    assert(t);
    assert(t->act < instance.n_actions());
    p -= 1;
    rel[p].assign_copy(rel[p + 1]);
    if (instance.actions[t->act].add.have_common_element(rel[p])) {
      for (HSPS::index_type i = 0; i < instance.actions[t->act].add.size(); i++)
	rel[p][instance.actions[t->act].add[i]] = false;
      for (HSPS::index_type i = 0; i < instance.actions[t->act].pre.size(); i++)
	rel[p][instance.actions[t->act].pre[i]] = true;
    }
  }
  HSPS::bool_vec rem(true, step_goals.size());
  p = 0;
  while ((p < path.length()) && (rem.first(true) != HSPS::no_such_index)) {
    HSPS::SeqProgState* sp = (HSPS::SeqProgState*)path[p]->state;
    assert(sp);
    if (p == 0) {
      end_state.assign_copy(sp->atom_set());
    }
    if (path[p]->bp_trans) {
      HSPS::SeqProgTrans* t = (HSPS::SeqProgTrans*)path[p]->bp_trans;
      if (instance.actions[t->act].add.have_common_element(rel[p])) {
	path[p]->bp_trans->insert(global_plan);
	std::cerr << "--> " << instance.actions[t->act].name << std::endl;
	end_state.subtract(instance.actions[t->act].del);
	end_state.insert(instance.actions[t->act].add);
      }
      else {
	std::cerr << "irrelevant: " << instance.actions[t->act].name
		  << std::endl;
      }
    }
    for (HSPS::index_type i = 0; i < step_goals.size(); i++)
      if (rem[i])
	if (sp->atom_set()[step_goals[i]])
	  rem[i] = false;
    if (rem.first(true) != HSPS::no_such_index)
      p += 1;
  }
  if (rem.first(true) == HSPS::no_such_index) {
    std::cerr << "step goals achieved in " << p << " steps" << std::endl;
  }
  else {
    std::cerr << "error: step goals not achieved!" << std::endl;
    exit(255);
  }
};

void StepSearchResult::solution(HSPS::State& s, HSPS::Transition* p, NTYPE cost)
{
  solution(s, cost);
}


int main(int argc, char *argv[]) {
  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  bool        opt_preprocess = true;
  bool        opt_preprocess_2 = false;
  bool        opt_rm_irrelevant = true;
  bool        opt_print_plan = true;
  HSPS::index_type opt_look_ahead = 0;
  bool        opt_pseudo_optimal = true;
  NTYPE       opt_scale_zero_cost = 1;
  NTYPE       opt_weight = 1;
  long        time_limit = 0;
  long        memory_limit = 0;
  HSPS::count_type node_limit = 0;
  int         verbose_level = 1;

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
    else if (strcmp(argv[k],"-rm") == 0) {
      opt_rm_irrelevant = true;
    }
    else if ((strcmp(argv[k],"-scale-zero-cost") == 0) && (k < argc - 1)) {
      opt_scale_zero_cost = A_TO_N(argv[++k]);
    }

    // limit-setting options
    else if ((strcmp(argv[k],"-t") == 0) && (k < argc - 1)) {
      time_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-y") == 0) && (k < argc - 1)) {
      memory_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-n") == 0) && (k < argc - 1)) {
      node_limit = atoi(argv[++k]);
    }

    // planner-specific options
    else if ((strcmp(argv[k],"-l") == 0) && (k < argc - 1)) {
      opt_look_ahead = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-w") == 0) && (k < argc - 1)) {
      opt_weight = A_TO_N(argv[++k]);
    }
    else if (strcmp(argv[k],"-f") == 0) {
      opt_pseudo_optimal = false;
    }

    // output/formatting options
    else if (strcmp(argv[k],"-no-plan") == 0) {
      opt_print_plan = false;
    }

    // misc options
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

  HSPS::Instance instance;
  HSPS::Preprocessor prep(instance, prep_stats);

  stats.enable_interrupt();
  if (time_limit > 0) stats.enable_time_out(time_limit);
  if (memory_limit > 0) stats.enable_memory_limit(memory_limit);
  if (node_limit > 0) stats.enable_node_expansion_limit(node_limit);

  prep_stats.start();
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
  prep_stats.stop();

  std::cerr << "instance " << instance.name << " built in "
	    << prep_stats.total_time() << " seconds" << std::endl;
  std::cerr << instance.n_atoms() << " atoms ("
	    << instance.goal_atoms.length() << " goals), "
	    << instance.n_actions() << " actions, "
	    << instance.n_invariants() << " invariants"
	    << std::endl;

  std::cerr << "computing goal order..." << std::endl;
  prep_stats.start();
  HSPS::graph lmg;
  HSPS::bool_vec goals(instance.goal_atoms, instance.n_atoms());
  HSPS::index_set fragile_goals;
  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++)
    if (goals[i]) {
      if (!instance.atoms[i].del_by.empty()) {
	goals[i] = false;
	fragile_goals.insert(i);
      }
    }
  prep.quick_landmark_graph(lmg, &goals, NULL, false);
  HSPS::bool_vec rem_goals(goals);
  HSPS::index_set_vec goal_strata;
  while (rem_goals.first(true) != HSPS::no_such_index) {
    HSPS::index_set& next_set = goal_strata.append();
    for (HSPS::index_type i = 0; i < instance.n_atoms(); i++)
      if (rem_goals[i]) {
	if (!lmg.predecessors(i).have_common_element(rem_goals))
	  next_set.insert(i);
      }
    assert(!next_set.empty());
    rem_goals.subtract(next_set);
  }
  HSPS::index_set stratified_goals;
  goal_strata.union_set(stratified_goals);
  prep_stats.stop();

  std::cerr << "finished in " << prep_stats.time() << " seconds" << std::endl;
  std::cerr << goal_strata.length() << " ordered goal sets" << std::endl;

  HSPS::CostACF action_cost(instance);
  HSPS::AnyACF modified_cost(instance.n_actions(), action_cost);
  if (opt_scale_zero_cost > 1) {
    for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
      if (action_cost(k) > 0)
	modified_cost.set_cost(k, action_cost(k) * opt_scale_zero_cost);
      else
	modified_cost.set_cost(k, 1);
    }
  }

  HSPS::ActionSequence the_plan;
  StepSearchResult result(instance, the_plan);
  HSPS::index_type next = 0;
  bool failed = false;
  while ((next < goal_strata.length()) &&
	 !stats.break_signal_raised() &&
	 !failed) {
    HSPS::index_set g(goal_strata[next]);
    for (HSPS::index_type k = next + 1;
	 (k <= (next + opt_look_ahead)) && (k < goal_strata.length()); k++)
      g.insert(goal_strata[k]);
    g.insert(fragile_goals);
    instance.set_goal(g);

    HSPS::bool_vec r_atoms(false, instance.n_atoms());
    HSPS::bool_vec r_acts(false, instance.n_actions());
    prep.compute_relevance(g, r_atoms, r_acts);
    for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
      instance.actions[k].sel = r_acts[k];

    std::cerr << "step " << next << ": g = ";
    instance.write_atom_set(std::cerr, g);
    std::cerr << ", " << r_acts.count(true)
	      << " of " << instance.n_actions()
	      << " actions relevant" << std::endl;

    HSPS::BFS search(search_stats, result, 1000007);
    result.setup(goal_strata[next], g, &(search.state_space()));
    if (opt_pseudo_optimal) {
      HSPS::ForwardLMCut2 h(instance, g, modified_cost, h_stats);
      HSPS::SeqProgState root(instance, h, modified_cost, result.end_state);
      search.set_weight(opt_weight);
      search.start(root);
    }
    else {
      HSPS::ForwardFF h(instance, modified_cost, g, h_stats);
      HSPS::SeqProgState root(instance, h, modified_cost, result.end_state);
      search.greedy = true;
      search.start(root);
    }
    if (search.solved()) {
      next += 1;
      NTYPE total_cost = 0;
      for (HSPS::index_type k = 0; k < the_plan.length(); k++)
	total_cost += action_cost(the_plan[k]);
      std::cerr << "total cost = " << total_cost << ", " << stats << std::endl;
    }
    else if (!stats.break_signal_raised()) {
      failed = true;
    }
  }

  bool solved = (!failed && !stats.break_signal_raised());

  std::cout << ";; stats: " << stats << std::endl;
  if (solved) {
    NTYPE total_cost = 0;
    for (HSPS::index_type k = 0; k < the_plan.length(); k++)
      total_cost += action_cost(the_plan[k]);
    std::cout << ";; plan cost = " << PRINT_NTYPE(total_cost) << std::endl;
    std::cout << ";; plan length = " << the_plan.length() << std::endl;
    if (opt_print_plan) {
      HSPS::PrintIPC print_plan(instance, std::cout);
      the_plan.output(print_plan);
    }
  }
  else {
    std::cout << ";; not solved" << std::endl;
  }

  return 0;
}
