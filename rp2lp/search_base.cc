
#include "search_base.h"
#include "heuristic.h"

BEGIN_HSPS_NAMESPACE

SearchResult::~SearchResult()
{
  // done
}

void SearchResult::write(std::ostream& to)
{
  // default:
  to << "<SEARCH_RESULT>";
}

count_type Result::solution_count()
{
  return n_sol;
}

bool Result::search_space_exhausted()
{
  return (max_ex == POS_INF);
}

void Result::reset()
{
  n_sol = 0;
  min_cost = POS_INF;
  max_lb = NEG_INF;
  max_ex = NEG_INF;
}

void Result::set_stop_condition(stop_condition c)
{
  sc = c;
}

void Result::set_n_to_find(count_type n)
{
  n_to_find = n;
  sc = stop_at_nth;
}

void Result::set_plan_set(PlanSet* ps)
{
  plans = ps;
}

void Result::lower_bound(NTYPE cost)
{
  max_lb = MAX(max_lb, cost);
}

void Result::solution(State& s, NTYPE cost)
{
  n_sol += 1;
  min_cost = MIN(min_cost, cost);
  if (plans) {
    Plan* p = plans->new_plan();
    if (p) {
      s.insert_path(*p);
      p->end();
    }
  }
}

void Result::solution(State& s, Transition* p, NTYPE cost)
{
  n_sol += 1;
  min_cost = MIN(min_cost, cost);
  if (plans) {
    Plan* pl = plans->new_plan();
    if (pl) {
      if (p)
	p->insert_path(*pl);
      pl->end();
    }
  }
}

void Result::no_more_solutions(NTYPE cost)
{
  max_ex = MAX(max_ex, cost);
}

bool Result::more()
{
  switch (sc) {
  case stop_at_first:
    return false;
  case stop_at_optimal:
    return (max_lb < min_cost);
  case stop_at_nth:
    return (n_sol < n_to_find);
  case stop_at_all_optimal:
    return (max_ex < min_cost);
  case stop_at_all:
    return true;
  default:
    assert(0);
  }
}

void Result::write(std::ostream& to)
{
  to << "(stop = " << sc
     << ", nsol = " << n_sol
     << ", min_cost = " << min_cost
     << ", max_lb = " << max_lb
     << ", max_ex = " << max_ex
     << ")";
}

#ifdef SEARCH_EXTRA_STATS
count_type rmaxx_count = 0;
count_type rmaxx_size = 0;
count_type rmaxx_succ = 0;
count_type rmaxc_count = 0;
count_type rminc_count = 0;
count_type rminc_size = 0;
count_type rminc_succ_size = 0;
double rminc_succ_size_ratio = 0.0;
count_type rminx_count = 0;
count_type rminx_size = 0;
count_type rminx_succ = 0;
count_type trie_count = 0;
count_type tries_applicable = 0;
count_type tries_within_bound = 0;
#endif

bool Statistics::running_print_max = false;

void* Statistics::downcast(long id)
{
  if (id == STATISTICS_FLAGS)
    return this;
  else
    return 0;
}

void Statistics::create_node(State& s)
{
  if (s.is_max()) max_nodes_created += 1;
  else min_nodes_created += 1;
  if (enabled_breaks(FLAG_NODES_CRT))
    if ((min_nodes_created + max_nodes_created +
	 total_min_nodes_created + total_max_nodes_created)
	> node_creation_limit) {
      std::cerr << "node creation limit reached: "
		<< (min_nodes_created + max_nodes_created +
		    total_min_nodes_created + total_max_nodes_created)
		<< std::endl;
      set_flags(FLAG_NODES_CRT);
    }
  index_type d = s.depth();
  if (d > max_depth) {
    max_depth = d;
    if (running_print_max) {
      ::std::cerr << "max depth: " << d << " (" << *this << ")" << ::std::endl;
    }
    if (enabled_breaks(FLAG_DEPTH))
      if (max_depth > depth_limit)
	set_flags(FLAG_DEPTH);
  }
  Statistics* p = (Statistics*)parent_as(STATISTICS_FLAGS);
  if (p) {
    p->create_node(s);
    update_parent_flag();
  }
}

void Statistics::expand_node(State& s)
{
  if (s.is_max()) max_nodes_expanded += 1;
  else min_nodes_expanded += 1;
  if (enabled_breaks(FLAG_NODES_EXP))
    if (total_nodes() > node_expansion_limit) {
      std::cerr << "node expansion limit reached: "
		<< total_nodes() << std::endl;
      set_flags(FLAG_NODES_EXP);
    }
  Statistics* p = (Statistics*)parent_as(STATISTICS_FLAGS);
  if (p) {
    p->expand_node(s);
    update_parent_flag();
  }
}

void Statistics::current_lower_bound(NTYPE b)
{
  if (b > max_lb) {
    max_lb = b;
    nodes_to_prove_lb = total_nodes();
#ifdef PRINT_EVOLUTION
    ::std::cerr << "EVO " << max_lb
	      << " " << nodes()
	      << " " << time()
	      << ::std::endl;
#endif
  }
  Statistics* p = (Statistics*)parent_as(STATISTICS_FLAGS);
  if (p) {
    p->current_lower_bound(b);
  }
}

void Statistics::begin_iteration()
{
  iterations_started += 1;
  if (enabled_breaks(FLAG_ITERATIONS))
    if ((total_iterations_started + iterations_started) > iteration_limit) {
      std::cerr << "iteration limit reached: "
		<< (total_iterations_started + iterations_started)
		<< " / " << (total_iterations_finished + iterations_finished)
		<< std::endl;
      set_flags(FLAG_ITERATIONS);
    }
  Statistics* p = (Statistics*)parent_as(STATISTICS_FLAGS);
  if (p) {
    p->begin_iteration();
    update_parent_flag();
  }
}

void Statistics::end_iteration()
{
  iterations_finished += 1;
  if (enabled_breaks(FLAG_ITERATIONS))
    if ((total_iterations_finished + iterations_finished) >= iteration_limit) {
      std::cerr << "iteration limit reached: "
		<< (total_iterations_started + iterations_started)
		<< " / " << (total_iterations_finished + iterations_finished)
		<< std::endl;
      set_flags(FLAG_ITERATIONS);
    }
  Statistics* p = (Statistics*)parent_as(STATISTICS_FLAGS);
  if (p) {
    p->end_iteration();
    update_parent_flag();
  }
}

void Statistics::enable_node_expansion_limit(count_type l)
{
  clear_flags(FLAG_NODES_EXP);
  set_breaks(FLAG_NODES_EXP);
  node_expansion_limit = l;
}

void Statistics::disable_node_expansion_limit()
{
  clear_breaks(FLAG_NODES_EXP);
}

void Statistics::enable_node_creation_limit(count_type l)
{
  clear_flags(FLAG_NODES_CRT);
  set_breaks(FLAG_NODES_CRT);
  node_creation_limit = l;
}

void Statistics::disable_node_creation_limit()
{
  clear_breaks(FLAG_NODES_CRT);
}

void Statistics::enable_node_limit(count_type l)
{
  clear_flags(FLAG_NODES_EXP | FLAG_NODES_CRT);
  set_breaks(FLAG_NODES_EXP | FLAG_NODES_CRT);
  node_expansion_limit = l;
  node_creation_limit = l;
}

void Statistics::disable_node_limit()
{
  clear_breaks(FLAG_NODES_EXP | FLAG_NODES_CRT);
}

void Statistics::enable_depth_limit(index_type l)
{
  clear_flags(FLAG_DEPTH);
  set_breaks(FLAG_DEPTH);
  depth_limit = l;
}

void Statistics::disable_depth_limit()
{
  clear_breaks(FLAG_DEPTH);
}

void Statistics::enable_iteration_limit(count_type l)
{
  clear_flags(FLAG_ITERATIONS);
  set_breaks(FLAG_ITERATIONS);
  iteration_limit = l;
}

void Statistics::disable_iteration_limit()
{
  clear_breaks(FLAG_ITERATIONS);
}

void Statistics::enable_eval_limit(count_type l)
{
  clear_flags(FLAG_EVAL);
  set_breaks(FLAG_EVAL);
  eval_limit = l;
}

void Statistics::disable_eval_limit()
{
  clear_breaks(FLAG_EVAL);
}

count_type Statistics::last_eval_count = 0;

void Statistics::check_for_update()
{
  unsigned long changed = 0;
  if (Heuristic::eval_count > last_eval_count) {
    last_eval_count = Heuristic::eval_count;
    changed = (changed | FLAG_EVAL);
  }
  if (changed)
    update_flags(changed);
  Stopwatch::check_for_update();
}

void Statistics::update_flags(unsigned long changed)
{
  if (enabled_breaks(FLAG_EVAL))
    if (evaluations() > eval_limit)
      set_flags(FLAG_EVAL);
  Stopwatch::update_flags(changed);
}

void Statistics::on_start()
{
  min_nodes_created = 0;
  max_nodes_created = 0;
  min_nodes_expanded = 0;
  max_nodes_expanded = 0;
  iterations_started = 0;
  iterations_finished = 0;
  start_eval_count = Heuristic::eval_count;
}

void Statistics::on_stop()
{
  total_min_nodes_created += min_nodes_created;
  total_max_nodes_created += max_nodes_created;
  total_min_nodes_expanded += min_nodes_expanded;
  total_max_nodes_expanded += max_nodes_expanded;
  total_iterations_started += iterations_started;
  total_iterations_finished += iterations_finished;
  total_eval_count += (Heuristic::eval_count - start_eval_count);
  min_nodes_created = 0;
  max_nodes_created = 0;
  min_nodes_expanded = 0;
  max_nodes_expanded = 0;
  iterations_started = 0;
  iterations_finished = 0;
}

void Statistics::on_reset()
{
  min_nodes_created = 0;
  max_nodes_created = 0;
  min_nodes_expanded = 0;
  max_nodes_expanded = 0;
  iterations_started = 0;
  iterations_finished = 0;
  total_min_nodes_created = 0;
  total_max_nodes_created = 0;
  total_min_nodes_expanded = 0;
  total_max_nodes_expanded = 0;
  total_iterations_started = 0;
  total_iterations_finished = 0;
  max_depth = 0;
  max_lb = 0;
  nodes_to_prove_lb = 0;
  total_eval_count = 0;
  start_eval_count = Heuristic::eval_count;
}

count_type Statistics::evaluations() const
{
  if (running())
    return (total_eval_count +
	    (Heuristic::eval_count - start_eval_count));
  else
    return total_eval_count;
}

void Statistics::print(::std::ostream& s) const
{
  s << total_nodes() << " nodes, peak depth: " << peak_depth()
    << ", BF: " << branching_factor() << ", ";
  Stopwatch::print(s);
}

void Statistics::print_brief(::std::ostream& s, const char* p)
{
  if (p) s << p;
  s << time() << " sec, "
    << min_nodes_created << "/" << max_nodes_created << " n.c., "
    << min_nodes_expanded << "/" << max_nodes_expanded << " n.x., "
    << nodes()/time() << " nodes/sec." << ::std::endl;
}

void Statistics::print_total(::std::ostream& s, const char* p)
{
  if (p) s << p;
  s << "total time: " << total_time() << " seconds" << ::std::endl;
  if (running() > 0) {
    if (p) s << p;
    s << "total nodes created: "
      << total_min_nodes_created + min_nodes_created << " min / "
      << total_max_nodes_created + max_nodes_created << " max" << ::std::endl;
    if (p) s << p;
    s << "total nodes expanded: "
      << total_min_nodes_expanded + min_nodes_expanded << " min / "
      << total_max_nodes_expanded + max_nodes_expanded << " max ("
      << total_nodes()/total_time() << " nodes/sec.)" << ::std::endl;
  }
  else {
    if (p) s << p;
    s << "total nodes created: " << total_min_nodes_created << " min / "
      << total_max_nodes_created << " max" << ::std::endl;
    if (p) s << p;
    s << "total nodes expanded: " << total_min_nodes_expanded << " min / "
      << total_max_nodes_expanded << " max ("
      << total_nodes()/total_time() << " nodes/sec.)" << ::std::endl;
  }
  s << "peak memory use: " << peak_memory() << "k heap, "
    << peak_stack_size() << "k stack"
#ifdef RSS_FROM_PSINFO
    << ", " << stats.peak_total_size() << "k total"
#endif
    << std::endl;
  s << "flags: " << flags() << std::endl;
}


int SearchAlgorithm::default_trace_level = 0;

SearchAlgorithm::SearchAlgorithm
(Statistics& s, SearchResult& r)
  : stats(s),
    result(r),
    is_solved(false),
    is_optimal(false),
    cost_limit(POS_INF),
    trace_level(default_trace_level)
{
  // done
}

SearchAlgorithm::SearchAlgorithm
(Statistics& s, SearchResult& r, NTYPE limit)
  : stats(s),
    result(r),
    is_solved(false),
    is_optimal(false),
    cost_limit(limit),
    trace_level(default_trace_level)
{
  // done
}

SearchAlgorithm::~SearchAlgorithm()
{
  // done
}

void SearchAlgorithm::set_problem_name(const Name* n)
{
  problem_name = n;
}

void SearchAlgorithm::set_trace_level(int level)
{
  trace_level = level;
}

void SearchAlgorithm::set_cost_limit(NTYPE c)
{
  cost_limit = c;
}

NTYPE SearchAlgorithm::get_cost_limit() const
{
  return cost_limit;
}

void SearchAlgorithm::reset()
{
  stats.reset();
  set_solved(false, false);
}

bool SearchAlgorithm::solved() const
{
  return is_solved;
}

bool SearchAlgorithm::optimal() const
{
  return (is_solved && is_optimal);
}

bool SearchAlgorithm::done() const
{
  return ((solved() && !result.more()) || stats.break_signal_raised());
}

void SearchAlgorithm::set_solved(bool s, bool o)
{
  is_solved = s;
  is_optimal = o;
}

void SearchAlgorithm::set_solved(bool s)
{
  is_solved = s;
  is_optimal = true;
}

NTYPE MultiSearchAlgorithm::resume(State& s)
{
  s.reevaluate();
  return resume(s, s.est_cost());
}

END_HSPS_NAMESPACE
