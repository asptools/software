#ifndef SEARCH_BASE_H
#define SEARCH_BASE_H

#include "config.h"
#include "search.h"
#include "stats.h"
#include "name.h"

BEGIN_HSPS_NAMESPACE

class SearchResult {
 public:
  SearchResult() { };
  virtual ~SearchResult();

  // non-strict lower bound proven (i.e., if any
  // solution exists, its cost is at least 'cost')
  virtual void lower_bound(NTYPE cost) { };

  // solution found; the state holds the solution path
  virtual void solution(State& s, NTYPE cost) { };

  // solution found; the state does not hold the correct
  // solution path, which is instead provided by the
  // transition (path) argument
  virtual void solution(State& s, Transition* p, NTYPE cost) { };

  // there are no more solutions (i.e., 'cost' is a
  // strict lower bound)
  virtual void no_more_solutions(NTYPE cost) { };

  // return true iff more solutions are desired
  virtual bool  more() = 0;

  virtual void write(std::ostream& to);
};

class Result : public SearchResult {
 public:
  enum stop_condition { stop_at_first,
			stop_at_optimal,
			stop_at_nth,
			stop_at_all_optimal,
			stop_at_all };

 private:
  stop_condition sc;
  count_type n_to_find;

  count_type n_sol;
  NTYPE      min_cost;
  NTYPE      max_lb;
  NTYPE      max_ex;

  PlanSet*   plans;

 public:
  Result() : sc(stop_at_optimal), n_to_find(1), n_sol(0),
    min_cost(POS_INF), max_lb(NEG_INF), max_ex(NEG_INF), plans(0) { };
  Result(PlanSet* ps) : sc(stop_at_optimal), n_to_find(1), n_sol(0),
    min_cost(POS_INF), max_lb(NEG_INF), max_ex(NEG_INF), plans(ps) { };
  virtual ~Result() { };

  void set_stop_condition(stop_condition c);
  void set_n_to_find(count_type n);
  void set_plan_set(PlanSet* ps);

  count_type solution_count();
  bool       search_space_exhausted();

  void reset();

  virtual void lower_bound(NTYPE cost);
  virtual void solution(State& s, NTYPE cost);
  virtual void solution(State& s, Transition* p, NTYPE cost);
  virtual void no_more_solutions(NTYPE cost);

  // return true iff more solutions are desired
  virtual bool  more();

  virtual void write(std::ostream& to);
};

class Statistics : public Stopwatch
{
  count_type min_nodes_created;
  count_type max_nodes_created;
  count_type min_nodes_expanded;
  count_type max_nodes_expanded;
  count_type iterations_started;
  count_type iterations_finished;
  count_type total_min_nodes_created;
  count_type total_max_nodes_created;
  count_type total_min_nodes_expanded;
  count_type total_max_nodes_expanded;
  count_type total_iterations_started;
  count_type total_iterations_finished;
  index_type max_depth;
  NTYPE      max_lb;
  count_type nodes_to_prove_lb;
  count_type total_eval_count;
  count_type start_eval_count;

  count_type node_expansion_limit;
  count_type node_creation_limit;
  index_type depth_limit;
  count_type iteration_limit;
  count_type eval_limit;

  static count_type last_eval_count;

 protected:
  virtual void* downcast(long id);
  virtual void check_for_update();
  virtual void update_flags(unsigned long changed);

  virtual void on_start();
  virtual void on_stop();
  virtual void on_reset();

 public:
  static bool running_print_max;

  Statistics() : min_nodes_created(0), max_nodes_created(0),
    min_nodes_expanded(0), max_nodes_expanded(0),
    iterations_started(0), iterations_finished(0),
    total_min_nodes_created(0), total_max_nodes_created(0),
    total_min_nodes_expanded(0), total_max_nodes_expanded(0),
    total_iterations_started(0), total_iterations_finished(0),
    max_depth(0), max_lb(0), nodes_to_prove_lb(0),
    total_eval_count(0), start_eval_count(0),
    node_expansion_limit(0), node_creation_limit(0),
    depth_limit(0), iteration_limit(0), eval_limit(0)
    { };

  Statistics(Stopwatch* p) : Stopwatch(p),
    min_nodes_created(0), max_nodes_created(0),
    min_nodes_expanded(0), max_nodes_expanded(0),
    iterations_started(0), iterations_finished(0),
    total_min_nodes_created(0), total_max_nodes_created(0),
    total_min_nodes_expanded(0), total_max_nodes_expanded(0),
    total_iterations_started(0), total_iterations_finished(0),
    max_depth(0), max_lb(0), nodes_to_prove_lb(0),
    total_eval_count(0), start_eval_count(0),
    node_expansion_limit(0), node_creation_limit(0),
    depth_limit(0), iteration_limit(0), eval_limit(0)
    { };

  void create_node(State& s);
  void expand_node(State& s);
  void current_lower_bound(NTYPE b);
  void begin_iteration();
  void end_iteration();

  static const long FLAG_NODES_EXP = 64;
  static const long FLAG_NODES_CRT = 128;
  static const long FLAG_DEPTH = 256;
  static const long FLAG_ITERATIONS = 512;
  static const long FLAG_EVAL = 1024;

  static const long STATISTICS_FLAGS = (FLAG_NODES_EXP | FLAG_NODES_CRT
					| FLAG_DEPTH | FLAG_ITERATIONS
					| FLAG_EVAL);

  void enable_node_expansion_limit(count_type l);
  void disable_node_expansion_limit();
  void enable_node_creation_limit(count_type l);
  void disable_node_creation_limit();
  // set/clear both node limits (to same value)
  void enable_node_limit(count_type l);
  void disable_node_limit();
  void enable_depth_limit(index_type l);
  void disable_depth_limit();
  // note: iteration limit is > on started
  void enable_iteration_limit(count_type l);
  void disable_iteration_limit();
  void enable_eval_limit(count_type l);
  void disable_eval_limit();

  virtual void print(::std::ostream& s) const;


  // methods to get additional information

  double branching_factor() const {
    return ((min_nodes_created + max_nodes_created)/
	    ((double)(min_nodes_expanded + max_nodes_expanded)));
  };

  count_type nodes() const {
    return (min_nodes_expanded + max_nodes_expanded);
  };

  index_type peak_depth() const {
    return max_depth;
  };

  count_type total_nodes() const {
    if (running())
      return (total_min_nodes_expanded + total_max_nodes_expanded +
	      min_nodes_expanded + max_nodes_expanded);
    else
      return (total_min_nodes_expanded + total_max_nodes_expanded);
  };

  count_type total_nodes_created() const {
    if (running())
      return (total_min_nodes_created + total_max_nodes_created +
	      min_nodes_created + max_nodes_created);
    else
      return (total_min_nodes_created + total_max_nodes_created);
  };

  count_type total_min_nodes() const {
    if (running()) return (total_min_nodes_expanded + min_nodes_expanded);
    else return total_min_nodes_expanded;
  };

  count_type total_max_nodes() const {
    if (running()) return (total_max_nodes_expanded + max_nodes_expanded);
    else return total_max_nodes_expanded;
  };

  NTYPE max_lower_bound() const {
    return max_lb;
  };

  count_type nodes_at_max_lower_bound() const {
    return nodes_to_prove_lb;
  };

  count_type iterations() {
    return iterations_started;
  };

  count_type complete_iterations() {
    return iterations_finished;
  };

  count_type evaluations() const;

  count_type total_iterations() {
    if (running()) return iterations_started + total_iterations_started;
    else return total_iterations_started;
  };

  count_type total_complete_iterations() {
    if (running()) return iterations_finished + total_iterations_finished;
    else return total_iterations_finished;
  };

  void print_brief(::std::ostream& s, const char* p = 0);
  void print_total(::std::ostream& s, const char* p = 0);
};

#ifdef SEARCH_EXTRA_STATS
extern count_type rmaxx_count;
extern count_type rmaxx_size;
extern count_type rmaxx_succ;
extern count_type rmaxc_count;
extern count_type rminc_count;
extern count_type rminc_size;
extern count_type rminc_succ_size;
extern double rminc_succ_size_ratio;
extern count_type rminx_count;
extern count_type rminx_size;
extern count_type rminx_succ;
extern count_type trie_count;
extern count_type tries_applicable;
extern count_type tries_within_bound;
#endif

class SearchAlgorithm : public Search {
  bool is_solved;
  bool is_optimal;
  NTYPE cost_limit;
 protected:
  Statistics& stats;
  SearchResult& result;
  const Name* problem_name;
  int   trace_level;

  void set_solved(bool s, bool o);
  void set_solved(bool s); // default: o = true
  // stats.reset + set_solved(false, false)
  void reset();

 public:
  static int default_trace_level;

  void set_problem_name(const Name* n);
  void set_trace_level(int level);

  void set_cost_limit(NTYPE c);
  NTYPE get_cost_limit() const;
  // returns solution cost or lower bound; implementation
  // is specific to each search algorithm
  virtual NTYPE cost() const = 0;

  SearchAlgorithm(Statistics& s,
		  SearchResult& r);
  SearchAlgorithm(Statistics& s,
		  SearchResult& r,
		  NTYPE limit);
  virtual ~SearchAlgorithm();

  virtual NTYPE start(State& s, NTYPE b) = 0;
  virtual NTYPE start(State& s) = 0;

  virtual bool solved() const;
  virtual bool optimal() const;
  virtual bool done() const;
};

class SingleSearchAlgorithm : public SearchAlgorithm {
 public:
  SingleSearchAlgorithm(Statistics& s,
			SearchResult& r)
    : SearchAlgorithm(s, r) { };
  SingleSearchAlgorithm(Statistics& s,
			SearchResult& r,
			NTYPE limit)
    : SearchAlgorithm(s, r, limit) { };
  virtual ~SingleSearchAlgorithm() { };

  virtual NTYPE resume() = 0;
};

class MultiSearchAlgorithm : public SearchAlgorithm {
 public:
  MultiSearchAlgorithm(Statistics& s,
		       SearchResult& r)
    : SearchAlgorithm(s, r) { };
  MultiSearchAlgorithm(Statistics& s,
		       SearchResult& r,
		       NTYPE limit)
    : SearchAlgorithm(s, r, limit) { };
  virtual ~MultiSearchAlgorithm() { };

  virtual NTYPE resume(State& s, NTYPE b) = 0;
  virtual NTYPE resume(State& s);
};

END_HSPS_NAMESPACE

#endif
