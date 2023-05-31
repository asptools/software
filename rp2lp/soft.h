#ifndef SOFT_H
#define SOFT_H

#include "config.h"
#include "problem.h"
#include "heuristic.h"
#include "enumerators.h"
#include "exec.h"
#include "search_base.h"
#include "resource.h"
#include <list>

BEGIN_HSPS_NAMESPACE

class SoftInstance : public Instance {
 public:

  struct SoftGoal {
    const Name* name;
    index_set atoms;
    NTYPE     weight;
    void*     src;

    SoftGoal() : name(0), atoms(EMPTYSET), weight(0), src(0) { };

    SoftGoal& operator=(const SoftGoal& g) {
      name = g.name;
      atoms = g.atoms;
      weight = g.weight;
      src = g.src;
      return *this;
    };

    bool operator==(const SoftGoal& g) const {
      return ((atoms == g.atoms) &&
	      (weight == g.weight));
    };

    bool is_sat(const bool* s) const {
      for (index_type k = 0; k < atoms.length(); k++)
	if (!s[atoms[k]]) return false;
      return true;
    };

    bool is_sat_init(const Instance& ins) const {
      for (index_type k = 0; k < atoms.length(); k++)
	if (!ins.atoms[atoms[k]].init) return false;
      return true;
    };

    bool is_sat(const index_set& s) const {
      for (index_type k = 0; k < atoms.length(); k++)
	if (!s.contains(atoms[k])) return false;
      return true;
    };
  };

  typedef lvector<SoftGoal> soft_goal_vec;

  soft_goal_vec soft;
  index_set     hard;
  NTYPE         null_value;

  SoftInstance();
  SoftInstance(const Name* n);
  ~SoftInstance() { };

  // build (add to) instance
  SoftGoal& new_soft_goal();

  // access instance information
  index_type n_soft() const { return soft.length(); };
  index_type n_hard() const { return hard.length(); };

  // test if empty plan is valid - true iff there are no hard goals,
  // or all hard goals are satisfied in initial state
  bool empty_plan_valid();

  // return the value of the empty plan, equal to null_value + value
  // of any soft goals satisfied in the initial state; if the empty
  // plan is not valid, return value is NEG_INF.
  NTYPE empty_plan_value();

  // remapping needs to be done after preprocessing
  void remap_hard_goals(const index_vec& atom_map);
  void remap_soft_goals(const index_vec& atom_map);

  long integrify_weights();

  // Keyder & Geffner's compilation
  void compile_KG(Instance& ins);

  // compilation with direct costs; works only for monotonic,
  // single-atom soft goals
  bool compile_direct_cost(Instance& ins);

  void create_decision_problem(const bool_vec& sel, Instance& ins);
  void create_decision_problem(const bool_vec& sel, NTYPE b, Instance& ins);
  NTYPE compute_epsilon();

  NTYPE eval_goal_state(const index_set& s);
  NTYPE eval_goal_state(const index_set& s, index_set& g);

  NTYPE eval_plan(Schedule& s);
  void  eval_plan_set(ScheduleSet& s, cost_vec& v);

  virtual void write_problem_goal(std::ostream& s) const;
  virtual void write_problem_metric(std::ostream& s) const;
  void write_goal_value_expression(std::ostream& s) const;
  void write_soft_goal_set(std::ostream& s, const index_set& set) const;
  void write_soft_goal_set(std::ostream& s, const bool_vec& set) const;

  virtual void print(std::ostream& s) const;
};

class DecisionProblemEnumerator : public IterativeEnumerator {
  SoftInstance&    instance;
  SubsetEnumerator selected;
  Heuristic&       h_cost;
  NTYPE            nb_min;

  index_set        g_sel;
  index_set        a_sel;
  NTYPE            v_sel;
  NTYPE            c_sel;
  NTYPE            nb_sel;

  bool find_next(bool more);

 public:
  DecisionProblemEnumerator(SoftInstance& ins, Heuristic& h, NTYPE b);
  virtual ~DecisionProblemEnumerator();

  virtual bool first();
  virtual bool next();

  NTYPE current_value() const;
  NTYPE current_min_cost() const;
  NTYPE current_max_cost() const;
  NTYPE current_min_nb() const;
  NTYPE current_max_nb() const;
  const index_set& current_soft_goals() const;
  const index_set& current_goal_atoms() const;
  void  create_decision_problem(Instance& ins);
};

class MaxValueSearch {
 protected:

  struct option {
    NTYPE  goal_value;
    NTYPE  est_cost;
    index_set goals;
    State* root;

    option() : goal_value(0), est_cost(0), goals(EMPTYSET), root(0) { };
    index_type n_goals() const { return goals.length(); };
    NTYPE est_value() const { return (goal_value - est_cost); };
  };

  typedef std::list<option> option_list;
  typedef option_list::iterator option_p;

  SoftInstance& instance;
  ACF& cost;
  Statistics& stats;
  Result& res;
  MultiSearchAlgorithm* search;

  option_list options;
  NTYPE       lb;
  bool        solved_flag;
  int         trace_level;

  void insert_option_in_list(const option& o);
  void init_option_list();

  void make_empty_plan();

  virtual void init_option(const index_set& selected, option& o) = 0;
  virtual NTYPE explore_next_option() = 0;

 public:
  MaxValueSearch(SoftInstance& i, ACF& c, Statistics& s, Result& r);
  MaxValueSearch(SoftInstance& i, ACF& c, Statistics& s, Result& r,
		 index_type tt_size, bool use_cc);
  ~MaxValueSearch();

  static index_type print_options_max;
  void print_option_list(std::ostream& s);

  index_type n_options() const { return options.size(); };
  NTYPE best_option_estimated_value();
  index_type best_option_size();

  void  init();
  NTYPE main();
  bool solved() { return solved_flag; }
};

class MaxNetBenefit : public MaxValueSearch {
  Heuristic& heuristic;
  RegressionResourceState* root_rs;

 protected:
  virtual void init_option(const index_set& selected, option& o);
  virtual NTYPE explore_next_option();

 public:
  MaxNetBenefit(SoftInstance& i, ACF& c, Statistics& s, Result& r,
		Heuristic& h)
    : MaxValueSearch(i, c, s, r), heuristic(h), root_rs(0)
    { };
  MaxNetBenefit(SoftInstance& i, ACF& c, Statistics& s, Result& r,
		Heuristic& h, index_type tt_size, bool use_cc)
    : MaxValueSearch(i, c, s, r, tt_size, use_cc), heuristic(h), root_rs(0)
    { };
  MaxNetBenefit(SoftInstance& i, ACF& c, Statistics& s, Result& r,
		Heuristic& h, RegressionResourceState* rs)
    : MaxValueSearch(i, c, s, r), heuristic(h), root_rs(rs)
    { };
  MaxNetBenefit(SoftInstance& i, ACF& c, Statistics& s, Result& r,
		Heuristic& h, RegressionResourceState* rs,
		index_type tt_size, bool use_cc)
    : MaxValueSearch(i, c, s, r, tt_size, use_cc), heuristic(h), root_rs(rs)
    { };
  ~MaxNetBenefit() { };
};

END_HSPS_NAMESPACE

#endif
