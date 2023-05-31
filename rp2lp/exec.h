#ifndef PLAN_EXEC_H
#define PLAN_EXEC_H

#include "config.h"
#include "search.h"
#include "problem.h"
// #include "stats.h"
#include "heuristic.h"

BEGIN_HSPS_NAMESPACE

class ExecError {
 public:
  enum ErrorType { error_unknown,
		   error_unsatisfied_precondition,
		   error_incompatible_actions,
		   error_resource_conflict,
		   error_resource_shortage,
		   error_unachieved_goal,
		   warning_redundant_action
  };

  enum ErrorSeverity { severity_none = 0,
		       severity_warning = 1,
		       severity_plan_failure = 2,
		       severity_execution_failure = 3
  };

 protected:
  ErrorType  toe;
  NTYPE      at;
  index_type step;

 public:
  static const char* error_type_string(ErrorType t);
  static const char* error_severity_string(ErrorSeverity s);

  ExecError() : toe(error_unknown), at(0), step(no_such_index) { };
  ExecError(ErrorType e, NTYPE t) : toe(e), at(t), step(no_such_index) { };
  ExecError(ErrorType e, NTYPE t, index_type s) : toe(e), at(t), step(s) { };
  virtual ~ExecError();

  void remap_step(const index_vec& map);

  ErrorType   type_of_error() const { return toe; };
  ErrorSeverity severity_of_error() const;
  NTYPE       time_of_error() const { return at; };
  index_type  step_of_error() const { return step; };

  virtual ExecError* copy() const;
  virtual void write(::std::ostream& s) const;
};

typedef lvector<ExecError*> exec_error_vec;
typedef svector<ExecError::ErrorType> error_type_set;
typedef svector<ExecError::ErrorSeverity> error_severity_set;

class ExecErrorSet : public exec_error_vec {
  index_type current_step;
  error_type_set ignored_error_types;
 public:
  ExecErrorSet()
    : exec_error_vec((ExecError*)0, 0),
    current_step(no_such_index),
    ignored_error_types() { };
  ~ExecErrorSet();

  void ignore_error_type(ExecError::ErrorType t);
  void ignore_error_severity(ExecError::ErrorSeverity s);
  void clear_ignored_error_types();

  // construction methods
  void new_error(ExecError* e);
  bool ignore(ExecError::ErrorType t);
  void remap_steps(const index_vec& map);

  NTYPE                    earliest_time_of_error();
  ExecError::ErrorSeverity greatest_error_severity();
  bool                     executable();
  bool                     valid();

  index_type count_of_type(ExecError::ErrorType t);

  ExecErrorSet* earliest();
  ExecErrorSet* all_of_type(ExecError::ErrorType t);
  ExecErrorSet* all_of_severity(ExecError::ErrorSeverity s);
  ExecErrorSet* earliest_of_type(ExecError::ErrorType t);

  virtual void write(::std::ostream& s) const;
};

class UnsatisfiedPreconditionError : public ExecError {
  Instance&  ins;
  index_type act;
  index_type pre;
  index_set  holds;
 public:
  UnsatisfiedPreconditionError(Instance& i,
			       index_type a,
			       index_type p,
			       NTYPE at)
    : ExecError(ExecError::error_unsatisfied_precondition, at),
    ins(i), act(a), pre(p) { };
  UnsatisfiedPreconditionError(Instance& i,
			       index_type a,
			       index_type p,
			       const index_set& h,
			       NTYPE at)
    : ExecError(ExecError::error_unsatisfied_precondition, at),
    ins(i), act(a), pre(p), holds(h) { };
  UnsatisfiedPreconditionError(Instance& i,
			       index_type a,
			       index_type p,
			       NTYPE at,
			       index_type s)
    : ExecError(ExecError::error_unsatisfied_precondition, at, s),
    ins(i), act(a), pre(p) { };
  UnsatisfiedPreconditionError(Instance& i,
			       index_type a,
			       index_type p,
			       const index_set& h,
			       NTYPE at,
			       index_type s)
    : ExecError(ExecError::error_unsatisfied_precondition, at, s),
    ins(i), act(a), pre(p), holds(h) { };
  virtual ~UnsatisfiedPreconditionError() { };

  const Instance::Action& action()       { return ins.actions[act]; };
  const Instance::Atom&   precondition() { return ins.atoms[pre]; };
  const index_set&        true_atoms()   { return holds; };

  virtual ExecError* copy();
  virtual void write(::std::ostream& s) const;
};

class IncompatibleActionError : public ExecError {
  Instance&  ins;
  index_type act0;
  index_type act1;
 public:
  IncompatibleActionError(Instance& i, index_type a0, index_type a1, NTYPE at)
    : ExecError(ExecError::error_incompatible_actions, at),
    ins(i), act0(a0), act1(a1) { };
  IncompatibleActionError(Instance& i, index_type a0, index_type a1,
			  NTYPE at, index_type s)
    : ExecError(ExecError::error_incompatible_actions, at, s),
    ins(i), act0(a0), act1(a1) { };
  virtual ~IncompatibleActionError() { };

  const Instance::Action& action()              { return ins.actions[act0]; };
  const Instance::Action& incompatible_action() { return ins.actions[act1]; };

  virtual ExecError* copy();
  virtual void write(::std::ostream& s) const;
};

class ResourceConflictError : public ExecError {
  Instance&  ins;
  index_type res;
  index_set  c_acts;
  index_set  c_steps;
 public:
  ResourceConflictError(Instance& i,
			index_type r,
			const index_set& as,
			const index_set& ss,
			NTYPE at) :
    ExecError(ExecError::error_resource_conflict, at), ins(i), res(r),
    c_acts(as), c_steps(ss) { };
  ResourceConflictError(Instance& i,
			index_type r,
			NTYPE at) :
    ExecError(ExecError::error_resource_conflict, at), ins(i), res(r) { };
  ResourceConflictError(Instance& i,
			index_type r,
			const index_set& as,
			const index_set& ss,
			NTYPE at,
			index_type s) :
    ExecError(ExecError::error_resource_conflict, at, s), ins(i), res(r),
    c_acts(as), c_steps(ss) { };
  ResourceConflictError(Instance& i,
			index_type r,
			NTYPE at,
			index_type s) :
    ExecError(ExecError::error_resource_conflict, at, s), ins(i), res(r) { };
  void add_action(index_type a)
    { c_acts.insert(a); };
  void add_action(index_type a, index_type s)
    { c_acts.insert(a); c_steps.insert(s); };
  virtual ~ResourceConflictError() { };

  const Instance::Resource& resource() { return ins.resources[res]; };
  const index_set& conflict_actions()  { return c_acts; };
  const index_set& conflict_steps()    { return c_steps; };

  virtual ExecError* copy();
  virtual void write(::std::ostream& s) const;
};

class ResourceShortageError : public ExecError {
  Instance&  ins;
  index_type res;
  NTYPE      avail;
  index_type act;
 public:
  ResourceShortageError(Instance& i,
			index_type r,
			NTYPE v,
			index_type a,
			NTYPE at) :
    ExecError(ExecError::error_resource_shortage, at), ins(i), res(r),
    avail(v), act(a) { };
  ResourceShortageError(Instance& i,
			index_type r,
			NTYPE v,
			index_type a,
			NTYPE at,
			index_type s) :
    ExecError(ExecError::error_resource_shortage, at, s), ins(i), res(r),
    avail(v), act(a) { };
  virtual ~ResourceShortageError() { };

  const Instance::Resource& resource()  { return ins.resources[res]; };
  NTYPE                     available() { return avail; };
  const Instance::Action&   action()    { return ins.actions[act]; };

  virtual ExecError* copy();
  virtual void write(::std::ostream& s) const;
};

class UnachievedGoalError : public ExecError {
  Instance&  ins;
  index_type atom;
 public:
  UnachievedGoalError(Instance& i, index_type g, NTYPE at)
    : ExecError(ExecError::error_unachieved_goal, at), ins(i), atom(g) { };
  UnachievedGoalError(Instance& i, index_type g, NTYPE at, index_type s)
    : ExecError(ExecError::error_unachieved_goal, at, s), ins(i), atom(g) { };
  virtual ~UnachievedGoalError() { };

  const Instance::Atom&   goal() { return ins.atoms[atom]; };

  virtual ExecError* copy();
  virtual void write(::std::ostream& s) const;
};

class RedundantActionWarning : public ExecError {
  Instance&  ins;
  index_type act;
 public:
  RedundantActionWarning(Instance& i, index_type a, NTYPE at) :
    ExecError(ExecError::warning_redundant_action, at), ins(i), act(a) { };
  RedundantActionWarning(Instance& i, index_type a, NTYPE at, index_type s) :
    ExecError(ExecError::warning_redundant_action, at, s), ins(i), act(a) { };
  virtual ~RedundantActionWarning() { };

  const Instance::Action& action()       { return ins.actions[act]; };

  virtual ExecError* copy();
  virtual void write(::std::ostream& s) const;
};

class ExecState : public ProgressionState {
 protected:
  Instance&  instance;

  // atoms true in the state
  // - atoms deleted by currently executing actions ARE set to false
  // - atoms added by currently executing actions ARE NOT set to true
  bool_vec   atoms;

  // remaining resources
  amt_vec    res;

 public:
  struct exec_act {
    index_type act;
    NTYPE      rem;
    index_type step;
    void*      esi; // ext. state info

    exec_act() :
      act(no_such_index), rem(ZERO), step(no_such_index), esi(0) { };
    exec_act(index_type a, NTYPE r) :
      act(a), rem(r), step(no_such_index), esi(0) { };
    exec_act(index_type a, NTYPE r, index_type s) :
      act(a), rem(r), step(s), esi(0) { };
    exec_act(const exec_act& e) :
      act(e.act), rem(e.rem), step(e.step), esi(e.esi) { };
    exec_act& operator=(const exec_act& e) {
      act = e.act;
      rem = e.rem;
      step = e.step;
      esi = e.esi;
    };
  };
  typedef lvector<exec_act> exec_act_vec;

 protected:
  // actions executing in the state
  exec_act_vec actions;

  // absolute time and dt w.r.t. previous state
  NTYPE      abs_t;
  NTYPE      delta_t;

  // duration of this state
  NTYPE      dur;

  // extended state hooks
  virtual bool check_extended_applicability
    (Instance::Action& a, ExecErrorSet* errors, index_type step);
  virtual void apply_extended_action_start_effects
    (Instance::Action& a, exec_act& e, bool_vec& s);
  virtual void apply_extended_action_end_effects
    (Instance::Action& a, exec_act& e, bool_vec& s);
  virtual void save_extended_state_info(exec_act& e);
  virtual void free_extended_state_info(exec_act& e);

  class EST : public ProgressionPath {
    index_set starting;
    NTYPE     delta;
    EST(const EST& from) {
      starting = from.starting;
      delta = from.delta;
    };
  public:
    EST(ExecState& s) {
      s.starting_actions(starting);
      if (s.is_final())
	delta = s.max_delta();
      else
	delta = s.min_delta();
    };
    int compare(const Transition& t) {
      const EST& est = (const EST&)t;
      if (starting < est.starting)
	return -1;
      else if (starting > est.starting)
	return 1;
      else if (delta < est.delta)
	return -1;
      else if (delta > est.delta)
	return 1;
      else
	return 0;
    };
    Transition* copy() {
      return new EST(*this);
    };
    void insert(Plan& p) {
      p.advance(delta);
      for (index_type k = 0; k < starting.size(); k++)
	p.insert(starting[k]);
    };
  };

 public:
  ExecState(Instance& i);
  ExecState(Instance& i, index_set g);
  ExecState(Instance& i, const bool_vec& g);
  ExecState(const ExecState& s);
  virtual ~ExecState();

  int trace_level;

  virtual Transition* transition();
  virtual NTYPE delta_cost();
  virtual void insert(Plan& p);
  virtual void write_plan(::std::ostream& s);
  virtual NTYPE est_cost();
  virtual bool is_final();
  virtual bool is_max();
  virtual NTYPE expand(Search& s, NTYPE bound);
  virtual void store(NTYPE cost, bool opt);
  virtual void reevaluate();
  virtual int compare(const State& s);
  virtual index_type hash();
  virtual State* new_state(index_set& s);
  virtual State* copy();
  virtual void write(::std::ostream& s);

  // query information about the state
  NTYPE current_time() const;
  NTYPE end_time() const;
  void current_atoms(index_set& atms) const;
  const bool_vec& current_atoms() const;
  void current_actions(index_set& acts) const;
  index_type n_current_actions() const;
  const exec_act_vec& current_actions() const;
  bool is_stable() const;
  void starting_actions(index_set& acts) const;
  void finishing_actions(index_set& acts) const;
  NTYPE min_delta() const;
  NTYPE max_delta() const;
  void current_resource_levels
    (amt_vec& avail, amt_vec& in_use) const;
  void current_resource_levels
    (index_type r, NTYPE& avail, NTYPE& in_use) const;

  // note: this method returns values for reusable resources only!
  void current_resource_use(amt_vec& res) const;

  bool test_conjunction(const index_set& set) const;
  bool test_disjunction(const index_set& set) const;
  bool check_atoms(const index_set& set) const;
  bool check_atoms(const index_set& set, index_set& holds) const;

  // construction methods: these should only be used by Schedule/ExecTrace
  // while simulating plans/constructing traces.
  bool is_final(ExecErrorSet* errors);
  bool applicable(index_type act,
		  ExecErrorSet* errors,
		  index_type step);
  void apply(index_type act,
	     ExecErrorSet* errors,
	     index_type step);
  void advance(NTYPE dt,
	       ExecErrorSet* errors);
  void finish(ExecErrorSet* errors);
  void clip(NTYPE at_t);
  void intersect(const bool_vec& atms);
};

class Timeline {
 public:
  virtual ~Timeline();

  virtual index_type n_intervals() = 0;
  virtual index_type n_points() = 0;

  // may return no_such_index
  virtual index_type interval_start_point(index_type i) = 0;
  virtual index_type interval_end_point(index_type i) = 0;

  // may cause error
  virtual NTYPE point_time(index_type i) = 0;
  virtual NTYPE interval_start_time(index_type i) = 0;
  virtual NTYPE interval_end_time(index_type i) = 0;

  NTYPE total_time();
};

typedef comparable_pair<index_set> index_set_pair;

enum plan_constraint_type {
  // the basic PDDL3 plan constraint types
  pc_at_end,
  pc_always,
  pc_sometime,
  pc_at_most_once,
  pc_sometime_before,
  pc_sometime_after,
  pc_within,
  pc_always_within,
  pc_hold_after,
  pc_hold_during,
  // my specialisations
  pc_at_end_disjunction,
  pc_always_disjunction,
  pc_sometime_disjunction,
  pc_at_most_once_disjunction,
  pc_never,
  pc_negated_sometime,
  pc_never_disjunction,
  pc_negated_sometime_disjunction,
  pc_disallowed,
  pc_contradiction
};

class ExecTrace : public Timeline, public lvector<ExecState*> {
  Instance& instance;
 public:
  ExecTrace(Instance& ins);
  ~ExecTrace();

  const Instance& get_instance() const { return instance; };

  // timeline implementations
  virtual index_type n_intervals();
  virtual index_type n_points();
  virtual index_type interval_start_point(index_type i);
  virtual index_type interval_end_point(index_type i);
  virtual NTYPE point_time(index_type i);
  virtual NTYPE interval_start_time(index_type i);
  virtual NTYPE interval_end_time(index_type i);

  ExecState* final_state();
  ExecTrace* copy();
  ExecTrace* necessary_trace();
  ExecTrace* stable_trace();

  // note: this method returns values for reusable resources only!
  void peak_resource_use(amt_vec& res);

  // test/extract trajectory constraints
  bool test_always(index_type p);
  bool test_always_conjunction(const index_set& c);
  bool test_always_disjunction(const index_set& d);
  bool test_sometime(index_type p);
  bool test_sometime_conjunction(const index_set& c);
  bool test_sometime_disjunction(const index_set& d);
  bool test_at_most_once(index_type p);
  bool test_at_most_once_conjunction(const index_set& c);
  bool test_at_most_once_disjunction(const index_set& c);
  bool test_sometime_after(index_type p, index_type q);
  bool test_sometime_before(index_type p, index_type q);
  void find_always(Heuristic* inc,
		   index_type lim,
		   const index_set_vec& c_in,
		   index_set_vec& c_out);
  void find_always(Heuristic* inc,
		   index_type lim,
		   index_set_vec& c_out);
  void find_always(Heuristic* inc,
		   index_set_vec& c_out);
  void extract_always_within(bool_matrix& c, cost_matrix& t);

  // extract sometime, sometime-before and sometime-after constraints
  void find_order_constraints(graph* lms,
			      index_set& e,
			      pair_set& sb,
			      pair_set& sa);

  // atoms true at any point (i.e., true sometime constraints)
  void true_at_any_time(bool_vec& atoms);

  // internal method (used by Schedule::simulate)
  void clip_last_state(NTYPE at_t);

  void write(::std::ostream& s);

 private:
  bool find_always(index_type k,
		   bool_vec& p,
		   Heuristic* inc,
		   index_type lim,
		   index_set_vec& c);
};

class BasicTimeline : public Timeline {
 protected:
  cost_vec points;
  bool     open_start;
  bool     open_end;

  void set_point(NTYPE t);
  void clip_start(NTYPE t);
  void clip_end(NTYPE t);

 public:
  BasicTimeline();
  virtual ~BasicTimeline();

  virtual index_type n_intervals();
  virtual index_type n_points();
  virtual index_type interval_start_point(index_type i);
  virtual index_type interval_end_point(index_type i);
  virtual NTYPE point_time(index_type i);
  virtual NTYPE interval_start_time(index_type i);
  virtual NTYPE interval_end_time(index_type i);

  void write(::std::ostream& s);
};

class ResourceProfile : public BasicTimeline {
 protected:
  Instance&     instance;
  index_type    res;
  amt_vec       avail;
  amt_vec       in_use;
  index_set_vec a_start;
  index_set_vec a_finish;
  NTYPE         max_req;
 public:
  ResourceProfile(Instance& ins, index_type r, ExecTrace& trace);
  virtual ~ResourceProfile();

  void set_makespan(NTYPE t);

  const Name* resource_name() const { return instance.resources[res].name; };
  NTYPE amount_available(index_type i);
  NTYPE amount_in_use(index_type i);
  NTYPE amount_free(index_type i);

  index_type first_use_interval(index_type i);
  index_type first_min_free_interval(index_type i);
  NTYPE possible_unexpected_loss_to(index_type i);
  NTYPE min_free_from(index_type i);

  NTYPE min_free();
  NTYPE peak_use();
  NTYPE min_peak_use();
  NTYPE tolerable_unexpected_loss();
  NTYPE total_consumption();

  void writeGantt(::std::ostream& s);
  static double GANTT_EXTRA_WIDTH;
  static double GANTT_EXTRA_HEIGHT;

  void write(::std::ostream& s);
};

class Schedule : public Plan {
 public:

  struct step {
    index_type act;
    NTYPE      at;
    index_type track;

    step() : act(no_such_index), at(0), track(no_such_index) { };
    step(index_type a, NTYPE t) : act(a), at(t), track(no_such_index) { };

    step& operator=(const step& s) {
      act = s.act;
      at = s.at;
      return *this;
    };

    bool operator==(const step& s) {
      return ((act == s.act) && (at == s.at));
    };
  };

  typedef lvector<step> step_vec;

 private:
  Instance& instance;
  step_vec  steps;
  NTYPE     end_t;

  index_set  action_set;
  index_vec  action_vec;
  index_type n_tracks;

  // plan annotations
  const Name* ann_name;
  bool        ann_optimal;

  // used during construction
  NTYPE     current_t;
  bool      finished;

  // internal util methods
  NTYPE next_start_time(NTYPE t) const;
  NTYPE next_finish_time(NTYPE t) const;
  NTYPE last_finish_time() const;
  void compute_action_set_and_vec();
  void assign_tracks();

  enum s_status { s_pending, s_ready, s_executing, s_finished };

  // internal recursive method
  bool construct_minimal_makespan(const index_vec& acts,
				  graph& prec,
				  const index_set& c,
				  NTYPE& best,
				  index_vec& sindex);

 public:
  static bool write_traits;
  int       trace_level;

  Schedule(Instance& i);
  Schedule(const Schedule& s);
  virtual ~Schedule();

  index_type length() const;
  NTYPE      makespan() const;
  NTYPE      cost() const;

  index_type n_steps() const { return steps.length(); };
  const step_vec&  plan_steps() const { return steps; };
  const index_set& plan_actions() const { return action_set; };
  const index_vec& step_actions() const { return action_vec; };
  void step_action_names(name_vec& nv);
  bool step_in_interval(index_type s, NTYPE i_start, NTYPE i_end) const;
  index_type step_action(index_type s) const;

  const Name* plan_name() const;
  bool        plan_is_optimal() const;

  // schedule traits (public only for convenience, treat as const).
  plan_trait_vec traits;

  // find first trait of given class
  const PlanTrait* find_trait(const char* cn);

  // plan input methods
  virtual void protect(index_type atom);
  virtual void insert(index_type act);
  virtual void advance(NTYPE delta);
  virtual void end();

  // explicit plan construction method: use with caution!
  void insert_step(NTYPE at, index_type act);

  // input schedule is assumed to be relative: start time is
  // set to zero by default, but can be changed after the
  // schedule is completed.
  void set_start_time(NTYPE at);

  // remove steps indicated as redundant by the set of warnings
  void reduce(ExecErrorSet* warnings);

  // clear the schedule (allows new schedule to be input/constructed)
  void clear();

  // set plan annotations
  virtual void set_name(const Name* n);
  virtual void set_optimal(bool o);
  virtual void add_trait(PlanTrait* t);

  // output and output with remapping of actions
  virtual void output(Plan& plan, bool mark_end = true) const;
  virtual void output(Plan& plan, const index_vec& act_map,
		      bool mark_end = true) const;

  // write the plan (in PDDL format)
  void write(::std::ostream& s, unsigned int c = Name::NC_DEFAULT) const;

  // write (a subset of) schedule steps
  void write_step_set(::std::ostream& s, const index_set& set) const;
  void write_steps(::std::ostream& s) const;

  // write the plan in XML format (errors, trace and precedence graph are
  // optional - if provided, the information they contain is also written)
  void writeXML(::std::ostream& s,
		ExecErrorSet* errors = 0,
		ExecTrace* trace = 0,
		graph* prec = 0,
		index_type id = 0) const;

  // write a Gantt chart of the schedule in LaTeX
  void writeGantt(::std::ostream& s) const;
  // parameters for (LaTeX) Gantt chart writing
  static double GANTT_UNIT_WIDTH;
  static double GANTT_UNIT_HEIGHT;
  static double GANTT_TEXT_XOFF;
  static double GANTT_TEXT_YOFF;
  static NTYPE  GANTT_TIME_MARK_INTERVAL;
  static bool   GANTT_ACTION_NAMES_ON_CHART;

  // basic simulation methods, high and low resolution
  bool simulate(ExecState* cs,
		ExecTrace* trace = 0,
		ExecErrorSet* errors = 0,
		bool finish = false) const;
  bool simulate(ExecTrace* trace = 0,
		ExecErrorSet* errors = 0,
		bool finish = false) const;
  bool simulate_low_resolution(ExecTrace* trace = 0,
			       ExecErrorSet* errors = 0,
			       bool finish = false) const;

  // short-cut methods: simulate and extract info from trace
  bool simulate(index_set& achieved,
		ExecErrorSet* errors = 0) const;
  bool simulate(amt_vec& rtl) const;

  // deorder the schedule (PRF algorithm)
  void deorder(graph& prec) const;
  void deorder_non_temporal(graph& prec) const;
  void deorder(weighted_graph& prec) const;

  // extract causal links (Chrpa's algorithm)
  void causal_links(index_set_graph& cl) const;

  // simple method of constructing a partially ordered plan: first
  // extract causal links, then order unordered interfering actions
  // as they are ordered in the input plan.
  void make_pop(index_set_graph& clt) const;

  // get the "raw" precedence graph (without deordering)
  void base_precedence_graph(graph& prec) const;

  index_type action_multiset_difference(const Schedule& s) const;

  // test if the given schedule is "essentially equivalent" to this;
  // NOTE: the test is done by deordering both plans, then checking if
  // there is a (action-preserving) mapping between the steps under
  // which their precedence graphs are equal; this test does not take
  // resources into account in any way;
  // the "c" vector, if non-nil, on return contains a correspondance
  // (1-1-mapping) between the steps of the two schedules if the return
  // value is true
  bool equivalent(const Schedule& s, index_vec* c = 0) const;

  // construct the schedule by scheduling the given set of actions subject
  // to given precedence constraints - returns false if scheduling fails, but
  // note that this does not necessarily mean the plan is unschedulable, if
  // use+consume resources are present (in this case, it may be necessary to
  // search for a working schedule). if map is non-nil, it is set (on success)
  // to a mapping from the input steps to the step indices in the final
  // schedule, i.e. map[i] = j <-> the i:th element of acts/prec appears as
  // step with index j in the schedule
  bool schedule(const index_vec& acts,
		const graph& prec,
		index_vec* map = 0);

  // construct a schedule from actions (acts) and precedence graph (prec),
  // returning true iff such a schedule can be constructed without any
  // resource conflict (if this is the case, the schedule is optimal w.r.t.
  // acts/prec); if construction fails, "cs" contains a set of conflicting
  // steps (indices into acts/prec) on return; map maps input steps to
  // schedule steps (as in method schedule)
  bool construct_conflict_free(const index_vec& acts,
			       const graph& prec,
			       index_set& cs,
			       index_vec& map);

  // find a scheduling of acts/prec with minimal makespan (through repeated
  // calls to the above method)
  bool construct_minimal_makespan(const index_vec& acts,
				  const graph& prec,
				  index_vec& map);


  // create a random sequential plan
  bool random_sequence(index_type ln_max,
		       index_type ln_avg,
		       bool continue_from_goal,
		       ExecTrace* trace,
		       RNG& rnd);
};

void compute_plan_dependency_graph
(Instance& instance,
 const index_vec& acts,
 const graph& prec,
 const index_set& goal,
 index_graph& dg,
 name_vec* oc_names);

void compute_plan_exec_graph
(Instance& instance,
 Schedule& schedule,
 index_graph& g,
 name_vec& node_name);

void compute_plan_transition_graph
(Instance& instance,
 Schedule& schedule,
 index_graph& g,
 name_vec& node_name);

// note: TraceSet owns the traces, but not the schedules; schedule
// pointer is only used to keep track of the plan-trace correspondance,
// and can be nil.
class TraceSet : public lvector< std::pair<ExecTrace*, Schedule*> > {
  Instance& instance;
 public:
  TraceSet(Instance& ins);
  ~TraceSet();

  void add_trace(ExecTrace* t, Schedule* s);
  void add_plan(Schedule* s, bool untimed = false);
  void add_plan_if_executable(Schedule* s, bool untimed = false);
  void add_plan_if_valid(Schedule* s, bool untimed = false);

  // test constraints against all traces, return vector of values
  void test_always_conjunction(const index_set& s, bool_vec& v);
  void test_always_disjunction(const index_set& s, bool_vec& v);
  void test_sometime_conjunction(const index_set& s, bool_vec& v);
  void test_sometime_disjunction(const index_set& s, bool_vec& v);
  void test_at_most_once_conjunction(const index_set& s, bool_vec& v);
  void test_at_most_once_disjunction(const index_set& s, bool_vec& v);

  // find constraints of type (always (or atom1 ... atomN)) that hold
  // in every trace in the set
  void find_common_always(Heuristic* inc, index_type lim, index_set_vec& c);

  void find_common_order_constraints(graph* lms,
				     index_set& e,
				     pair_set& sb,
				     pair_set& sa);

  // find a set of at-most-once constraints, of size at most lim,
  // separating trace p from other traces. both conjunctions and
  // disjunctions are considered: vector is_d indicates which each
  // output constraint (set of atoms) is.
  void find_separating_always(index_type p,
			      index_type lim,
			      index_set_vec& c,
			      bool_vec& is_neg,
			      bool_matrix& ex);
  void find_separating_sometime(index_type p,
				index_type lim,
				index_set_vec& c,
				bool_vec& is_neg,
				bool_matrix& ex);
  void find_separating_at_most_once(index_type p,
				    index_type lim,
				    index_set_vec& c,
				    bool_vec& is_d,
				    bool_matrix& ex);

  struct PlanConstraint {
    plan_constraint_type pc_type;
    index_set            c_atoms; // constraint
    index_set            t_atoms; // trigger
    bool_vec             ex;
    PlanConstraint() : pc_type(pc_at_end) { };
    PlanConstraint(const PlanConstraint& pc) : pc_type(pc.pc_type),
	 c_atoms(pc.c_atoms), t_atoms(pc.t_atoms), ex(pc.ex) { };
    PlanConstraint(plan_constraint_type t, const index_set& s, index_type n) :
      pc_type(t), c_atoms(s), t_atoms(EMPTYSET), ex(false, n) { };
    PlanConstraint(plan_constraint_type t, const index_set& s, const bool_vec& x) :
      pc_type(t), c_atoms(s), t_atoms(EMPTYSET), ex(x) { };

    index_type ltype() const; // position in type lattice
    bool below(index_type lt) const; // strict-less-than in type lattice
    bool dominates(const PlanConstraint& pc) const;
    bool is_dominated_by(const PlanConstraint& pc) const;
  };

  typedef lvector<PlanConstraint> plan_constraint_vec;

  void find_separating_constraints(index_type p,
				   index_type s_min,
				   index_type s_max,
				   plan_constraint_vec& pc);

  void extract_sc_candidates(weighted_vec<index_type,NTYPE>& sc_always,
			     weighted_vec<index_type,NTYPE>& sc_sometime,
			     weighted_vec<index_type,NTYPE>& sc_once,
			     weighted_vec<index_pair,NTYPE>& sc_before);
  void calculate_sc_weights(weighted_vec<index_type,NTYPE>& sc_always,
			    weighted_vec<index_type,NTYPE>& sc_sometime,
			    weighted_vec<index_type,NTYPE>& sc_once,
			    weighted_vec<index_pair,NTYPE>& sc_before);

 private:
  bool add_constraint_if_not_dominated(const index_set& nc,
				       bool nisx,
				       const bool_vec& nex,
				       index_set_vec& c,
				       bool_vec& isx,
				       bool_matrix& ex,
				       index_type& max_ex) const;
  bool add_non_dominated(const PlanConstraint& pc,
			 plan_constraint_vec& c) const;
};

typedef lvector<Schedule*> plan_vec;

class ScheduleSet : public plan_vec, public PlanSet {
 protected:
  struct ScheduleProperties {
    bool    valid;
    ExecTrace* trace;
    index_type n_resources;
    amt_vec total_consumption;
    amt_vec peak_use;
    amt_vec tolerable_loss;
    NTYPE   makespan;

    ScheduleProperties()
      : valid(true), trace(0), n_resources(0), total_consumption(0, 0),
	peak_use(0, 0), tolerable_loss(POS_INF, 0), makespan(0) { };
    ScheduleProperties(index_type n)
      : valid(true), trace(0), n_resources(n), total_consumption(0, n),
	peak_use(0, n), tolerable_loss(POS_INF, n), makespan(0) { };

    bool dominates(const ScheduleProperties& p);
  };

  typedef lvector<ScheduleProperties*> prop_vec;

  Instance& instance;
  prop_vec  props;

  ScheduleProperties* compute_properties(Schedule* s);
  void cache_properties(index_type i, ScheduleProperties* p);

  bool dominated(const ScheduleProperties& p, prop_vec& pv);
  void dominated(prop_vec& pv, bool_vec& dom);
  bool dominated(const ScheduleProperties& p);

  void replace_schedule_with_properties
    (index_type i, Schedule* s,	ScheduleProperties* p);
  void add_schedule_with_properties
    (Schedule* s, ScheduleProperties* p);

  void find_tweaks(Schedule* src_plan,
		   ScheduleProperties* src_plan_props,
		   const index_vec& acts,
		   const graph& prec,
		   index_type res,
		   ResourceProfile* rp,
		   pair_vec& tweaks,
		   plan_vec& new_scheds,
		   prop_vec& new_props);

  void explore_options(Schedule* s,
		       ScheduleProperties* p,
		       index_type res,
		       plan_vec& new_scheds,
		       prop_vec& new_props);

  void sequential_variations(const Name* src_name,
			     index_type& n,
			     const index_vec acts,
			     const graph& p,
			     const graph& p0,
			     const pair_set& v,
			     index_type d,
			     index_type d_min,
			     index_type d_max);

 public:
  ScheduleSet(Instance& i);
  ScheduleSet(ScheduleSet& s, const bool_vec& sel);
  ~ScheduleSet();

  int trace_level;

  // note: schedule set (this) assumes ownership of the schedule, EVEN in
  // case the schedule is not different (it is then immediately deleted)
  void add_schedule(Schedule* s);
  void add_schedule_if_different(Schedule* s);

  TraceSet* trace_set(bool valid_only, bool exec_only, bool untimed = false);

  void reduce_plans();
  void filter_invalid_plans();
  void filter_unschedulable_plans();
  void filter_equivalent_plans();
  void add_distinguishing_traits_1();
  void add_distinguishing_traits_2();

  bool common_precedence_constraints(graph& prec);
  bool separating_precedence_constraints
    (ScheduleSet& s, pair_set& d0, pair_set& d1);

  void explore_options(bool cross_dominance_check, bool src_dominance_check);
  void sequential_variations(Schedule* s, index_type d_min, index_type d_max);
  index_type random_sequential_variations(Schedule* s,
					  index_type d_min,
					  index_type d_max,
					  index_type n,
					  RNG& rng);

  void write_deordered_graphs(::std::ostream& s, bool w_names = true);
  void writeXML(::std::ostream& s);

  virtual Plan* new_plan();
  virtual void output(PlanSet& to);
  virtual void output(PlanSet& to, const bool_vec& s);

  void remove(bool_vec& set);
  void clear();
  // clear the schedule set without deleteing the objects (plans) in it
  void clear_no_delete();
};

class ScheduleTrait : public PlanTrait {
 protected:
  Schedule*    plan;
  ScheduleSet* plan_set;
  bool is_min;
  bool is_max;
  bool is_unique;

  void write_meta_short(::std::ostream& s) const;
  void write_meta_attributes(::std::ostream& s) const;
 public:
  ScheduleTrait(Schedule* p)
    : plan(p), plan_set(0), is_min(false), is_max(false), is_unique(false) { };
  ScheduleTrait(Schedule* p, ScheduleSet* s)
    : plan(p), plan_set(s), is_min(false), is_max(false), is_unique(false) { };
  virtual ~ScheduleTrait();

  void set_min()    { is_min = true; };
  void set_max()    { is_max = true; };
  void set_unique() { is_unique = true; };

  virtual const PlanTrait* cast_to(const char* class_name) const;
  virtual void write_short(::std::ostream& s) const = 0;
  virtual void write_detail(::std::ostream& s) const;
  virtual void writeXML(::std::ostream& s) const = 0;
};

class EquivalentTo : public ScheduleTrait {
  Schedule* s_eq;
  index_vec cor;
 public:
  EquivalentTo(Schedule* p, ScheduleSet* ss, Schedule* s, const index_vec& c)
    : ScheduleTrait(p, ss), s_eq(s), cor(c) { };
  virtual ~EquivalentTo() { };
  virtual const PlanTrait* cast_to(const char* class_name) const;
  virtual void write_short(::std::ostream& s) const;
  virtual void write_detail(::std::ostream& s) const;
  virtual void writeXML(::std::ostream& s) const;
};

class DerivedFrom : public ScheduleTrait {
  Schedule* s_src;
  pair_set  e_prec;
 public:
  DerivedFrom(Schedule* p, ScheduleSet* ss, Schedule* s, const pair_set& e)
    : ScheduleTrait(p, ss), s_src(s), e_prec(e) { };
  DerivedFrom(Schedule* p, ScheduleSet* ss, Schedule* s, const index_pair& e)
    : ScheduleTrait(p, ss), s_src(s) { e_prec.assign_singleton(e); };
  DerivedFrom(Schedule* p, ScheduleSet* ss,
	      DerivedFrom* a, const index_vec& m,
	      const index_pair& e);

  virtual ~DerivedFrom() { };
  virtual const PlanTrait* cast_to(const char* class_name) const;
  virtual void write_short(::std::ostream& s) const;
  virtual void write_detail(::std::ostream& s) const;
  virtual void writeXML(::std::ostream& s) const;
};

class PlanPrecedenceRelation : public ScheduleTrait {
  graph prec;
 public:
  PlanPrecedenceRelation(Schedule* s, const graph& p)
    : ScheduleTrait(s), prec(p) { };
  PlanPrecedenceRelation(Schedule* s, const graph& p, const index_vec& m)
    : ScheduleTrait(s) { prec.copy_and_rename(p, m); };
  virtual ~PlanPrecedenceRelation() { };
  virtual const PlanTrait* cast_to(const char* class_name) const;
  virtual void write_short(::std::ostream& s) const;
  virtual void writeXML(::std::ostream& s) const;

  const graph& precedence_relation() const { return prec; };
};

class PlanActionOccurs : public ScheduleTrait {
  Instance& instance;
  index_type act;
  index_type n_of_times;
 public:
  PlanActionOccurs(Schedule* p, Instance& i, index_type a, index_type n)
    : ScheduleTrait(p), instance(i), act(a), n_of_times(n) { };
  virtual ~PlanActionOccurs() { };
  virtual const PlanTrait* cast_to(const char* class_name) const;
  virtual void write_short(::std::ostream& s) const;
  virtual void writeXML(::std::ostream& s) const;
};

class PlanStepOrder : public ScheduleTrait {
  Instance& instance;
  index_pair order;
 public:
  PlanStepOrder(Schedule* p, Instance& i, const index_pair& o)
    : ScheduleTrait(p), instance(i), order(o) { };
  virtual ~PlanStepOrder() { };
  virtual const PlanTrait* cast_to(const char* class_name) const;
  virtual void write_short(::std::ostream& s) const;
  virtual void writeXML(::std::ostream& s) const;

  const index_pair& precedence() const { return order; };
};

class PlanFeatureValue : public ScheduleTrait {
 public:
  enum plan_feature_type { makespan,
			   cost,
			   resource_peak_use,
			   resource_total_consumption,
			   resource_tolerable_loss };
 private:
  plan_feature_type ftype;
  Instance&      instance;
  index_type        index;
  NTYPE             value;
 public:
  PlanFeatureValue(Schedule* p,
		   plan_feature_type t,
		   Instance& i,
		   index_type x,
		   NTYPE v) :
    ScheduleTrait(p), ftype(t), instance(i), index(x), value(v) { };
  virtual ~PlanFeatureValue() { };
  virtual const PlanTrait* cast_to(const char* class_name) const;
  virtual void write_short(::std::ostream& s) const;
  virtual void writeXML(::std::ostream& s) const;
};

// util method: compute all extensions of a precedence order on a set of
// actions that ensure feasibility w.r.t. resource consumption and use.
// returns true iff at least one feasible order exists; if rfps != NIL, all
// feasible solutions are generated stored in the vector.
bool feasible(Instance& ins,
	      const index_vec& acts,
	      const graph& prec,
	      graph_vec* rfps = 0,
	      index_type* rff = 0);

// sub-methods internal to feasible(...)
bool feasible(Instance& ins,
	      const index_vec& acts,
	      index_type r,
	      graph& uc,
	      graph_vec* rfps = 0,
	      index_type* rff = 0);
bool feasible(Instance& ins,
	      const index_vec& acts,
	      index_type r,
	      const index_set& a,
	      index_type i,
	      graph& uc,
	      graph_vec* rfps = 0,
	      index_type* rff = 0);
bool feasible(Instance& ins,
	      const index_vec& acts,
	      index_type r,
	      const index_set& a,
	      index_type i_a,
	      const index_set& b,
	      index_type i_b,
	      NTYPE c_max,
	      graph& uc,
	      graph_vec* rfps = 0,
	      index_type* rff = 0);

class PlanName : public Name {
  const Name* src;
  const char* desc;
  index_type  index;
 public:
  PlanName(const char* s, index_type i)
    : src(0), desc(s), index(i) { };
  PlanName(const Name* n, const char* s, index_type i)
    : src(n), desc(s), index(i) { };
  virtual ~PlanName() { };

  virtual void write(::std::ostream& s, unsigned int c = NC_DEFAULT) const;
};

::std::ostream& operator<<(::std::ostream& s, ExecError::ErrorType t);
::std::ostream& operator<<(::std::ostream& s, ExecError::ErrorSeverity sv);

END_HSPS_NAMESPACE

#endif
