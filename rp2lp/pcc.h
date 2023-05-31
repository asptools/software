#ifndef PLAN_CONSTRAINT_CALCULUS_H
#define PLAN_CONSTRAINT_CALCULUS_H

#include "config.h"
//#include "problem.h"
#include "exec.h"
#include "base.h"
#include "mutex.h"
#include "stats.h"
#include "preprocess.h"

BEGIN_HSPS_NAMESPACE

/// ppc representation and some basic utils

// A ppc struct represents a Preference over a Plan Constraint.

struct ppc {
  const Name* name;
  NTYPE     weight;
  plan_constraint_type pct;
  index_set s_c; // condition set
  bool c_dis;    // condition set is disjunctive
  index_set s_t; // trigger set
  bool t_dis;    // trigger set is disjunctive
  PDDL_Base::Preference* src;

  ppc()
    : name(0), weight(0), pct(pc_at_end), c_dis(false), t_dis(false),
      src(0) { };
  ppc(const Name* n, NTYPE w, plan_constraint_type t,
      const index_set& c, bool d, PDDL_Base::Preference* s)
    : name(n), weight(w), pct(t), s_c(c), c_dis(d), s_t(EMPTYSET), t_dis(false),
      src(s) { };
  ppc(const Name* n, NTYPE w, plan_constraint_type t,
      const index_set& a_c, bool d_c, const index_set& a_t, bool d_t,
      PDDL_Base::Preference* s)
    : name(n), weight(w), pct(t), s_c(a_c), c_dis(d_c), s_t(a_t), t_dis(d_t),
      src(s) { };

  ppc& operator=(const ppc& p) {
    name = p.name;
    weight = p.weight;
    pct = p.pct;
    s_c.assign_copy(p.s_c);
    c_dis = p.c_dis;
    s_t.assign_copy(p.s_t);
    t_dis = p.t_dis;
    src = p.src;
  };

  // write preference in PDDL syntax
  void write(std::ostream& s, const Instance& ins) const;
  // write constraint only (ignoring name and weight) in PDDL syntax
  void write_PDDL_constraint(std::ostream& s, const Instance& ins) const;

  // print constraint in "logic-like" syntax
  void print(std::ostream& s, const Instance& ins) const;

  // enforce, compile and test the constraint
  void enforce(Instance& ins, index_vec& map) const;
  index_type compile(Instance& ins, index_vec& map) const;
  bool test(ExecTrace* t) const;
};

typedef lvector<ppc> ppc_vec;

bool extract_ppcs(PDDL_Base* b, Instance& ins, ppc_vec& ppcv);


/// new PCC implementation

class InferenceTracker;

class Propagator {
  friend class InferenceTracker;

 protected:
  const Instance& ins;
  Stopwatch& stats;

  // State formula representation (must be implemented by concrete subclass)

  // State formulas are indexed 0 .. n_state_formulas - 1; this field must
  // be set by subclass constructor.
  index_type n_state_formulas;

  // Logical relations between state formulas, atoms and actions.
  // These must be initialised by subclass constructor.

  //  g \in sf_implies[f] iff f implies g:
  index_set_vec sf_implies;
  //  f \in sf_implied_by[g] iff f implies g:
  index_set_vec sf_implied_by;
  //  p \in sf_implies_atom[f] iff f implies p (note that p is an atom):
  index_set_vec sf_implies_atom;
  //  sf_initial[f] == true iff f holds in init state:
  bool_vec sf_initial;

  //  a \in act_implies[f] iff f is true before or after applying a:
  index_set_vec act_implies;
  //  a \in act_implies_not[f] iff f is false before or after applying a:
  index_set_vec act_implies_not;

  //  a \in act_change[2*f] iff a changes f to true;
  //  a \in act_change[(2*f)+1] iff a changes f to false;
  index_set_vec act_change;

  Propagator(const Instance& i, Stopwatch& s,
	     const graph* lm, const graph* na);
  virtual ~Propagator();

  index_type add_constraint(plan_constraint_type t, index_type f);
  index_type add_constraint(plan_constraint_type t,
			    index_type ft, index_type fc);

  // Check mutex(f, g), where f and g are state formula indices.
  virtual bool mutex(index_type f, index_type g) const = 0;

  virtual void print_state_formula(index_type f, std::ostream& s) const = 0;

  // Check for triggering of new conditional constraints, given (updated)
  // set of allowed actions. Triggered constraints must be asserted by
  // calling the assert_<constraint-type> methods below. The method must
  // return true iff any of the assert calls do (indicating a contradiction),
  // and false otherwise.
  virtual bool check_conditional_constraints
    (const bool_vec& allowed_actions,
     const bool_vec& current_never,
     const adjacency_list_graph& current_prec,
     InferenceTracker* it);

  // Asserting constraints:
  // note: these return true iff a contradiction has been detected; only
  // sometime/never contradictions are checked.
  bool assert_never(index_type f);
  bool assert_sometime(index_type f);
  // note: assert_sometime_before/after follow the PDDL3 convention:
  // ft is the triggering formula and fc the consequence.
  bool assert_sometime_before(index_type ft, index_type fc);
  bool assert_sometime_after(index_type ft, index_type fc);

  // A hook for subclasses to add steps at reset.
  virtual void reset_hook();

  // Optional problem knowledge: landmarks (sb's), never-after's and
  // triggered (conditional) landmarks. These must all be expressed in
  // terms of state formula indices, but may cover only a subset of
  // sf indices, as long as that subset is [0 .. m'-1] for some m' <= m
  // (m being the number of state formulas)
  const graph* landmarks;
  const graph* neverafter;

 private:
  // Internal representation of constraints. We only need to keep
  // track of the constraint type and the index of the (one or two)
  // state formulas involved: if only one, then pcf.first is used,
  // if two, the pcf.first is the trigger and pcf.second the
  // consequence.
  std::vector<plan_constraint_type> pct;
  pair_vec                          pcf;

  // State of inference:
  // never[f] == true iff never(f).
  // sometime[f] == true iff sometime(f).
  // allowed[a] == true iff action a is still allowed.
  // sometime_before: (adjacency-list-)graph with edges (f,g) such that
  //  sometime-before(f, g) (i.e., f is the trigger).
  // sometime_after: (adjacency-list-)graph with edges (f,g) such that
  //  sometime-after(f, g) (i.e., f is the trigger).
  // prec holds the sometime-before relation + implications:
  //  prec.adjacent(f, g) == true iff sometime-before(f, g) (trigger ->
  //  consequence) or if f implies g.
  // prec0 holds the prec graph with implications and landmark sb's only;
  //  it is initialised on the first call to reset, and copies into prec
  //  on subsequent reset calls.
  // wprec holds the weak precedence relation:
  //  wprec.adjacent(f, g) == true iff sometime-after(f, g) or f implies g.
  bool_vec never;
  bool_vec sometime;
  bool_vec allowed;
  bool allowed_actions_changed;
  adjacency_list_graph sometime_before;
  adjacency_list_graph sometime_after;
  adjacency_list_graph prec;
  adjacency_list_graph wprec;
  adjacency_list_graph prec0;
  bool sometime_before_changed;
  bool sometime_after_changed;

  InferenceTracker* it;

  // Steps of the inference algorithm:

  bool run1(const index_set& s);
  void reset();
  bool checkAMO(const index_set& s);

 public:
  void set_tracker(InferenceTracker* t);
  InferenceTracker* tracker() { return it; };

  // Check consistency of a given subset of the specified constraints.
  // Returns true iff a contradiction is proven.
  bool run(const index_set& s);
};

class InferenceTracker {
 public:

  struct assertion {
    plan_constraint_type pc;
    index_type a;
    index_type b;

    assertion()
      : pc(pc_at_end), a(no_such_index), b(no_such_index) { };
    assertion(plan_constraint_type c)
      : pc(c), a(no_such_index), b(no_such_index) { };
    assertion(plan_constraint_type c, index_type a1)
      : pc(c), a(a1), b(no_such_index) { };
    assertion(plan_constraint_type c, index_type a1, index_type b1)
      : pc(c), a(a1), b(b1) { };

    bool match(plan_constraint_type t, index_type f1, index_type f2) const;
  };

  struct step {
    assertion ass;
    const char* rule;
    index_vec pre;
    bool valid;

    step() : rule(NULL), valid(false) { };
    step(assertion a, const char* r) : ass(a), rule(r), valid(true) { };
  };

 private:

  struct compare_assertions {
    bool operator()(const assertion& a1, const assertion& a2) const
    {
      return ((a1.pc < a2.pc) ||
	      ((a1.pc == a2.pc) && (a1.a < a2.a)) ||
	      ((a1.pc == a2.pc) && (a1.a == a2.a) && (a1.b < a2.b)));
    }
  };

  Propagator& propagator;
  std::vector<step> steps;
  std::map<assertion, index_type, compare_assertions> index;
  bool verbose;
  bool debug_mode;

  void mark_steps(index_type i, bool_vec& m) const;
  void print_step(index_type i, std::ostream& to) const;

  index_type find_transitive(index_type alpha,
			     index_type beta,
			     const adjacency_list_graph& g,
			     plan_constraint_type pc,
			     const char* rule);

 public:
  InferenceTracker(Propagator& p);
  ~InferenceTracker();

  void set_verbose(bool yes) { verbose = yes; };
  void set_debug_mode(bool yes) { debug_mode = yes; };

  void print_assertion(const assertion& a, std::ostream& to) const;
  void print_proof(const assertion& a, std::ostream& to) const;
  void print_proof(std::ostream& to) const;
  void print_all(std::ostream& to) const;

  void extract_proof(const assertion& a, std::vector<step>& p) const;
  bool extract_proof_premises(const assertion& a, index_set& p) const;

  // check and print proofs of all constraints currently asserted
  // in the propagator's internal data structures.
  void validate() const;

  // premise: A alpha
  void premise_A(index_type alpha, const char* rule);
  // premise: E alpha
  void premise_E(index_type alpha, const char* rule);
  // premise: sometime-before(alpha, beta) (or: beta SB alpha)
  void premise_SB(index_type alpha, index_type beta, const char* rule);
  // premise: sometime-after(alpha, beta) (or: beta SA alpha)
  void premise_SA(index_type alpha, index_type beta, const char* rule);
  // premise: AMO alpha
  void premise_AMO(index_type alpha, const char* rule);

  // A alpha => E alpha (because there is always a state)
  void infer_E_from_A(index_type alpha, const char* rule);

  // A alpha & mutex(alpha, beta) => N beta
  void infer_N_from_A_mutex(index_type alpha, index_type beta,
			    const char* rule);

  // A alpha & (del(act) -> !alpha) => D act
  void infer_D_from_A(index_type act, index_type alpha, const char* rule);

  // N alpha & (pre(act)|add(act) -> alpha) => D act
  void infer_D_from_N(index_type act, index_type alpha, const char* rule);

  // AMO alpha & (act in ActChT(alpha)) => D act
  void infer_D_from_AMO(index_type act, index_type alpha, const char* rule);

  // // sometime-before(alpha, gamma) | (alpha -> gamma) &
  // // sometime-before(gamma, beta) | (gamma -> beta)
  // //  => sometime-before(alpha, beta)
  // void infer_transitive_SB(index_type alpha, index_type beta, const char* rule);
  // // sometime-after(alpha, gamma) | (alpha -> gamma) &
  // // sometime-after(gamma, beta) | (gamma -> beta)
  // //  => sometime-after(alpha, beta)
  // void infer_transitive_SA(index_type alpha, index_type beta, const char* rule);

  // N beta & beta SB alpha => N alpha
  void infer_N_from_N_SB(index_type alpha, index_type beta, const char* rule);
  // E beta & alpha SB beta => E alpha
  void infer_E_from_E_SB(index_type alpha, index_type beta, const char* rule);

  // N beta & beta SA alpha => N alpha
  void infer_N_from_N_SA(index_type alpha, index_type beta, const char* rule);
  // E beta & alpha SA beta => E alpha
  void infer_E_from_E_SA(index_type alpha, index_type beta, const char* rule);

  // N beta & (alpha -> beta) => N alpha
  void infer_N_from_N_implies(index_type alpha, index_type beta,
			      const char* rule);
  // E beta & (beta -> alpha) => E alpha
  void infer_E_from_E_implies(index_type alpha, index_type beta,
			      const char* rule);

  // alpha SB beta & beta SB alpha => N alpha
  void infer_N_from_SB_cycle(index_type alpha, index_type beta,
			     const char* rule);

  // beta SB alpha & alpha NA beta => N alpha
  void infer_N_from_SB_NA(index_type alpha, index_type beta,
			  const char* rule);

  // (alpha SA beta) & (beta SA alpha) & mutex(alpha, beta) => N alpha
  void infer_N_from_SA_cycle(index_type alpha, index_type beta,
			     const char* rule);

  // beta SA alpha & beta NA alpha => N alpha
  void infer_N_from_SA_NA(index_type alpha, index_type beta,
			  const char* rule);

  // N alpha & E alpha => CONTRADICTION
  void contradiction_from_N_E(index_type alpha, const char* rule);

  // E alpha & E beta & (alpha NA beta) & (beta NA alpha) => CONTRADICTION
  void contradiction_from_NA_cycle(index_type alpha, index_type beta,
				   const char* rule);

  // |R[I]| > |S| => CONTRADICTION
  // R is the set of atoms (implied E's)
  // I is an independent subset of R (indexes R)
  // S is a set AMO-action sets (indices into act_change) that covers
  // the still allowed adders of each atom in R;
  //  AMO f and f not initial => (2*f) may be in S
  //  AMO f => (2*f)+1 may be in S
  void contradiction_from_AMO(const index_set& R, const index_set& I,
			      const index_set& S, const char* rule);

  // { Da | a in acts } => N alpha
  void infer_N_from_D(index_type alpha, const index_set& acts,
		      const char* rule);

  // { Da | a in acts } => beta SB alpha
  void infer_SB_from_D(index_type alpha, index_type beta,
		       const index_set& acts, const char* rule);
};

class PCC : public Propagator {
  // Representation of state formulas: The set of state formulas appearing
  // in constraints is fixed (doesn't change during propagation), and
  // formulas are indexed 0 .. N-1. The set always includes all atoms in
  // the instance, which appear first in the order. They are followed by
  // a set of atom sets (stored in the vector "atomsets"), which represent
  // either conjunctions or disjunctions (indicated by the vector
  // "disjunctive")
  index_set_vec atomsets;
  bool_vec disjunctive;

  //virtual index_type n_state_formulas() const
  //{ return ins.n_atoms() + atomsets.size(); };

  bool is_atomic(index_type f) const {
    assert(f < n_state_formulas);
    return (f < ins.n_atoms());
  }

  bool is_disjunctive(index_type f) const {
    assert(f < n_state_formulas);
    if (f >= ins.n_atoms())
      return disjunctive[f - ins.n_atoms()];
    else
      return false;
  }

  // Map a con/dis atom set to a state formula index, adding it
  // to atomsets if necessary. (Should only be called during
  // initialisation.)
  index_type map_state_formula(const index_set& s, bool dis);

 protected:
  // Implementation of virtuals from the Propagator class:
  virtual bool mutex(index_type f, index_type g) const;
  virtual void print_state_formula(index_type f, std::ostream& s) const;
  virtual bool check_conditional_constraints
    (const bool_vec& allowed_actions,
     const bool_vec& current_never,
     const adjacency_list_graph& current_prec,
     InferenceTracker* it);
  virtual void reset_hook();

  // Internal helper methods:
 private:
  // Check implication f -> g, where f is a conjunction (set of atoms)
  // and g is a state formula index.
  bool conjunction_implies(const index_set& f, index_type g) const;

  // Check implication f -> g, where f is a disjunction (set of atoms)
  // and g is a state formula index.
  bool disjunction_implies(const index_set& f, index_type g) const;

  // Check if conjunction (set of atoms) f is mutex with g (if some
  // element in f is mutex with g):
  bool conjunction_is_mutex(const index_set& f, index_type g) const;

  // Check if disjunction (set of atoms) f is mutex with g (only if
  // every element in f is mutex with g):
  bool disjunction_is_mutex(const index_set& f, index_type g) const;

  // Check if del(a) makes f false.
  bool action_deletes_formula(index_type a, index_type f) const;

  // Initialise internal representation with a vector of ppcs.
  void init(const ppc_vec& ppcs);

  // helper methods used in init.
  bool check_initial(index_type f) const;
  bool check_implies(index_type f, index_type g) const;
  bool check_action_change_to_true(index_type a, index_type f) const;
  bool check_action_change_to_false(index_type a, index_type f) const;

  // Additional problem knowledge:
  StaticMutex* mx;
  set_edge_vec* trlm;
  bool_vec rem_trlm;

  // Expanded never-after graph: The graph is enriched with edges between
  // non-atomic formulas, derived by implications.
  graph* xnag;

  // Goals must be added as at-end constraints; but then when run gets
  // called, it's with a set that doesn't include the goals, so we have
  // to remember them and add them.
  index_set goals;

public:
  PCC(Instance& i, Stopwatch& s, const ppc_vec& v);
  PCC(Instance& i, Stopwatch& s, const ppc_vec& v, StaticMutex* m,
      const graph* lmg, set_edge_vec* t, const graph* nag);
  ~PCC();

  bool run(const index_set& s);
};

// The ConsistencyTest class provides a uniform interface to several
// consistency test methods.
class ConsistencyTest {
 protected:
  const ppc_vec& ppcs;
  index_set conflict;

  // Greedy conflict minimisation, implemented by calling the virtual
  // test method.
  void minimise_conflict();

 public:
  ConsistencyTest(const ppc_vec& constraints);
  virtual ~ConsistencyTest();

  // Test a given set (of trajectory constraints), return true iff
  // a contradiction is proven. If find_conflict is true, a conflict
  // (inconsistent subset of s) will be generated and stored; it is
  // retrieved with the last_conflict method below. If min_conflict
  // is true, the returned conflict should be subset-minimal.
  virtual bool test(const index_set& s,
		    bool find_conflict = false,
		    bool min_conflict = false) = 0;

  // Returns the conflict generated by the last failed test.
  // This may be garbage, if the most recent test did not find a
  // contradiction.
  const index_set& last_conflict() { return conflict; };
};

class PCCTest : public ConsistencyTest {
  Stopwatch& pcc_test_stats;
  count_type& n_pcc_tests;
  set_edge_vec trlm;
  PCC* tester;

 public:
  bool verbose1;
  bool verbose2;

  PCCTest(const ppc_vec& constraints,
	  Instance& ins,
	  Preprocessor& prep,
	  Stopwatch& prep_stats,
	  Stopwatch& test_stats,
	  count_type& n_tests);
  ~PCCTest();

  virtual bool test(const index_set& s, bool find_conflict, bool min_conflict);
};

class hmTest : public ConsistencyTest {
  Instance& instance;
  Stopwatch& hm_test_stats;
  count_type& n_hm_built;
  count_type& n_hm_tests;

  bool generate_conflict(const index_set& s);

 public:
  bool verbose;
  bool opt_H2;

  hmTest(const ppc_vec& constraints,
	 Instance& ins,
	 Stopwatch& test_stats,
	 count_type& n_built,
	 count_type& n_tests);
  ~hmTest();

  virtual bool test(const index_set& s, bool find_conflict, bool min_conflict);
};

class SerialTest : public ConsistencyTest {
  ConsistencyTest& _test1;
  ConsistencyTest& _test2;

 public:
  SerialTest(const ppc_vec& constraints,
	     ConsistencyTest& test1,
	     ConsistencyTest& test2);
  ~SerialTest();

  virtual bool test(const index_set& s, bool find_conflict, bool min_conflict); 
};

END_HSPS_NAMESPACE

#endif
