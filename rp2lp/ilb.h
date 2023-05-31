
#ifndef ILB_H
#define ILB_H

#include "config.h"
#include "mutex.h"
#include "search.h"
#include "stats.h"

#include "cost_table.h"
#include "lmcut.h"

//#define USE_CACHE
//#define RANDOMIZE
#define CPLEX_INCREMENTAL

#ifdef USE_CACHE
#include <ext/hash_map>
#endif

#ifdef RANDOMIZE
#include "rng.h"
#endif

#ifdef CPLEX_INCREMENTAL
#include <ilcplex/ilocplex.h>
#endif

BEGIN_HSPS_NAMESPACE

class EvalWithMetaAtoms : public Heuristic {
  Heuristic& base_h;
  const rule_set& meta_atom_map;

 public:
  EvalWithMetaAtoms(Heuristic& h, const rule_set& map)
    : Heuristic(h.get_instance()), base_h(h), meta_atom_map(map) { };
  ~EvalWithMetaAtoms() { };

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
};

class ILB {
  Instance& ins;
  bool is_sub; // flag this ILB object as subroutine of another ILB.
  index_type n_real_atoms;
  index_type n_real_actions;
  const ACF& cost;
  Heuristic* inc;
  Stopwatch& stats;

 public:
 
	int probnum;
   
   int **SstepOrder;
   int **adjorddes;
   int *nadjorddes;
   int *sadjorddes;
   int **adjorddesrev;
   int *nadjorddesrev;
   int *sadjorddesrev;
   int **triangles;
   int ntriangles = 0;
   int striangles = 1000;
   int betterhplus = 0;
   int betterlmcut = 0;
   

  Stopwatch hs_stats;
  Stopwatch hss_stats;
  Stopwatch hsa_stats;
  Stopwatch lm_stats;
  Stopwatch h1_stats;
  Stopwatch ra_stats;
  Stopwatch rpn_stats;
  Stopwatch ce_stats;
  Stopwatch sce_stats;
  count_type calls_to_hplus;
  count_type calls_to_hplus_ce;
  count_type calls_to_newlm;
  count_type calls_to_hs_opt;
  double     relevant_actions_ratio;
  double     sum_density;
  count_type calls_to_hs_apx;
  count_type hs_apx_improve;
  count_type sum_width;
  count_type n_width;

  count_type hs_nodes;
  count_type hs_cache_hits;
  count_type hs_cache_miss;
  count_type hs_lb_calls;
  count_type hs_lb1_max;
  count_type hs_lb2_max;
  count_type hs_lb3_max;
  count_type hs_lb4_max;
  count_type hs_split1;
  count_type hs_splits;
  count_type hs_branch;

  count_type hpp_iterations;
  count_type conflict_set_size;
  double     pc_actions_ratio;

  // h+ with cond effs:
  count_type n_ce_compiled; // total
  count_type max_ce_compiled; // per action
  double     cce_actions_ratio;

 private:
#ifdef NTYPE_RATIONAL
  long ac_div;
#endif

  class set_compare {
  public:
    bool operator()(const HSPS::bool_vec& s0, const HSPS::bool_vec& s1) const;
  };

  class set_hash {
  public:
    static set_hash_function f;
    unsigned int operator()(const HSPS::bool_vec& s0) const;
    unsigned int operator()(const HSPS::index_set& s0) const;
  };

#ifdef USE_CACHE
  typedef __gnu_cxx::hash_map<bool_vec, NTYPE, set_hash, set_compare> set_map;
#endif

  rule_set        meta_atom_map;

  bool_vec        relevant_atoms;
  bool_vec        relevant_actions;
  bool            FA_and_ALM_avail;
  index_set_vec   first_achievers;
  index_set_vec   atom_landmarks;

  index_set_vec   landmarks;
  index_set       lm_actions;
  index_set_vec   lm_with;
  SBR             lm_conflict;
  index_set_vec   lm_conflict_set;
  cost_vec        lm_min_cost;
  cost_vec        lm_wd;
  lvector<index_vec> lm_branch_order;
  index_set       dominated;
#ifdef USE_CACHE
  set_map         store_lb;
#endif
  index_vec       score;
  index_set       rp_actions;
  NTYPE           rp_cost;
  index_set_graph rp_dg;
  Reachability    reach;

  CostTable*      h1;

#ifdef RANDOMIZE
  LC_RNG rng;
#endif

#ifdef CPLEX_INCREMENTAL
  IloEnv cplex_env;
  IloModel cplex_model;
  IloBoolVarArray cplex_vars;
  index_type cplex_model_actions;
  index_type cplex_model_atoms;
  index_type cplex_model_landmarks;
  IloCplex cplex_solver;
#endif

  // utils for copying stuff to/from subinstances
  void copy_options(const ILB& ilb, int vl_adjustment = 0);
  void add_hplus_stats(const ILB& ilb);

  void init_preprocessing();
  // compute first achievers and atom-atom landmark relation; on first
  // call after init (i.e., if FA_and_ALM_avail flag is false), atom
  // landmarks are re-initialised; on subsequent calls, only new atom
  // landmarks are added, and the method returns true if any new
  // landmark was found.
  bool compute_FA_and_ALM(const index_set& init);

  // compute L1C1 relevance; this is valid for delete-relaxed problem
  void compute_relevant(const index_set& init, const index_set& goal);
  void L1C1(const index_set& goal, bool_vec& ratm, bool_vec& ract);

  // remove (init-dependent) relaxed dominated actions from the set
  // of relevant actions
  void remove_dominated_actions(const index_set& init);

  // find an optimal hitting set over current set of landmarks
  NTYPE hitting_set(const bool_vec& rem,
		    bool_vec& out,
		    bool_vec& set,
		    NTYPE acc,
		    index_set& best,
		    NTYPE& ub,
		    NTYPE lb);

  // alt. implementation of optimal hitting set, using splitting into
  // independent subproblems
  NTYPE hitting_set_special_case_1(index_type rem, const bool_vec& out,
				   index_type& opt);
  NTYPE hitting_set_special_case_2(index_type rem1, index_type rem2,
				   const bool_vec& out, index_set& opt);
  NTYPE hitting_set_split(const bool_vec& rem, bool_vec& out,
			  bool_vec& set, NTYPE acc, index_set& best,
			  NTYPE& ub, NTYPE lb, index_type depth);
  NTYPE hitting_set_branch(const bool_vec& rem, bool_vec& out,
			   bool_vec& set, NTYPE acc, NTYPE est,
			   index_set& best, NTYPE& ub, NTYPE lb,
			   index_type depth);

  NTYPE weighted_degree(index_type i);
  void update_weighted_degree();
  NTYPE hitting_set_lb_1(const bool_vec& rem, const bool_vec& out);
  NTYPE hitting_set_lb_2(const bool_vec& rem, const bool_vec& out);
  NTYPE hitting_set_lb_3(const bool_vec& rem, const bool_vec& out);
  NTYPE hitting_set_lb_4(const bool_vec& rem, const bool_vec& out);
  NTYPE hitting_set_lb(const bool_vec& rem, const bool_vec& out);

  // forward checking: returns true if rem + out is inconsistent
  // (called only with USE_OUT).
  bool hitting_set_fc(const bool_vec& rem, const bool_vec& out);

  // find a not-necessarily-optimal hitting set over current landmarks
  NTYPE apx_hitting_set(index_set& hs, NTYPE bound = POS_INF);

  // optimal hitting set using extern solver (compile option)
  NTYPE hitting_set_extern(index_set& best, NTYPE& ub, NTYPE lb);

  void dump_hitting_set_problem_wcnf(NTYPE lb, NTYPE ub, std::ostream& to);
  void dump_hitting_set_problem_lp(NTYPE lb, NTYPE ub, std::ostream& to);

  // admin methods
  bool is_dominated(index_type a, const index_set& s) const;
  index_type add_new_landmark(const index_set& new_lm);

  // different ways of generating a new landmark (note: rp_actions
  // is an implicit parameters to all of them)

  void add_redundant_actions_to_sets(bool_vec& a, bool_vec& b);
  index_type make_new_landmark_ST(const index_set& init,
				  const index_set& goal,
				  const index_set& rp_set);

  void select_pcf(const index_set& g, bool_vec& rel);
  index_type make_new_landmark_BC(const index_set& init,
				  const index_set& goal,
				  const index_set& rp_set);

  // initialise ILA with lmcut loop
  NTYPE initialise_with_lmcut(const index_set& init,
			      const index_set& goal);

  // compute relaxed plan using the iterative landmark algorithm
  NTYPE compute_relaxed_plan_ILA(const index_set& init,
				 const index_set& goal,
				 NTYPE bound);

  // compute relaxed plan using extern (MIP) solver
  NTYPE compute_relaxed_plan_extern(const index_set& init,
				    const index_set& goal);

  //// obsolete
  // // standard FF-style non-optimal relaxed plan extraction.
  // // h1c is the h^1 computed with action costs; h1u is an h^1
  // // computed using any strictly positive cost (this is needed
  // // for tie-breaking to avoid cycles).
  // void extract_rp(Heuristic* h1c, Heuristic* h1u,
  // 		  const index_set& g, bool_vec& holds, index_vec& stack);

  // compute a non-optimal relaxed plan.
  NTYPE compute_relaxed_plan_FF(const index_set& init,
				const index_set& goal);

  // recursive subroutine of compute_relaxed_plan_BB.
  void rp_bb(const index_set& goal, bool_vec& holds, index_type next,
	     NTYPE lb, const index_vec& action_choice_order,
	     const index_set& zero_cost_actions, bool_vec& rem_acts,
	     bool_vec& plan, NTYPE acc, index_set& best, NTYPE& ub);

  // recursive subroutine of compute_relaxed_plan_BDGBT.
  NTYPE rp_regress(bool_vec& goals, bool_vec& supported, graph& prec,
		   bool_vec& plan, NTYPE acc, NTYPE& ub, index_set& best);

  bool write_hplus_wcnf(const index_set& init,
			const index_set& goal,
			std::ostream& to,
			count_type climit = 0);

 public:
  // branch-and-bound forward search for relaxed plan.
  NTYPE compute_relaxed_plan_BB(const index_set& init, const index_set& goal);

  // regression branch-and-bound relaxed planning procedure from
  // Bartak, Dvorak, Gemrot, Brom & Toropila (FLAIRS 2012).
  NTYPE compute_relaxed_plan_BDGBT(const index_set& init,
				   const index_set& goal);

 private:

  void filter_flaws(index_set_vec& cs);

  // note: if the rp is optimal, only zero cost actions can be
  // redundant, so if is_optimal == true, only those are checked.
  void compute_non_redundant_rp(const index_set& init,
				const index_set& goal,
				bool is_optimal);

  void compute_dependency_graph(const index_set& init,
				const index_set& goal,
				index_set_graph& dg);

  void compute_dependency_closure(const index_set& rp,
				  const index_set_graph& dg,
				  index_type n_from,
				  index_type n_to,
				  index_set& dc_labels);

  void partition_applicable(const index_vec& rp,
			    const bool_vec& rem_rp,
			    const bool_vec& rstate,
			    index_set_vec& apps);

  typedef std::pair<index_pair, index_type> node_conflict;
  typedef svector<node_conflict> node_conflict_set;

  index_type estimate_conflict_weight(const index_set_graph& dg,
				      const index_pair& c,
				      index_type& ncd);
  void choose_atoms_on_path(const index_set_graph& dg,
			    const index_type n_from,
			    const index_type n_to,
			    index_set& as);
  void atom_sets_on_path(const index_set_graph& dg,
			 const index_type n_from,
			 const index_type n_to,
			 index_set_vec& ps);
  void node_conflict_to_atom_conflicts(const index_set_graph& dg,
				       const node_conflict& nc,
				       index_set_vec& acs);
  void node_conflicts_to_atom_conflicts(const index_set_graph& dg,
					const node_conflict_set& ncs,
					index_set_vec& acs);
  void node_conflict_to_atom_conflicts_with_dc(const index_set& rp,
					       const index_set_graph& dg,
					       const node_conflict& nc,
					       index_set_vec& acs);

  void meta_atom_to_atom_set(index_type p,
			     const rule_set& map,
			     index_set& set);

  void conflicts2(index_type n_fail,
		  const index_set& c_fail,
		  lvector<bool_vec>& state,
		  lvector<bool_vec>& rstate,
		  const index_vec& rp,
		  index_type s_fail,
		  const index_set_vec& seq,
		  bool check_ce,
		  node_conflict_set& cs);

  bool exec_rp2(bool_vec& rem_rp,
		lvector<bool_vec>& state,
		lvector<bool_vec>& rstate,
		const index_set& goal,
		index_type step,
		index_set_vec& seq,
		lvector<node_conflict_set>& cs);

  bool check_rp2(const index_set& init,
		 const index_set& goal,
		 index_set_vec& seq,
		 index_set_vec& cs);

  bool simple_exec_rp(bool_vec& rem_rp,
		      index_type step,
		      lvector<bool_vec>& state,
		      const index_set& goal,
		      const index_vec& rp,
		      index_set_vec& seq);

  bool rp_is_plan(const index_set& init,
		  const index_set& goal,
		  const index_vec& rp,
		  index_set_vec& seq);


  struct PotentialNodeConflict {
    index_pair pair;
    index_type dp; // deleted precondition
    index_set  rescue;
    index_type weight;
    PotentialNodeConflict()
      : pair(index_pair(no_such_index)), dp(no_such_index), weight(0) { };
    PotentialNodeConflict(index_type d, index_type f, index_type p)
      : pair(index_pair(d, f)), dp(p), weight(0) { };
  };

  class PNCWeightIncreasing : public lvector<PotentialNodeConflict>::order {
  public:
    virtual bool operator()(const PotentialNodeConflict& v0,
			    const PotentialNodeConflict& v1) const;
  };

  void choose_next_conflict(const lvector<PotentialNodeConflict>& pc,
			    bool_vec& rem,
			    graph& prec,
			    index_set& current,
			    index_type w_current,
			    index_set& best,
			    index_type w_best);

  void remove_self_pairs(index_set_vec& cs);

  // check if cs contains any atom set not represented in the current
  // meta-atom map.
  bool contains_new_meta_atom(const index_set_vec& cs) const;

  void conflicts3(const index_set& init,
		  const index_set& goal,
		  index_set_vec& cs);

  bool check_rp3(const index_set& init,
		 const index_set& goal,
		 index_set_vec& seq,
		 index_set_vec& cs);

  // debug method; runs both check_rp2 and check_rp3 to verify
  // that they agree on plan validity.
  //bool check_rp_both(const index_set& init,
  //		     const index_set& goal,
  //		     index_set_vec& seq,
  //		     index_set_vec& cs);

  // switching method: calls some check_rpX method, depending
  // on the check_rp_strategy option.
  bool check_rp(const index_set& init,
		const index_set& goal,
		index_set_vec& seq,
		index_set_vec& cs);

  // helper methods for hplus_with_ce:

  index_type get_ce_atom(const Instance& irce,
			 index_type first_rce_atom,
			 index_type act);

  bool sequence_dg(const Instance& irce,
		   index_type first_rce_atom,
		   index_set_graph& dg,
		   index_vec& seq);

  bool schedule_ce(const Instance& irce,
		   index_type first_rce_atom,
		   index_type first_rce_action,
		   index_set_graph& dg,
		   index_set& failed);

  // schedule ce's non-optimally; this may insert new non-ce actions
  // into the dg.
  void schedule_ce_nonopt(const Instance& irce,
			  index_type first_rce_atom,
			  index_type first_rce_action,
			  index_set_graph& dg,
			  index_set& failed);

  // schedule_ce_rec is the internal, recursive part of schedule_ce
  bool schedule_ce_rec(const Instance& irce,
		       index_type first_rce_atom,
		       index_type first_rce_action,
		       index_set_graph& dg,
		       index_set& failed);

  // internal version of hplus_with_ce; returns the plan in two
  // components: rp is the sequence of actions; rp_ce is a
  // corresponding sequence of sets of triggered ce's (indexed
  // into the set of cond-adds of the action).
  NTYPE hplus_with_ce(const index_set& init, const index_set& goal,
		      NTYPE bound, index_vec& rp, index_set_vec& rp_ce);


  // helper methods for hplusplus_with_ce

  void partition_applicable_with_ce(const index_vec& rp,
				    const index_set_vec& rp_ce,
				    const bool_vec& rem_rp,
				    const bool_vec& rstate,
				    index_set_vec& apps);

  void choose_from_node_conflict_set(const index_set_graph& dg,
				     const node_conflict_set& ncs,
				     index_set_vec& cs);

  // choose-methods using full dc
  void compute_dependency_closure_in_split_dg
    (const index_vec& rp,
     const index_set_vec& rp_ce,
     const pair_vec& split_index,
     const index_set_graph& dg,
     index_type n_from,
     index_type n_to,
     index_set& dc_labels);

  void node_conflict_to_atom_conflicts_with_dc_in_split_dg
    (const index_vec& rp,
     const index_set_vec& rp_ce,
     const pair_vec& split_index,
     const index_set_graph& dg1,
     const index_set_graph& dg2,
     const node_conflict& nc,
     index_set_vec& acs);

  void choose_from_node_conflict_set_with_dc
    (const index_vec& rp,
     const index_set_vec& rp_ce,
     const pair_vec& split_index,
     const index_set_graph& dg1,
     const index_set_graph& dg2,
     const node_conflict_set& ncs,
     index_set_vec& cs);

  void split_deleters(const pair_vec& split_index,
		      const index_set_graph& split_dg,
		      node_conflict_set& ncs);

  bool exec_rp_with_ce(bool_vec& rem_rp,
		       lvector<bool_vec>& state,
		       lvector<bool_vec>& rstate,
		       const index_set& goal,
		       const index_vec& rp,
		       const index_set_vec& rp_ce,
		       const pair_vec& split_index,
		       const index_set_graph& split_dg,
		       const index_set_graph& dg2,
		       bool with_full_dc,
		       index_type step,
		       index_set_vec& seq,
		       index_set_vec& cs);

  void make_rp_instance(const index_vec& rp,
			const index_set_vec& rp_ce,
			Instance& rp_ins,
			pair_vec& split_index);

  void normalise_rp_with_ce(index_vec& rp,
			    index_set_vec& rp_ce);

  void compute_split_dependency_graph(const index_vec& rp,
				      index_set_vec& rp_ce,
				      pair_vec& split_index,
				      index_set_graph& dg);
  void add_back_edges(const pair_vec& split_index,
		      index_set_graph& split_dg);

  bool check_rp_with_ce(const index_set& init,
			const index_set& goal,
			const index_vec& rp,
			const index_set_vec& rp_ce,
			bool with_full_dc,
			index_set_vec& seq,
			index_set_vec& cs);

  // some general utilities

  void print_rp(std::ostream& s,
		const index_set& rp);

  void print_rpdg(std::ostream& s,
		  const index_set& rp,
		  const index_set_graph& dg);

  void print_rp_with_ce(std::ostream& s,
			const index_vec& rp,
			const index_set_vec& rp_ce);

  bool validate_rp_with_ce(const index_vec& rp,
			   const index_set_vec& rp_ce);

  static void new_lb(NTYPE lb);
  static void reset_hlb();

 public:
  int   verbose_level;
  int   check_rp_strategy;
  bool  remove_dominated_conditions;
  bool  prune_relaxed_irrelevant;
  bool  prune_relaxed_dominated;
  bool  iterated_preprocessing;
  bool  ILA_use_approximate;
  bool  ILA_use_saturation;
  bool  ILA_use_lmcut;
  bool  zero_cost_fill;
  bool  hitting_set_use_split;
  bool  hitting_set_use_dominance;

  class ConflictModifier {
  public:
    ConflictModifier();
    virtual ~ConflictModifier();
    virtual void apply(ILB& ilb, index_set_vec& cs);
  };
  friend class ConflictModifier;
  ConflictModifier* conflict_mod;

  static NTYPE hlb;

  ILB(Instance& ins, const ACF& cost, Heuristic* inc, Stopwatch& stats);
  ~ILB();

  const Instance& instance() const { return ins; };

  NTYPE hplus(const index_set& init,
	      const index_set& goal,
	      NTYPE bound,
	      Plan** sol);

  NTYPE hplus_with_selected_actions(const index_set& init,
				    const index_set& goal,
				    const index_set& acts,
				    NTYPE bound,
				    Plan** sol);

  NTYPE hplus_with_ce(const index_set& init, const index_set& goal,
		      NTYPE bound, Plan** sol);

  NTYPE hplus_with_ce_time_indexed(const index_set& init,
				   const index_set& goal,
				   NTYPE bound,
				   index_vec& rp, index_set_vec& rp_ce);

  NTYPE hplusplus(NTYPE bound, Plan** sol);

  // hplusplus1 does a single h++ iteration, starting with compiling
  // in an input set of conflicts (cs_in), then computing a relaxed
  // plan, checking it, and returning an updated bound, a new set of
  // conflict (cs_out) or a plan (sol, if given); on return, the
  // completed flag is set to true iff the iteration was completed,
  // i.e., not interrupted by a stat limit. For the first call, set
  // the initialise flag to true.
  NTYPE hplusplus1(const index_set_vec& cs_in,
		   index_set_vec& cs_out,
		   Plan** sol,
		   bool& completed,
		   bool initialise = false);

  NTYPE hplusplus_with_ce(NTYPE bound, Plan** sol);

  // access some internal state variables from the last call
  const index_set_vec& action_landmarks() const { return landmarks; };
  const index_set& relaxed_plan() const { return rp_actions; };
  const rule_set& get_meta_atom_map() const { return meta_atom_map; };

  bool validate_rp(const index_set& rp);

  void print_hplus_stats(std::ostream& s) const;
  void print_hplusplus_stats(std::ostream& s) const;

  void save_hplus_wcnf(const index_set& init,
		       const index_set& goal,
		       const char* pname,
		       count_type climit = 0);

  // direct access to internal methods -- calling these will mess up the
  // ILB object state -- do not call unless you know what you're doing!!!
  bool is_relaxed_plan(const bool_vec& set);
  const index_set& make_new_landmark(const index_set& set);
  void compute_relevant_actions();

  // debugging method
  void test(const index_vec& plan);

 private:
  ///
  /// Garbage and obsolete stuff
  ///


  // yet another implementation of optimal hitting set
  // NTYPE hitting_set_rc(index_set& hs);

  // void conflicts
  //   (index_type failed, lvector<bool_vec>& state, lvector<bool_vec>& rstate,
  //    const index_set& goal, index_type s_fail, index_set_vec& seq,
  //    pair_set& cs);
  //
  // void conflicts_simple_min
  //   (index_type failed, lvector<bool_vec>& state, lvector<bool_vec>& rstate,
  //    const index_set& goal, index_type s_fail, index_set_vec& seq,
  //    pair_set& cs);
  //
  // bool exec_rp(bool_vec& rem_rp,
  //	       lvector<bool_vec>& state,
  //	       lvector<bool_vec>& rstate,
  //	       const index_set& goal,
  //	       index_type step,
  //	       index_set_vec& seq,
  //	       pair_set& failed);
  //
  // bool check_rp(const index_set& init,
  //		const index_set& goal,
  //		index_set_vec& seq,
  //		pair_set& cs);

  // struct RPState : public State {
  //   Instance&  ins;
  //   const ACF& cost;
  //   RegressionLMCut& hlmc;
  //   CostTable& h1;
  //   index_type d;
  //   index_pair e_new;
  //   NTYPE      delta;
  //   NTYPE      est;
  //   index_set  plan;
  //   graph      prec;
  //   pair_set   oc;
  // 
  //   RPState(Instance& i,
  // 	    const ACF& c,
  // 	    RegressionLMCut& hlmc,
  // 	    CostTable& h1,
  // 	    const index_set& goal);
  //   RPState(const RPState& s);
  //   virtual ~RPState();
  // 
  //   virtual Transition* transition();
  //   virtual NTYPE delta_cost();
  //   virtual NTYPE acc_cost();
  //   virtual index_type depth() const;
  //   virtual NTYPE est_cost();
  //   virtual bool  is_final();
  //   virtual bool  is_max();
  //   virtual NTYPE expand(Search& s, NTYPE bound);
  //   virtual void reevaluate();
  //   virtual int compare(const State& s);
  //   virtual index_type hash();
  //   virtual State* copy();
  //   virtual void write(::std::ostream& s);
  // };
  // 
  // struct RPTransition : public Transition {
  //   index_pair edge;
  //   RPTransition(const index_pair& e);
  //   RPTransition(const Transition& t);
  // 
  //   int compare(const Transition& t);
  //   void insert(Plan& p);
  //   void insert_path(Plan& p);
  //   void write(::std::ostream& s) const;
  // };
  // 
  // struct RP_BFS_SearchResult : public SearchResult {
  //   index_set plan;
  //   bool solved;
  // 
  //   RP_BFS_SearchResult();
  //   virtual ~RP_BFS_SearchResult();
  // 
  //   virtual void solution(State& s, NTYPE cost);
  //   virtual void solution(State& s, Transition* p, NTYPE cost);
  //   virtual bool more();
  // };
  // 
  // struct RPState2 : public State {
  //   Instance&        ins;
  //   const ACF&       cost;
  //   const index_vec& action_choice_order;
  //   const index_vec& zero_cost_actions;
  //   ForwardLMCut&    hlmc;
  //   const index_set& goal;
  //   bool_vec         holds;
  //   index_type       last_choice;
  //   NTYPE            est;
  // 
  //   RPState2(Instance& i,
  // 	     const ACF& c,
  // 	     const index_vec& o,
  // 	     const index_vec& z,
  // 	     ForwardLMCut& h,
  // 	     const index_set& s0,
  // 	     const index_set& g);
  //   RPState2(const RPState2& s);
  //   virtual ~RPState2();
  // 
  //   RPState2* apply_choice(index_type k);
  //   void apply_fixpoint();
  //   void reevaluate();
  // 
  //   virtual Transition* transition();
  //   virtual NTYPE delta_cost();
  //   virtual NTYPE est_cost();
  //   virtual bool  is_final();
  //   virtual bool  is_max();
  //   virtual NTYPE expand(Search& s, NTYPE bound);
  //   virtual int compare(const State& s);
  //   virtual index_type hash();
  //   virtual State* copy();
  //   virtual void write(::std::ostream& s);
  // };
  // 
  // NTYPE compute_relaxed_plan_by_search(const index_set& init,
  // 				       const index_set& goal);

};

class ForwardHPlus : public Heuristic {
  ILB ilb;
  index_set goal;
 public:
  ForwardHPlus(Instance& i, const ACF& c, const index_set& g, Stopwatch& s)
    : Heuristic(i), ilb(i, c, 0, s), goal(g) 
  {
    ilb.verbose_level = 0;
    ilb.prune_relaxed_irrelevant = true;
    ilb.prune_relaxed_dominated = false;
    ilb.ILA_use_approximate = true;
    ilb.ILA_use_saturation = true;
    ilb.zero_cost_fill = true;
    ilb.hitting_set_use_split = true;
    ilb.hitting_set_use_dominance = true;
  };
  ~ForwardHPlus() { };

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);

  void print_stats(std::ostream& s) const;
};

END_HSPS_NAMESPACE

#endif
