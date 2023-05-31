#ifndef ADDITIVE_H
#define ADDITIVE_H

#include "config.h"
#include "cost_table.h"
#include "rng.h"

BEGIN_HSPS_NAMESPACE

class AH : public Heuristic {
  void max_cost_mset(index_type m, const index_set& g, index_set& s);
  void goal_relevant_remaining_actions(const index_set& g,
				       index_set& a,
				       bool_vec& rem);

 protected:
  CostTable* Hmax;
  lvector<CostTable*> H_vec;
  Stopwatch& stats;

 public:
  static count_type Hmax_wins;
  static count_type Hsum_wins;
  static count_type draws;

  static bool use_linear_scan_eval;

  AH(Instance& i, Stopwatch& s);
  virtual ~AH();

  index_type n_additive_components() { return H_vec.length(); };
  CostTable* additive_component(index_type i) { return H_vec[i]; };
  CostTable* max_component() { return Hmax; };

  // compute additive H1/H2 with given action partitioning
  void compute_additive(const ACF& cost, const index_set_vec& p, bool useH2);

  // compute additive H1/H2 with fractional costs distributed over given
  // action sets (same as additive when p is a partitioning).
  void compute_fractional(const ACF& cost, index_set_vec& p, bool useH2);

  // compute standard H1/H2 (enables eval = max(std H1/H2, additive H1/H2))
  void compute_max(const ACF& cost, bool useH2);

  // disable standard H1/H2
  void disable_max();

  void compute_with_relevance_partitioning
    (const ACF& cost, const index_set& g);
  // void compute_with_random_relevance_partitioning
  //   (const ACF& cost, const index_set& g, RNG& rnd, bool useH2);
  void compute_with_iterative_assignment
    (const ACF& cost, const index_set& g, index_set_vec& app,
     bool useH2, bool fractional, const index_set& g_limit);

  void compute_with_iterative_assignment_1
    (const ACF& cost, const index_set& g,
     bool useH2, bool fractional, bool optimal,
     const index_set& g_limit);
  void compute_with_iterative_assignment_2
    (const ACF& cost, const index_set& g, bool useH2, bool fractional,
     const index_set& g_limit);

  void compute_with_layered_partitioning
    (const ACF& cost, const index_set& g);
  void compute_layered_action_partition
    (index_type m, const index_set& g, index_set_vec& p, index_type pg,
     bool_vec& rem);

  void compute_with_new_decomposition
    (const ACF& cost, const index_set& g, bool useH2);
  void compute_bottom_up_2
    (const ACF& cost, const index_set& g, index_type p_max, bool useH2);
  void compute_bottom_up
    (const ACF& cost, const index_set& g, bool do_glue, bool useH2);
  void compute_with_k_cuts
    (const ACF& cost, const index_set& g, index_type k, bool useH2);

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);

  virtual void write_eval(const index_set& s, std::ostream& st,
			  char* p = 0, bool e = true);
  virtual void write_eval(const bool_vec& s, std::ostream& st,
			  char* p = 0, bool e = true);

  virtual NTYPE incremental_eval(const index_set& s, index_type i_new);
  virtual NTYPE incremental_eval(const bool_vec& s, index_type i_new);
  virtual NTYPE eval_to_bound(const index_set& s, NTYPE bound);
  virtual NTYPE eval_to_bound(const bool_vec& s, NTYPE bound);

  virtual void store(const index_set& s, NTYPE v, bool opt);
  virtual void store(const bool_vec& s, NTYPE v, bool opt);
};

END_HSPS_NAMESPACE

#endif
