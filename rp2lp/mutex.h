#ifndef STATIC_MUTEX_H
#define STATIC_MUTEX_H

#include "config.h"
#include "problem.h"
#include "heuristic.h"

BEGIN_HSPS_NAMESPACE

// implementation of h^2 optimised for static mutex
// calculation (i.e., all zero-cost actions)

class StaticMutex : public Heuristic {
  s2index   pairix;
  bool_vec  tab;

  void compute(const index_set& init, const bool_vec* aa = 0);

 public:
  StaticMutex(Instance& ins, bool compute_on_construction = true);
  StaticMutex(Instance& ins, const index_set& init);
  StaticMutex(Instance& ins, const index_set& init, const bool_vec& aa);
  StaticMutex(Instance& ins, const index_set& init, const bool_vec* aa);
  StaticMutex(Instance& ins, const StaticMutex& mx);

  index_type size() const { return pairix.size(); };

  void recompute();
  void recompute(const index_set& init);
  void recompute(const bool_vec& aa);
  void recompute(const index_set& init, const bool_vec& aa);

  bool unreachable(index_type i) const;
  bool mutex(index_type i, index_type j) const;
  bool mutex(const index_set& s) const;
  bool mutex(const index_set& s, index_type i) const;
  bool mutex(const index_set& s1, const index_set& s2) const;

  void write(std::ostream& s) const;

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
  virtual NTYPE incremental_eval(const index_set& s, index_type i_new);
  virtual NTYPE incremental_eval(const bool_vec& s, index_type i_new);
  virtual NTYPE eval(index_type atom);
  virtual NTYPE eval(index_type atom1, index_type atom2);
};

// similarly optimised h^1 for reachability

class Reachability : public Heuristic {
  bool_vec  tab;
  //bool_vec  atab;

  index_vec rem_pre;
  index_vec q;
  index_type q_head;
  index_type q_tail;

  void compute_main(read_only_vector_with_default<bool>& aa);
  void compute(const index_set& init, const bool_vec* aa = 0);
  void compute(const bool_vec& init, const bool_vec* aa = 0);

 public:
  Reachability(Instance& ins, bool compute_on_construction = true);
  Reachability(Instance& ins, const index_set& init);
  Reachability(Instance& ins, const index_set& init, const bool_vec& aa);
  Reachability(Instance& ins, const index_set& init, const bool_vec* aa);

  void init_structs();

  struct reachability_state {
    bool_vec tab;
    //bool_vec atab;
    index_vec rem_pre;
    index_type q_head;
  };

  void save_state(reachability_state& s) const;
  void restore_state(const reachability_state& s);

  void recompute();
  void recompute(const index_set& init);
  void recompute(const bool_vec& aa);
  void recompute(const index_set& init, const bool_vec& aa);
  void recompute(const index_set& init, const index_set& aa);
  void recompute_bv_init(const bool_vec& init);
  void update(const bool_vec& aa);

  // special form of update: new_aa is the only newly allowed
  // action whose preconditions are reachable.
  void update(const bool_vec& aa, index_type new_aa);

  // compute reachability until goal is true, and output the order
  // in which actions were applied; returns false if goal never reached.
  bool order_relaxed_plan(const index_set& init, const index_set& goal,
			  const index_set& aa, index_vec& ao);

  const bool_vec& reachable() const { return tab; };

  bool unreachable(index_type i) const;
  bool unreachable(const index_set& s) const;
  bool unreachable(const bool_vec& s) const;

  // check (un)reachability of an actions precondition
  bool unreachable_action(index_type a) const { return (rem_pre[a] > 0); };

  // find the subset of s that is unreachable
  void unreachable_subset(const index_set& s, index_set& u) const;

  void write(std::ostream& s) const;

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
  virtual NTYPE incremental_eval(const index_set& s, index_type i_new);
  virtual NTYPE incremental_eval(const bool_vec& s, index_type i_new);
  virtual NTYPE eval(index_type atom);
  virtual NTYPE eval(index_type atom1, index_type atom2);
};

// compute atom landmark graph via P^2 + Reachability
// atms is the set of atoms to consider (on either side of the
// landmark relation), acts is the set of "allowed" actions (i.e.,
// all actions not in this set are not considered when computing
// reachability)
void landmark_graph_viaP2(Instance& ins,
			  const bool_vec& atms,
			  const bool_vec& acts,
			  graph& g);

// versions of the landmark_graph_viaP2 procedure omitting
// initialisation steps (computation of initial h^2,
// construction of P^2)
void landmark_graph_viaP2(Instance& ins,
			  const bool_vec& atms,
			  const bool_vec& acts,
			  Heuristic* inc,
			  graph& g);

void landmark_graph_viaP2(const bool_vec& atms,
			  const bool_vec& acts,
			  Instance& insP2,
			  s2index& pair_map,
			  index_vec& act_map,
			  graph& g);

// generate a list of landmark relations that do not hold (are not
// in the graph) but would hold iff the associated set of actions
// are removed.
typedef lvector< std::pair<index_set, index_pair> > set_edge_vec;

void landmark_graph_triggered_edges(Instance& ins,
				    graph& lmg,
				    set_edge_vec& tev);

class ForwardReachabilityCheck : public Reachability {
  index_set goals;

 public:
  ForwardReachabilityCheck(Instance& i, const index_set& g);
  virtual ~ForwardReachabilityCheck();

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
};

END_HSPS_NAMESPACE

#endif
