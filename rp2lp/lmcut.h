#ifndef LM_CUT_H
#define LM_CUT_H

#include "config.h"
#include "cost_table.h"

BEGIN_HSPS_NAMESPACE

class LMCutBase : public Heuristic {
 protected:
  const ACF& init_cost;
  Stopwatch& stats;
  CostTable* table;

  void extend_goal_set(const index_set& sgset,
		       const AnyACF& costs,
		       bool_vec& ext_goal_set);
  void find_cut(const bool_vec& ext_goal_set,
		const AnyACF& costs,
		bool_vec& cut);

  void select_relevant(const index_set& g, bool_vec& rel);
  NTYPE compute2(const bool_vec& s, const index_set& g);

 public:
  LMCutBase(Instance& i, const ACF& c, Stopwatch& s);
  virtual ~LMCutBase();

  NTYPE compute(const bool_vec& s, const index_set& g,
		const ACF& c, const bool_vec* a = 0);
  NTYPE compute(const index_set& s, const index_set& g,
		const ACF& c, const bool_vec* a = 0);
  NTYPE compute(const bool_vec& s, const index_set& g);
  NTYPE compute(const index_set& s, const index_set& g);

  NTYPE get_landmark(const bool_vec& s, const index_set& g, bool_vec& lm);
};

class ForwardLMCut : public LMCutBase {
 protected:
  index_set goals;

 public:
  ForwardLMCut(Instance& i, const index_set& g, const ACF& c, Stopwatch& s);
  virtual ~ForwardLMCut() { };

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
};

class ForwardLMCut2 : public ForwardLMCut {
 public:
  ForwardLMCut2(Instance& i, const index_set& g, const ACF& c, Stopwatch& s)
    : ForwardLMCut(i, g, c, s) { };
  virtual ~ForwardLMCut2() { };

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
};

class RegressionLMCut : public LMCutBase {
  bool_vec inits;

 public:
  RegressionLMCut(Instance& i, const ACF& c, Stopwatch& s);
  RegressionLMCut(Instance& i, const index_set& s0,
		  const ACF& c, Stopwatch& s);
  virtual ~RegressionLMCut() { };

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);

  NTYPE compute(const index_set& g, const ACF& c);
};

class RegressionLMCutP2 : public LMCutBase {
  Instance base_ins; // the real (non-P2) instance
  const s2index& atom_map;
  const index_vec& action_map;
  bool_vec inits;

 protected:
  NTYPE compute(const bool_vec& s, const index_set& g);

 public:
  // constructor args:
  // p2_ins: the P^2 instance
  // c: an action cost function for the original instance
  // base_ins: the original instance
  // pm: atom index map (s2index base_ins -> p2_ins)
  // am: action index map (p2_ins -> base_ins)
  RegressionLMCutP2(Instance& p2_ins, const ACF& c, Instance& base_ins,
		    const s2index& pm, const index_vec& am, Stopwatch& s);
  virtual ~RegressionLMCutP2() { };

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
};

END_HSPS_NAMESPACE

#endif
