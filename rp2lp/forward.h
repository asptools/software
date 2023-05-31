#ifndef SEQ_PROGRESSION_STATE_H
#define SEQ_PROGRESSION_STATE_H

#include "config.h"
#include "atomset.h"
#include "resource.h"
#include "cost_table.h"
#include "preprocess.h"

BEGIN_HSPS_NAMESPACE

struct SeqProgTrans : public ProgressionPath, public SequentialTransition
{
  SeqProgTrans(index_type a, NTYPE d) : SequentialTransition(a, d) { };
  virtual int compare(const Transition& t) {
    return SequentialTransition::compare((const SeqProgTrans&)t);
  };
  virtual Transition* copy() {
    return new SeqProgTrans(act, delta);
  };
  virtual void insert(Plan& p) {
    SequentialTransition::insert(p);
  };
  virtual void write(std::ostream& s) const {
    SequentialTransition::write(s);
  };
  virtual void insert_path(Plan& p) {
    ProgressionPath::insert_path(p);
  };
};

class SeqProgState : public AtomSetState, public StateFactory
{
 protected:
  const ACF& cost;

  // action leading from predecessor
  index_type     act;

  // consumable resource state (NIL if not planning with resources)
  BasicResourceState* res;

  virtual bool applicable(Instance::Action& a);

  // return (a copy of) the state that results from applying action a to
  // the present state (present state becomes predecessor of the new state).
  virtual SeqProgState* apply(Instance::Action& a);

  bool is_update(Instance::Action& a) {
    return (a.del.empty() && IS_ZERO(cost(a.index)) && (res == 0));
  };

  void apply_update_actions();

 public:
  SeqProgState(Instance& i, Heuristic& h, const ACF& c);
  SeqProgState(Instance& i, Heuristic& h, const ACF& c, const index_set& s);
  SeqProgState(Instance& i, Heuristic& h, const ACF& c, const bool_vec& s);

  SeqProgState(Instance& i, Heuristic& h, const ACF& c,
	       BasicResourceState* r);
  SeqProgState(Instance& i, Heuristic& h, const ACF& c, const index_set& s,
	       BasicResourceState* r);
  SeqProgState(Instance& i, Heuristic& h, const ACF& c, const bool_vec& s,
	       BasicResourceState* r);

  SeqProgState(const SeqProgState& s);
  virtual ~SeqProgState();

  static  bool separate_update_actions;

  virtual bool is_final();
  virtual Transition* transition();
  virtual NTYPE delta_cost();
  virtual NTYPE expand(Search& s, NTYPE bound);
  virtual void store(NTYPE cost, bool opt);
  virtual void set_predecessor(State* p);

  virtual void insert(Plan& p);
  virtual void insert_path(Plan& p);
  virtual void write_plan(std::ostream& s);
  virtual State* new_state(const index_set& s, State* p);
  virtual State* new_state(const bool_vec& s, State* p);
  virtual State* copy();

  virtual int compare(const State& s);
  index_type hash();
  virtual void write(std::ostream& s);
};

class SeqCProgState : public SeqProgState
{
 public:
  SeqCProgState(Instance& i, Heuristic& h, const ACF& c)
    : SeqProgState(i, h, c) { };
  SeqCProgState(Instance& i, Heuristic& h, const ACF& c, const index_set& s)
    : SeqProgState(i, h, c, s) { };
  SeqCProgState(Instance& i, Heuristic& h, const ACF& c, const bool_vec& s)
    : SeqProgState(i, h, c, s) { };

  SeqCProgState(Instance& i, Heuristic& h, const ACF& c, BasicResourceState* r)
    : SeqProgState(i, h, c, r) { };
  SeqCProgState(Instance& i, Heuristic& h, const ACF& c,
	       const index_set& s, BasicResourceState* r)
    : SeqProgState(i, h, c, s, r) { };
  SeqCProgState(Instance& i, Heuristic& h, const ACF& c,
	       const bool_vec& s, BasicResourceState* r)
    : SeqProgState(i, h, c, s, r) { };

  SeqCProgState(const SeqProgState& s)
    : SeqProgState(s) { };

  virtual NTYPE expand(Search& s, NTYPE bound);

  virtual State* new_state(const index_set& s, State* p);
  virtual State* new_state(const bool_vec& s, State* p);
  virtual State* copy();

  virtual int compare(const State& s);
  index_type hash();
  virtual void write(std::ostream& s);
};


class RedSeqProgState : public SeqProgState {
  ReducedActionSet& ras;
  index_type d_max;
 public:
  RedSeqProgState(Instance& i, Heuristic& h,
		  const ACF& c, ReducedActionSet& r, index_type d = 0)
    : SeqProgState(i, h, c), ras(r), d_max(d) { };
  RedSeqProgState(Instance& i, Heuristic& h,
		  const ACF& c, const index_set& s,
		  ReducedActionSet& r, index_type d = 0)
    : SeqProgState(i, h, c, s), ras(r), d_max(d) { };
  RedSeqProgState(Instance& i, Heuristic& h,
		  const ACF& c, const bool_vec& s,
		  ReducedActionSet& r, index_type d = 0)
    : SeqProgState(i, h, c, s), ras(r), d_max(d) { };

  RedSeqProgState(const RedSeqProgState& s)
    : SeqProgState(s), ras(s.ras), d_max(s.d_max) { };

  virtual NTYPE expand(Search& s, NTYPE bound);

  virtual State* new_state(const index_set& s, State* p);
  virtual State* new_state(const bool_vec& s, State* p);
  virtual State* copy();
};


class RestrictedSeqProgState : public SeqProgState
{
 protected:
  const bool_vec& p_atoms;
  bool_vec& f_atoms;

 public:
  RestrictedSeqProgState(Instance& i, const bool_vec& pa, bool_vec& fa,
			 Heuristic& h, const ACF& c, const index_set& s)
    : SeqProgState(i, h, c, s), p_atoms(pa), f_atoms(fa) { };
  RestrictedSeqProgState(Instance& i, const bool_vec& pa, bool_vec& fa,
			 Heuristic& h, const ACF& c, const bool_vec& s)
    : SeqProgState(i, h, c, s), p_atoms(pa), f_atoms(fa) { };

  RestrictedSeqProgState(Instance& i, const bool_vec& pa, bool_vec& fa,
			 Heuristic& h, const ACF& c, const index_set& s,
			 BasicResourceState* r)
    : SeqProgState(i, h, c, s, r), p_atoms(pa), f_atoms(fa) { };
  RestrictedSeqProgState(Instance& i, const bool_vec& pa, bool_vec& fa,
			 Heuristic& h, const ACF& c, const bool_vec& s,
			 BasicResourceState* r)
    : SeqProgState(i, h, c, s, r), p_atoms(pa), f_atoms(fa) { };

  RestrictedSeqProgState(const RestrictedSeqProgState& s)
    : SeqProgState(*this), p_atoms(s.p_atoms), f_atoms(s.f_atoms) { };
  virtual ~RestrictedSeqProgState() { };

  virtual NTYPE expand(Search& s, NTYPE bound);

  virtual State* new_state(const index_set& s, State* p);
  virtual State* new_state(const bool_vec& s, State* p);
  virtual State* copy();
};

class RelaxedSeqProgState : public SeqProgState
{
 protected:
  index_set x_atoms;
  bool_vec  x_action;

  void compute_x_actions();
  void closure();

  virtual bool applicable(Instance::Action& a);
  virtual SeqProgState* apply(Instance::Action& a);

 public:
  RelaxedSeqProgState(Instance& i, const index_set& x, Heuristic& h,
		      const ACF& c, const index_set& s)
    : SeqProgState(i, h, c, s), x_atoms(x)
  {
    compute_x_actions();
    closure();
  };

  RelaxedSeqProgState(Instance& i, const index_set& x, Heuristic& h,
		      const ACF& c, const bool_vec& g)
    : SeqProgState(i, h, c, g), x_atoms(x)
  {
    compute_x_actions();
    closure();
  };

  RelaxedSeqProgState(const RelaxedSeqProgState& s);
  virtual ~RelaxedSeqProgState() { };

  virtual State* new_state(const index_set& s, State* p);
  virtual State* new_state(const bool_vec& s, State* p);
  virtual State* copy();
};

class FwdUnitHeuristic : public Heuristic
{
 public:
  FwdUnitHeuristic(Instance& ins) : Heuristic(ins) { };
  virtual ~FwdUnitHeuristic() { };

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
};


// FF heuristic can also be used as a ReducedActionSet object
// (returns "helpful" actions)

class ForwardFF : public Heuristic, public ReducedActionSet
{
  const ACF& cost;
  index_set goals;
  CostTable* table;
  weighted_vec<index_type, NTYPE> queue;
  bool_vec props;
  bool_vec acts;

 public:
  ForwardFF(Instance& i, const ACF& c, const index_set& g, Stopwatch& s);
  virtual ~ForwardFF();

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);

  // compute reduced action set (i.e., helpful actions)
  virtual void compute(const bool_vec& s);
};

END_HSPS_NAMESPACE

#endif
