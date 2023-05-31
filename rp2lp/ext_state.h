#ifndef EXTENDED_STATE_H
#define EXTENDED_STATE_H

#include "config.h"
#include "exec.h"
#include "forward.h"
#include "base.h"

BEGIN_HSPS_NAMESPACE

class ExtendedExecState : public ExecState {
 protected:

  // extended state hooks
  virtual bool check_extended_applicability
    (Instance::Action& a, ExecErrorSet* errors, index_type step);
  virtual void apply_extended_action_start_effects
    (Instance::Action& a, exec_act& e, bool_vec& s);
  virtual void apply_extended_action_end_effects
    (Instance::Action& a, exec_act& e, bool_vec& s);
  virtual void save_extended_state_info(exec_act& e);
  virtual void free_extended_state_info(exec_act& e);

 public:
  ExtendedExecState(Instance& i) : ExecState(i) { };
  ExtendedExecState(Instance& i, index_set g) : ExecState(i, g) { };
  ExtendedExecState(Instance& i, const bool_vec& g) : ExecState(i, g) { };
  ExtendedExecState(const ExtendedExecState& s) : ExecState(s) { };
  virtual ~ExtendedExecState() { };

  virtual State* new_state(index_set& s);
  virtual State* copy();
};


class ExtendedSeqProgState : public SeqProgState {
 protected:
  virtual bool applicable(Instance::Action& a);
  virtual SeqProgState* apply(Instance::Action& a);

 public:
  ExtendedSeqProgState(Instance& i, Heuristic& h, const ACF& c)
    : SeqProgState(i, h, c) { };
  ExtendedSeqProgState(Instance& i, Heuristic& h, const ACF& c,
		       const index_set& s)
    : SeqProgState(i, h, c, s) { };
  ExtendedSeqProgState(Instance& i, Heuristic& h, const ACF& c,
			const bool_vec& s)
    : SeqProgState(i, h, c, s) { };

  ExtendedSeqProgState(const ExtendedSeqProgState& s)
    : SeqProgState(s) { };

  virtual State* new_state(const index_set& s, State* p);
  virtual State* new_state(const bool_vec& s, State* p);
  virtual State* copy();
};

END_HSPS_NAMESPACE

#endif
