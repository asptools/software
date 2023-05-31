
#include "ext_state.h"
#include "base.h"

BEGIN_HSPS_NAMESPACE

bool ExtendedExecState::check_extended_applicability
(Instance::Action& a, ExecErrorSet* errors, index_type step)
{
  bool pre_ok = true;
  assert(a.src);
  ptr_pair* p = (ptr_pair*)a.src;
  PDDL_Base::ActionSymbol* act = (PDDL_Base::ActionSymbol*)p->first;
  // if (act->n_dis_pre > 0) {
  //   ptr_table* ins = (ptr_table*)p->second;
  //   ptr_table::key_vec* args = ins->key_sequence();
  //   for (index_type k = 0; k < act->n_param; k++)
  //     act->param[k]->value = (PDDL_Base::Symbol*)((*args)[k + 1]);
  //   delete args;
  //   for (index_type k = 0; k < act->n_dset_pre; k++) {
  //     index_set s;
  //     act->dset_pre[k]->instantiate(instance, s);
  //     bool d_sat = false;
  //     for (index_type i = 0; (i < s.length()) && !d_sat; i++) {
  // 	bool a_sat = atoms[s[i]];
  // 	for (index_type j = 0; (j < instance.n_actions()) && !a_sat; j++) {
  // 	  if (actions[j] && (remains[j] <= dt)) {
  // 	    a_sat = instance.actions[j].add.contains(s[i]);
  // 	  }
  // 	}
  // 	if (a_sat) d_sat = true;
  //     }
  //     if (!d_sat) {
  // 	if (trace_level > 0) {
  // 	  ::std::cerr << "action " << a.name << ": precondition ";
  // 	  act->dset_pre[k]->print(::std::cerr);
  // 	  ::std::cerr << " not satisfied" << ::std::endl;
  // 	}
  // 	pre_ok = false;
  //     }
  //   }
  // }
  return pre_ok;
}

void ExtendedExecState::apply_extended_action_start_effects
(Instance::Action& a, exec_act& e, bool_vec& s)
{
  assert(a.src);
  ExtendedExecState* start_state = (ExtendedExecState*)e.esi;
  assert(start_state);
  PDDL_Base::ActionSymbol* act =
    PDDL_Base::src_action_symbol((ptr_pair*)a.src);
  if (act->cond_eff.length() == 0) return;
  for (index_type k = 0; k < act->cond_eff.length(); k++) {
    rule_set pe;
    rule_set ne;
    act->cond_eff[k]->instantiate_conditional(instance, pe, ne);
    for (index_type i = 0; i < ne.length(); i++) {
      bool app = true;
      for (index_type j = 0; (j < ne[i].antecedent.length()) && app; j++)
	if (!start_state->atoms[ne[i].antecedent[j]]) app = false;
      if (app) {
	s[ne[i].consequent] = false;
      }
    }
  }
}

void ExtendedExecState::apply_extended_action_end_effects
(Instance::Action& a, exec_act& e, bool_vec& s)
{
  assert(a.src);
  ExtendedExecState* start_state = (ExtendedExecState*)e.esi;
  assert(start_state);
  PDDL_Base::ActionSymbol* act =
    PDDL_Base::src_action_symbol((ptr_pair*)a.src);
  if (act->cond_eff.length() == 0) return;
  for (index_type k = 0; k < act->cond_eff.length(); k++) {
    rule_set pe;
    rule_set ne;
    act->cond_eff[k]->instantiate_conditional(instance, pe, ne);
    for (index_type i = 0; i < pe.length(); i++) {
      bool app = true;
      for (index_type j = 0; (j < pe[i].antecedent.length()) && app; j++)
	if (!start_state->atoms[pe[i].antecedent[j]]) app = false;
      if (app) {
	s[pe[i].consequent] = true;
      }
    }
  }
}

void ExtendedExecState::save_extended_state_info(exec_act& e)
{
  e.esi = copy();
}

void ExtendedExecState::free_extended_state_info(exec_act& e)
{
  ExtendedExecState* s = (ExtendedExecState*)e.esi;
  delete s;
  e.esi = 0;
}

State* ExtendedExecState::new_state(index_set& s)
{
  return new ExtendedExecState(instance, s);
}

State* ExtendedExecState::copy()
{
  return new ExtendedExecState(*this);
}


bool ExtendedSeqProgState::applicable(Instance::Action& a)
{
  if (!SeqProgState::applicable(a)) return false;
  if (PDDL_Base::compile_away_disjunctive_preconditions) return true;
  ptr_pair* p = (ptr_pair*)a.src;
  PDDL_Base::ActionSymbol* act = (PDDL_Base::ActionSymbol*)p->first;
//   if (act->n_dset_pre > 0) {
//     ptr_table* ins = (ptr_table*)p->second;
//     ptr_table::key_vec* args = ins->key_sequence();
//     for (index_type k = 0; k < act->n_param; k++)
//       act->param[k]->value = (PDDL_Base::Symbol*)((*args)[k + 1]);
//     delete args;
//     for (index_type k = 0; k < act->n_dset_pre; k++) {
//       index_set s;
//       act->dset_pre[k]->instantiate(instance, s);
//       bool is_sat = false;
//       for (index_type i = 0; (i < s.length()) && !is_sat; i++)
// 	if (set[s[i]]) is_sat = true;
//       if (!is_sat) return false;
//     }
//   }
  return true;
}

SeqProgState* ExtendedSeqProgState::apply(Instance::Action& a)
{
  ExtendedSeqProgState* s = (ExtendedSeqProgState*)copy();
  for (index_type k = 0; k < a.add.length(); k++) {
    s->set[a.add[k]] = true;
  }
  for (index_type k = 0; k < a.del.length(); k++) {
    s->set[a.del[k]] = false;
  }
  if (!PDDL_Base::compile_away_conditional_effects) {
    ptr_pair* p = (ptr_pair*)a.src;
    PDDL_Base::ActionSymbol* act = (PDDL_Base::ActionSymbol*)p->first;
    if (act->cond_eff.length() > 0) {
      ptr_table* ins = (ptr_table*)p->second;
      ptr_table::key_vec* args = ins->key_sequence();
      for (index_type k = 0; k < act->param.length(); k++)
	act->param[k]->value = (PDDL_Base::Symbol*)((*args)[k + 1]);
      delete args;
      for (index_type k = 0; k < act->cond_eff.length(); k++) {
	rule_set ne;
	rule_set pe;
	act->cond_eff[k]->instantiate_conditional(instance, pe, ne);
	for (index_type i = 0; i < ne.length(); i++) {
	  bool app = true;
	  for (index_type j = 0; (j < ne[i].antecedent.length()) && app; j++)
	    if (!set[ne[i].antecedent[j]]) app = false;
	  if (app) {
	    s->set[ne[i].consequent] = false;
	  }
	}
	for (index_type i = 0; i < pe.length(); i++) {
	  bool app = true;
	  for (index_type j = 0; (j < pe[i].antecedent.length()) && app; j++)
	    if (!set[pe[i].antecedent[j]]) app = false;
	  if (app) {
	    s->set[pe[i].consequent] = true;
	  }
	}
      }
    }
  }
  s->count();
  s->State::set_predecessor(this);
  s->act = a.index;
  s->eval();
  return s;
}

State* ExtendedSeqProgState::new_state(const index_set& s, State* p)
{
  ExtendedSeqProgState* new_s =
    new ExtendedSeqProgState(instance, heuristic, cost, s);
  new_s->State::set_predecessor(p);
  return new_s;
}

State* ExtendedSeqProgState::new_state(const bool_vec& s, State* p)
{
  ExtendedSeqProgState* new_s =
    new ExtendedSeqProgState(instance, heuristic, cost, s);
  new_s->State::set_predecessor(p);
  return new_s;
}

State* ExtendedSeqProgState::copy()
{
  return new ExtendedSeqProgState(*this);
}

END_HSPS_NAMESPACE
