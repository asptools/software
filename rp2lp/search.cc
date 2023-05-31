
#include "search.h"

BEGIN_HSPS_NAMESPACE

Transition::~Transition()
{
  // done
}

void Transition::write(::std::ostream& s) const
{
  s << "<TRANSITION>";
}

void Transition::write_path(::std::ostream& s) const
{
  s << "<TRANSITION PATH>";
}

void ProgressionPath::insert_path(Plan& p)
{
  if (predecessor()) {
    predecessor()->insert_path(p);
  }
  insert(p);
}

void ProgressionPath::write_path(::std::ostream& s) const
{
  if (predecessor()) {
    predecessor()->write_path(s);
    s << ", ";
  }
  write(s);
}

void RegressionPath::insert_path(Plan& p)
{
  Transition* t = this;
  while (t != 0) {
    t->insert(p);
    t = t->predecessor();
  }
}

void RegressionPath::write_path(::std::ostream& s) const
{
  const Transition* t = this;
  while (t != 0) {
    if (t != this) s << ", ";
    t->write(s);
    t = t->predecessor();
  }
}

int SequentialTransition::compare(const SequentialTransition& t)
{
  // std::cerr << "comparing ";
  // write(std::cerr);
  // std::cerr << " and ";
  // st.write(std::cerr);
  // std::cerr << std::endl;
  if (act < t.act)
    return -1;
  else if (act > t.act)
    return 1;
  else if (delta < t.delta)
    return -1;
  else if (delta > t.delta)
    return 1;
  else
    return 0;
}

void SequentialTransition::insert(Plan& p)
{
  p.insert(act);
  p.advance(delta);
}

void SequentialTransition::write(::std::ostream& s) const
{
  s << act << ":" << delta;
}

int ParallelTransition::compare(const ParallelTransition& t)
{
  //const ParallelTransition& pt = (const ParallelTransition&)t;
  if (acts < t.acts)
    return -1;
  else if (acts > t.acts)
    return 1;
  else if (noops < t.noops)
    return -1;
  else if (noops > t.noops)
    return 1;
  else if (delta < t.delta)
    return -1;
  else if (delta > t.delta)
    return 1;
  else
    return 0;
}

void ParallelTransition::insert(Plan& p)
{
  for (index_type k = 0; k < acts.length(); k++)
    p.insert(acts[k]);
  for (index_type k = 0; k < noops.length(); k++)
    p.protect(noops[k]);
  p.advance(delta);
}

void ParallelTransition::write(::std::ostream& s) const
{
  s << acts << ":" << noops << ":" << delta;
}

State::~State()
{
  // done
}

const State* State::predecessor() const
{
  return pre;
}

State* State::predecessor()
{
  return pre;
}

void State::set_predecessor(State* p)
{
  pre = p;
}

Transition* State::transition_path()
{
  Transition* t0 = transition();
  Transition* t = t0;
  State* s = this;
  while (s->predecessor()) {
    s = s->predecessor();
    Transition* t1 = s->transition();
    t->set_predecessor(t1);
    t = t1;
  }
  return t0;
}

NTYPE State::acc_cost()
{
  if (predecessor()) {
    return delta_cost() + predecessor()->acc_cost();
  }
  else {
    return 0;
  }
}

index_type State::depth() const
{
  if (predecessor()) {
    return 1 + predecessor()->depth();
  }
  else {
    return 0;
  }
}

void State::store(NTYPE cost, bool opt)
{
  std::cerr << "error: store() operation not supported by state" << std::endl;
  exit(255);
}

void State::reevaluate()
{
  // default implement: does nothing
}

void State::insert(Plan& p)
{
  // default implement: does nothing
}

void State::insert_path(Plan& p)
{
  Transition* tp = transition_path();
  if (tp != NULL)
    tp->insert_path(p);
}

void State::write_plan(::std::ostream& s)
{
  // default implement: does nothing
}

void State::write_eval(::std::ostream& s, const char* p, bool e)
{
  if (p) s << p << " ";
  s << est_cost();
  if (e) s << ::std::endl;
}

State* State::copy_path()
{
  State* s = copy();
  if (predecessor()) {
    s->set_predecessor(predecessor()->copy_path());
  }
  return s;
}

void State::delete_path()
{
  if (predecessor()) {
    predecessor()->delete_path();
  }
  delete this;
}

int State::compare_path(const State* s)
{
  if (s == 0) return 1;
  int c = compare(*s);
  if (c == 0) {
    if (predecessor() == 0) return -1;
    return predecessor()->compare_path(s->predecessor());
  }
  return c;
}

void State::write_path(::std::ostream& s)
{
  State* p = predecessor();
  if (p) {
    p->write_path(s);
  }
  write(s);
  s << ::std::endl;
}

void State::write_path_as_graph(::std::ostream& s)
{
  State* p = predecessor();
  if (p) {
    p->write_path_as_graph(s);
  }
  s << "S" << depth() << " [label=\"";
  write(s);
  if (is_final()) {
    s << ",style=\"bold\"";
  }
  s << "\"];" << ::std::endl;
  if (p) {
    s << "S" << depth() << " -> S" << p->depth()
      << " [label=\"";
    Transition* t = p->transition();
    if (t) {
      t->write(s);
      delete t;
    }
    else {
      s << "?";
    }
    s << "\"];" << ::std::endl;
  }
}

void ProgressionState::insert_path(Plan& p)
{
  if (predecessor()) {
    predecessor()->insert_path(p);
  }
  insert(p);
}

void RegressionState::insert_path(Plan& p)
{
  State* sp = this;
  while (sp != 0) {
    sp->insert(p);
    sp = sp->predecessor();
  }
}

PlanTrait::~PlanTrait()
{
  // done
}

const PlanTrait* PlanTrait::cast_to(const char* n) const
{
  return 0;
}

Plan::~Plan()
{
  // done
}

void Plan::output(Plan& to, bool mark_end)
{
  // default: do nothing, excep call end() on 'to' if requested
  if (mark_end) {
    to.end();
  }
}

void Plan::set_name(const Name* n)
{
  // default: ignore
}

void Plan::set_optimal(bool o)
{
  // default: ignore
}

void Plan::add_trait(PlanTrait* t)
{
  delete t;
}

Search::~Search()
{
  // done
}

NoSearch::~NoSearch()
{
  // done
}

void NoSearch::reset()
{
  _solved = false;
}

NTYPE NoSearch::new_state(State& s, NTYPE bound)
{
  if (s.is_final()) _solved = true;
  return s.est_cost();
}

bool NoSearch::solved() const
{
  return _solved;
}

bool NoSearch::optimal() const
{
  return false;
}

bool NoSearch::done() const
{
  return false;
}

// Transitions::Transitions()
//   : state_vec((State*)0, 0), target_state(0)
// {
//   // done
// }
// 
// Transitions::Transitions(State* from, State* to, NTYPE d)
//   : state_vec((State*)0, 0), target_state(0)
// {
//   find(from, to, d, true);
// }
// 
// Transitions::~Transitions()
// {
//   clear();
// }
// 
// void Transitions::clear()
// {
//   for (index_type k = 0; k < length(); k++)
//     delete (*this)[k];
//   state_vec::clear();
// }
// 
// bool Transitions::find(State* from, State* to, NTYPE d, bool x)
// {
//   target_state = to;
//   delta_bound = d;
//   bound_is_exact = x;
//   from->expand(*this, delta_bound + to->est_cost());
//   target_state = 0;
//   return (length() > 0);
// }
// 
// NTYPE Transitions::new_state(State& s, NTYPE bound)
// {
//   if (!target_state) {
//     ::std::cerr << "error: Transitions::new_state called with state = ";
//     s.write(::std::cerr);
//     ::std::cerr << " and target_state == nil" << ::std::endl;
//     exit(255);
//   }
//   if (((s.delta_cost() == delta_bound) || !bound_is_exact) &&
//       (target_state->compare(s) == 0)) {
//     append(s.copy());
//   }
//   return s.est_cost();
// }
// 
// bool Transitions::solved() const
// {
//   return (length() > 0);
// }
// 
// bool Transitions::optimal() const
// {
//   return false;
// }
// 
// bool Transitions::done() const
// {
//   return false;
// }

StateFactory::~StateFactory()
{
  // done
}

PlanSet::~PlanSet()
{
  // done
}

void PlanSet::output(PlanSet& to)
{
  // done
}

void PlanSet::output(PlanSet& to, const bool_vec& s)
{
  // done
}

END_HSPS_NAMESPACE
