
#include "forward.h"

BEGIN_HSPS_NAMESPACE

// #define TRACE_PRINT_LOTS

bool SeqProgState::separate_update_actions = false;

SeqProgState::SeqProgState
(Instance& i, Heuristic& h, const ACF& c)
  : AtomSetState(i, h), cost(c), act(no_such_index), res(0)
{
  if (separate_update_actions) {
    apply_update_actions();
  }
}

SeqProgState::SeqProgState
(Instance& i, Heuristic& h, const ACF& c, const index_set& s)
  : AtomSetState(i, h, s), cost(c), act(no_such_index), res(0)
{
  if (separate_update_actions) {
    apply_update_actions();
  }
}

SeqProgState::SeqProgState
(Instance& i, Heuristic& h, const ACF& c, const bool_vec& s)
  : AtomSetState(i, h, s), cost(c), act(no_such_index), res(0)
{
  if (separate_update_actions) {
    apply_update_actions();
  }
}

SeqProgState::SeqProgState
(Instance& i, Heuristic& h, const ACF& c, BasicResourceState* r)
  : AtomSetState(i, h), cost(c), act(no_such_index), res(r)
{
  if (separate_update_actions) {
    apply_update_actions();
  }
}

SeqProgState::SeqProgState
(Instance& i, Heuristic& h, const ACF& c, const index_set& s,
 BasicResourceState* r)
  : AtomSetState(i, h, s), cost(c), act(no_such_index), res(r)
{
  if (separate_update_actions) {
    apply_update_actions();
  }
}

SeqProgState::SeqProgState
(Instance& i, Heuristic& h, const ACF& c, const bool_vec& s,
 BasicResourceState* r)
  : AtomSetState(i, h, s), cost(c), act(no_such_index), res(r)
{
  if (separate_update_actions) {
    apply_update_actions();
  }
}

SeqProgState::SeqProgState(const SeqProgState& s)
  : AtomSetState(s), cost(s.cost), act(s.act), res(0)
{
  if (s.res) res = s.res->copy();
}

SeqProgState::~SeqProgState()
{
  if (res) delete res;
}

bool SeqProgState::applicable(Instance::Action& a)
{
  bool app = true;
  for (index_type i = 0; (i < a.pre.length()) && app; i++)
    if (!set[a.pre[i]]) app = false;
  if (res) {
    if (app) {
      if (!res->applicable(a)) app = false;
    }
  }
  return app;
}

SeqProgState* SeqProgState::apply(Instance::Action& a)
{
  SeqProgState* s = (SeqProgState*)copy();
  for (index_type k = 0; k < a.add.length(); k++) {
    if (!s->set[a.add[k]]) s->size += 1;
    s->set[a.add[k]] = true;
  }
  for (index_type k = 0; k < a.del.length(); k++) {
    if (s->set[a.del[k]]) s->size -= 1;
    s->set[a.del[k]] = false;
  }
  if (separate_update_actions) {
    s->apply_update_actions();
  }
  s->State::set_predecessor(this);
  s->act = a.index;
  s->eval();
  if (s->res) {
    s->res->apply(a);
  }
  return s;
}

void SeqProgState::apply_update_actions()
{
#ifdef TRACE_PRINT_LOTS
  std::cerr << "state BEFORE applying update actions: " << *this << std::endl;
#endif
  bool done = false;
  while (!done) {
    done = true;
    for (index_type k = 0; k < instance.n_actions(); k++) {
      Instance::Action& a = instance.actions[k];
      if (is_update(a) && applicable(a)) {
	for (index_type i = 0; i < instance.actions[k].add.length(); i++) {
	  if (!set[instance.actions[k].add[i]]) done = false;
	  set[instance.actions[k].add[i]] = true;
	}
      }
    }
  }
  count();
#ifdef TRACE_PRINT_LOTS
  std::cerr << "state AFTER applying update actions: " << *this << std::endl;
#endif
}

void SeqProgState::set_predecessor(State* p)
{
  if (p == 0) {
    act = no_such_index;
    pre = 0;
    return;
  }
  SeqProgState* s = (SeqProgState*)p;
  for (index_type k = 0; k < instance.n_actions(); k++) {
    Instance::Action& a = instance.actions[k];
    if (a.sel && s->applicable(a)) {
      SeqProgState* new_s = s->apply(a);
      int d = compare(*new_s);
      delete new_s;
      if (d == 0) {
	act = k;
	pre = s;
	return;
      }
    }
  }
  std::cerr << "error in SeqProgState::set_predecessor: state " << *this
	    << " can not be a successor of " << *s << std::endl;
  exit(255);
}

bool SeqProgState::is_final()
{
  for (index_type k = 0; k < instance.goal_atoms.length(); k++)
    if (!set[instance.goal_atoms[k]]) return false;
  return true;
}

Transition* SeqProgState::transition()
{
  if (act != no_such_index) {
    return new SeqProgTrans(act, cost(act));
  }
  else {
    return 0;
  }
}

NTYPE SeqProgState::delta_cost() {
  if (act == no_such_index)
    return 0;
  else
    return cost(act);
}

NTYPE SeqProgState::expand(Search& s, NTYPE bound) {
  NTYPE c_min = POS_INF;
  for (index_type k = 0; k < instance.n_actions(); k++) {
    Instance::Action& a = instance.actions[k];
    if (a.sel && (!separate_update_actions || !is_update(a))) {
      if (applicable(a)) {
#ifdef TRACE_PRINT_LOTS
	std::cerr << depth() << ". action " << a.name << " is applicable"
		  << std::endl;
#endif
	SeqProgState* new_s = apply(a);
#ifdef TRACE_PRINT_LOTS
	std::cerr << depth() << ". successor state: " << *new_s << std::endl;
#endif
	if (FINITE(new_s->est) && (cost(a.index) + new_s->est <= bound)) {
	  NTYPE c_new =
	    cost(a.index) + s.new_state(*new_s, bound - cost(a.index));
	  if (s.solved()) {
	    if (size > max_set_size_encountered)
	      max_set_size_encountered = size;
	  }
	  if (s.done()) {
	    delete new_s;
	    return c_new;
	  }
	  else {
	    c_min = MIN(c_min, c_new);
	  }
	}
	else {
	  c_min = MIN(c_min, cost(a.index) + new_s->est);
	}
	delete new_s;
      }
    }
#ifdef TRACE_PRINT_LOTS
    else {
      std::cerr << depth() << ". action " << a.name
		<< " is non-selectable | update" << std::endl;
    }
#endif
  }
  return c_min;
}

void SeqProgState::store(NTYPE cost, bool opt) {
  if (res) {
    if (!res->is_root()) return;
  }
  heuristic.store(set, cost, opt);
}

void SeqProgState::insert(Plan& p) {
  if (act != no_such_index) {
    p.insert(act);
    for (index_type k = 0; k < instance.n_atoms(); k++)
      if (set[k] && !instance.actions[act].add.contains(k))
	p.protect(k);
    p.advance(instance.actions[act].dur);
  }
}

void SeqProgState::insert_path(Plan& p)
{
  if (predecessor()) {
    predecessor()->insert_path(p);
  }
  insert(p);
}

void SeqProgState::write_plan(std::ostream& s)
{
  if (act != no_such_index) {
    s << instance.actions[act].name;
  }
}

State* SeqProgState::copy() {
  return new SeqProgState(*this);
}

State* SeqProgState::new_state(const index_set& s, State* p) {
  SeqProgState* new_s =
    (res ? new SeqProgState(instance, heuristic, cost, s, res->new_state())
         : new SeqProgState(instance, heuristic, cost, s));
  new_s->State::set_predecessor(p);
  return new_s;
}

State* SeqProgState::new_state(const bool_vec& s, State* p) {
  SeqProgState* new_s =
    (res ? new SeqProgState(instance, heuristic, cost, s, res->new_state())
         : new SeqProgState(instance, heuristic, cost, s));
  new_s->State::set_predecessor(p);
  return new_s;
}

int SeqProgState::compare(const State& s)
{
  int c = AtomSetState::compare(s);
  if (c != 0) return c;
  SeqProgState& ss = (SeqProgState&)s;
  if (res) {
    if (ss.res)
      return res->compare(*ss.res);
    else
      return -1;
  }
  else {
    if (ss.res) return 1;
  }
  return 0;
}

index_type SeqProgState::hash()
{
  index_type h = AtomSetState::hash();
  if (res) h = (h + res->hash());
  return h;
}

void SeqProgState::write(std::ostream& s) {
  instance.write_atom_set(s, set);
  if (res) {
    s << " (";
    res->write(s);
    s << ")";
  }
}

NTYPE SeqCProgState::expand(Search& s, NTYPE bound) {
  NTYPE c_min = POS_INF;
  for (index_type k = 0; k < instance.n_actions(); k++) {
    Instance::Action& a = instance.actions[k];
    if (a.sel && (!separate_update_actions || !is_update(a))) {
      bool app = applicable(a);
      if (app) {
#ifdef TRACE_PRINT_LOTS
	std::cerr << depth() << ". action " << instance.actions[k].name
		  << " is applicable" << std::endl;
#endif
	if ((act != no_such_index) && (act < a.index) &&
	    instance.commutative(act, a.index)) app = false;
#ifdef TRACE_PRINT_LOTS
	if (!app) {
	  std::cerr << depth() << ". commutativity cut: action "
		    << instance.actions[k].name << ", predecessor "
		    << instance.actions[act].name << std::endl;
	}
#endif
      }
      if (app) {
	SeqCProgState* new_s = (SeqCProgState*)apply(a);
#ifdef TRACE_PRINT_LOTS
	std::cerr << depth() << ". successor state: " << *new_s << std::endl;
#endif
	if (FINITE(new_s->est) && (cost(a.index) + new_s->est <= bound)) {
	  NTYPE c_new =
	    cost(a.index) + s.new_state(*new_s, bound - cost(a.index));
	  if (s.solved()) {
	    if (size > max_set_size_encountered)
	      max_set_size_encountered = size;
	  }
	  if (s.done()) {
	    delete new_s;
	    return c_new;
	  }
	  else {
	    c_min = MIN(c_min, c_new);
	  }
	}
	else {
	  c_min = MIN(c_min, cost(a.index) + new_s->est);
	}
	delete new_s;
      }
    }
  }
  return c_min;
}

State* SeqCProgState::copy() {
  return new SeqCProgState(*this);
}

State* SeqCProgState::new_state(const index_set& s, State* p) {
  SeqCProgState* new_s =
    (res ? new SeqCProgState(instance, heuristic, cost, s, res)
         : new SeqCProgState(instance, heuristic, cost, s));
  new_s->State::set_predecessor(p);
  return new_s;
}

State* SeqCProgState::new_state(const bool_vec& s, State* p) {
  SeqCProgState* new_s =
    (res ? new SeqCProgState(instance, heuristic, cost, s, res)
         : new SeqCProgState(instance, heuristic, cost, s));
  new_s->State::set_predecessor(p);
  return new_s;
}

int SeqCProgState::compare(const State& s)
{
  SeqCProgState& ss = (SeqCProgState&)s;
  if (act < ss.act) return -1;
  if (act > ss.act) return 1;
  int c = AtomSetState::compare(s);
  if (c != 0) return c;
  if (res) {
    if (ss.res)
      return res->compare(*ss.res);
    else
      return -1;
  }
  else {
    if (ss.res) return 1;
  }
  return 0;
}

index_type SeqCProgState::hash()
{
  index_type h = (AtomSetState::hash() + act);
  if (res) h = (h + res->hash());
  return h;
}

void SeqCProgState::write(std::ostream& s) {
  if (act != no_such_index)
    s << "(" << instance.actions[act].name << ") ";
  else
    s << "() ";
  instance.write_atom_set(s, set);
  if (res) {
    s << " (";
    res->write(s);
    s << ")";
  }
}


NTYPE RedSeqProgState::expand(Search& s, NTYPE bound)
{
  NTYPE c_min = POS_INF;
  ras.compute(set);
  for (index_type k = 0; k < instance.n_actions(); k++) {
    Instance::Action& a = instance.actions[k];
    if (a.sel && (ras[k] || (d_max > 0)))
      if (applicable(a)) {
	RedSeqProgState* new_s = (RedSeqProgState*)apply(a);
	if (!ras[k]) new_s->d_max = (d_max - 1);
	if (FINITE(new_s->est) && (cost(a.index) + new_s->est <= bound)) {
	  NTYPE c_new =
	    cost(a.index) + s.new_state(*new_s, bound - cost(a.index));
	  if (s.solved()) {
	    if (size > max_set_size_encountered)
	      max_set_size_encountered = size;
	  }
	  if (s.done()) {
	    delete new_s;
	    return c_new;
	  }
	  else {
	    c_min = MIN(c_min, c_new);
	  }
	}
	else {
	  c_min = MIN(c_min, cost(a.index) + new_s->est);
	}
	delete new_s;
      }
  }
  return c_min;
}

State* RedSeqProgState::new_state(const index_set& s, State* p)
{
  RedSeqProgState* new_s =
    new RedSeqProgState(instance, heuristic, cost, s, ras);
  new_s->State::set_predecessor(p);
  return new_s;
}

State* RedSeqProgState::new_state(const bool_vec& s, State* p)
{
  RedSeqProgState* new_s =
    new RedSeqProgState(instance, heuristic, cost, s, ras);
  new_s->State::set_predecessor(p);
  return new_s;
}

State* RedSeqProgState::copy()
{
  return new RedSeqProgState(*this);
}


NTYPE RestrictedSeqProgState::expand(Search& s, NTYPE bound)
{
  NTYPE c_min = POS_INF;
  for (index_type k = 0; k < instance.n_actions(); k++) if (instance.actions[k].sel) {
    Instance::Action& a = instance.actions[k];
    if (applicable(a)) {
      RestrictedSeqProgState* new_s = (RestrictedSeqProgState*)apply(a);
      if (FINITE(new_s->est) && (cost(a.index) + new_s->est <= bound)) {
	bool allowed = true;
	for (index_type i = 0; i < a.del.length(); i++)
	  if (p_atoms[a.del[i]]) {
	    allowed = false;
	    f_atoms[a.del[i]] = true;
	  }
	if (allowed) {
	  NTYPE c_new =
	    cost(a.index) + s.new_state(*new_s, bound - cost(a.index));
	  if (s.done()) {
	    delete new_s;
	    return c_new;
	  }
	  else {
	    c_min = MIN(c_min, c_new);
	  }
	}
      }
      else {
	c_min = MIN(c_min, cost(a.index) + new_s->est);
      }
      delete new_s;
    }
  }
  return c_min;
}

State* RestrictedSeqProgState::new_state(const index_set& s, State* p)
{
  RestrictedSeqProgState* new_s =
    (res ? new RestrictedSeqProgState(instance, p_atoms, f_atoms,
				      heuristic, cost, s, res)
         : new RestrictedSeqProgState(instance, p_atoms, f_atoms,
				      heuristic, cost, s));
  new_s->State::set_predecessor(p);
  return new_s;
}

State* RestrictedSeqProgState::new_state(const bool_vec& s, State* p)
{
  RestrictedSeqProgState* new_s =
    (res ? new RestrictedSeqProgState(instance, p_atoms, f_atoms,
				      heuristic, cost, s, res)
         : new RestrictedSeqProgState(instance, p_atoms, f_atoms,
				      heuristic, cost, s));
  new_s->State::set_predecessor(p);
  return new_s;
}

State* RestrictedSeqProgState::copy()
{
  return new RestrictedSeqProgState(*this);
}


RelaxedSeqProgState::RelaxedSeqProgState(const RelaxedSeqProgState& s)
  : SeqProgState(s), x_atoms(s.x_atoms), x_action(s.x_action)
{
  // done
}

void RelaxedSeqProgState::RelaxedSeqProgState::compute_x_actions()
{
  x_action.assign_value(false, instance.n_actions());
  for (index_type k = 0; k < instance.n_actions(); k++) {
    if (instance.actions[k].del.first_common_element(x_atoms) != no_such_index)
      x_action[k] = true;
  }
}

void RelaxedSeqProgState::closure()
{
  bool done = false;
  while (!done) {
    done = true;
    for (index_type k = 0; k < instance.n_actions(); k++) if (!x_action[k])
      if (SeqProgState::applicable(instance.actions[k]))
	for (index_type i = 0; i < instance.actions[k].add.length(); i++) {
	  if (!set[instance.actions[k].add[i]]) done = false;
	  set[instance.actions[k].add[i]] = true;
	}
  }
  count();
}

bool RelaxedSeqProgState::applicable(Instance::Action& a)
{
  if (!x_action[a.index]) return false;
  return SeqProgState::applicable(a);
}

SeqProgState* RelaxedSeqProgState::apply(Instance::Action& a)
{
  RelaxedSeqProgState* s = new RelaxedSeqProgState(*this);
  for (index_type k = 0; k < a.add.length(); k++)
    s->set[a.add[k]] = true;
  for (index_type k = 0; k < a.del.length(); k++)
    s->set[a.del[k]] = false;
  s->closure();
  s->State::set_predecessor(this);
  s->act = a.index;
  s->eval();
  return s;
}

State* RelaxedSeqProgState::new_state(const index_set& s, State* p)
{
  RelaxedSeqProgState* new_s =
    new RelaxedSeqProgState(instance, x_atoms, heuristic, cost, s);
  new_s->State::set_predecessor(p);
  return new_s;
}

State* RelaxedSeqProgState::new_state(const bool_vec& s, State* p)
{
  RelaxedSeqProgState* new_s =
    new RelaxedSeqProgState(instance, x_atoms, heuristic, cost, s);
  new_s->State::set_predecessor(p);
  return new_s;
}

State* RelaxedSeqProgState::copy()
{
  return new RelaxedSeqProgState(*this);
}

NTYPE FwdUnitHeuristic::eval(const index_set& s)
{
  if (s.contains(instance.goal_atoms))
    return 0;
  else
    return 1;
}

NTYPE FwdUnitHeuristic::eval(const bool_vec& s)
{
  if (s.contains(instance.goal_atoms))
    return 0;
  else
    return 1;
}

ForwardFF::ForwardFF
(Instance& i, const ACF& c, const index_set& g, Stopwatch& s)
  : Heuristic(i), ReducedActionSet(i), cost(c), goals(g), table(0)
{
  table = new CostTable(i, s);
}

ForwardFF::~ForwardFF()
{
  delete table;
}

NTYPE ForwardFF::eval(const index_set& s)
{
  bool_vec s1(s, Heuristic::instance.n_atoms());
  return eval(s1);
}

NTYPE ForwardFF::eval(const bool_vec& s)
{
  table->compute_H1(cost, s);
  props.assign_value(false, Heuristic::instance.n_atoms());
  queue.clear();
  for (index_type k = 0; k < goals.size(); k++) {
    NTYPE v = table->eval(goals[k]);
    if (INFINITE(v)) return POS_INF;
    queue.insert_decreasing(goals[k], v);
    props[goals[k]] = true;
  }
  acts.assign_value(false, Heuristic::instance.n_actions());
  while (queue.length() > 0) {
    index_type p = queue[0].value;
    queue.remove(0);
    if (table->eval(p) > 0) {
      NTYPE c_min = POS_INF;
      index_type a_min = no_such_index;
      for (index_type i = 0; i < Heuristic::instance.atoms[p].add_by.size(); i++) {
	index_type a = Heuristic::instance.atoms[p].add_by[i];
	if (acts[a]) {
	  c_min = 0;
	  a_min = a;
	}
	else {
	  NTYPE v = table->eval(Heuristic::instance.actions[a].pre);
	  if (v < c_min) {
	    c_min = v;
	    a_min = a;
	  }
	}
      }
      assert(FINITE(c_min));
      assert(a_min < Heuristic::instance.n_actions());
      acts[a_min] = true;
      for (index_type i = 0; i < Heuristic::instance.actions[a_min].pre.size(); i++) {
	NTYPE v = table->eval(Heuristic::instance.actions[a_min].pre[i]);
	if (!props[Heuristic::instance.actions[a_min].pre[i]]) {
	  queue.insert_decreasing(Heuristic::instance.actions[a_min].pre[i], v);
	  props[Heuristic::instance.actions[a_min].pre[i]] = true;
	}
      }
    }
  }
  return acts.count(true);
}

void ForwardFF::compute(const bool_vec& s)
{
  NTYPE v = eval(s);
  for (index_type k = 0; k < Heuristic::instance.n_actions(); k++) {
    (*this)[k] = false;
    if (acts[k])
      if (s.contains(Heuristic::instance.actions[k].pre))
	(*this)[k] = true;
  }
}

END_HSPS_NAMESPACE
