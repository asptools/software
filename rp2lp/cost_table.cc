
#include "cost_table.h"
#include "graph.h"
#include "enumerators.h"

#ifdef BUILD_WITH_BOOSTING
#include "idao.h"
#include "collector.h"
#include "plans.h"
#endif

BEGIN_HSPS_NAMESPACE

PairIndexFunction::PairIndexFunction(index_type n)
  : index_vec(0, n)
{
  for (index_type k = 1; k < n; k++)
    (*this)[k] = (*this)[k - 1] + (n - k) + 1;
}

PairIndexFunction::~PairIndexFunction()
{
  // done
}

index_type PairIndexFunction::operator()(index_type i, index_type j) const
{
  sort2(i, j);
  return ((*this)[i] + (j - i));
}

index_type PairIndexFunction::max_in() const
{
  return length() - 1;
}

index_type PairIndexFunction::max_out() const
{
  return (*this)(length() - 1, length() - 1);
}

bool CostNode::eval_set_mark = false;
#ifdef EVAL_EXTRA_STATS
count_type CostNode::eval_rec_count = 0;
#endif

#ifdef BUILD_WITH_BOOSTING
bool CostTableBoost::strong_conflict_detection = true;
bool CostTableBoost::ultra_weak_conflict_detection = false;
bool CostTableBoost::boost_new_entries = true;
count_type CostTableBoost::n_entries_boosted = 0;
count_type CostTableBoost::n_entries_solved = 0;
count_type CostTableBoost::n_entries_discarded = 0;
count_type CostTableBoost::n_entries_created = 0;
count_type CostTableBoost::n_boost_searches = 0;
count_type CostTableBoost::n_boost_searches_with_cd = 0;
#endif

CostNode::CostNode(index_type size)
  : _size(size),
    _depth(0),
    _first(0),
    _store(0),
    _default(0),
    _max(POS_INF),
    _prev(0)
{
  _store = new Value[_size];
  for (index_type k = 0; k < _size; k++) {
    _store[k].val = _default;
    _store[k].next = 0;
  }
}

CostNode::CostNode(index_type size, NTYPE v_default)
  : _size(size),
    _depth(0),
    _first(0),
    _store(0),
    _default(v_default),
    _max(POS_INF),
    _prev(0)
{
  _store = new Value[_size];
  for (index_type k = 0; k < _size; k++) {
    _store[k].val = _default;
    _store[k].next = 0;
  }
}

CostNode::CostNode(CostNode* p, index_type f)
  : _size(p->_size),
    _depth(p->_depth + 1),
    _first(f),
    _store(0),
    _default(p->_default),
    _max(POS_INF),
    _prev(p)
{
  _store = new Value[_size];
  for (index_type k = 0; k < _size; k++) {
    _store[k].val = _default;
    _store[k].next = 0;
  }
}

CostNode::~CostNode()
{
  for (index_type k = 0; k < _size; k++)
    if (_store[k].next) delete _store[k].next;
  delete [] _store;
}

void CostNode::clear()
{
  for (index_type k = 0; k < _size; k++) {
    if (_store[k].next) delete _store[k].next;
    _store[k].clear(_default);
  }
  _max = POS_INF;
}

CostNode::Value* CostNode::find(const index_set& s)
{
  Value* v = 0;
  CostNode* n = this;
  for (index_type k = 0; k < s.length(); k++) {
    if (!n) return 0;
    v = n->val_p(s[k]);
    n = v->next;
  }
  return v;
}

CostNode::Value* CostNode::find(const bool_vec& s)
{
  Value* v = 0;
  CostNode* n = this;
  for (index_type k = 0; k < _size; k++) if (s[k]) {
    if (!n) return 0;
    v = n->val_p(k);
    n = v->next;
  }
  return v;
}

void CostNode::find_all
(const index_set& s, index_set_vec& a)
{
  find_all(s, 0, EMPTYSET, a);
}

void CostNode::find_all
(const index_set& s, index_type i, const index_set& c, index_set_vec& a)
{
  for (index_type k = i; k < s.length(); k++) {
    index_set nc(c);
    nc.insert(s[k]);
    a.append_if_new(nc);
    if (val(s[k]).next) {
      val(s[k]).next->find_all(s, k + 1, nc, a);
    }
  }
}

CostNode::Value& CostNode::insert(const index_set& s)
{
  assert(s.length() > 0);
  CostNode* n = this;
  for (index_type k = 0; k < s.length() - 1; k++) {
    n = n->next_p(s[k]);
  }
  return n->val(s[s.length() - 1]);
}

CostNode::Value& CostNode::insert(const bool_vec& s)
{
  index_type last = no_such_index;
  for (index_type k = _size; (k > 0) && (last == no_such_index); k--)
    if (s[k - 1]) last = k - 1;
  assert(last != no_such_index);
  CostNode* n = this;
  for (index_type k = 0; k < last; k++) if (s[k]) {
    n = n->next_p(k);
  }
  return n->val(last);
}

void CostNode::store(index_type p, NTYPE v, bool opt)
{
  val(p) = Value(v, opt);
}

void CostNode::store(index_type p, NTYPE v)
{
  val(p) = v;
}

void CostNode::set_max(NTYPE v)
{
  if (v > _max) {
    _max = v;
    if (_prev) _prev->set_max(v);
  }
}

void CostNode::store(const index_set& s, NTYPE v, bool opt)
{
  if (s.length() > 0) {
    CostNode* n = this;
    for (index_type i = 0; i < s.length() - 1; i++) n = n->next_p(s[i]);
    n->store(s[s.length() - 1], v, opt);
  }
}

void CostNode::store(const bool_vec& s, NTYPE v, bool opt)
{
  CostNode* p = 0;
  CostNode* n = this;
  index_type l = no_such_index;
  for (index_type k = 0; k < _size; k++) if (s[k]) {
    p = n;
    n = n->next_p(k);
    l = k;
  }
  if (p) p->store(l, v, opt);
}

void CostNode::store(const index_set& s, NTYPE v)
{
  if (s.length() > 0) {
    CostNode* n = this;
    for (index_type i = 0; i < s.length() - 1; i++) n = n->next_p(s[i]);
    n->store(s[s.length() - 1], v);
  }
}

void CostNode::store(const bool_vec& s, NTYPE v)
{
  CostNode* p = 0;
  CostNode* n = this;
  index_type l = no_such_index;
  for (index_type k = 0; k < _size; k++) if (s[k]) {
    p = n;
    n = n->next_p(k);
    l = k;
  }
  if (p) p->store(l, v);
}

NTYPE CostNode::eval(const index_set& s)
{
  return eval(s, 0);
}

NTYPE CostNode::eval(const bool_vec& s)
{
  return eval(s, 0);
}

extended_cost CostNode::extended_eval(const index_set& s)
{
  Value* v = find(s);
  if (v) {
    if (v->opt) {
      return extended_cost(v->val, v->opt);
    }
    else {
      return extended_cost(eval(s, 0), false);
    }
  }
  else {
    return extended_cost(eval(s, 0), false);
  }
}

extended_cost CostNode::extended_eval(const bool_vec& s)
{
  Value* v = find(s);
  if (v) {
    if (v->opt) {
      return extended_cost(v->val, v->opt);
    }
    else {
      return extended_cost(eval(s, 0), false);
    }
  }
  else {
    return extended_cost(eval(s, 0), false);
  }
}

NTYPE CostNode::eval_to_bound(const index_set& s, NTYPE bound)
{
  return eval_to_bound(s, 0, bound);
}

NTYPE CostNode::eval_to_bound(const bool_vec& s, NTYPE bound)
{
  return eval_to_bound(s, 0, bound);
}

NTYPE CostNode::eval(const index_set& s, index_type i)
{
  Heuristic::eval_count += 1;
  NTYPE v_new = 0;
  for (index_type k = i; k < s.length(); k++) {
    if (INFINITE(val(s[k]).val)) return POS_INF;
    v_new = MAX(v_new, val(s[k]).val);
    if (val(s[k]).next) {
#ifdef EVAL_EXTRA_STATS
      eval_rec_count += 1;
#endif
      v_new = MAX(v_new, val(s[k]).next->eval(s, k + 1));
    }
    val(s[k]).mark = true;
  }
  return v_new;
}

NTYPE CostNode::eval(const bool_vec& s, index_type i)
{
  Heuristic::eval_count += 1;
  NTYPE v_new = 0;
  for (index_type k = i; k < _size; k++) if (s[k]) {
    if (INFINITE(val(k).val)) return POS_INF;
    v_new = MAX(v_new, val(k).val);
    if (val(k).next) {
#ifdef EVAL_EXTRA_STATS
      eval_rec_count += 1;
#endif
      v_new = MAX(v_new, val(k).next->eval(s, k + 1));
    }
    val(k).mark = true;
  }
  return v_new;
}

NTYPE CostNode::eval_to_bound(const index_set& s, index_type i, NTYPE bound)
{
  NTYPE v = 0;
  for (index_type k = i; k < s.length(); k++) {
    v = MAX(v, val(s[k]).val);
    val(s[k]).mark = true;
    if (v > bound) return v;
    if (val(s[k]).next)
      v = MAX(v, val(s[k]).next->eval_to_bound(s, k + 1, bound));
  }
  return v;
}

NTYPE CostNode::eval_to_bound(const bool_vec& s, index_type i, NTYPE bound)
{
  NTYPE v = 0;
  for (index_type k = i; k < _size; k++) if (s[k]) {
    v = MAX(v, val(k).val);
    val(k).mark = true;
    if (v > bound) return v;
    if (val(k).next)
      v = MAX(v, val(k).next->eval_to_bound(s, k + 1, bound));
  }
  return v;
}

NTYPE CostNode::incremental_eval(const bool_vec& s, index_type i, index_type i_new)
{
  Heuristic::eval_count += 1;
  NTYPE v = val(i_new).val;
  if (val(i_new).next) {
#ifdef EVAL_EXTRA_STATS
    eval_rec_count += 1;
#endif
    v = MAX(v, val(i_new).next->eval(s, i_new + 1));
  }
  for (index_type k = i; k < i_new; k++) if (s[k]) {
    if (val(k).next) {
#ifdef EVAL_EXTRA_STATS
      eval_rec_count += 1;
#endif
      v = MAX(v, val(k).next->incremental_eval(s, k + 1, i_new));
    }
  }
  return v;
}

NTYPE CostNode::incremental_eval(const bool_vec& s, index_type i_new)
{
  return incremental_eval(s, 0, i_new);
}

NTYPE CostNode::incremental_eval(const index_set& s, index_type i_new)
{
  index_set s_new(s);
  s_new.insert(i_new);
  return eval(s_new);
}

NTYPE CostNode::eval_min(const index_set& s)
{
#ifdef PRINT_DEBUG
  std::cerr << this << ": " << s << std::endl;
#endif

  NTYPE v_min = POS_INF;
  CostNode* n = this;
  index_type l = no_such_index;
  for (index_type k = 0; (k < s.length()) && (n != 0); k++) {
#ifdef PRINT_DEBUG
    std::cerr << "k = " << k
	      << ", s[k] = " << s[k]
	      << std::endl;
#endif
    if (s[k] >= _first) {
      if (l == no_such_index) l = s[k];
      Value& v = n->val(s[k]);
#ifdef PRINT_DEBUG
      std::cerr << "n = " << n
		<< ", val = " << v.val
		<< ", next = " << v.next
		<< std::endl;
#endif
      if (k == (s.length() - 1)) v_min = v.val;
      n = v.next;
    }
  }

  // v_min = value of {atoms on path to this node} U s
  // l = smallest index greater than _first that is in s
#ifdef PRINT_DEBUG
  std::cerr << "v_min = " << v_min << ", l = " << l << std::endl;
#endif
  if (l != no_such_index) {
    for (index_type k = _first; k <= l; k++) {
      CostNode* n1 = val(k).next;
      if (n1) {
	NTYPE v1 = n1->eval_min(s);
	v_min = MIN(v_min, v1);
      }
    }
  }
  else { // s is a subset of {atoms on path to this node}
    for (index_type k = _first; k < _size; k++) {
      v_min = MIN(v_min, val(k).val);
      CostNode* n1 = val(k).next;
      if (n1) {
	NTYPE v1 = n1->eval_min(s);
	v_min = MIN(v_min, v1);
      }
    }
  }

  return v_min;
}

NTYPE CostNode::max_finite() const
{
  NTYPE v = 0;
  for (index_type k = _first; k < _size; k++)
    if (FINITE(_store[k].val)) {
      v = MAX(v, _store[k].val);
      if (_store[k].next) {
	v = MAX(v, _store[k].next->max_finite());
      }
    }
  return v;
}

void CostNode::compute_max()
{
  NTYPE v = 0;
  for (index_type k = _first; k < _size; k++) {
    v = MAX(v, _store[k].val);
    if (_store[k].next) {
      _store[k].next->compute_max();
      v = MAX(v, _store[k].next->_max);
    }
  }
  _max = v;
}

void CostNode::clear_marks()
{
  for (index_type k = _first; k < _size; k++) {
    _store[k].mark = false;
    if (_store[k].next) _store[k].next->clear_marks();
  }
}

void CostNode::write(std::ostream& s, index_set q) const
{
  index_type l = q.length();
  q.inc_length();
  for (index_type k = _first; k < _size; k++) {
    q[l] = k;
    if (true || (_store[k].val != _default) || (_store[k].opt)) {
#ifdef PRINT_DEBUG
      s << "(" << this << ") ";
#endif
      s << q << ": " << _store[k] << std::endl;
    }
    if (_store[k].next) _store[k].next->write(s, q);
  }
  q.dec_length();
}

void CostNode::write_pddl(std::ostream& s, Instance& ins, index_set q) const
{
  assert(ins.n_atoms() >= _size);
  index_type l = q.length();
  q.inc_length();
  for (index_type k = _first; k < _size; k++) {
    q[l] = k;
    if ((_store[k].val != _default) || (_store[k].opt)) {
      s << std::endl << '(';
      for (index_type i = 0; i < q.length(); i++) {
	if (i > 0) s << ' ';
	ins.atoms[q[i]].name->write(s, Name::NC_PDDL);
      }
      s << ')';
      if (_store[k].opt) s << " :opt";
      if (INFINITE(_store[k].val)) s << " :inf";
      else s << ' ' << _store[k].val;
    }
    if (_store[k].next) _store[k].next->write_pddl(s, ins, q);
  }
  q.dec_length();
}

void CostNode::write(std::ostream& s) const
{
  index_set q;
  write(s, q);
}

void CostNode::write_pddl(std::ostream& s, Instance& ins) const
{
  index_set q;
  s << "(:heuristic";
  write_pddl(s, ins, q);
  s << ")" << std::endl;
}

CostTable::CostTable(const Instance& i, Stopwatch& s)
  : Heuristic(i),
    CostNode(i.n_atoms()),
    stats(s),
    pre_cost(0),
    per_cost(0)
{
  // done
}

CostTable::CostTable(const Instance& i, Stopwatch& s, NTYPE v_default)
  : Heuristic(i),
    CostNode(i.n_atoms(), v_default),
    stats(s),
    pre_cost(0),
    per_cost(0)
{
  // done
}

CostTable::~CostTable()
{
  if (pre_cost) {
    delete pre_cost;
    pre_cost = 0;
  }
  if (per_cost) {
    delete per_cost;
    per_cost = 0;
  }
}

CostTable::Entry* CostTable::Entry::first()
{
  Entry* e = this;
  while (e->prev) e = e->prev;
  return e;
}

CostTable::Entry* CostTable::Entry::last()
{
  Entry* e = this;
  while (e->next) e = e->next;
  return e;
}

CostTable::Entry* CostTable::Entry::move_down()
{
  if (next) if (val.val > next->val.val) {
    Entry* p = next;
    while ((val.val > p->val.val) && p->next) p = p->next;
    if (p->val.val <= val.val) {
      unlink();
      place_after(p);
    }
    else {
      unlink();
      place_before(p);
    }
  }
  return first();
}

CostTable::Entry* CostTable::Entry::move_up()
{
  if (prev) if (val.val < prev->val.val) {
    Entry* p = prev;
    while ((val.val < p->val.val) && p->prev) p = p->prev;
    if (p->val.val <= val.val) {
      unlink();
      place_after(p);
    }
    else {
      unlink();
      place_before(p);
    }
  }
  return first();
}

void CostTable::Entry::unlink()
{
  if (prev) prev->next = next;
  if (next) next->prev = prev;
  prev = 0;
  next = 0;
}

void CostTable::Entry::place_before(Entry* p)
{
  prev = p->prev;
  if (p->prev) p->prev->next = this;
  next = p;
  p->prev = this;
}

void CostTable::Entry::place_after(Entry* p)
{
  next = p->next;
  if (p->next) p->next->prev = this;
  p->next = this;
  prev = p;
}

void CostTable::Entry::delete_list()
{
  Entry* l = this;
  while (l) {
    Entry* e = l;
    l = l->next;
    delete e;
  }
}

count_type CostTable::Entry::list_length()
{
  count_type n = 0;
  for (Entry* e = this; e; e = e->next) n += 1;
  return n;
}

CostTable::Entry* CostTable::insert_entry
(CostTable::Entry* e, CostTable::Entry* l)
{
  if (!l) {
    return e;
  }
  if (e->val.val < l->val.val) {
    e->next = l;
    l->prev = e;
    return e;
  }
  Entry* p = l;
  while ((p->val.val <= e->val.val) && p->next) {
    if (p->set == e->set) return l;
    p = p->next;
  }
  if (p->val.val <= e->val.val) {
    e->place_after(p);
  }
  else {
    e->place_before(p);
  }
  return l;
}

CostTable::Entry* CostTable::prepend_entry
(CostTable::Entry* e, CostTable::Entry* l)
{
  if (!l) {
    return e;
  }
  e->next = l;
  l->prev = e;
  return e;
}

CostTable::Entry* CostTable::remove_entry
(CostTable::Entry* e, CostTable::Entry* l)
{
  if (e == l) {
    if (e->next) e->next->prev = 0;
    Entry* r = e->next;
    e->next = 0;
    e->prev = 0;
    return r;
  }
  else {
    if (e->next) e->next->prev = e->prev;
    if (e->prev) e->prev->next = e->next;
    e->next = 0;
    e->prev = 0;
    return l;
  }
}

CostTable::Entry* CostTable::copy_list(Entry* l)
{
  Entry* r = 0;
  while (l) {
    Entry* e = new Entry(l->set, l->val);
    r = insert_entry(e, r);
    l = l->next;
  }
  return r;
}

CostTable::Entry* CostTable::find_entry(index_set& s, Entry* l)
{
  for (Entry* e = l; e; e = e->next)
    if (e->set == s) return e;
  return 0;
}

count_type CostTable::write_list(std::ostream& s, Entry* l)
{
  count_type n = 0;
  for (Entry* e = l; e; e = e->next) {
    instance.write_atom_set(s, e->set);
    s << ": " << e->val << " (" << eval(e->set) << ")" << std::endl;
    n += 1;
  }
  return n;
}

void CostTable::store_list(Entry* l)
{
  for (Entry* e = l; e; e = e->next)
    store(e->set, e->val.val, e->val.opt);
}

CostTable::Entry* CostTable::atoms()
{
  Entry* l = 0;
  index_set s;
  s.set_length(1);
  for (index_type i = 0; i < instance.n_atoms(); i++) {
    s[0] = i;
    Entry* e = new Entry(s, val(i));
    l = insert_entry(e, l);
  }
  return l;
}

CostTable::Entry* CostTable::pairs()
{
  Entry* l = atoms();
  index_set s;
  s.set_length(2);
  for (index_type i = 0; i < instance.n_atoms(); i++) {
    for (index_type j = i + 1; j < instance.n_atoms(); j++) {
      CostNode::Value v = next(i).val(j);
      if (!v.opt && FINITE(v.val)) {
	s[0] = i;
	s[1] = j;
	Entry* e = new Entry(s, v);
	l = insert_entry(e, l);
      }
    }
  }
  return l;
}

bool CostTable::interesting(const index_set& s)
{
  NTYPE v = eval(s);
  for (index_type k = 0; k < s.length(); k++) {
    index_set sub_s(s);
    sub_s.remove(k);
    if (eval(sub_s) >= v) return false;
  }
  return true;
}

CostTable::Entry* CostTable::entries
(CostNode* n, index_set s, Entry* l, const index_set* filter,
 bool ign_zero, bool ign_opt, bool ign_inf, bool req_mark, bool sort)
{
  index_type i = s.length();
  s.inc_length();
  for (index_type k = n->first(); k < instance.n_atoms(); k++) {
    s[i] = k;
    bool accept = true;
    if (accept && req_mark) {
      if (!n->val(k).mark) accept = false;
    }
    if (accept && (filter != 0)) {
      if (!filter->contains(k)) accept = false;
    }
    if (accept && ign_zero) {
      if (n->val(k).val == 0) accept = false;
    }
    if (accept && ign_opt) {
      if (n->val(k).opt) accept = false;
    }
    if (accept && ign_inf) {
      if (INFINITE(n->val(k).val)) accept = false;
    }
    if (accept) {
      if (sort) {
	l = insert_entry(new Entry(s, n->val(k)), l);
      }
      else {
	l = prepend_entry(new Entry(s, n->val(k)), l);
      }
    }
    if (n->val(k).next)
      l = entries(n->val(k).next, s, l, filter, ign_zero, ign_opt, ign_inf,
		  req_mark, sort);
  }
  s.dec_length();
  return l;
}

count_type CostTable::count_entries
(CostNode* n, index_set s, bool ign_zero, bool ign_opt, bool ign_inf)
{
  count_type l = 0;
  index_type i = s.length();
  s.inc_length();
  for (index_type k = n->first(); k < instance.n_atoms(); k++) {
    s[i] = k;
    if (((n->val(k).val > 0) || !ign_zero) &&
	(!n->val(k).opt || !ign_opt) &&
	(FINITE(n->val(k).val) || !ign_inf))
    {
      l += 1;
    }
    if (n->val(k).next)
      l += count_entries(n->val(k).next, s, ign_zero, ign_opt, ign_inf);
  }
  s.dec_length();
  return l;
}

CostTable::Entry* CostTable::entries()
{
  index_set s;
  return entries(this, s, 0, 0, false, false, false, false, true);
}

CostTable::Entry* CostTable::entries
(bool ign_zero, bool ign_opt, bool ign_inf, bool req_mark)
{
  index_set s;
  return entries(this, s, 0, 0, ign_zero, ign_opt, ign_inf, req_mark, true);
}

CostTable::Entry* CostTable::unsorted_entries
(bool ign_zero, bool ign_opt, bool ign_inf, bool req_mark)
{
  index_set s;
  return entries(this, s, 0, 0, ign_zero, ign_opt, ign_inf, req_mark, false);
}

CostTable::Entry* CostTable::boostable_entries()
{
  index_set s;
  return entries(this, s, 0, 0, true, true, true, false, true);
}

CostTable::Entry* CostTable::boostable_entries(const index_set& f)
{
  index_set s;
  return entries(this, s, 0, &f, true, true, true, false, true);
}

CostTable::Entry* CostTable::boost_cd_entries()
{
  index_set s;
  return entries(this, s, 0, 0, true, false, true, false, true);
}

CostTable::Entry* CostTable::marked_entries(bool ign_nb)
{
  index_set s;
  return entries(this, s, 0, 0, ign_nb, ign_nb, ign_nb, true, true);
}

CostTable::Entry* CostTable::entries(const index_set& f)
{
  index_set s;
  return entries(this, s, 0, &f, false, false, false, false, true);
}

count_type CostTable::count_entries(bool ign_zero, bool ign_opt, bool ign_inf)
{
  index_set s;
  return count_entries(this, s, ign_zero, ign_opt, ign_inf);
}

void CostTable::fill(index_set& s, index_type d)
{
  index_type l = s.length();
  s.inc_length();
  for (index_type k = (l > 0 ? s[l-1] + 1 : 0); k < instance.n_atoms(); k++) {
    s[l] = k;
    store(s, eval(s));
    if (d > s.length()) fill(s, d);
  }
  s.dec_length();
}

void CostTable::fill(index_type to_depth)
{
  index_set s;
  fill(s, to_depth);
}

void CostTable::clear()
{
  CostNode::clear();
  if (pre_cost) {
    delete pre_cost;
    pre_cost = 0;
  }
  if (per_cost) {
    delete per_cost;
    per_cost = 0;
  }
}

NTYPE CostTable::eval(index_type i)
{
  return val(i).val;
}

NTYPE CostTable::eval(index_type i, index_type j)
{
  if (i == j) {
    return val(i).val;
  }
  else {
    sort2(i, j);
    NTYPE v = MAX(val(i).val, val(j).val);
    if (val(i).next) {
      v = MAX(v, next(i).val(j).val);
    }
    return v;
  }
}

NTYPE CostTable::eval(index_type i, index_type j, index_type k)
{
  if (i == j) {
    return eval(i, k);
  }
  else if (i == k) {
    return eval(i, j);
  }
  else if (j == k) {
    return eval(i, j);
  }
  else {
    sort3(i, j, k);
    NTYPE v = MAX(MAX(val(i).val, val(j).val), val(k).val);
    if (val(i).next) {
      v = MAX(v, MAX(next(i).val(j).val, next(i).val(k).val));
      if (next(i).val(j).next) {
	v = MAX(v, next(i).next(j).val(k).val);
      }
    }
    if (val(j).next) {
      v = MAX(v, next(j).val(k).val);
    }
    return v;
  }
}

NTYPE CostTable::action_pre_cost(index_type i)
{
  if (pre_cost) return pre_cost[i];
  else return eval(instance.actions[i].pre);
}

NTYPE CostTable::action_per_cost(index_type i)
{
  if (per_cost) return per_cost[i];
  else {
    index_set per_set(instance.actions[i].pre);
    per_set.subtract(instance.actions[i].del);
    return eval(per_set);
  }
}

bool CostTable::update(index_type i, NTYPE v, bool_vec& f)
{
#ifdef TRACE_PRINT_LOTS
  std::cerr << "update: " << i << "." << instance.atoms[i].name
	    << " = " << v << " (" << val(i) << ")" << std::endl;
#endif
  if (SAFE_GT(val(i).val, v)) {
    val(i) = CostNode::Value(v, false);
    for (index_type k = 0; k < instance.atoms[i].req_by.length(); k++)
      f[instance.atoms[i].req_by[k]] = true;
    return true;
  }
  return false;
}

bool CostTable::update(index_type i, NTYPE v)
{
  if (SAFE_GT(val(i).val, v)) {
#ifdef TRACE_PRINT_LOTS
    std::cerr << "update: " << i << "." << instance.atoms[i].name
	      << " = " << v << " (" << val(i) << ")" << std::endl;
#endif
    val(i) = CostNode::Value(v, false);
    return true;
  }
  return false;
}

bool CostTable::update
(index_type i, index_type j, NTYPE v)
{
  if (i == j) {
    CostNode::Value& v_store = val(i);
    if (SAFE_GT(v_store.val, v)) {
#ifdef TRACE_PRINT_LOTS
      std::cerr << "update: " << i << "." << instance.atoms[i].name
		<< ", " << j << "." << instance.atoms[j].name
		<< " = " << v << " (" << v_store << ")"
		<< std::endl;
#endif
      v_store = CostNode::Value(v, false);
      return true;
    }
    return false;
  }
  else {
    sort2(i, j);
    CostNode::Value& v_store = next(i).val(j);
    if (SAFE_GT(v_store.val, v)) {
#ifdef TRACE_PRINT_LOTS
      std::cerr << "update: " << i << "." << instance.atoms[i].name
		<< ", " << j << "." << instance.atoms[j].name
		<< " = " << v << " (" << v_store << ")"
		<< std::endl;
#endif
      v_store = CostNode::Value(v, false);
      return true;
    }
    return false;
  }
}

bool CostTable::update
(index_type i, index_type j, index_type l, NTYPE v)
{
  if (i == j) {
    return update(i, l, v);
  }
  else if (i == l) {
    return update(i, j, v);
  }
  else if (j == l) {
    return update(i, j, v);
  }
  else {
    sort3(i, j, l);
    CostNode::Value& v_store = next(i).next(j).val(l);
#ifdef TRACE_PRINT_LOTS
    std::cerr << "update: " << i << "." << instance.atoms[i].name
	      << ", " << j << "." << instance.atoms[j].name
	      << ", " << l << "." << instance.atoms[l].name
	      << " = " << v << " (" << v_store << ")"
	      << std::endl;
#endif
    if (SAFE_GT(v_store.val, v)) {
      v_store = CostNode::Value(v, false);
      return true;
    }
    return false;
  }
}

// void CostTable::closure(const ACF& cost)
// {
//   bool done = false;
//   while (!done) {
//     done = true;
//     for (index_type k = 0; k < instance.n_actions(); k++) {
//       if (stats.break_signal_raised()) return;
//       NTYPE c_pre = eval(instance.actions[k].pre);
//       if (FINITE(c_pre)) {
// 	for (index_type i = 0; i < instance.actions[k].add.length(); i++)
// 	  if (update(instance.actions[k].add[i], c_pre + cost(k)))
// 	    done = false;
//       }
//     }
//   }
// }
//
// void CostTable::closure(const ACF& cost, const index_set& nd)
// {
//   bool r_sel[instance.n_actions()];
//   for (index_type k = 0; k < instance.n_actions(); k++) {
//     r_sel[k] = instance.actions[k].sel;
//     if (instance.actions[k].del.first_common_element(nd) != no_such_index)
//       instance.actions[k].sel = false;
//   }
//   closure(cost);
//   for (index_type k = 0; k < instance.n_actions(); k++)
//     instance.actions[k].sel = r_sel[k];
// }

bool CostTable::compatible
(index_type a, index_type b, bool opt_resources) const
{
  // check non-interference (i.e. non-strict commutativity)
  if (!instance.non_interfering(a, b)) return false;
  // check locks
  if (!instance.lock_compatible(a, b)) return false;
  // check resources
  if (opt_resources) {
    if (!instance.resource_compatible(a, b)) return false;
  }
  return true;
}

void CostTable::compute_H1(const ACF& cost)
{
  bool_vec init(instance.init_atoms, instance.n_atoms());
  compute_H1(cost, init);
}

void CostTable::compute_H1
(const ACF& cost, const bool_vec& init, const bool_vec* aa)
{
  stats.start();
  bool_vec f(false, instance.n_actions());
  if (!pre_cost) pre_cost = new NTYPE[instance.n_actions()];

  read_only_vector_with_default<bool> aa1(aa, true);

  // init pre_cost array
  for (index_type k = 0; k < instance.n_actions(); k++) {
    pre_cost[k] = POS_INF;
    if ((instance.actions[k].pre.length() == 0) && aa1[k])
      f[k] = true;
  }

  // init table
  for (index_type k = 0; k < instance.n_atoms(); k++) {
    if (init[k]) {
      val(k) = CostNode::Value(0, true);
      for (index_type i = 0; i < instance.atoms[k].req_by.length(); i++)
	if (aa1[instance.atoms[k].req_by[i]])
	  f[instance.atoms[k].req_by[i]] = true;
    }
    else {
      val(k) = CostNode::Value(POS_INF, false);
    }
  }

  bool done = false;
  while (!done) {
    done = true;
    for (index_type k = 0; k < instance.n_actions(); k++)
      if (f[k] && aa1[k]) {
	if (stats.break_signal_raised()) return;
	NTYPE c_pre = eval(instance.actions[k].pre);
	if (FINITE(c_pre) && (c_pre < pre_cost[k])) {
	  for (index_type i = 0; i < instance.actions[k].add.length(); i++)
	    update(instance.actions[k].add[i], c_pre + cost(k), f);
	  pre_cost[k] = c_pre;
	  done = false;
	}
	f[k] = false;
      }
  }
  stats.stop();
}


void CostTable::extend_goal_set
(const index_set& sgset, const AnyACF& costs, bool_vec& ext_goal_set)
{
  if (stats.break_signal_raised()) return;
  // find all h^1-maximisers in the subgoal set, and check if any of them
  // already belongs to the extended goal set
  index_type a_max = no_such_index;
  NTYPE c_max = 0;
  bool  max_in_set = false;
  for (index_type k = 0; k < sgset.size(); k++) {
    if (eval(sgset[k]) > c_max) {
      c_max = eval(sgset[k]);
      max_in_set = false;
      if (ext_goal_set[sgset[k]])
	max_in_set = true;
      a_max = sgset[k];
    }
    else if (eval(sgset[k]) == c_max) {
      if (ext_goal_set[sgset[k]])
	max_in_set = true;
    }
  }
  // if no maximiser is in the extended goal set, we have to add one
  // (a_max), and recurse on all actions that add it and have zero cost.
  // if a_max == no_such_index, all atoms in sgset have h^1 value zero.
  if (!max_in_set && (a_max != no_such_index)) {
    ext_goal_set[a_max] = true;
    for (index_type k = 0; k < instance.atoms[a_max].add_by.size(); k++) {
      index_type i = instance.atoms[a_max].add_by[k];
      if (IS_ZERO(costs(i))) {
	extend_goal_set(instance.actions[i].pre, costs, ext_goal_set);
      }
    }
  }
}

void CostTable::find_cut
(const bool_vec& ext_goal_set, const AnyACF& costs, bool_vec& cut)
{
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (!IS_ZERO(costs(k)))
      if (instance.actions[k].add.have_common_element(ext_goal_set))
	if (FINITE(eval(instance.actions[k].pre)))
	  cut[k] = true;
}

NTYPE CostTable::compute_lmcut
(const ACF& cost, const bool_vec& s,
 const index_set& g, const bool_vec* a)
{
  stats.start();
  AnyACF costs(instance.n_actions(), cost);
  compute_H1(costs, s, a);
  if (INFINITE(eval(g)))
    return POS_INF;
  NTYPE total_cost = 0;
  while (eval(g) > 0) {
    // std::cerr << "h^1 table:" << std::endl;
    // table->write_pddl(std::cerr, instance);
    if (stats.break_signal_raised()) return 0;
    bool_vec ext_goal_set(false, instance.n_atoms());
    extend_goal_set(g, costs, ext_goal_set);
    // std::cerr << "extended goal set: ";
    // instance.write_atom_set(std::cerr, ext_goal_set);
    // std::cerr << std::endl;
    bool_vec allowed_acts(true, instance.n_actions());
    if (a) {
      for (index_type k = 0; k < instance.n_actions(); k++)
	if (!(*a)[k]) allowed_acts[k] = false;
    }
    for (index_type k = 0; k < instance.n_actions(); k++)
      if (allowed_acts[k] &&
	  instance.actions[k].pre.have_common_element(ext_goal_set))
	allowed_acts[k] = false;
    compute_H1(costs, s, &allowed_acts);
    // std::cerr << "revised h^1 table:" << std::endl;
    // table->write_pddl(std::cerr, instance);
    bool_vec cut(false, instance.n_actions());
    find_cut(ext_goal_set, costs, cut);
    cut.intersect(allowed_acts);
    NTYPE c_cut = costs.min_cost(cut);
    // std::cerr << "cut: ";
    // instance.write_action_set(std::cerr, cut);
    // std::cerr << ", cost = " << c_cut << std::endl;
    assert(c_cut > 0);
    total_cost += c_cut;
    costs.decrease(cut, c_cut);
    compute_H1(costs, s);
  }
  // std::cerr << "total cost = " << total_cost << std::endl;
  stats.stop();
  return total_cost;
}

NTYPE CostTable::compute_lmcut
(const ACF& cost, const index_set& s,
 const index_set& g, const bool_vec* a)
{
  bool_vec s1(s, instance.n_atoms());
  return compute_lmcut(cost, s1, g, a);
}

void CostTable::compute_H1max
(const ACF& cost, const bool_vec& init, const bool_vec* aa)
{
  stats.start();
  bool_vec f(false, instance.n_actions());
  if (!pre_cost) pre_cost = new NTYPE[instance.n_actions()];

  read_only_vector_with_default<bool> aa1(aa, true);

  // init pre_cost array
  for (index_type k = 0; k < instance.n_actions(); k++) {
    pre_cost[k] = POS_INF;
    if ((instance.actions[k].pre.length() == 0) && aa1[k])
      f[k] = true;
  }

  // init table
  for (index_type k = 0; k < instance.n_atoms(); k++) {
    if (init[k]) {
      val(k) = CostNode::Value(0, true);
      for (index_type i = 0; i < instance.atoms[k].req_by.length(); i++)
	if (aa1[instance.atoms[k].req_by[i]])
	  f[instance.atoms[k].req_by[i]] = true;
    }
    else {
      val(k) = CostNode::Value(POS_INF, false);
    }
  }

  bool done = false;
  while (!done) {
    done = true;
    for (index_type k = 0; k < instance.n_actions(); k++) if (f[k]) {
      if (stats.break_signal_raised()) return;
      NTYPE c_pre = eval(instance.actions[k].pre);
      if (FINITE(c_pre) && (c_pre < pre_cost[k])) {
	for (index_type i = 0; i < instance.actions[k].add.length(); i++)
	  update(instance.actions[k].add[i], MAX(c_pre, cost(k)), f);
	pre_cost[k] = c_pre;
	done = false;
      }
      f[k] = false;
    }
  }
  stats.stop();
}

void CostTable::compute_H2(const ACF& cost)
{
  bool_vec init(instance.init_atoms, instance.n_atoms());
  compute_H2(cost, init);
}

void CostTable::compute_H2
(const ACF& cost, const bool_vec& init, const bool_vec* aa)
{
  if (trace_level > 2) {
    std::cerr << "computing H2..." << std::endl;
  }
  stats.start();

  // init table
  for (index_type i = 0; i < instance.n_atoms(); i++) {
    if (init[i])
      val(i) = CostNode::Value(0, true);
    else
      val(i) = CostNode::Value(POS_INF, false);
    for (index_type j = i + 1; j < instance.n_atoms(); j++) {
      if (init[i] && init[j])
	next(i).val(j) = CostNode::Value(0, true);
      else
	next(i).val(j) = CostNode::Value(POS_INF, false);
    }
  }

  read_only_vector_with_default<bool> aa1(aa, true);

  bool done = false;
  while (!done) {
    if (trace_level > 2) {
      std::cerr << "begin iteration..." << std::endl;
    }
    done = true;
    for (index_type a = 0; a < instance.n_actions(); a++) if (aa1[a]) {
      if (stats.break_signal_raised()) return;
      const Instance::Action& act = instance.actions[a];
      NTYPE c_pre = eval(act.pre);
      if (FINITE(c_pre)) {
#ifdef TRACE_PRINT_LOTS
	std::cerr << "updating over " << a << "." << act.name
		  << " (pre = " << act.pre << ", add = " << act.add
		  << ", del = " << act.del << ", cost = " << cost(a)
		  << ")..." << std::endl;
#endif
	for (index_type i = 0; i < instance.n_atoms(); i++) {
	  // if i in act.add...
	  if (act.add.contains(i)) {
	    // update cost of i
	    if (update(i, c_pre + cost(a))) done = false;
	    for (index_type j = 0; j < act.add.length(); j++)
	      // update cost of {i,q} for every q in act.add
	      if (act.add[j] > i) // q < i have already been updated (as i)
		if (update(act.add[j], i, c_pre + cost(a))) done = false;
	  }
	  // else if i not in act.del (i.e. case act + noop(i))...
	  else if (!act.del.contains(i)) {
	    // update cost of {i,q} for every q in act.add
	    for (index_type j = 0; j < act.add.length(); j++) {
	      // compute H(act.pre U {i})
	      NTYPE c_join = MAX(c_pre, eval(i));
	      if (!act.pre.contains(i)) {
		for (index_type k = 0; k < act.pre.length(); k++)
		  c_join = MAX(eval(act.pre[k], i), c_join);
	      }
	      if (FINITE(c_join))
		if (update(act.add[j], i, c_join + cost(a))) done = false;
	    }
	  }
	} // for each atom i
      } // if FINITE(c_pre)
    } // for each action
    if (trace_level > 2) {
      std::cerr << "end iteration (done = " << done << ")" << std::endl;
    }
  }
  stats.stop();
}

void CostTable::compute_H3(const ACF& cost, const bool_vec* aa)
{
  stats.start();

  // init table
  for (index_type i = 0; i < instance.n_atoms(); i++) {
    if (instance.atoms[i].init)
      val(i) = CostNode::Value(0, true);
    else
      val(i) = CostNode::Value(POS_INF, false);
    for (index_type j = i + 1; j < instance.n_atoms(); j++) {
      if (instance.atoms[i].init && instance.atoms[j].init)
	next(i).val(j) = CostNode::Value(0, true);
      else
	next(i).val(j) = CostNode::Value(POS_INF, false);
      for (index_type l = j + 1; l < instance.n_atoms(); l++) {
	if (instance.atoms[i].init &&
	    instance.atoms[j].init &&
	    instance.atoms[l].init)
	  next(i).next(j).val(l) = CostNode::Value(0, true);
	else
	  next(i).next(j).val(l) = CostNode::Value(POS_INF, false);
      }
    }
  }

  read_only_vector_with_default<bool> aa1(aa, true);

  bool done = false;
  while (!done) {
    done = true;
    for (index_type a = 0; a < instance.n_actions(); a++) if (aa1[a]) {
      if (stats.break_signal_raised()) return;
      const Instance::Action& act = instance.actions[a];
      NTYPE c_pre = eval(act.pre);
      if (FINITE(c_pre)) {
#ifdef TRACE_PRINT_LOTS
	std::cerr << "updating over " << a << "." << act.name
		  << " (pre = " << act.pre << ", add = " << act.add
		  << ", del = " << act.del << ", cost = " << cost(a)
		  << ")..." << std::endl;
#endif
	for (index_type i = 0; i < instance.n_atoms(); i++) {
	  // if i in act.add...
	  if (act.add.contains(i)) {
	    // update cost of i
	    if (update(i, c_pre + cost(a))) done = false;
	    for (index_type j = 0; j < act.add.length(); j++)
	      // update cost of {i,q} for every q in act.add
	      if (act.add[j] > i) { // q < i have already been updated (as i)
		if (update(i, act.add[j], c_pre + cost(a))) done = false;
		// update cost of {i,q,q'} for every q' in act.add
		for (index_type l = j + 1; l < act.add.length(); l++)
		  if (act.add[l] > i)
		    if (update(i, act.add[j], act.add[l], c_pre + cost(a)))
		      done = false;
	      }
	  }
	  // else if i not in act.del (i.e. case act + noop(i))...
	  else if (!act.del.contains(i)) {
	    // compute H(act.pre U {i})
	    NTYPE c_join = MAX(c_pre, eval(i));
	    if (!act.pre.contains(i)) {
	      for (index_type k = 0; k < act.pre.length(); k++)
		c_join = MAX(eval(act.pre[k], i), c_join);
	    }
	    if (FINITE(c_join)) {
	      // update cost of {i,q} and {i,q,q'} for every q,q' in act.add
	      for (index_type j = 0; j < act.add.length(); j++) {
		if (update(i, act.add[j], c_join + cost(a))) done = false;
		for (index_type l = j + 1; l < act.add.length(); l++)
		  if (update(i, act.add[j], act.add[l], c_join + cost(a)))
		    done = false;
	      }
	      // finally, check case act + noop(i) + noop(j)
	      for (index_type j = i + 1; j < instance.n_atoms(); j++)
		if (!act.del.contains(j)) {
		  // compute H(act.pre U {i,j})
		  NTYPE c_jj = MAX(c_join, eval(i, j));
		  if (!act.pre.contains(j)) {
		    for (index_type k = 0; k < act.pre.length(); k++)
		      c_jj = MAX(eval(act.pre[k], j), c_jj);
		  }
		  if (FINITE(c_jj)) {
		    // update cost of {i,j,q} for every q in act.add
		    for (index_type l = 0; l < act.add.length(); l++)
		      if (update(i, j, act.add[l], c_jj + cost(a)))
			done = false;
		  }
		}
	    } // FINITE(c_join)
	  }
	} // for each atom i
      } // if FINITE(c_pre)
    } // for each action
  }
  stats.stop();
}

void CostTable::compute_H2C(const ACF& cost, bool opt_resources)
{
  stats.start();

  // init table
  for (index_type i = 0; i < instance.n_atoms(); i++) {
    if (instance.atoms[i].init)
      val(i) = CostNode::Value(0, true);
    else
      val(i) = CostNode::Value(POS_INF, false);
    for (index_type j = i + 1; j < instance.n_atoms(); j++) {
      if (instance.atoms[i].init && instance.atoms[j].init)
	next(i).val(j) = CostNode::Value(0, true);
      else
	next(i).val(j) = CostNode::Value(POS_INF, false);
    }
  }

  bool done = false;
  while (!done) {
    done = true;
    for (index_type a = 0; a < instance.n_actions(); a++) {
      if (stats.break_signal_raised()) return;
      const Instance::Action& act = instance.actions[a];
      NTYPE c_pre = eval(act.pre);
      if (FINITE(c_pre)) {
#ifdef TRACE_PRINT_LOTS
	std::cerr << "updating over " << a << "." << act.name
		  << "..." << std::endl;
#endif
	for (index_type i = 0; i < instance.n_atoms(); i++) {
	  // if i in act.add...
	  if (act.add.contains(i)) {
	    // update cost of i
	    if (update(i, c_pre + cost(a))) done = false;
	    for (index_type j = 0; j < act.add.length(); j++)
	      // update cost of {i,q} for every q in act.add
	      if (act.add[j] > i) // q < i have already been updated (as i)
		if (update(act.add[j], i, c_pre + cost(a))) done = false;
	  }
	  // else if i not in act.del (i.e. case act + noop(i))...
#ifdef SUPPORT_VOLATILE_ATOMS
	  else if (!act.del.contains(i) && !instance.atoms[i].volatile) {
#else
          else if (!act.del.contains(i)) {
#endif
	    // update cost of {i,q} for every q in act.add
	    for (index_type j = 0; j < act.add.length(); j++) {
	      // compute H(act.pre U {i})
	      NTYPE c_join = MAX(c_pre, eval(i));
	      if (!act.pre.contains(i)) {
		for (index_type k = 0; k < act.pre.length(); k++)
		  c_join = MAX(eval(act.pre[k], i), c_join);
	      }
	      if (FINITE(c_join))
		if (update(act.add[j], i, c_join + cost(a))) done = false;
	    }
	  }
	}
	// update over act + act_b for every act_b compatible with act
	// (only need to update over b > a, since {a,b} and {b,a} symmetric
	for (index_type b = a + 1; b < instance.n_actions(); b++) {
	  const Instance::Action& act_b = instance.actions[b];
	  if (compatible(a, b, opt_resources)) {
	    // compute H(act.pre U act_b.pre)
	    NTYPE c_b = eval(act_b.pre);
	    if (FINITE(c_b)) {
	      NTYPE c_join = MAX(c_pre, c_b);
	      for (index_type i = 0; i < act.pre.length(); i++)
		for (index_type j = 0; j < act_b.pre.length(); j++)
		  c_join = MAX(eval(act.pre[i], act_b.pre[j]), c_join);
	      // c_conc is H of act + act_b concurrently
	      NTYPE c_conc = MAX(MAX(c_pre + cost(a), c_b + cost(b)),
				 c_join + MIN(cost(a), cost(b)));
	      if (FINITE(c_conc)) {
		// only need to update pairs across act/act_b (pairs added
		// by act alone are handled above, pairs added by act_b alone
		// are updated when the outer loop reaches b)
		for (index_type i = 0; i < act.add.length(); i++)
		  for (index_type j = 0; j < act_b.add.length(); j++)
		    if (update(act.add[i], act_b.add[j], c_conc))
		      done = false;
	      } // if FINTITE(c_conc)
	    } // if FINITE(c_b)
	  }
	}
      } // if FINITE(c_pre)
    } // for each action
  }
  stats.stop();
}

void CostTable::verify_H2C
(const ACF& cost, bool opt_resources, bool opt_verbose)
{
  for (index_type i = 0; i < instance.n_atoms(); i++) {
    const Instance::Atom& atom1 = instance.atoms[i];
    if (atom1.init) {
      if (opt_verbose)
	std::cerr << "atom " << atom1.name << " is initial" << std::endl;
    }
    else {
      if (opt_verbose) {
	std::cerr << "checking atom " << atom1.name << "..." << std::endl;
      }
      NTYPE c_tab = eval(atom1.index);
      NTYPE c_min = POS_INF;
      for (index_type k = 0; k < atom1.add_by.length(); k++) {
	const Instance::Action& act = instance.actions[atom1.add_by[k]];
	NTYPE c_act = eval(act.pre) + cost(act.index);
	if (opt_verbose) {
	  std::cerr << act.name << ": " << c_act << std::endl;
	}
	c_min = MIN(c_act, c_min);
      }
      if (opt_verbose) {
	std::cerr << "min " << atom1.name << " = " << c_min << ", ";
	if (c_min != c_tab) {
	  std::cerr << "error: H(" << atom1.name << ") = " << c_tab << std::endl;
	}
	else {
	  std::cerr << "ok" << std::endl;
	}
      }
      else {
	if (c_min > c_tab) {
	  std::cerr << "improvement: H(" << atom1.name << ") = " << c_min
		    << " (stored: " << c_tab << ")" << std::endl;
	}
      }
    }
    for (index_type j = i+1; j < instance.n_atoms(); j++) {
      const Instance::Atom& atom2 = instance.atoms[j];
      if (atom1.init && atom2.init) {
	if (opt_verbose)
	  std::cerr << "pair {" << atom1.name << ", " << atom2.name
		    << "} is initial" << std::endl;
      }
      else {
	if (opt_verbose) {
	  std::cerr << "checking pair {" << atom1.name
		    << ", " << atom2.name << "}..." << std::endl;
	}
	NTYPE c_tab = eval(atom1.index, atom2.index);
	NTYPE c_min = POS_INF;
	for (index_type k = 0; k < atom1.add_by.length(); k++) {
	  const Instance::Action& act1 = instance.actions[atom1.add_by[k]];
	  if (!act1.del.contains(atom2.index)) {
	    // act1 adds both atoms
	    if (act1.add.contains(atom2.index)) {
	      NTYPE c_both = eval(act1.pre) + cost(act1.index);
	      if (opt_verbose) {
		std::cerr << act1.name << ": " << c_both << std::endl;
	      }
	      c_min = MIN(c_both, c_min);
	    }
	    else {
	      // act1 + NOOP(atom2)
	      index_set s0(act1.pre);
	      s0.insert(atom2.index);
	      NTYPE c_1n = eval(s0) + cost(act1.index);
	      if (opt_verbose) {
		std::cerr << act1.name << ", NOOP(" << atom2.name << "): "
			  << c_1n << std::endl;
	      }
	      c_min = MIN(c_1n, c_min);
	      for (index_type l = 0; l < atom2.add_by.length(); l++) {
		const Instance::Action& act2 =
		  instance.actions[atom2.add_by[l]];
		if (compatible(act1.index, act2.index, opt_resources) &&
		    !act2.add.contains(atom1.index)) {
		  index_set s(act1.pre);
		  s.insert(act2.pre);
		  NTYPE c_conc;
		  if (cost(act1.index) > cost(act2.index)) {
		    NTYPE c_1 = eval(act1.pre) + cost(act1.index);
		    NTYPE c_join = eval(s) + cost(act2.index);
		    c_conc = MAX(c_1, c_join);
		  }
		  else {
		    NTYPE c_2 = eval(act2.pre) + cost(act2.index);
		    NTYPE c_join = eval(s) + cost(act1.index);
		    c_conc = MAX(c_2, c_join);
		  }
		  if (opt_verbose) {
		    std::cerr << act1.name << ", " << act2.name << ": "
			      << c_conc << std::endl;
		  }
		  c_min = MIN(c_conc, c_min);
		}
	      }
	    }
	  }
	}
	for (index_type k = 0; k < atom2.add_by.length(); k++) {
	  const Instance::Action& act = instance.actions[atom2.add_by[k]];
	  if (!act.del.contains(atom1.index) &&
	      !act.add.contains(atom1.index)) {
	    index_set s(act.pre);
	    s.insert(atom1.index);
	    NTYPE c_2n = eval(s) + cost(act.index);
	    if (opt_verbose) {
	      std::cerr << "NOOP(" << atom1.name << "), " << act.name << ": "
			<< c_2n << std::endl;
	    }
	    c_min = MIN(c_2n, c_min);
	  }
	}
	if (opt_verbose) {
	  std::cerr << "min {" << atom1.name << ", " << atom2.name << "} = "
		    << c_min << ", ";
	  if (c_min != c_tab) {
	    std::cerr << "error: H(" << atom1.name << ", " << atom2.name
		      << ") = " << c_tab << std::endl;
	  }
	  else {
	    std::cerr << "ok" << std::endl;
	  }
	}
	else {
	  if (c_min > c_tab) {
	    std::cerr << "improvement: H(" << atom1.name << "," << atom2.name
		      << ") = " << c_min << " (stored: " << c_tab << ")"
		      << std::endl;
	  }
	}
      }
    }
  }
}

#ifdef BUILD_WITH_BOOSTING

CostTable::Entry* CostTableBoost::boost_cd
(Entry* list, StateFactory& search_space, HashTable* solved_tab,
 index_type cd_size_limit, const index_set& cd_atoms_ignore,
 NTYPE cost_limit, const index_set& cost_limit_set,
 count_type wps_limit, bool scale_wps)
{
  Conflicts boost_search_res(instance);
  Statistics boost_search_stats(&bstats);
  IDAO boost_search(boost_search_stats, boost_search_res, solved_tab);
  boost_search.set_trace_level(trace_level - 2);
  boost_search.set_store_cost(false);

  if (trace_level > 2) {
    std::cerr << "entries in: ";
    write_list(std::cerr, list);
    std::cerr << std::endl;
  }

  bstats.start();
  NTYPE current_bound = 0;
  count_type n_entries_in = (list ? list->list_length() : 0);
  count_type n_before = n_entries_created;

  while (list) {
    Entry* e = list;

    if (INFINITE(e->val.val)) {
      bstats.stop();
      count_type n_entries_out = (list ? list->list_length() : 0);
      n_entries_boosted += (n_entries_in - n_entries_out);
      if (boost_new_entries)
	n_entries_boosted += (n_entries_created - n_before);
      if (trace_level > 0) {
	std::cerr << "done at INF (" << n_entries_solved << " solved, "
		  << n_entries_created << " created, "
		  << n_entries_discarded << " discarded, "
		  << (list ? list->list_length() : 0) << " in list, "
		  << bstats.nodes() << " nodes, "
		  << bstats.time() << " sec.)" << std::endl;
      }
      if (trace_level > 2) {
	std::cerr << "entries out: ";
	write_list(std::cerr, list);
	std::cerr << std::endl;
      }
      return list;
    }

    if (e->val.val > cost_limit) {
      cost_limit = eval(cost_limit_set);
      if (e->val.val > cost_limit) {
	bstats.stop();
	count_type n_entries_out = (list ? list->list_length() : 0);
	n_entries_boosted += (n_entries_in - n_entries_out);
	if (boost_new_entries)
	  n_entries_boosted += (n_entries_created - n_before);
	if (trace_level > 0) {
	  std::cerr << "done at " << e->val.val
		    << " (" << n_entries_solved << " solved, "
		    << n_entries_created << " created, "
		    << n_entries_discarded << " discarded, "
		    << (list ? list->list_length() : 0) << " in list, "
		    << bstats.nodes() << " nodes, "
		    << bstats.time() << " sec.)" << std::endl;
	}
	if (trace_level > 2) {
	  std::cerr << "entries out: ";
	  write_list(std::cerr, list);
	  std::cerr << std::endl;
	}
	return list;
      }
      else if (trace_level > 0) {
	std::cerr << "limit is " << cost_limit << std::endl;
      }
    }

    NTYPE boost_limit = cost_limit;
    if (e->next) {
      boost_limit = e->next->val.val;
    }

    list = e->next;
    e->unlink();

    if (trace_level > 0) {
      if (e->val.val > current_bound) {
	std::cerr << "at " << e->val.val
		  << " (" << n_entries_solved << " solved, "
		  << n_entries_created << " created, "
		  << n_entries_discarded << " discarded, "
		  << (list ? list->list_length() : 0) << " in list, "
		  << bstats.nodes() << " nodes, "
		  << bstats.time() << " sec.)" << std::endl;
	current_bound = e->val.val;
      }
    }

    count_type w_lim = wps_limit;
    bool wps_break_flag = false;
    if (scale_wps) {
      if (2*(e->val.work) > w_lim) w_lim = 2*(e->val.work);
    }
    else {
      if (e->val.work > w_lim) wps_break_flag = true;
    }

    if (!wps_break_flag) {
      State* boost_state = search_space.new_state(e->set, 0);
      NTYPE c_est = boost_state->est_cost();

      NTYPE c_new = 0;
      if (e->set.length() < cd_size_limit) {
	if (trace_level > 2) {
	  std::cerr << "boosting " << *boost_state
		    << " (cost = " << c_est
		    << ", limit = " << boost_limit
		    << ", cd = "
		    << (strong_conflict_detection ? "strong" :
			(ultra_weak_conflict_detection ? "ultra weak" : "weak"))
		    << ", wps = " << w_lim
		    << ")... " << std::endl;
	}
	if (strong_conflict_detection || ultra_weak_conflict_detection) {
	  boost_search_res.set_stop_condition(Result::stop_at_all_optimal);
	}
	else {
	  boost_search_res.set_stop_condition(Result::stop_at_first);
	}
	boost_search_res.reset();
	boost_search.set_cost_limit(boost_limit);
	boost_search_stats.enable_eval_limit(boost_search_stats.evaluations() + w_lim);
	c_new = boost_search.start(*boost_state);
	n_boost_searches_with_cd += 1;
      }
      else {
	if (trace_level > 2) {
	  std::cerr << "boosting " << *boost_state
		    << " (cost = " << c_est
		    << ", limit = " << boost_limit
		    << ", cd = off, wps = " << w_lim
		    << ")... " << std::endl;
	}
	boost_search_res.reset();
	boost_search_res.set_stop_condition(Result::stop_at_first);
	boost_search.set_cost_limit(boost_limit);
	boost_search_stats.enable_eval_limit(boost_search_stats.evaluations() + w_lim);
	c_new = boost_search.start(*boost_state);
      }

      wps_break_flag = ((c_new == e->val.val) && !boost_search.solved());
      if (wps_break_flag) e->val.work = w_lim;

      e->val = CostNode::Value(c_new, boost_search.solved());
      n_boost_searches += 1;

      if (trace_level > 1) {
	if (wps_break_flag) {
	  std::cerr << *boost_state << " (cost = " << c_est
		    << ") discarded (wps = " << w_lim
		    << ", " << bstats << ")" << std::endl;
	}
	else {
	  std::cerr << *boost_state << ": cost " << c_est
		    << ", max " << boost_limit
		    << ", new " << e->val << " (" << bstats << ")"
		    << std::endl;
	}
      }

      if (boost_search.solved()) {
	if (e->set.length() < cd_size_limit) {
	  list = create_new_entries(e->set,
				    (ultra_weak_conflict_detection ?
				     boost_search_res.pos_deleted() :
				     boost_search_res.nec_deleted()),
				    cd_atoms_ignore, list);
	}
	e->val.opt = true;
	delete e;
	n_entries_solved += 1;
      }
      else if (wps_break_flag) {
	n_entries_discarded += 1;
	delete e;
      }
      else {
	list = insert_entry(e, list);
      }

      delete boost_state;
    }

    // scale_wps is false && e->work > wps_limit
    else {
      n_entries_discarded += 1;
      delete e;
    }

    if (bstats.break_signal_raised()) {
      bstats.stop();
      count_type n_entries_out = (list ? list->list_length() : 0);
      n_entries_boosted += (n_entries_in - n_entries_out);
      if (boost_new_entries)
	n_entries_boosted += (n_entries_created - n_before);
      if (trace_level > 2) {
	std::cerr << "entries out: ";
	write_list(std::cerr, list);
	std::cerr << std::endl;
      }
      return list;
    }
  }

  bstats.stop();
  count_type n_entries_out = (list ? list->list_length() : 0);
  n_entries_boosted += (n_entries_in - n_entries_out);
  if (boost_new_entries)
    n_entries_boosted += (n_entries_created - n_before);
  if (trace_level > 0) {
    std::cerr << "done (" << n_entries_solved << " solved, "
	      << n_entries_created << " created, "
	      << n_entries_discarded << " discarded, "
	      << (list ? list->list_length() : 0) << " in list, "
	      << bstats.nodes() << " nodes, "
	      << bstats.time() << " sec.)" << std::endl;
  }
  if (trace_level > 2) {
    std::cerr << "entries out: ";
    write_list(std::cerr, list);
    std::cerr << std::endl;
  }
  return list;
}

CostTable::Entry* CostTableBoost::create_new_entries
(const index_set& set, const index_set& conflict_set,
 const index_set& cd_ignore, CostTable::Entry* list)
{
  for (index_type k = 0; k < conflict_set.length(); k++) {
    if (!set.contains(conflict_set[k]) &&
	!cd_ignore.contains(conflict_set[k])) {
      index_set new_set(set);
      new_set.insert(conflict_set[k]);
#ifdef ASSUME_UNIT_COST
      NTYPE v = eval(new_set) + (strong_conflict_detection ? 1 : 0);
#else
      NTYPE v = eval(new_set);
#endif
      if (FINITE(v)) {
	CostNode::Value& new_set_val = insert(new_set);
	if (!new_set_val.opt) {
	  if (new_set_val.val == 0) n_entries_created += 1;
	  if (v > new_set_val.val) new_set_val.val = v;
	  if (boost_new_entries && (find_entry(new_set, list) == 0)) {
	    Entry* new_entry = new Entry(new_set, new_set_val);
	    if (trace_level > 2) {
	      instance.write_atom_set(std::cerr, new_set);
	      std::cerr << " new at " << new_set_val << std::endl;
	    }
	    list = insert_entry(new_entry, list);
	  }
	}
      }
    }
  }
  return list;
}

CostTable::Entry* CostTableBoost::boost
(Entry* list, StateFactory& search_space, HashTable* solved_tab,
 NTYPE cost_limit, const index_set& cost_limit_set,
 count_type wps_limit, bool scale_wps)
{
  Result boost_search_res;
  Statistics boost_search_stats(&bstats);
  IDAO boost_search(boost_search_stats, boost_search_res, solved_tab);
  boost_search.set_trace_level(trace_level - 2);
  boost_search_res.set_n_to_find(1);
  boost_search.set_store_cost(false);

  bstats.start();
  NTYPE current_bound = 0;
  count_type n_entries_in = (list ? list->list_length() : 0);

  while (list) {
    Entry* e = list;

    if (INFINITE(e->val.val)) {
      bstats.stop();
      count_type n_entries_out = (list ? list->list_length() : 0);
      n_entries_boosted += (n_entries_in - n_entries_out);
      if (trace_level > 0) {
	std::cerr << "done at INF (" << n_entries_solved << " solved, "
		  << n_entries_discarded << " discarded, "
		  << (list ? list->list_length() : 0) << " in list, "
		  << bstats.nodes() << " nodes, "
		  << bstats.time() << " sec.)" << std::endl;
      }
      return list;
    }

    if (e->val.val > cost_limit) {
      cost_limit = eval(cost_limit_set);
      if (e->val.val > cost_limit) {
	bstats.stop();
	count_type n_entries_out = (list ? list->list_length() : 0);
	n_entries_boosted += (n_entries_in - n_entries_out);
	if (trace_level > 0) {
	  std::cerr << "done at " << e->val.val
		    << " (" << n_entries_solved << " solved, "
		    << n_entries_discarded << " discarded, "
		    << (list ? list->list_length() : 0) << " in list, "
		    << bstats.nodes() << " nodes, "
		    << bstats.time() << " sec.)" << std::endl;
	}
	return list;
      }
      else if (trace_level > 0) {
	std::cerr << "limit is " << cost_limit << std::endl;
      }
    }

    NTYPE boost_limit = cost_limit;
    if (e->next) {
      boost_limit = e->next->val.val;
    }

    list = e->next;
    e->unlink();

    if (trace_level > 0) {
      if (e->val.val > current_bound) {
	std::cerr << "at " << e->val.val
		  << " (" << n_entries_solved << " solved, "
		  << n_entries_discarded << " discarded, "
		  << (list ? list->list_length() : 0) << " in list, "
		  << bstats.nodes() << " nodes, "
		  << bstats.time() << " sec.)" << std::endl;
	current_bound = e->val.val;
      }
    }

    if (e->val.opt) {
      delete e;
    }
    else {
      count_type w_lim = wps_limit;
      bool wps_break_flag = false;
      if (scale_wps) {
	if (2*(e->val.work) > w_lim) w_lim = 2*(e->val.work);
      }
      else {
	if (e->val.work > w_lim) wps_break_flag = true;
      }

      if (!wps_break_flag) {
	State* boost_state = search_space.new_state(e->set, 0);
	NTYPE c_est = boost_state->est_cost();

	if (trace_level > 2) {
	  std::cerr << "boosting " << *boost_state << " at " << c_est
		    << "..." << std::endl;
	}

	boost_search.set_cost_limit(boost_limit);
	boost_search_stats.enable_eval_limit(boost_search_stats.evaluations() + w_lim);
	NTYPE c_new = boost_search.start2(*boost_state);
	wps_break_flag = ((c_new == e->val.val) && !boost_search.solved());
	e->val = CostNode::Value(c_new, boost_search.solved());
	if (wps_break_flag) e->val.work = w_lim;
	n_boost_searches += 1;

	if (trace_level > 1) {
	  if (wps_break_flag) {
	    std::cerr << *boost_state << " (cost = " << c_est
		      << ") discarded (wps = " << w_lim
		      << ", " << bstats << ")" << std::endl;
	  }
	  else {
	    std::cerr << *boost_state << ": cost " << c_est
		      << ", max " << boost_limit
		      << ", new " << e->val << " (" << bstats << ")"
		      << std::endl;
	  }
	}

//  	if (!boost_search.solved()) {
//  	  std::cerr << "cost improved, running quick boost pass..." << std::endl;
//  	  quick_boost_pass(list, search_space, cost_limit, cost_limit_set);
//  	}

	if (boost_search.solved()) {
	  e->val.opt = true;
	  delete e;
	  n_entries_solved += 1;
	}
	else if (wps_break_flag) {
	  n_entries_discarded += 1;
	  delete e;
	}
	else {
	  list = insert_entry(e, list);
	}
	delete boost_state;
      }

      else {
	n_entries_discarded += 1;
	delete e;
      }

      if (bstats.break_signal_raised()) {
	bstats.stop();
	count_type n_entries_out = (list ? list->list_length() : 0);
	n_entries_boosted += (n_entries_in - n_entries_out);
	return list;
      }
    }
  }

  bstats.stop();
  count_type n_entries_out = (list ? list->list_length() : 0);
  n_entries_boosted += (n_entries_in - n_entries_out);
  if (trace_level > 0) {
    std::cerr << "done (" << n_entries_solved << " solved, "
	      << n_entries_discarded << " discarded, "
	      << (list ? list->list_length() : 0) << " in list, "
	      << bstats.nodes() << " nodes, "
	      << bstats.time() << " sec.)" << std::endl;
  }
  return list;
}


void CostTableBoost::quick_boost_pass
(Entry* list, StateFactory& search_space,
 NTYPE cost_limit, const index_set& cost_limit_set)
{
  NoSearch ns;

  bstats.start();
  NTYPE current_bound = 0;
  count_type n_pushed = 0;

  Entry* e = list;
  while (e) {

    if (INFINITE(e->val.val)) {
      bstats.stop();
      if (trace_level > 0) {
	std::cerr << "done at INF (" << n_pushed << " pushed, "
		  << bstats.time() << " sec.)" << std::endl;
      }
      return;
    }

    if (e->val.val > cost_limit) {
      cost_limit = eval(cost_limit_set);
      if (e->val.val > cost_limit) {
	bstats.stop();
	if (trace_level > 0) {
	  std::cerr << "done at " << e->val.val
		    << " (" << n_pushed << " pushed, "
		    << bstats.time() << " sec.)" << std::endl;
	}
	return;
      }
      else if (trace_level > 0) {
	std::cerr << "limit is " << cost_limit << std::endl;
      }
    }

    if (trace_level > 0) {
      if (e->val.val > current_bound) {
	std::cerr << "at " << e->val.val
		  << " (" << n_pushed << " pushed, "
		  << bstats.time() << " sec.)" << std::endl;
	current_bound = e->val.val;
      }
    }

    if (!e->val.opt) {
      State* boost_state = search_space.new_state(e->set, 0);
      NTYPE c_est = boost_state->est_cost();
      NTYPE c_new = boost_state->expand(ns, c_est);
      if (c_new > c_est) {
	if (trace_level > 1) {
	  std::cerr << *boost_state << " pushed from "
		    << c_est << " to " << c_new << std::endl;
	}
	e->val = CostNode::Value(c_new, ns.solved());
	n_pushed += 1;
      }
      ns.reset();
      delete boost_state;
    }

    if (bstats.break_signal_raised()) {
      bstats.stop();
      return;
    }

    e = e->next;
  }

  bstats.stop();
  if (trace_level > 0) {
    std::cerr << "done (" << n_pushed << " pushed, "
	      << bstats.time() << " sec.)" << std::endl;
  }
}

#endif // BUILD_WITH_BOOSTING

void CostTable::compute_action_cost()
{
  if (!pre_cost) pre_cost = new NTYPE[instance.n_actions()];
  if (!per_cost) per_cost = new NTYPE[instance.n_actions()];
  for (index_type k = 0; k < instance.n_actions(); k++) {
    pre_cost[k] = eval(instance.actions[k].pre);
    index_set per_set(instance.actions[k].pre);
    per_set.subtract(instance.actions[k].del);
    per_cost[k] = eval(per_set);
  }
}

void CostTable::compute_hm_graph_min
(const ACF& cost,
 const index_set& s,
 bool strict,
 index_set_graph& g)
{
  index_type n = g.node_with_label(s);
  assert(n != no_such_index);
  std::cerr << "establishers for " << s << ":" << std::endl;
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (instance.actions[k].add.have_common_element(s) &&
	!instance.actions[k].del.have_common_element(s)) {
      index_set r(s);
      r.subtract(instance.actions[k].add);
      r.insert(instance.actions[k].pre);
      bool rel = true;
      if (strict) {
	if ((eval(r) + cost(k)) > eval(s)) rel = false;
      }
      std::cerr << k << ": " << r << " (" << rel << ")" << std::endl;
      if (rel) {
	index_type rn = g.node_with_label(r);
	if (rn == no_such_index) {
	  rn = g.add_node(r);
	  CostNode::Value* v = find(r);
	  if (v)
	    compute_hm_graph_min(cost, r, strict, g);
	  else
	    compute_hm_graph_max(cost, r, strict, g);
	}
	// std::cerr << "s = " << s << ", n = " << n << ", a = " << k
	// 		<< ", r = " << r << ", rn = " << rn << std::endl;
	if (!g.adjacent(rn, n)) {
	  g.add_edge(rn, n, EMPTYSET);
	}
	g.edge_label(rn, n).insert(k);
      }
    }
}

void CostTable::compute_hm_graph_max
(const ACF& cost,
 const index_set& s,
 bool strict,
 index_set_graph& g)
{
  index_type n = g.node_with_label(s);
  assert(n != no_such_index);
  index_set_vec ss;
  find_all(s, ss);
  // std::cerr << "find_all(" << s << ") = " << ss << std::endl;
  for (index_type k = 0; k < ss.length(); k++) {
    index_type sn = g.node_with_label(ss[k]);
    if (sn == no_such_index) {
      sn = g.add_node(ss[k]);
      compute_hm_graph_min(cost, ss[k], strict, g);
    }
    if (!g.adjacent(sn, n))
      g.add_edge(sn, n);
  }
}

void CostTable::compute_hm_graph
(const ACF& cost, const index_set& s, bool strict, index_set_graph& g)
{
  g.init(1);
  g.node_label(0) = s;
  CostNode::Value* v = find(s);
  if (v) {
    compute_hm_graph_min(cost, s, strict, g);
  }
  else {
    compute_hm_graph_max(cost, s, strict, g);
  }
}

void CostTable::compute_complete_hm_graph
(const ACF& cost, index_type m, const index_set& s, bool strict,
 index_set_graph& g)
{
  compute_hm_graph(cost, s, strict, g);
  mSubsetEnumerator e(instance.n_atoms(), m);
  bool more = e.first();
  while (more) {
    index_set ss(e.current_set());
    index_type ssn = g.node_with_label(ss);
    if (ssn == no_such_index) {
      ssn = g.add_node(ss);
      CostNode::Value* v = find(ss);
      if (v)
	compute_hm_graph_min(cost, ss, strict, g);
      else
	compute_hm_graph_max(cost, ss, strict, g);
    }
    more = e.next();
  }
}

bool CostTable::simple_cas
(const index_set& g, const ACF& cost, index_type m, index_set& cas)
{
  index_set_vec stk;
  if (g.length() <= m)
    return simple_cas_min(g, eval(g), cost, m, EMPTYSET, stk, cas);
  else
    return simple_cas_max(g, eval(g), cost, m, EMPTYSET, stk, cas);
}

bool CostTable::simple_cas_max
(const index_set& g, NTYPE c, const ACF& cost, index_type m,
 const index_set& cas_in, index_set_vec& stk, index_set& cas_out)
{
  assert(g.length() > m);
  std::cerr << stk << "/max: " << g << ", " << c << ", " << cas_in
	    << std::endl;
  if (c <= 0) {
    cas_out.assign_copy(cas_in);
    return true;
  }
  bool solved = false;
  mSubsetEnumerator e(g.length(), m);
  bool more = e.first();
  while (more) {
    index_set sg;
    e.current_set(g, sg);
    if ((eval(sg) >= c) && (stk.first(sg) == no_such_index)) {
      index_set sg_cas;
      bool ok = simple_cas_min(sg, c, cost, m, cas_in, stk, sg_cas);
      if (ok) {
	if (!solved) { // first solution
	  cas_out.assign_copy(sg_cas);
	  solved = true;
	}
	else if (sg_cas.length() < cas_out.length()) { // better solution
	  cas_out.assign_copy(sg_cas);
	}
      }
    }
    more = e.next();
  }
  std::cerr << stk << "/max: return " << solved << ", " << cas_out
	    << std::endl;
  return solved;
}

bool CostTable::simple_cas_min
(const index_set& g, NTYPE c, const ACF& cost, index_type m,
 const index_set& cas_in, index_set_vec& stk, index_set& cas_out)
{
  assert(g.length() <= m);
  std::cerr << stk << "/min: " << g << ", " << c << ", " << cas_in
	    << std::endl;
  cas_out.assign_copy(cas_in);
  if (c <= 0) {
    return true;
  }
  if (stk.first(g) != no_such_index) {
    std::cerr << stk << "/min: cycle, returning " << cas_out << std::endl;
    return true;
  }
  stk.append(g);
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (!instance.actions[k].del.have_common_element(g) &&
	instance.actions[k].add.have_common_element(g)) {
      std::cerr << stk << "/min: checking action " << k << std::endl;
      index_set rg(g);
      rg.subtract(instance.actions[k].add);
      rg.insert(instance.actions[k].pre);
      if (FINITE(eval(rg))) {
	if ((eval(rg) + cost(k)) < c) {
	  stk.dec_length();
	  std::cerr << stk << "/min: fail (eval)" << std::endl;
	  return false;
	}
	bool ok = ((c - cost(k)) <= 0);
	if (!ok) {
	  if (cost(k) > 0) cas_out.insert(k);
	  index_set new_cas;
	  if (rg.length() <= m)
	    ok = simple_cas_min(rg, c - cost(k), cost, m, cas_out, stk, new_cas);
	  else
	    ok = simple_cas_max(rg, c - cost(k), cost, m, cas_out, stk, new_cas);
	  if (ok) {
	    cas_out.assign_copy(new_cas);
	  }
	  else {
	    stk.dec_length();
	    std::cerr << stk << "/min: fail (recursive)" << std::endl;
	    return false;
	  }
	}
      }
    }
  stk.dec_length();
  std::cerr << stk << "/min: return " << cas_out << std::endl;
  return true;
}

bool CostTable::critical_tree
(const index_set& g, const ACF& cost, index_type m, index_set_vec& cta)
{
  index_set_vec ps;
  bool_vec ra(false, instance.n_actions());
  if (g.length() <= m)
    return critical_tree_min(g, eval(g), cost, m, ps, ra, cta);
  else
    return critical_tree_max(g, eval(g), cost, m, ps, ra, cta);
}

bool CostTable::critical_tree_max
(const index_set& g, NTYPE c, const ACF& cost, index_type m,
 index_set_vec& ps, bool_vec& ra, index_set_vec& cta)
{
  std::cerr << "enter ct-max (" << g << " at " << c << ")..." << std::endl;
  assert(g.length() > m);
  cta.clear();
  if (c <= 0) {
    std::cerr << "ct-max (" << g << " at " << c << "): solved by empty set"
	      << std::endl;
    cta.append(EMPTYSET);
    return true;
  }
  mSubsetEnumerator e(g.length(), m);
  bool found_max = false;
  bool more = e.first();
  while (more) {
    index_set sg;
    e.current_set(g, sg);
    if (eval(sg) >= c) {
      if (ps.first(sg) == no_such_index) {
	std::cerr << "ct-max: trying " << sg << " >= " << c << std::endl;
	index_set_vec sg_cta;
	bool ok = critical_tree_min(sg, c, cost, m, ps, ra, sg_cta);
	if (ok) {
	  for (index_type i = 0; i < sg_cta.length(); i++)
	    cta.append_if_new(sg_cta[i]);
	  found_max = true;
	}
      }
      else {
	std::cerr << "ct-max: skipping " << sg << " as its a cycle"
		  << std::endl;
      }
    }
    more = e.next();
  }
  std::cerr << "ct-max (" << g << " at " << c << "): solved = "
	    << found_max << ", " << cta.length() << " solutions: "
	    << cta << std::endl;
  return found_max;
}

bool CostTable::critical_tree_min
(const index_set& g, NTYPE c, const ACF& cost, index_type m,
 index_set_vec& ps, bool_vec& ra, index_set_vec& cta)
{
  std::cerr << "enter ct-min: " << g << " at " << c << std::endl;
  assert(g.length() <= m);
  cta.clear();
  cta.append(EMPTYSET);
  if (c <= 0) {
    std::cerr << "ct-min (" << g << " at " << c << "): solved by emptyset"
	      << std::endl;
    return true;
  }
  if (ps.first(g) != no_such_index) {
    std::cerr << "ct-min (" << g << " at " << c << "): solved by cycle"
	      << std::endl;
    return true;
  }
  ps.append(g);
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (!instance.actions[k].del.have_common_element(g) &&
	instance.actions[k].add.have_common_element(g)) {
      NTYPE cak = (ra[k] ? 0 : cost(k));
      index_set rg(g);
      rg.subtract(instance.actions[k].add);
      rg.insert(instance.actions[k].pre);
      std::cerr << "ct-min: action " << k << " -> " << rg << " >= "
		<< c - cak << std::endl;
      std::cerr << "ct-min: currently " << cta.length() << " solutions"
		<< std::endl;
      if (eval(rg) < (c - cak)) {
	std::cerr << "ct-min (" << g << " at " << c << "): "
		  << rg << " >= " << c - cak << " has no solution"
		  << std::endl;
	cta.clear();
	ps.dec_length(1);
	return false;
      }
      index_set_vec new_cta;
      index_set_vec rg_cta;
      bool ok;
      if (rg.length() <= m)
	ok = critical_tree_min(rg, c - cak, cost, m, ps, ra, rg_cta);
      else
	ok = critical_tree_max(rg, c - cak, cost, m, ps, ra, rg_cta);
      if (ok) {
	rg_cta.insert_in_all(k);
	new_cta.combinations_by_union(cta, rg_cta);
	new_cta.remove_duplicate_elements();
      }
      std::cerr << "ct-min: " << rg << " >= " << c - cak << ": "
		<< new_cta.length() << " solutions not relaxing action "
		<< k << std::endl;
      // can try to relax action k
      if ((eval(rg) - cost(k)) >= c) {
	bool nec_k = true;
	for (index_type i = 0; (i < cta.length()) && nec_k; i++)
	  if (!cta[i].contains(k)) nec_k = false;
	if (!nec_k) {
	  std::cerr << "ct-min: trying to relax action " << k << std::endl;
	  // ra[k] = true;
	  rg_cta.clear();
	  if (rg.length() <= m)
	    ok = critical_tree_min(rg, c, cost, m, ps, ra, rg_cta);
	  else
	    ok = critical_tree_max(rg, c, cost, m, ps, ra, rg_cta);
	  if (ok) {
	    for (index_type i = 0; i < cta.length(); i++)
	      if (!cta[i].contains(k))
		for (index_type j = 0; j < rg_cta.length(); j++) {
		  new_cta.append(cta[i]);
		  new_cta[new_cta.length() - 1].insert(rg_cta[j]);
		}
	  }
	  new_cta.remove_duplicate_elements();
	  // ra[k] = false;
	  std::cerr << "ct-min " << rg << " >= " << c - cak << ": "
		    << new_cta.length() << " solutions relaxing action "
		    << k << std::endl;
	}
      }
      cta.assign_copy(new_cta);
      if (cta.length() == 0) {
	ps.dec_length();
	return false;
      }
    }
  ps.dec_length(1);
  std::cerr << "ct-min (" << g << " at " << c << "): returning "
	    << cta.length() << " solutions: " << cta << std::endl;
  return true;
}

LinearScanH2Eval::LinearScanH2Eval(const Instance& i, Stopwatch& s)
  : CostTable(i, s), n_pairs(0), pairs(0), values(0), complete(false)
{
  // done
}

LinearScanH2Eval::~LinearScanH2Eval()
{
  if (pairs)
    delete[] pairs;
  if (values)
    delete[] values;
}

void LinearScanH2Eval::compile_finite()
{
  if (pairs)
    delete[] pairs;
  if (values)
    delete[] values;
  n_pairs = 0;
  complete = true;
  Entry* list = entries(true, false, true, false);
  if (list == 0) return;
  for (Entry* e = list; e != 0; e = e->next)
    if (e->set.length() <= 2)
      n_pairs += 1;
    else
      complete = false;
  if (n_pairs == 0) return;
  pairs = new index_pair[n_pairs];
  values = new NTYPE[n_pairs];
  index_type i = 0;
  Entry* last = list->last();
  while (last != 0) {
    assert(i < n_pairs);
    if (last->set.length() == 1) {
      pairs[i].first = last->set[0];
      pairs[i].second = last->set[0];
      values[i] = last->val.val;
      i += 1;
      last = last->prev;
    }
    else if (last->set.length() == 2) {
      pairs[i].first = last->set[0];
      pairs[i].second = last->set[1];
      values[i] = last->val.val;
      i += 1;
      last = last->prev;
    }
    else {
      Entry* p = last->prev;
      last->unlink();
      delete last;
      last = p;
    }
  }
  assert(i == n_pairs);
  list->delete_list();
  std::cerr << "h^2 compiled: " << n_pairs << " pairs, complete = "
	    << complete << std::endl;
}

NTYPE LinearScanH2Eval::eval(const bool_vec& s)
{
  for (index_type k = 0; k < n_pairs; k++)
    if (s[pairs[k].first] && s[pairs[k].second])
      return values[k];
  if (complete)
    return 0;
  else
    return CostNode::eval(s);
}


ForwardH1::ForwardH1
(const Instance& i, const index_set& g, const ACF& c, Stopwatch& s)
  : Heuristic(i), goals(g), cost(c), table(0)
{
  table = new CostTable(i, s);
}

ForwardH1::~ForwardH1()
{
  delete table;
}

NTYPE ForwardH1::eval(const index_set& s)
{
  bool_vec s1(s, instance.n_atoms());
  table->compute_H1(cost, s1);
  return table->eval(goals);
}

NTYPE ForwardH1::eval(const bool_vec& s)
{
  table->compute_H1(cost, s);
  return table->eval(goals);
}

ForwardH2::ForwardH2
(const Instance& i, const index_set& g, const ACF& c, Stopwatch& s)
  : Heuristic(i), goals(g), cost(c), table(0)
{
  table = new CostTable(instance, s);
}

ForwardH2::~ForwardH2()
{
  delete table;
}

NTYPE ForwardH2::eval(const index_set& s)
{
  bool_vec s1(s, instance.n_atoms());
  table->compute_H2(cost, s1);
  return table->eval(goals);
}

NTYPE ForwardH2::eval(const bool_vec& s)
{
  table->compute_H2(cost, s);
  return table->eval(goals);
}

END_HSPS_NAMESPACE
