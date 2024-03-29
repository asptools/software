
#define PRINT_COMPARE_EVAL

#include "heuristic.h"

BEGIN_HSPS_NAMESPACE

count_type Heuristic::eval_count = 0;
int Heuristic::default_trace_level = 1;

Heuristic::~Heuristic()
{
  // done
}

void Heuristic::set_trace_level(int level)
{
  trace_level = level;
}

extended_cost Heuristic::extended_eval(const index_set& s)
{
  return extended_cost(eval(s), false);
}

extended_cost Heuristic::extended_eval(const bool_vec& s)
{
  return extended_cost(eval(s), false);
}

void Heuristic::write_eval
(const index_set& s, ::std::ostream& st, char* p, bool e)
{
  if (p) st << p << ' ';
  st << eval(s);
  if (e) st << ::std::endl;
}

void Heuristic::write_eval
(const bool_vec& s, ::std::ostream& st, char* p, bool e)
{
  if (p) st << p << ' ';
  st << eval(s);
  if (e) st << ::std::endl;
}

NTYPE Heuristic::eval_precondition(const Instance::Action& a)
{
  return eval(a.pre);
}

NTYPE Heuristic::incremental_eval(const index_set& s, index_type i_new)
{
  index_set s_new(s);
  s_new.insert(i_new);
  return eval(s_new);
}

NTYPE Heuristic::incremental_eval(const bool_vec& s, index_type i_new)
{
  bool_vec s_new(s);
  for (index_type k = 0; k < instance.n_atoms(); k++) s_new[k] = s[k];
  s_new[i_new] = true;
  return eval(s_new);
}

NTYPE Heuristic::eval_to_bound(const index_set& s, NTYPE bound)
{
  return eval(s);
}

NTYPE Heuristic::eval_to_bound(const bool_vec& s, NTYPE bound)
{
  return eval(s);
}

void Heuristic::store(const index_set& s, NTYPE v, bool opt)
{
  // does nothing
}

void Heuristic::store(const bool_vec& s, NTYPE v, bool opt)
{
  // does nothing
}

NTYPE Heuristic::eval(index_type atom)
{
  index_set s;
  s.assign_singleton(atom);
  return eval(s);
}

NTYPE Heuristic::eval(index_type atom1, index_type atom2)
{
  index_set s;
  s.assign_singleton(atom1);
  s.insert(atom2);
  return eval(s);
}

void Heuristic::compute_heuristic_graph(const ACF& cost, graph& g)
{
  g.init(instance.n_atoms());
  for (index_type k = 0; k < instance.n_actions(); k++) {
    NTYPE c = eval(instance.actions[k].pre) + cost(k);
    for (index_type i = 0; i < instance.actions[k].add.length(); i++)
      if (c <= eval(instance.actions[k].add[i]))
	g.add_edge(instance.actions[k].pre, instance.actions[k].add[i]);
  }
}

NTYPE ZeroHeuristic::eval(const index_set& s)
{
  eval_count += 1;
  return 0;
}

NTYPE ZeroHeuristic::eval(const bool_vec& s)
{
  eval_count += 1;
  return 0;
}

EvalActionCache::EvalActionCache(Instance& ins, Heuristic& h)
  : Heuristic(ins), base_h(h), cache(0, instance.n_actions())
{
  for (index_type k = 0; k < instance.n_actions(); k++)
    cache[k] = base_h.eval(instance.actions[k].pre);
}

NTYPE EvalActionCache::eval(const index_set& s)
{
  return base_h.eval(s);
}

NTYPE EvalActionCache::eval(const bool_vec& s)
{
  return base_h.eval(s);
}

NTYPE EvalActionCache::eval_precondition(const Instance::Action& a)
{
  return cache[a.index];
}

NTYPE EvalActionCache::incremental_eval(const index_set& s, index_type i_new)
{
  return base_h.incremental_eval(s, i_new);
}

NTYPE EvalActionCache::incremental_eval(const bool_vec& s, index_type i_new)
{
  return base_h.incremental_eval(s, i_new);
}

NTYPE EvalActionCache::eval_to_bound(const index_set& s, NTYPE bound)
{
  return base_h.eval_to_bound(s, bound);
}

NTYPE EvalActionCache::eval_to_bound(const bool_vec& s, NTYPE bound)
{
  return base_h.eval_to_bound(s, bound);
}

NTYPE RegressionInvariantCheck::eval(const index_set& s)
{
  NTYPE val = base_h.eval(s);
  if (FINITE(val)) {
    for (index_type k = 0; (k < instance.n_invariants()) && FINITE(val); k++)
      if (instance.invariants[k].verified || !verified_invariants_only)
	if (!instance.eval_invariant_in_partial_state(s, instance.invariants[k])) val = POS_INF;
  }
  return val;
}

NTYPE RegressionInvariantCheck::eval(const bool_vec& s)
{
  NTYPE val = base_h.eval(s);
  if (FINITE(val)) {
    for (index_type k = 0; (k < instance.n_invariants()) && FINITE(val); k++)
      if (instance.invariants[k].verified || !verified_invariants_only)
	if (!instance.eval_invariant_in_partial_state(s, instance.invariants[k])) val = POS_INF;
  }
  return val;
}

NTYPE RegressionInvariantCheck::incremental_eval
(const index_set& s, index_type i_new)
{
  NTYPE val = base_h.incremental_eval(s, i_new);
  if (FINITE(val)) {
    index_set s_new(s);
    s_new.insert(i_new);
    for (index_type k = 0; (k < instance.n_invariants()) && FINITE(val); k++)
      if (instance.invariants[k].verified || !verified_invariants_only)
	if (!instance.eval_invariant_in_partial_state(s_new, instance.invariants[k])) val = POS_INF;
  }
  return val;
}

NTYPE RegressionInvariantCheck::incremental_eval
(const bool_vec& s, index_type i_new)
{
  NTYPE val = base_h.incremental_eval(s, i_new);
  if (FINITE(val)) {
    bool_vec s_new(s);
    s_new[i_new] = true;
    for (index_type k = 0; (k < instance.n_invariants()) && FINITE(val); k++)
      if (instance.invariants[k].verified || !verified_invariants_only)
	if (!instance.eval_invariant_in_partial_state(s_new, instance.invariants[k])) val = POS_INF;
  }
  return val;
}

NTYPE RegressionInvariantCheck::eval_to_bound
(const index_set& s, NTYPE bound)
{
  NTYPE val = base_h.eval(s);
  if (val < bound) {
    for (index_type k = 0; (k < instance.n_invariants()) && FINITE(val); k++)
      if (instance.invariants[k].verified || !verified_invariants_only)
	if (!instance.eval_invariant_in_partial_state(s, instance.invariants[k])) val = POS_INF;
  }
  return val;
}

NTYPE RegressionInvariantCheck::eval_to_bound
(const bool_vec& s, NTYPE bound)
{
  NTYPE val = base_h.eval(s);
  if (val < bound) {
    for (index_type k = 0; (k < instance.n_invariants()) && FINITE(val); k++)
      if (instance.invariants[k].verified || !verified_invariants_only)
	if (!instance.eval_invariant_in_partial_state(s, instance.invariants[k])) val = POS_INF;
  }
  return val;
}

WeightedHeuristic::WeightedHeuristic(Instance& ins, Heuristic& h, NTYPE w)
  : Heuristic(ins), base_h(h), weight(w)
{
  // done
}

void WeightedHeuristic::set_weight(NTYPE w)
{
  weight = w;
}

NTYPE WeightedHeuristic::eval(const index_set& s)
{
  return (weight * base_h.eval(s));
}

NTYPE WeightedHeuristic::eval(const bool_vec& s)
{
  return (weight * base_h.eval(s));
}

NTYPE WeightedHeuristic::eval_precondition(const Instance::Action& a)
{
  return (weight * base_h.eval_precondition(a));
}

NTYPE WeightedHeuristic::incremental_eval(const index_set& s, index_type i_new)
{
  return (weight * base_h.incremental_eval(s, i_new));
}

NTYPE WeightedHeuristic::incremental_eval(const bool_vec& s, index_type i_new)
{
  return (weight * base_h.incremental_eval(s, i_new));
}

NTYPE WeightedHeuristic::eval_to_bound(const index_set& s, NTYPE bound)
{
  return (weight * base_h.eval_to_bound(s, bound));
}

NTYPE WeightedHeuristic::eval_to_bound(const bool_vec& s, NTYPE bound)
{
  return (weight * base_h.eval_to_bound(s, bound));
}

NTYPE RoundUp::eval(const index_set& s)
{
  NTYPE v = h.eval(s);
  return CEIL_TO(v, d);
  // if (FRAC(v) > 0)
  //   return FLOOR(v) + 1;
  // else
  //   return v;
}

NTYPE RoundUp::eval(const bool_vec& s)
{
  NTYPE v = h.eval(s);
  return CEIL_TO(v, d);
  // if (FRAC(v) > 0)
  //   return FLOOR(v) + 1;
  // else
  //   return v;
}

NTYPE RoundUp::incremental_eval(const index_set& s, index_type i_new)
{
  NTYPE v = h.incremental_eval(s, i_new);
  return CEIL_TO(v, d);
  // if (FRAC(v) > 0)
  //   return FLOOR(v) + 1;
  // else
  //   return v;
}

NTYPE RoundUp::incremental_eval(const bool_vec& s, index_type i_new)
{
  NTYPE v = h.incremental_eval(s, i_new);
  return CEIL_TO(v, d);
  // if (FRAC(v) > 0)
  //   return FLOOR(v) + 1;
  // else
  //   return v;
}

NTYPE RoundUp::eval_to_bound(const index_set& s, NTYPE bound)
{
  NTYPE v = h.eval_to_bound(s, bound);
  return CEIL_TO(v, d);
  // if (FRAC(v) > 0)
  //   return FLOOR(v) + 1;
  // else
  //   return v;
}

NTYPE RoundUp::eval_to_bound(const bool_vec& s, NTYPE bound)
{
  NTYPE v = h.eval_to_bound(s, bound);
  return CEIL_TO(v, d);
  // if (FRAC(v) > 0)
  //   return FLOOR(v) + 1;
  // else
  //   return v;
}

NTYPE Combine2ByMax::eval(const index_set& s)
{
  NTYPE v0 = h0.eval(s);
  NTYPE v1 = h1.eval(s);
  if (trace_level > 1) {
    if (v0 != v1) {
      ::std::cerr << "Combine2ByMax: ";
      instance.write_atom_set(::std::cerr, s);
      ::std::cerr << ", v0 = " << v0 << ", v1 = " << v1 << ::std::endl;
    }
  }
  return MAX(v0, v1);
}

NTYPE Combine2ByMax::eval(const bool_vec& s)
{
  NTYPE v0 = h0.eval(s);
  NTYPE v1 = h1.eval(s);
  if (trace_level > 1) {
    if (v0 != v1) {
      ::std::cerr << "Combine2ByMax: ";
      instance.write_atom_set(::std::cerr, s);
      ::std::cerr << ", v0 = " << v0 << ", v1 = " << v1 << ::std::endl;
    }
  }
  return MAX(v0, v1);
}

NTYPE Combine2ByMax::incremental_eval(const index_set& s, index_type i_new)
{
  NTYPE v0 = h0.incremental_eval(s, i_new);
  NTYPE v1 = h1.incremental_eval(s, i_new);
  return MAX(v0, v1);
}

NTYPE Combine2ByMax::incremental_eval(const bool_vec& s, index_type i_new)
{
  NTYPE v0 = h0.incremental_eval(s, i_new);
  NTYPE v1 = h1.incremental_eval(s, i_new);
  return MAX(v0, v1);
}

NTYPE Combine2ByMax::eval_to_bound(const index_set& s, NTYPE bound)
{
  NTYPE v0 = h0.eval_to_bound(s, bound);
  if (v0 > bound) return v0;
  NTYPE v1 = h1.eval_to_bound(s, bound);
  return MAX(v0, v1);
}

NTYPE Combine2ByMax::eval_to_bound(const bool_vec& s, NTYPE bound)
{
  NTYPE v0 = h0.eval_to_bound(s, bound);
  if (v0 > bound) return v0;
  NTYPE v1 = h1.eval_to_bound(s, bound);
  return MAX(v0, v1);
}

NTYPE CombineNByMax::eval(const index_set& s)
{
  NTYPE v_max = 0;
  for (index_type k = 0; k < h_vec.length(); k++)
    v_max = MAX(v_max, h_vec[k]->eval(s));
  return v_max;
}

NTYPE CombineNByMax::eval(const bool_vec& s)
{
  NTYPE v_max = 0;
  for (index_type k = 0; k < h_vec.length(); k++)
    v_max = MAX(v_max, h_vec[k]->eval(s));
  return v_max;
}

NTYPE CombineNByMax::incremental_eval(const index_set& s, index_type i_new)
{
  NTYPE v_max = 0;
  for (index_type k = 0; k < h_vec.length(); k++)
    v_max = MAX(v_max, h_vec[k]->incremental_eval(s, i_new));
  return v_max;
}

NTYPE CombineNByMax::incremental_eval(const bool_vec& s, index_type i_new)
{
  NTYPE v_max = 0;
  for (index_type k = 0; k < h_vec.length(); k++)
    v_max = MAX(v_max, h_vec[k]->incremental_eval(s, i_new));
  return v_max;
}

NTYPE CombineNByMax::eval_to_bound(const index_set& s, NTYPE bound)
{
  NTYPE v_max = 0;
  for (index_type k = 0; k < h_vec.length(); k++) {
    v_max = MAX(v_max, h_vec[k]->eval(s));
    if (v_max > bound) return v_max;
  }
  return v_max;
}

NTYPE CombineNByMax::eval_to_bound(const bool_vec& s, NTYPE bound)
{
  NTYPE v_max = 0;
  for (index_type k = 0; k < h_vec.length(); k++) {
    v_max = MAX(v_max, h_vec[k]->eval(s));
    if (v_max > bound) return v_max;
  }
  return v_max;
}

NTYPE CombineNBySum::eval(const index_set& s)
{
  NTYPE v_sum = 0;
  for (index_type k = 0; k < h_vec.length(); k++)
    v_sum += h_vec[k]->eval(s);
  return v_sum;
}

NTYPE CombineNBySum::eval(const bool_vec& s)
{
  NTYPE v_sum = 0;
  for (index_type k = 0; k < h_vec.length(); k++)
    v_sum += h_vec[k]->eval(s);
  return v_sum;
}

NTYPE CombineNBySum::incremental_eval(const index_set& s, index_type i_new)
{
  NTYPE v_sum = 0;
  for (index_type k = 0; k < h_vec.length(); k++)
    v_sum += h_vec[k]->incremental_eval(s, i_new);
  return v_sum;
}

NTYPE CombineNBySum::incremental_eval(const bool_vec& s, index_type i_new)
{
  NTYPE v_sum = 0;
  for (index_type k = 0; k < h_vec.length(); k++)
    v_sum += h_vec[k]->incremental_eval(s, i_new);
  return v_sum;
}

NTYPE CombineNBySum::eval_to_bound(const index_set& s, NTYPE bound)
{
  NTYPE v_sum = 0;
  for (index_type k = 0; k < h_vec.length(); k++) {
    v_sum += h_vec[k]->eval(s);
    if (v_sum > bound) return v_sum;
  }
  return v_sum;
}

NTYPE CombineNBySum::eval_to_bound(const bool_vec& s, NTYPE bound)
{
  NTYPE v_sum = 0;
  for (index_type k = 0; k < h_vec.length(); k++) {
    v_sum += h_vec[k]->eval(s);
    if (v_sum > bound) return v_sum;
  }
  return v_sum;
}


// NTYPE Combine2ByRandomChoice::eval(const index_set& s)
// {
//   unsigned long i = rng.random_in_range(alpha.divisor());
//   if (i < alpha.numerator())
//     return h0.eval(s);
//   else
//     return h1.eval(s);
// }
// 
// NTYPE Combine2ByRandomChoice::eval(const bool_vec& s)
// {
//   unsigned long i = rng.random_in_range(alpha.divisor());
//   if (i < alpha.numerator())
//     return h0.eval(s);
//   else
//     return h1.eval(s);
// }

NTYPE HX::eval(const index_set& s)
{
  assert(!X.empty());
  index_type c = s.first_common_element(X);
  if (c == no_such_index) {
    bool_vec sx(s, instance.n_atoms());
    NTYPE v = POS_INF;
    for (index_type k = 0; k < X.length(); k++) {
      assert(!sx[X[k]]);
      sx[X[k]] = true;
      v = MIN(v, h0.eval(sx));
      sx[X[k]] = false;
    }
    return v;
  }
  else {
    return h0.eval(s);
  }
}

NTYPE HX::eval(const bool_vec& s)
{
  assert(!X.empty());
  index_type c = X.first_common_element(s);
  if (c == no_such_index) {
    bool_vec sx(s);
    NTYPE v = POS_INF;
    for (index_type k = 0; k < X.length(); k++) {
      assert(!sx[X[k]]);
      sx[X[k]] = true;
      v = MIN(v, h0.eval(sx));
      sx[X[k]] = false;
    }
    return v;
  }
  else {
    return h0.eval(s);
  }
}

NTYPE HX::incremental_eval(const index_set& s, index_type i_new)
{
  assert(!X.empty());
  if (X.contains(i_new)) {
    return h0.incremental_eval(s, i_new);
  }
  else {
    index_type c = s.first_common_element(X);
    if (c == no_such_index) {
      bool_vec sx(s, instance.n_atoms());
      sx[i_new] = true;
      NTYPE v = POS_INF;
      for (index_type k = 0; k < X.length(); k++) {
	assert(!sx[X[k]]);
	sx[X[k]] = true;
	v = MIN(v, h0.eval(sx));
	sx[X[k]] = false;
      }
      return v;
    }
    else {
      return h0.incremental_eval(s, i_new);
    }
  }
}

NTYPE HX::incremental_eval(const bool_vec& s, index_type i_new)
{
  if (X.contains(i_new)) {
    return h0.incremental_eval(s, i_new);
  }
  else {
    assert(!X.empty());
    index_type c = X.first_common_element(s);
    if (c == no_such_index) {
      bool_vec sx(s);
      bool i_new_set = sx[i_new];
      sx[i_new] = true;
      NTYPE v = POS_INF;
      for (index_type k = 0; k < X.length(); k++) {
	assert(!sx[X[k]]);
	sx[X[k]] = true;
	v = MIN(v, h0.eval(sx));
	sx[X[k]] = false;
      }
      sx[i_new] = i_new_set;
      return v;
    }
    else {
      return h0.incremental_eval(s, i_new);
    }
  }
}

void HX::write_eval(const index_set& s, std::ostream& st, char* p, bool e)
{
  assert(!X.empty());
  if (p) st << p << ' ';
  index_type c = s.first_common_element(X);
  if (c == no_such_index) {
    st << "max {";
    bool_vec sx(s, instance.n_atoms());
    NTYPE v = POS_INF;
    for (index_type k = 0; k < X.length(); k++) {
      assert(!sx[X[k]]);
      if (k > 0) st << ", ";
      st << "X=" << instance.atoms[X[k]].name << ": ";
      sx[X[k]] = true;
      h0.write_eval(s, st, 0, false);
      v = MIN(v, h0.eval(sx));
      sx[X[k]] = false;
    }
    st << "} = " << v;
  }
  else {
    st << "X=" << instance.atoms[c].name << ": ";
    h0.write_eval(s, st, 0, false);
  }
  if (e) st << ::std::endl;
}

void HX::write_eval(const bool_vec& s, std::ostream& st, char* p, bool e)
{
  assert(!X.empty());
  if (p) st << p << ' ';
  index_type c = X.first_common_element(s);
  if (c == no_such_index) {
    st << "max {";
    bool_vec sx(s);
    NTYPE v = POS_INF;
    for (index_type k = 0; k < X.length(); k++) {
      assert(!sx[X[k]]);
      if (k > 0) st << ", ";
      st << "X=" << instance.atoms[X[k]].name << ": ";
      sx[X[k]] = true;
      h0.write_eval(s, st, 0, false);
      v = MIN(v, h0.eval(sx));
      sx[X[k]] = false;
    }
    st << "} = " << v;
  }
  else {
    st << "X=" << instance.atoms[c].name << ": ";
    h0.write_eval(s, st, 0, false);
  }
  if (e) st << ::std::endl;
}

NTYPE AtomMapAdapter::eval(const index_set& s)
{
  index_set xs;
  for (index_type k = 0; k < s.length(); k++)
    if (map[s[k]] != no_such_index) xs.insert(map[s[k]]);
  return base_h.eval(xs);
}

NTYPE AtomMapAdapter::eval(const bool_vec& s)
{
  index_set xs;
  for (index_type k = 0; k < instance.n_atoms(); k++)
    if (s[k] && (map[k] != no_such_index)) xs.insert(map[k]);
  return base_h.eval(xs);
}

NTYPE AtomMapAdapter::incremental_eval(const index_set& s, index_type i_new)
{
  index_set xs;
  for (index_type k = 0; k < s.length(); k++)
    if (map[s[k]] != no_such_index) xs.insert(map[s[k]]);
  if (map[i_new] != no_such_index) 
    return base_h.incremental_eval(xs, map[i_new]);
  else
    return base_h.eval(xs);
}

NTYPE AtomMapAdapter::incremental_eval(const bool_vec& s, index_type i_new)
{
  index_set xs;
  for (index_type k = 0; k < instance.n_atoms(); k++)
    if (s[k] && (map[k] != no_such_index)) xs.insert(map[k]);
  if (map[i_new] != no_such_index) 
    return base_h.incremental_eval(xs, map[i_new]);
  else
    return base_h.eval(xs);
}

NTYPE ToP2Adapter::eval(const index_set& s)
{
  bool_vec s2(false, map.n_pairs());
  for (index_type i = 0; i < s.size(); i++) {
    s2[map(s[i],s[i])] = true;
    for (index_type j = i + 1; j < s.size(); j++)
      s2[map(s[i],s[j])] = true;
  }
  NTYPE v = base_h.eval(s2);
  return v;
}

NTYPE ToP2Adapter::eval(const bool_vec& s)
{
  bool_vec s2(false, map.n_pairs());
  for (index_type i = 0; i < instance.n_atoms(); i++)
    if (s[i]) {
      s2[map(i,i)] = true;
      for (index_type j = i + 1; j < instance.n_atoms(); j++)
	if (s[j])
	  s2[map(i,j)] = true;
    }
  NTYPE v = base_h.eval(s2);
  return v;
}

count_type CompareEval::lower = 0;
count_type CompareEval::equal = 0;
count_type CompareEval::higher = 0;

NTYPE CompareEval::eval(const index_set& s)
{
  NTYPE base_val = base_h.eval(s);
  NTYPE alt_val = alt_h.eval(s);
  if (alt_val < base_val) {
    lower += 1;
  }
  else if (alt_val > base_val) {
    higher += 1;
  }
  else {
    equal += 1;
  }
#ifdef PRINT_COMPARE_EVAL
  ::std::cout << "COMPARE ";
  instance.write_atom_set(::std::cout, s);
  ::std::cout << " : " << base_val << " / " << alt_val << ::std::endl;
#endif
  if ((trace_level > 0) && (alt_val != base_val)) {
    base_h.set_trace_level(trace_level);
    ::std::cerr << "evaluating base heuristic..." << ::std::endl;
    base_h.eval(s);
    base_h.set_trace_level(0);
    alt_h.set_trace_level(trace_level);
    ::std::cerr << "evaluating alternative heuristic..." << ::std::endl;
    alt_h.eval(s);
    alt_h.set_trace_level(0);
  }
  if (max_h_val) {
    return MAX(base_val, alt_val);
  }
  else {
    return base_val;
  }
}

NTYPE CompareEval::eval(const bool_vec& s)
{
  NTYPE base_val = base_h.eval(s);
  NTYPE alt_val = alt_h.eval(s);
  if (alt_val < base_val) {
    lower += 1;
  }
  else if (alt_val > base_val) {
    higher += 1;
  }
  else {
    equal += 1;
  }
#ifdef PRINT_COMPARE_EVAL
  ::std::cout << "COMPARE ";
  instance.write_atom_set(::std::cout, s);
  ::std::cout << " : " << base_val << " / " << alt_val << ::std::endl;
#endif
  if ((trace_level > 0) && (alt_val != base_val)) {
    base_h.set_trace_level(trace_level);
    ::std::cerr << "evaluating base heuristic..." << ::std::endl;
    base_h.eval(s);
    base_h.set_trace_level(0);
    alt_h.set_trace_level(trace_level);
    ::std::cerr << "evaluating alternative heuristic..." << ::std::endl;
    alt_h.eval(s);
    alt_h.set_trace_level(0);
  }
  if (max_h_val) {
    return MAX(base_val, alt_val);
  }
  else {
    return base_val;
  }
}

CompleteNegationAdapter::CompleteNegationAdapter
(Instance& ins, const pair_vec& p, Heuristic& h)
  : Heuristic(ins), h_base(h), pn_map(p), sc(false, p.length() * 2)
{
  // done
}

CompleteNegationAdapter::~CompleteNegationAdapter()
{
  // done
}

NTYPE CompleteNegationAdapter::eval(const index_set& s)
{
  bool_vec sv(s, instance.n_atoms());
  return eval(sv);
}

NTYPE CompleteNegationAdapter::eval(const bool_vec& s)
{
  for (index_type k = 0; k < pn_map.length(); k++) {
    sc[pn_map[k].first] = s[pn_map[k].first];
    sc[pn_map[k].second] = !s[pn_map[k].first];
  }
  return h_base.eval(sc);
}

NTYPE ACF::min_cost(index_type n) const
{
  NTYPE c_min = POS_INF;
  for (index_type k = 0; k < n; k++)
    c_min = MIN(c_min, (*this)(k));
  return c_min;
}

NTYPE ACF::max_cost(index_type n) const
{
  NTYPE c_max = NEG_INF;
  for (index_type k = 0; k < n; k++)
    c_max = MAX(c_max, (*this)(k));
  return c_max;
}

NTYPE ACF::avg_cost(index_type n) const
{
  // std::cerr << "computing average cost..." << std::endl;
  NTYPE c_sum = 0;
  for (index_type k = 0; k < n; k++) {
    NTYPE c_k = (*this)(k);
    // std::cerr << "cost(" << k << ") = " << c_k << std::endl;
    c_sum += c_k;
  }
  // std::cerr << "returning " << c_sum << " / " << n << " = "
  //	    << (c_sum / n) << std::endl;
  return c_sum / n;
}

NTYPE ACF::sum(const index_set& s) const
{
  NTYPE t = 0;
  for (index_type i = 0; i < s.size(); i++)
    t += (*this)(s[i]);
  return t;
}

NTYPE ACF::cost_gcd(index_type n) const
{
#ifdef NTYPE_RATIONAL
  bool first = true;
  NTYPE e = 0;
  for (index_type k = 0; k < n; k++) {
    NTYPE ck = (*this)(k);
    if (ck != 0) {
      if (first) {
	e = ck;
	first = false;
      }
      else {
	e = rational::rgcd(e, ck);
      }
    }
  }
  return e;
#else
  std::cerr << "error: gcd calculation not implemented for NTYPE_FLOAT"
	    << std::endl;
  exit(255);
#endif
}

NTYPE ACF::min_cost(const index_set& s) const
{
  NTYPE m = POS_INF;
  for (index_type k = 0; k < s.size(); k++)
    m = MIN(m, (*this)(s[k]));
  return m;
}

NTYPE ACF::max_cost(const index_set& s) const
{
  NTYPE m = NEG_INF;
  for (index_type k = 0; k < s.size(); k++)
    m = MAX(m, (*this)(s[k]));
  return m;
}

index_type ACF::min_cost_action(const index_set& s) const
{
  assert(!s.empty());
  index_type i_min = 0;
  NTYPE c_min = (*this)(s[0]);
  for (index_type k = 1; k < s.size(); k++)
    if ((*this)(s[k]) < c_min) {
      i_min = k;
      c_min = (*this)(s[k]);
    }
  return s[i_min];
}

index_type ACF::max_cost_action(const index_set& s) const
{
  assert(!s.empty());
  index_type i_max = 0;
  NTYPE c_max = (*this)(s[0]);
  for (index_type k = 1; k < s.size(); k++)
    if ((*this)(s[k]) > c_max) {
      i_max = k;
      c_max = (*this)(s[k]);
    }
  return s[i_max];
}

NTYPE UnitACF::operator()(index_type a) const
{
  return 1;
}

NTYPE UnitACF::min_cost(index_type n) const
{
  return 1;
}

NTYPE UnitACF::max_cost(index_type n) const
{
  return 1;
}

NTYPE UnitACF::avg_cost(index_type n) const
{
  return 1;
}

AnyACF::AnyACF(index_type n)
  : costs(0, n)
{
  // done
}

AnyACF::AnyACF(index_type n, NTYPE c)
  : costs(c, n)
{
  // done
}

AnyACF::AnyACF(index_type n, const ACF& c)
  : costs(0, n)
{
  for (index_type k = 0; k < n; k++)
    costs[k] = c(k);
}

NTYPE AnyACF::min_cost(const bool_vec& acts) const
{
  assert(acts.size() <= costs.size());
  NTYPE c_min = POS_INF;
  for (index_type k = 0; k < acts.size(); k++)
    if (acts[k])
      c_min = MIN(c_min, costs[k]);
  return c_min;
}

void AnyACF::decrease(const bool_vec& acts, NTYPE v)
{
  assert(acts.size() <= costs.size());
  for (index_type k = 0; k < acts.size(); k++)
    if (acts[k]) {
      assert(costs[k] >= v);
      costs[k] = costs[k] - v;
    }
}

void AnyACF::decrease(const index_set& acts, NTYPE v)
{
  for (index_type k = 0; k < acts.size(); k++) {
    assert(acts[k] < costs.size());
    assert(costs[acts[k]] >= v);
    costs[acts[k]] = costs[acts[k]] - v;
  }
}

void AnyACF::set_cost(index_type a, NTYPE v)
{
  assert(a < costs.length());
  costs[a] = v;
}

void AnyACF::extend_to(index_type n_new, NTYPE v)
{
  if (n_new < costs.length()) return;
  index_type n_old = costs.length();
  costs.set_length(n_new);
  for (index_type k = n_old; k < n_new; k++)
    costs[k] = v;
}

NTYPE AnyACF::operator()(index_type a) const
{
  assert(a < costs.size());
  return costs[a];
}

NTYPE AnyACF::min_cost(index_type n) const
{
  assert(n == costs.size());
  return cost_vec_util::min(costs);
}

NTYPE AnyACF::max_cost(index_type n) const
{
  assert(n == costs.size());
  return cost_vec_util::max(costs);
}

NTYPE AnyACF::avg_cost(index_type n) const
{
  assert(n == costs.size());
  return (cost_vec_util::sum(costs) / costs.size());
}

NTYPE ZeroACF::operator()(index_type a) const
{
  return 0;
}

NTYPE ZeroACF::min_cost(index_type n) const
{
  return 0;
}

NTYPE ZeroACF::max_cost(index_type n) const
{
  return 0;
}

NTYPE ZeroACF::avg_cost(index_type n) const
{
  return 0;
}

NTYPE CostACF::operator()(index_type a) const
{
  return instance.actions[a].cost;
}

NTYPE PlusOne::operator()(index_type a) const
{
  return baseACF(a) + 1;
}

FracACF::FracACF(const ACF& b, index_type l)
  : baseACF(b), df(1, l)
{
  // done
}

FracACF::FracACF(const ACF& b, index_type l, NTYPE f)
  : baseACF(b), df(f, l)
{
  // done
}

FracACF::~FracACF()
{
  // done
}

void FracACF::set(index_type a, NTYPE f)
{
  assert(a < df.length());
  df[a] = f;
}

void FracACF::set(const index_set& d, NTYPE f)
{
  for (index_type k = 0; k < d.length(); k++) {
    assert(d[k] < df.length());
    df[d[k]] = f;
  }
}

NTYPE FracACF::operator()(index_type a) const
{
  assert(a < df.length());
  return (baseACF(a) * df[a]);
}

void DiscountACF::count_only(const bool_vec& d)
{
  assert(d.length() == discounted.length());
  for (index_type k = 0; k < discounted.length(); k++)
    discounted[k] = !d[k];
}

NTYPE DiscountACF::operator()(index_type a) const
{
  if (discounted[a]) return 0;
  else return baseACF(a);
}

NTYPE MakespanACF::operator()(index_type a) const
{
  return instance.actions[a].dur;
}

ResourceConsACF::ResourceConsACF(Instance& i, index_type r)
  : instance(i), resource_id(r)
{
  assert(r < instance.n_resources());
}

NTYPE ResourceConsACF::operator()(index_type a) const
{
  return instance.actions[a].cons[resource_id];
}

ResourceReqACF::ResourceReqACF(Instance& i, index_type r)
  : instance(i), resource_id(r)
{
  assert(r < instance.n_resources());
}

NTYPE ResourceReqACF::operator()(index_type a) const
{
  return instance.actions[a].req(resource_id);
}

END_HSPS_NAMESPACE
