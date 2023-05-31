
#include "enumerators.h"

BEGIN_HSPS_NAMESPACE

IterativeEnumerator::~IterativeEnumerator()
{
  // done
}

SubsetEnumerator::SubsetEnumerator(index_type _n)
  : n(_n), in(false, _n)
{
  // done
}

bool SubsetEnumerator::first()
{
  in.assign_value(false, n);
  return true;
}

bool SubsetEnumerator::next()
{
  index_type p = 0;
  while ((p < n) && in[p]) p += 1;
  if (p < n) {
    for (index_type k = 0; k < p; k++) in[k] = false;
    in[p] = true;
    return true;
  }
  else {
    return false;
  }
}

bool ReverseSubsetEnumerator::first()
{
  in.assign_value(true, n);
  return true;
}

bool ReverseSubsetEnumerator::next()
{
  if (n == 0) return false; // stupid special case...
  index_type n_false = in.count(false);
  if (n_false == 0) {
    in[0] = false;
    return true;
  }
  index_type p = in.first(false);
  assert(p != no_such_index);
  // if the last n_false elements are out, next set will have
  // n_false + 1 elements out; we set those to be the n_false + 1
  // first elements
  if ((p + n_false) == n) {
    // unless the number of elements out are already the whole set,
    // in which case there is no next...
    if (n_false == n) return false;
    for (index_type i = 0; i < n_false + 1; i++) in[i] = false;
    for (index_type i = n_false + 1; i < n; i++) in[i] = true;
  }
  // otherwise, find the last false that has a true after it; this
  // will be moved up by one, and the falses that come after it
  // lined up immediately after
  else {
    index_type j = n - 1;
    index_type n_after = 0;
    while (!in[j]) {
      assert(j > p);
      n_after += 1;
      j -= 1;
    }
    while (in[j]) {
      assert(j > p);
      j -= 1;
    }
    assert(j >= p);
    assert(!in[j]);
    assert(in[j + 1]);
    //std::cerr << "j = " << j << ", n_after = " << n_after << std::endl;
    in[j] = true;
    for (index_type i = j + 1; i <= (j + 1 + n_after); i++) in[i] = false;
    for (index_type i = (j + 1 + n_after + 1); i < n; i++) in[i] = true;
  }
  return true;
}

const bool_vec& SubsetEnumerator::current_set() const
{
  return in;
}

void SubsetEnumerator::current_set(const index_set& elements, index_set& set)
{
  set.clear();
  for (index_type k = 0; k < n; k++)
    if (in[k]) set.append(elements[k]);
}

void SubsetEnumerator::current_set(index_set& set)
{
  set.clear();
  for (index_type k = 0; k < n; k++)
    if (in[k]) set.append(k);
}

index_type SubsetEnumerator::current_set_size()
{
  index_type s = 0;
  for (index_type k = 0; k < n; k++) if (in[k]) s += 1;
  return s;
}

void SubsetEnumerator::all_sets(index_set_vec& sets)
{
  sets.clear();
  bool ok = first();
  while (ok) {
    index_set& s = sets.append();
    current_set(s);
    ok = next();
  }
}

mSubsetEnumerator::mSubsetEnumerator(index_type _n, index_type _m)
  : SubsetEnumerator(_n), m(_m)
{
  assert(m > 0);
}

count_type mSubsetEnumerator::m_of_n(index_type n, index_type m)
{
  if (m > n) return 0;
  if (m > (n / 2))
    m = (n - m);
  count_type r = 1;
  // ::std::cerr << "(0) m = " << m << ", n = " << n << ::std::endl;
  for (index_type k = n; k > (n - m); k--) {
    r = (r * k);
    // ::std::cerr << "(1) k = " << k << ", r = " << r << ::std::endl;
  }
  for (index_type k = 2; k <= m; k++) {
    r = (r / k);
    // ::std::cerr << "(2) k = " << k << ", r = " << r << ::std::endl;
  }
  return r;
}

count_type mSubsetEnumerator::n_sets()
{
  if (m > n) return 1;
  return m_of_n(n, m);
}

bool mSubsetEnumerator::first()
{
  for (index_type k = 0; (k < m) && (k < n); k++) in[k] = true;
  for (index_type k = m; k < n; k++) in[k] = false;
  return true;
}

bool mSubsetEnumerator::next()
{
  if ((m >= n) || (n == 0)) return false;
  if (!in[n - 1]) {
    index_type p = n - 1;
    while (!in[p]) p -= 1;
    assert(in[p]);
    in[p] = false;
    in[p + 1] = true;
    return true;
  }
  else {
    index_type p_last_out = n - 1;
    while (in[p_last_out] && (p_last_out > 0)) p_last_out -= 1;
    assert(!in[p_last_out]);
    index_type p_next_in = p_last_out;
    while (!in[p_next_in] && (p_next_in > 0)) p_next_in -= 1;
    if (!in[p_next_in]) return false;
    index_type n_rem_in = n - p_last_out;
//     ::std::cerr << "in = " << in
// 	      << ", p_last_out = " << p_last_out
// 	      << ", p_next_in = " << p_next_in
// 	      << ", n_rem_in = " << n_rem_in
// 	      << ::std::endl;
    in[p_next_in] = false;
    for (index_type k = 0; k < n_rem_in; k++) in[p_next_in + k + 1] = true;
    for (index_type k = p_next_in + n_rem_in + 1; k < n; k++) in[k] = false;
    return true;
  }
}

kAssignmentEnumerator::kAssignmentEnumerator(index_type _n, index_type _k)
  : n(_n), k(_k), a(0, _k)
{
  assert(k > 0);
}

bool kAssignmentEnumerator::first()
{
  for (index_type i = 0; i < n; i++) a[i] = 0;
  return true;
}

bool kAssignmentEnumerator::next()
{
  index_type e = n - 1;
  while (e > 0) {
    a[e] = ((a[e] + 1) % k);
    if (a[e] > 0) return true;
    e -= 1;
  }
  a[0] = ((a[0] + 1) % k);
  if (a[0] > 0) return true;
  return false;
}

void kAssignmentEnumerator::current_assignment(index_set_vec& sets)
{
  sets.set_length(k);
  for (index_type i = 0; i < k; i++) sets[i].clear();
  for (index_type i = 0; i < n; i++) sets[a[i]].insert(i);
}

CorrespondanceEnumerator::CorrespondanceEnumerator
(const index_vec& v0, const index_vec& v1)
  : a(v0), b(v1), c(v0.length(), no_such_index, false), f(true, v0.length())
{
  // done
}

index_type CorrespondanceEnumerator::first_free
(index_type x, const index_vec& vec, const bool_vec& f_vec)
{
  for (index_type i = 0; i < vec.length(); i++)
    if ((vec[i] == x) && f_vec[i]) return i;
  return no_such_index;
}

index_type CorrespondanceEnumerator::next_free
(index_type x, const index_vec& vec, const bool_vec& f_vec, index_type k)
{
  for (index_type i = k + 1; i < vec.length(); i++)
    if ((vec[i] == x) && f_vec[i]) return i;
  return no_such_index;
}

bool CorrespondanceEnumerator::find(index_type p)
{
  if (p == a.length()) return true;
  index_type i = first_free(a[p], b, f);
  while (i != no_such_index) {
    c[p] = i;
    f[i] = false;
    bool ok = find(p + 1);
    if (ok) return true;
    c[p] = no_such_index;
    f[i] = true;
    i = next_free(a[p], b, f, i);
  }
  return false;
}

bool CorrespondanceEnumerator::first()
{
  if (a.length() != b.length()) return false;
  c.assign_value(no_such_index, a.length());
  f.assign_value(true, a.length());
  return find(0);
}

bool CorrespondanceEnumerator::next()
{
  index_type p = a.length();
  while (p > 0) {
    index_type i = c[p - 1];
    f[i] = true;
    c[p - 1] = no_such_index;
    i = next_free(a[p - 1], b, f, i);
    while (i != no_such_index) {
      c[p - 1] = i;
      f[i] = false;
      bool ok = find(p);
      if (ok) return true;
      c[p - 1] = no_such_index;
      f[i] = true;
      i = next_free(a[p - 1], b, f, i);
    }
    p -= 1;
  }
  return false;
}

void write_correspondance(::std::ostream& s, const index_vec& c)
{
  for (index_type k = 0; k < c.length(); k++) {
    if (k > 0) s << ", ";
    s << k << " <-> " << c[k];
  }
}

PermutationEnumerator::PermutationEnumerator(index_type n)
  : e(index_vec(0, n), index_vec(0, n))
{
  // done
}

bool PermutationEnumerator::first()
{
  return e.first();
}

bool PermutationEnumerator::next()
{
  return e.next();
}

void RecursivekPartitionEnumerator::partition(index_type cn, index_type ck)
{
  if (done) return;
  if (cn < ck) return;
  if (cn == ck) {
    for (index_type i = 0; i < ck; i++) ass[i] = i;
    solution();
  }
  else if (ck == 1) {
    for (index_type i = 0; i < cn; i++) ass[i] = 0;
    solution();
  }
  else {
    ass[cn - 1] = ck - 1;
    partition(cn - 1, ck - 1);
    for (index_type i = 0; i < ck; i++) {
      ass[cn - 1] = i;
      partition(cn - 1, ck);
    }
  }
}

void RecursivekPartitionEnumerator::construct()
{
  sets.assign_value(index_set(), k);
  for (index_type i = 0; i < n; i++)
    sets[ass[i]].insert(i);
}

void RecursivekPartitionEnumerator::construct(const index_set& set)
{
  sets.assign_value(index_set(), k);
  for (index_type i = 0; i < n; i++)
    sets[ass[i]].insert(set[i]);
}

void RecursivekPartitionEnumerator::solution()
{
  // does nothing
}

RecursivekPartitionEnumerator::RecursivekPartitionEnumerator
(index_type _n, index_type _k)
  : ass(0, _n), n(_n), k(_k), sets(EMPTYSET, _k), done(false)
{
  // done
}

RecursivekPartitionEnumerator::~RecursivekPartitionEnumerator()
{
  // done
}

void RecursivekPartitionEnumerator::partition()
{
  done = false;
  partition(n, k);
}


RecursivePartitionEnumerator::RecursivePartitionEnumerator(index_type _n)
  : RecursivekPartitionEnumerator(_n, _n)
{
  // done
}

void RecursivePartitionEnumerator::partition()
{
  for (k = 1; k <= n; k++) {
    RecursivekPartitionEnumerator::partition();
    if (done) return;
  }
}

void RecursivePartitionEnumerator::partition_bounded
(index_type min, index_type max)
{
  if (min < 1) min = 1;
  if (max > n) max = n;
  for (k = min; k <= max; k++) {
    RecursivekPartitionEnumerator::partition();
    if (done) return;
  }  
}

void CountPartitions::solution()
{
  c += 1;
}

index_type CountPartitions::count()
{
  c = 0;
  partition();
  return c;
}

void PrintPartitions::solution()
{
  construct();
  ::std::cerr << sets << ::std::endl;
}

END_HSPS_NAMESPACE
