
#include "index_type.h"

BEGIN_HSPS_NAMESPACE

const index_set EMPTYSET;
index_vec_util::decreasing_index_order index_vec_util::decreasing;
index_vec_util::increasing_index_order index_vec_util::increasing;

pair_vec_util::decreasing_first_order pair_vec_util::decreasing_on_first;
pair_vec_util::decreasing_second_order pair_vec_util::decreasing_on_second;
pair_vec_util::increasing_first_order pair_vec_util::increasing_on_first;
pair_vec_util::increasing_second_order pair_vec_util::increasing_on_second;

void index_vec_util::fill(index_vec& vec, index_type max)
{
  vec.set_length(max);
  for (index_type k = 0; k < max; k++) vec[k] = k;
}

index_type index_vec_util::min(const index_vec& vec, index_type def)
{
  index_type a_min = vec.arg_min();
  if (a_min != no_such_index)
    return vec[a_min];
  else
    return def;
}

index_type index_vec_util::min
(const index_vec& vec, const index_set& sel, index_type def)
{
  index_type a_min = vec.arg_min(sel);
  if (a_min != no_such_index)
    return vec[a_min];
  else
    return def;
}

index_type index_vec_util::min
(const index_vec& vec, const bool_vec& sel, index_type def)
{
  index_type a_min = vec.arg_min(sel);
  if (a_min != no_such_index)
    return vec[a_min];
  else
    return def;
}

index_type index_vec_util::max(const index_vec& vec, index_type def)
{
  index_type a_max = vec.arg_max();
  if (a_max != no_such_index)
    return vec[a_max];
  else
    return def;
}

index_type index_vec_util::max
(const index_vec& vec, const index_set& sel, index_type def)
{
  index_type a_max = vec.arg_max(sel);
  if (a_max != no_such_index)
    return vec[a_max];
  else
    return def;
}

index_type index_vec_util::max
(const index_vec& vec, const bool_vec& sel, index_type def)
{
  index_type a_max = vec.arg_max(sel);
  if (a_max != no_such_index)
    return vec[a_max];
  else
    return def;
}

index_type index_vec_util::sum(const index_vec& vec)
{
  index_type s = 0;
  for (index_type k = 0; k < vec.size(); k++)
    s += vec[k];
  return s;
}

int index_vec_util::compare(const index_vec& v0, const index_vec& v1)
{
  if (v0.length() < v1.length()) return -1;
  if (v0.length() > v1.length()) return 1;
  for (index_type k = 0; k < v0.length(); k++) {
    if (v0[k] < v1[k]) return -1;
    if (v0[k] > v1[k]) return 1;
  }
  return 0;
}

index_type index_vec_util::hash(const index_vec& vec)
{
  index_type v = 0;
  for (index_type k = 0; k < vec.length(); k++)
    v = ((v + (vec[k] * k)) % LARGE_PRIME);
  return v;
}

int index_vec_util::compare(const index_vec& v1) const
{
  return compare(*this, v1);
}

void index_vec_util::fill(index_type max)
{
  fill(*this, max);
}

index_type index_vec_util::hash() const
{
  return hash(*this);
}

index_set::index_set(const bool* _arr, index_type n)
{
  for (index_type k = 0; k < n; k++)
    if (_arr[k]) append(k);
}

index_set::index_set(const bool_vec& _vec)
{
  for (index_type k = 0; k < _vec.length(); k++)
    if (_vec[k]) append(k);
}

index_set::index_set(const index_set& s0, const index_set& s)
{
  for (index_type k = 0; k < s.length(); k++) {
    assert(s[k] < s0.length());
    append(s0[s[k]]);
  }
}

index_set::index_set(const index_set& s0, const bool_vec& s)
{
  assert(s.length() >= s0.length());
  for (index_type k = 0; k < s0.length(); k++) if (s[k]) {
    append(s0[k]);
  }
}

index_set::index_set(const index_set& s0, const index_vec& map)
{
  for (index_type k = 0; k < s0.length(); k++) {
    assert(s0[k] < map.length());
    if (map[s0[k]] != no_such_index)
      insert(map[s0[k]]);
  }
}

index_set::index_set(index_type n, const index_set& s)
{
  index_type k = 0;
  index_type p = 0;
  while (k < n) {
    if (p < s.size()) {
      if (k == s[p]) {
	p += 1;
      }
      else {
	assert(k < s[p]);
	append(k);
      }
    }
    else {
      append(k);
    }
    k += 1;
  }
}

index_set::index_set(index_type e)
{
  set_length(1);
  (*this)[0] = e;
}

index_set::index_set(index_type e1, index_type e2)
{
  if (e1 < e2) {
    set_length(2);
    (*this)[0] = e1;
    (*this)[1] = e2;
  }
  else if (e2 < e1) {
    set_length(2);
    (*this)[0] = e2;
    (*this)[1] = e1;
  }
  else {
    set_length(1);
    (*this)[0] = e1;
  }
}

void index_set::insert(const index_type& v)
{
  svector<index_type>::insert(v);
}

void index_set::insert(const index_set& set)
{
  svector<index_type>::insert(set);
}

void index_set::insert(const index_vec& vec)
{
  svector<index_type>::insert(vec);
}

void index_set::insert(const bool_vec& set)
{
  for (index_type k = 0; k < set.length(); k++)
    if (set[k]) insert(k);
}

void index_set::intersect(const index_set& set)
{
  svector<index_type>::intersect(set);
}

void index_set::intersect(const bool_vec& set)
{
  index_type k = 0;
  while (k < size()) {
    if ((*this)[k] < set.length()) {
      if (!set[(*this)[k]]) {
	remove(k);
      }
      else {
	k += 1;
      }
    }
    else {
      remove(k, length());
    }
  }
}

void index_set::subtract(const index_vec& vec)
{
  svector<index_type>::subtract(vec);
}

void index_set::subtract(const bool_vec& set)
{
  index_type k = 0;
  while (k < length()) {
    if ((*this)[k] < set.length()) {
      if (set[(*this)[k]]) {
	remove(k);
      }
      else {
	k += 1;
      }
    }
    else {
      return;
    }
  }
}

void index_set::subtract(const index_type& v)
{
  svector<index_type>::subtract(v);
}

bool index_set::contains(const bool_vec& set) const
{
  index_type k0 = 0;
  index_type k1 = 0;
  while ((k0 < length()) && (k1 < set.length())) {
    if (set[k1])
      if (k1 != (*this)[k0])
	return false;
    if (k1 == (*this)[k0])
      k0 += 1;
    k1 += 1;
  }
  while (k1 < set.length()) {
    if (set[k1])
      return false;
    k1 += 1;
  }
  return true;
}

bool index_set::contained_in(const lvector<bool>& set) const
{
  for (index_type k = 0; k < size(); k++) {
    if ((*this)[k] > set.length()) return false;
    else if (!set[(*this)[k]]) return false;
  }
  return true;
}

void index_set::fill(index_type max)
{
  index_vec_util::fill(*this, max);
}

void index_set::assign_remap
(const index_set& set, const index_vec& map)
{
  clear();
  for (index_type k = 0; k < set.length(); k++) {
    assert(set[k] < map.length());
    if (map[set[k]] != no_such_index)
      insert(map[set[k]]);
  }
}

void index_set::remap(const index_vec& map)
{
  index_set s0(*this);
  assign_remap(s0, map);
}

index_type index_set::first_common_element(const index_set& vec) const
{
  index_pair p = svector<index_type>::first_common(vec);
  if (p.first != no_such_index)
    return (*this)[p.first];
  else
    return no_such_index;
}

index_type index_set::first_common_element(const index_vec& vec) const
{
  index_pair p = svector<index_type>::first_common(vec);
  if (p.first != no_such_index)
    return (*this)[p.first];
  else
    return no_such_index;
}

bool index_set::have_common_element(const index_set& vec) const
{
  return (first_common_element(vec) != no_such_index);
}

bool index_set::have_common_element(const bool_vec& set) const
{
  return (first_common_element(set) != no_such_index);
}

index_type index_set::first_common_element(const bool_vec& vec) const
{
  for (index_type i = 0; i < length(); i++) {
    if ((*this)[i] < vec.length()) {
      if (vec[(*this)[i]])
	return (*this)[i];
    }
    else {
      return no_such_index;
    }
  }
  return no_such_index;
}

index_type index_set::first_common_element(const bool* vec, index_type n) const
{
  for (index_type i = 0; i < length(); i++) {
    if ((*this)[i] < n) {
      if (vec[(*this)[i]])
	return (*this)[i];
    }
    else {
      return no_such_index;
    }
  }
  return no_such_index;
}

index_type index_set::count_common(const index_set& vec) const
{
  return svector<index_type>::count_common(vec);
}

index_type index_set::count_common(const bool_vec& set) const
{
  index_type c = 0;
  for (index_type k = 0; k < size(); k++)
    if ((*this)[k] < set.length())
      if (set[(*this)[k]]) c += 1;
  return c;
}

bool* index_set::copy_to(bool* s, index_type n) const
{
  for (index_type k = 0; k < n; k++)
    s[k] = false;
  for (index_type k = 0; k < size(); k++)
    if ((*this)[k] < n) s[(*this)[k]] = true;
  return s;
}

bool_vec::bool_vec(const index_set& set, index_type l)
  : lvector<bool>(false, l)
{
  for (index_type k = 0; k < set.length(); k++) {
    if (set[k] < l) (*this)[set[k]] = true;
  }
}

void bool_vec::complement()
{
  for (index_type i = 0; i < length(); i++) (*this)[i] = !((*this)[i]);
}

void bool_vec::insert(const lvector<bool>& vec)
{
  for (index_type i = 0; i < length(); i++) if (vec[i]) (*this)[i] = true;
}

void bool_vec::insert(const index_set& set)
{
  for (index_type i = 0; i < set.length(); i++) {
    assert(set[i] < length());
    (*this)[set[i]] = true;
  }
}

void bool_vec::intersect(const lvector<bool>& vec)
{
  for (index_type i = 0; i < length(); i++) if (!vec[i]) (*this)[i] = false;
}

void bool_vec::intersect(const index_set& set)
{
  for (index_type i = 0; i < length(); i++)
    if (!set.contains(i)) (*this)[i] = false;
}

void bool_vec::subtract(const lvector<bool>& vec)
{
  for (index_type i = 0; i < length(); i++) if (vec[i]) (*this)[i] = false;
}

void bool_vec::subtract(const index_set& set)
{
  for (index_type i = 0; i < set.length(); i++) {
    assert(set[i] < length());
    (*this)[set[i]] = false;
  }
}

bool bool_vec::subset(const lvector<bool>& vec) const
{
  assert(length() == vec.length());
  for (index_type i = 0; i < length(); i++)
    if ((*this)[i] && !vec[i]) return false;
  return true;
}

bool bool_vec::strict_subset(const lvector<bool>& vec) const
{
  assert(length() == vec.length());
  bool strict = false;
  for (index_type i = 0; i < length(); i++) {
    if ((*this)[i] && !vec[i]) return false;
    if (vec[i] && !(*this)[i]) strict = true;
  }
  return strict;
}

bool bool_vec::superset(const lvector<bool>& vec) const
{
  assert(length() == vec.length());
  for (index_type i = 0; i < length(); i++)
    if (!(*this)[i] && vec[i]) return false;
  return true;
}

bool bool_vec::strict_superset(const lvector<bool>& vec) const
{
  assert(length() == vec.length());
  bool strict = false;
  for (index_type i = 0; i < length(); i++) {
    if (!(*this)[i] && vec[i]) return false;
    if (!vec[i] && (*this)[i]) strict = true;
  }
  return strict;
}

bool bool_vec::contains(const bool& v) const
{
  for (index_type k = 0; k < length(); k++)
    if ((*this)[k] == v) return true;
  return false;
}

bool bool_vec::contains(const lvector<bool>& set) const
{
  index_type n = (set.length() > length() ? set.length() : length());
  for (index_type k = 0; k < n; k++)
    if (set[k] && !(*this)[k]) return false;
  return true;
}

bool bool_vec::contains(const index_set& set) const
{
  for (index_type k = 0; k < set.length(); k++)
    if (!(*this)[set[k]]) return false;
  return true;
}

bool bool_vec::contains_any(const index_set& set) const
{
  for (index_type k = 0; k < set.length(); k++)
    if ((*this)[set[k]]) return true;
  return false;
}

index_type bool_vec::first_common_element(const index_set& vec) const
{
  return vec.first_common_element(*this);
}

index_type bool_vec::first_common_element(const lvector<bool>& vec) const
{
  for (index_type k = 0; (k < length()) && (k < vec.length()); k++)
    if ((*this)[k] && vec[k])
      return k;
  return no_such_index;
}

index_type bool_vec::count_common(const lvector<bool>& vec) const
{
  index_type c = 0;
  index_type l = (size() < vec.size() ? size() : vec.size());
  for (index_type k = 0; k < l; k++)
    if ((*this)[k] && vec[k]) c += 1;
  return c;
}

index_type bool_vec::count_common(const index_set& set) const
{
  return set.count_common(*this);
}

index_set& bool_vec::copy_to(index_set& set) const
{
  set.set_length(0);
  for (index_type i = 0; i < length(); i++) if ((*this)[i]) set.insert(i);
  return set;
}

index_set& bool_vec::insert_into(index_set& set) const
{
  for (index_type i = 0; i < length(); i++) if ((*this)[i]) set.insert(i);
  return set;
}

index_set& bool_vec::subtract_from(index_set& set) const
{
  for (index_type i = 0; i < length(); i++) if ((*this)[i]) set.subtract(i);
  return set;
}

bool* bool_vec::copy_to(bool* s, index_type n) const
{
  for (index_type i = 0; i < n; i++) s[i] = (*this)[i];
}

int bool_vec::compare(const lvector<bool>& vec) const
{
  if (length() < vec.length()) return -1;
  if (length() > vec.length()) return 1;
  for (index_type k = 0; k < length(); k++) {
    if (vec[k] && !(*this)[k]) return -1;
    if (!vec[k] && (*this)[k]) return 1;
  }
  return 0;
}

index_type bool_vec::hash() const
{
  index_type v = 0;
  index_type vi = 0;
  for (index_type k = 0; k < length(); k++) {
    if ((*this)[k]) vi += (1 << (k % (INDEX_TYPE_BITS - 1)));
    if ((k % (INDEX_TYPE_BITS - 1)) == 0) {
      v += vi;
      vi = 0;
    }
  }
}

void mapping::extend_domain(index_type n, index_vec& map)
{
  if (map.size() >= n) return;
  index_type m = map.size();
  map.set_length(n);
  for (index_type k = m; k < n; k++) map[k] = k;
}

void mapping::delete_index_map
(index_type n, index_type i, index_vec& map)
{
  assert(i < n);
  index_vec_util::fill(map, n);
  map[i] = no_such_index;
  for (index_type k = i + 1; k < n; k++)
    map[k] = (k - 1);
}

void mapping::subset_map
(index_type n, const index_set s, index_vec& map)
{
  map.assign_value(no_such_index, n);
  for (index_type k = 0; k < s.size(); k++) {
    assert(s[k] < n);
    map[s[k]] = k;
  }
}

void mapping::all_to_one_map
(index_type n, index_type i, index_vec& map)
{
  map.assign_value(i, n);
}

bool mapping::invert_map(const index_vec& map, index_vec& inv, index_type m)
{
  index_type n = m;
  for (index_type k = 0; k < map.length(); k++)
    if (map[k] != no_such_index)
      if (map[k] > n) n = map[k];
  inv.assign_value(no_such_index, n + 1);
  bool ok = true;
  for (index_type k = 0; k < map.length(); k++)
    if (map[k] != no_such_index) {
      if (inv[map[k]] == no_such_index)
	inv[map[k]] = k;
      else
	ok = false;
    }
  return ok;
}

void mapping::compose
(const index_vec& m0, const index_vec& m1, index_vec& cm)
{
  cm.set_length(m0.length());
  for (index_type k = 0; k < m0.length(); k++)
    if (m0[k] != no_such_index) {
      assert(m0[k] < m1.length());
      cm[k] = m1[m0[k]];
    }
    else {
      cm[k] = no_such_index;
    }
}

bool mapping::remap
(const index_vec& m0, const index_vec& m1, index_vec& rm)
{
  bool ok = invert_map(m1, rm);
  if (!ok) return false;
  compose(rm, m0, rm);
  compose(rm, m1, rm);
  return true;
}

void mapping::map_image
(const index_vec& map, const index_vec& vec, index_vec& img)
{
  img.clear();
  for (index_type k = 0; k < vec.length(); k++) {
    assert(vec[k] < map.length());
    if (map[vec[k]] != no_such_index)
      img.append(map[vec[k]]);
  }
}

void mapping::map_image
(const index_vec& map, const index_set& s, index_set& img)
{
  for (index_type k = 0; k < s.size(); k++) {
    assert(s[k] < map.length());
    if (map[s[k]] != no_such_index)
      img.insert(map[s[k]]);
  }
}

void mapping::inverse_map_image
(const index_vec& map, index_type v, index_set& img)
{
  img.clear();
  for (index_type k = 0; k < map.length(); k++)
    if (map[k] == v)
      img.insert(k);
}

void mapping::inverse_map_image
(const index_vec& map, const index_set& vs, index_set& img)
{
  img.clear();
  for (index_type k = 0; k < map.length(); k++)
    if (vs.contains(map[k]))
      img.insert(k);
}

index_type mapping::range(const index_vec& map, index_type d)
{
  if (d == 0) return 0;
  assert(map.length() >= d);
  index_type m = 0;
  for (index_type k = 0; k <= (d - 1); k++)
    if (map[k] != no_such_index)
      if (map[k] > m) m = map[k];
  return m + 1;
}

void mapping::range_set(const index_vec& map, index_set& s)
{
  s.clear();
  for (index_type k = 0; k < map.size(); k++)
    if (map[k] != no_such_index)
      s.insert(map[k]);
}

void mapping::compact(const index_vec& m0, index_vec& cm)
{
  index_set s;
  range_set(m0, s);
  cm.set_length(m0.size());
  for (index_type k = 0; k < m0.size(); k++)
    if (m0[k] != no_such_index) {
      index_type i = s.first(m0[k]);
      assert(i < s.size());
      cm[k] = i;
    }
    else {
      cm[k] = no_such_index;
    }
}

void sparse_mapping::dense_to_sparse
(const index_vec& dm, pair_vec sm)
{
  sm.clear();
  for (index_type k = 0; k < dm.length(); k++)
    if (dm[k] != no_such_index)
      sm.append(index_pair(k, dm[k]));
}

void sparse_mapping::sparse_to_dense
(const pair_vec& sm, index_vec dm)
{
  dm.clear();
  index_type m = no_such_index;
  for (index_type k = 0; k < sm.length(); k++)
    if ((m == no_such_index) || (sm[k].first > m))
      m = sm[k].first;
  dm.assign_value(no_such_index, m);
  for (index_type k = 0; k < sm.length(); k++)
    if (dm[sm[k].first] == no_such_index)
      dm[sm[k].first] = sm[k].second;
}

index_type sparse_mapping::map_image
(const pair_vec& map, index_type x)
{
  for (index_type k = 0; k < map.length(); k++)
    if (map[k].first == x)
      return map[k].second;
    else if (map[k].first > x)
      return no_such_index;
  return no_such_index;
}

void sparse_mapping::map_image
(const pair_vec& map, const index_vec& vec, index_vec& img)
{
  img.clear();
  for (index_type k = 0; k < vec.length(); k++) {
    index_type i = map_image(map, vec[k]);
    if (i != no_such_index)
      img.append(i);
  }
}

void sparse_mapping::inverse_map_image
(const pair_vec& map, index_type x, index_set& img)
{
  img.clear();
  for (index_type k = 0; k < map.length(); k++)
    if (map[k].second == x)
      img.insert(map[k].first);
}

void sparse_mapping::inverse_map_image
(const pair_vec& map, const index_set& x, index_set& img)
{
  img.clear();
  for (index_type k = 0; k < map.length(); k++)
    if (x.contains(map[k].second))
      img.insert(map[k].first);
}

index_type index_set_vec::minimum_cardinality() const
{
  if (length() == 0) return no_such_index;
  index_type lmin = (*this)[0].length();
  for (index_type k = 1; k < length(); k++)
    if ((*this)[k].length() < lmin)
      lmin = (*this)[k].length();
  return lmin;
}

index_type index_set_vec::maximum_cardinality() const
{
  if (length() == 0) return no_such_index;
  index_type lmax = (*this)[0].length();
  for (index_type k = 1; k < length(); k++)
    if ((*this)[k].length() > lmax)
      lmax = (*this)[k].length();
  return lmax;
}

index_type index_set_vec::selected_minimum_cardinality
(const index_set& sel) const
{
  if (sel.length() == 0) return no_such_index;
  assert(sel[0] < length());
  index_type lmin = (*this)[sel[0]].length();
  for (index_type k = 1; k < sel.length(); k++) {
    assert(sel[k] < length());
    if ((*this)[sel[k]].length() < lmin)
      lmin = (*this)[sel[k]].length();
  }
  return lmin;
}

index_type index_set_vec::selected_maximum_cardinality
(const index_set& sel) const
{
  if (sel.length() == 0) return no_such_index;
  assert(sel[0] < length());
  index_type lmax = (*this)[sel[0]].length();
  for (index_type k = 1; k < sel.length(); k++) {
    assert(sel[k] < length());
    if ((*this)[sel[k]].length() > lmax)
      lmax = (*this)[sel[k]].length();
  }
  return lmax;
}

index_type index_set_vec::first_minimum_cardinality_set() const
{
  if (length() == 0) return no_such_index;
  index_type lmin = (*this)[0].length();
  index_type imin = 0;
  for (index_type k = 1; k < length(); k++)
    if ((*this)[k].length() < lmin) {
      lmin = (*this)[k].length();
      imin = k;
    }
  return imin;
}

index_type index_set_vec::first_maximum_cardinality_set() const
{
  if (length() == 0) return no_such_index;
  index_type lmax = (*this)[0].length();
  index_type imax = 0;
  for (index_type k = 1; k < length(); k++)
    if ((*this)[k].length() > lmax) {
      lmax = (*this)[k].length();
      imax = k;
    }
  return imax;
}

index_type index_set_vec::first_superset(const index_set& set) const
{
  index_type k = 0;
  while (k < length()) {
    if ((*this)[k].contains(set))
      return k;
    k += 1;
  }
  return no_such_index;
}

index_type index_set_vec::first_superset(const bool_vec& set) const
{
  index_type k = 0;
  while (k < length()) {
    if ((*this)[k].contains(set))
      return k;
    k += 1;
  }
  return no_such_index;
}

index_type index_set_vec::first_strict_superset(const index_set& set) const
{
  index_type k = 0;
  while (k < length()) {
    if ((*this)[k].contains(set) && !set.contains((*this)[k]))
      return k;
    k += 1;
  }
  return no_such_index;
}

index_type index_set_vec::first_contains(index_type e) const
{
  index_type k = 0;
  while (k < length()) {
    if ((*this)[k].contains(e))
      return k;
    k += 1;
  }
  return no_such_index;
}

index_type index_set_vec::next_contains(index_type e, index_type p) const
{
  index_type k = p + 1;
  while (k < length()) {
    if ((*this)[k].contains(e))
      return k;
    k += 1;
  }
  return no_such_index;
}

index_type index_set_vec::first_subset(const index_set& set) const
{
  index_type k = 0;
  while (k < length()) {
    if (set.contains((*this)[k]))
      return k;
    k += 1;
  }
  return no_such_index;
}

index_type index_set_vec::first_subset(const bool_vec& set) const
{
  index_type k = 0;
  while (k < length()) {
    if (set.contains((*this)[k]))
      return k;
    k += 1;
  }
  return no_such_index;
}

index_type index_set_vec::first_strict_subset(const index_set& set) const
{
  index_type k = 0;
  while (k < length()) {
    if (set.contains((*this)[k]) && !(*this)[k].contains(set))
      return k;
    k += 1;
  }
  return no_such_index;
}

index_set& index_set_vec::union_set(index_set& set) const
{
  set.clear();
  for (index_type k = 0; k < length(); k++)
    set.insert((*this)[k]);
  return set;
}

index_set& index_set_vec::selected_union_set
(const index_set& sel, index_set& set) const
{
  set.clear();
  for (index_type k = 0; k < sel.length(); k++) {
    assert(sel[k] < length());
    set.insert((*this)[sel[k]]);
  }
  return set;
}

index_set& index_set_vec::intersection_set(index_set& set) const
{
  if (length() == 0) {
    set.clear();
    return set;
  }
  set.assign_copy((*this)[0]);
  for (index_type k = 1; k < length(); k++)
    set.intersect((*this)[k]);
  return set;
}

void index_set_vec::append_if_not_subset(const index_set& set)
{
  for (index_type k = 0; k < length(); k++)
    if ((*this)[k].contains(set)) return;
  append(set);
}

void index_set_vec::append_if_not_superset(const index_set& set)
{
  for (index_type k = 0; k < length(); k++)
    if (set.contains((*this)[k])) return;
  append(set);
}

void index_set_vec::append_if_new(const index_set& set)
{
  for (index_type k = 0; k < length(); k++)
    if (set == (*this)[k]) return;
  append(set);
}

void index_set_vec::insert_maximal(const index_set& set)
{
  for (index_type k = 0; k < length(); k++)
    if ((*this)[k].contains(set)) return;
  index_type r = 0;
  index_type w = 0;
  while (r < length()) {
    if (!set.contains((*this)[r])) {
      if (w < r) (*this)[w] = (*this)[r];
      w += 1;
    }
    r += 1;
  }
  set_length(w);
  append(set);
}

void index_set_vec::insert_minimal(const index_set& set)
{
  for (index_type k = 0; k < length(); k++)
    if (set.contains((*this)[k])) return;
  index_type r = 0;
  index_type w = 0;
  while (r < length()) {
    if (!(*this)[r].contains(set)) {
      if (w < r) (*this)[w] = (*this)[r];
      w += 1;
    }
    r += 1;
  }
  set_length(w);
  append(set);
}

void index_set_vec::maximal(bool_vec& is_max)
{
  is_max.assign_value(true, length());
  for (index_type i = 0; i < length(); i++)
    if (is_max[i])
      for (index_type j = 0; j < length(); j++)
	if ((j != i) && is_max[j])
	  if ((*this)[i].contains((*this)[j]))
	    is_max[j] = false;
}

void index_set_vec::minimal(bool_vec& is_min)
{
  is_min.assign_value(true, length());
  for (index_type i = 0; i < length(); i++)
    if (is_min[i])
      for (index_type j = 0; j < length(); j++)
	if ((j != i) && is_min[j])
	  if ((*this)[j].contains((*this)[i]))
	    is_min[j] = false;
}

void index_set_vec::reduce_to_maximal()
{
  bool_vec d(false, length());
  for (index_type i = 0; i < length(); i++) if (!d[i])
    for (index_type j = 0; j < length(); j++) if ((j != i) && !d[j])
      if ((*this)[i].contains((*this)[j]))
	d[j] = true;
  assert(d.count(false) > 0);
  remove(d);
}

void index_set_vec::reduce_to_minimal()
{
  bool_vec d(false, length());
  for (index_type i = 0; i < length(); i++)
    if (!d[i])
      for (index_type j = 0; j < length(); j++)
	if ((j != i) && !d[j])
	  if ((*this)[j].size() >= (*this)[i].size())
	    if ((*this)[j].contains((*this)[i]))
	      d[j] = true;
  assert(d.count(false) > 0);
  remove(d);
}

void index_set_vec::remove_sets_size_le(index_type l)
{
  index_type r = 0;
  index_type w = 0;
  while (r < length()) {
    if ((*this)[r].length() > l) {
      if (w < r) (*this)[w] = (*this)[r];
      w += 1;
    }
    r += 1;
  }
  set_length(w);
}

void index_set_vec::remove_empty_sets()
{
  index_type r = 0;
  index_type w = 0;
  while (r < length()) {
    if (!(*this)[r].empty()) {
      if (w < r) (*this)[w] = (*this)[r];
      w += 1;
    }
    r += 1;
  }
  set_length(w);
}

void index_set_vec::insert_in_all(index_type i)
{
  for (index_type k = 0; k < length(); k++)
    (*this)[k].insert(i);
}

void index_set_vec::insert_in_all(const index_set& set)
{
  for (index_type k = 0; k < length(); k++)
    (*this)[k].insert(set);
}

void index_set_vec::subtract_from_all(index_type i)
{
  for (index_type k = 0; k < length(); k++)
    (*this)[k].subtract(i);
}

void index_set_vec::subtract_from_all(const index_set& set)
{
  for (index_type k = 0; k < length(); k++)
    (*this)[k].subtract(set);
}

void index_set_vec::combinations_by_union(const index_set_vec& sv)
{
  index_type l = length();
  index_type k = 0;
  while (k < sv.length()) {
    // i_from .. i_to are the last l elements
    index_type i_from = length() - l;
    index_type i_to = length();
    if (k < (sv.length() - 1)) {
      // append copies of last l elements
      for (index_type i = i_from; i < i_to; i++)
	append((*this)[i]);
    }
    // add sv[k] to elements in [i_from .. i_to]
    for (index_type i = i_from; i < i_to; i++)
      (*this)[i].insert(sv[k]);
    k += 1;
  }
  // special case not handled by above algorithm
  if (sv.length() == 0) clear();
}

void index_set_vec::combinations_by_union
(const index_set_vec& sv1, const index_set_vec& sv2)
{
  for (index_type i = 0; i < sv1.length(); i++)
    for (index_type j = 0; j < sv2.length(); j++) {
      append(sv1[i]);
      (*this)[length() - 1].insert(sv2.length());
    }
}

void bool_matrix::complement()
{
  for (index_type i = 0; i < rows(); i++)
    for (index_type j = 0; j < columns(); j++)
      (*this)[i][j] = !((*this)[i][j]);
}

void bool_matrix::insert(const bool_matrix& m)
{
  assert(m.rows() == rows());
  assert(m.columns() == columns());
  for (index_type i = 0; i < rows(); i++)
    for (index_type j = 0; j < columns(); j++)
      if (m[i][j])
	(*this)[i][j] = true;
}

void bool_matrix::intersect(const bool_matrix& m)
{
  assert(m.rows() == rows());
  assert(m.columns() == columns());
  for (index_type i = 0; i < rows(); i++)
    for (index_type j = 0; j < columns(); j++)
      if (!m[i][j])
	(*this)[i][j] = false;
}

void bool_matrix::subtract(const bool_matrix& m)
{
  assert(m.rows() == rows());
  assert(m.columns() == columns());
  for (index_type i = 0; i < rows(); i++)
    for (index_type j = 0; j < columns(); j++)
      if (m[i][j])
	(*this)[i][j] = false;
}

void bool_matrix::multiply(const bool_matrix& m0, const bool_matrix& m1)
{
  assert(m0.columns() == m1.rows());
  set_size(m0.rows(), m1.columns());
  for (index_type i = 0; i < rows(); i++)
    for (index_type j = 0; j < columns(); j++) {
      bool p = false;
      for (index_type k = 0; (k < m0.columns()) && !p; k++)
	if (m0[i][k] && m1[k][j]) p = true;
      (*this)[i][j] = p;
    }
}

void bool_matrix::transitive_closure()
{
  assert(rows() == columns());
  for (index_type k = 0; k < rows(); k++)
    for (index_type i = 0; i < rows(); i++)
      if ((*this)[i][k])
	for (index_type j = 0; j < rows(); j++)
	  if ((*this)[k][j])
	    (*this)[i][j] = true;
}

void bool_matrix::and_of_rows(lvector<bool>& v) const
{
  v.assign_value(true, columns());
  for (index_type j = 0; j < columns(); j++)
    for (index_type i = 0; (i < rows()) && v[j]; i++)
      if (!(*this)[i][j]) v[j] = false;
}

void bool_matrix::and_of_rows(const bool_vec& s, lvector<bool>& v) const
{
  v.assign_value(true, columns());
  for (index_type j = 0; j < columns(); j++)
    for (index_type i = 0; (i < rows()) && v[j]; i++)
      if (s[i])
	if (!(*this)[i][j]) v[j] = false;
}

void bool_matrix::or_of_rows(lvector<bool>& v) const
{
  v.assign_value(false, columns());
  for (index_type j = 0; j < columns(); j++)
    for (index_type i = 0; (i < rows()) && !v[j]; i++)
      if ((*this)[i][j]) v[j] = true;
}

void bool_matrix::or_of_rows(const bool_vec& s, lvector<bool>& v) const
{
  v.assign_value(false, columns());
  for (index_type j = 0; j < columns(); j++)
    for (index_type i = 0; (i < rows()) && !v[j]; i++)
      if (s[i])
	if ((*this)[i][j]) v[j] = true;
}

void bool_matrix::and_of_columns(lvector<bool>& v) const
{
  v.assign_value(true, rows());
  for (index_type i = 0; i < rows(); i++)
    for (index_type j = 0; (j < columns()) && v[i]; j++)
      if (!(*this)[i][j]) v[i] = false;
}

void bool_matrix::and_of_columns(const bool_vec& s, lvector<bool>& v) const
{
  v.assign_value(true, rows());
  for (index_type i = 0; i < rows(); i++)
    for (index_type j = 0; (j < columns()) && v[i]; j++)
      if (s[j])
	if (!(*this)[i][j]) v[i] = false;
}

void bool_matrix::or_of_columns(lvector<bool>& v) const
{
  v.assign_value(false, rows());
  for (index_type i = 0; i < rows(); i++)
    for (index_type j = 0; (j < columns()) && !v[i]; j++)
      if ((*this)[i][j]) v[i] = true;
}

void bool_matrix::or_of_columns(const bool_vec& s, lvector<bool>& v) const
{
  v.assign_value(false, rows());
  for (index_type i = 0; i < rows(); i++)
    for (index_type j = 0; (j < columns()) && !v[i]; j++)
      if (s[j])
	if ((*this)[i][j]) v[i] = true;
}

bool equivalence::operator()(index_type a, index_type b) const
{
  return (canonical(a) == canonical(b));
}

index_type equivalence::canonical(index_type a) const
{
  //assert(a < length());
  index_type x = a;
  index_type y = (*this)[x];
  while (x != y) {
    //assert(y < length());
    x = y;
    y = (*this)[x];
  }
  // path compression
  if ((*this)[a] != y) (*((index_vec*)this))[a] = y;
  return y;
}

void equivalence::extend(index_type a)
{
  for (index_type k = length(); k < a; k++) append(k);
}

void equivalence::merge(index_type a, index_type b)
{
  index_type c_a = canonical(a);
  index_type c_b = canonical(b);
  (*this)[c_a] = c_b;
}

void equivalence::merge(const index_set& set)
{
  for (index_type i = 1; i < set.length(); i++)
    merge(set[0], set[i]);
}

void equivalence::merge(const index_set& sa, const index_set& sb)
{
  for (index_type i = 0; i < sa.length(); i++)
    for (index_type j = 0; j < sb.length(); j++)
      merge(sa[i], sb[j]);
}

void equivalence::merge(const equivalence& eq)
{
  for (index_type k = 0; k < eq.length(); k++) {
    index_type i = eq.canonical(k);
    if (i != k) merge(i, k);
  }
}

void equivalence::reset()
{
  for (index_type k = 0; k < length(); k++) ((*this)[k]) = k;
}

void equivalence::reset(index_type n)
{
  set_length(n);
  reset();
}

void equivalence::canonical_set(index_set& set) const
{
  index_set tmp(set);
  set.clear();
  for (index_type k = 0; k < tmp.length(); k++)
    set.insert(canonical(tmp[k]));
}

index_type equivalence::n_squeezed() const
{
  index_type n = 0;
  for (index_type k = 0; k < length(); k++) if (canonical(k) != k) n += 1;
  return n;
}

index_type equivalence::n_classes() const
{
  index_type n = 0;
  for (index_type k = 0; k < length(); k++) if (canonical(k) == k) n += 1;
  return n;
}

void equivalence::canonical_elements(index_set& set) const
{
  for (index_type k = 0; k < length(); k++)
    if (canonical(k) == k) set.insert(k);
}

void equivalence::class_elements(index_type c, index_set& set) const
{
  for (index_type k = 0; k < length(); k++)
    if (canonical(k) == c) set.insert(k);
}

void equivalence::classes(index_set_vec& sets) const
{
  index_set c;
  canonical_elements(c);
  sets.assign_value(EMPTYSET, c.length());
  for (index_type k = 0; k < c.length(); k++)
    class_elements(c[k], sets[k]);
}

index_type equivalence::class_number(index_vec& cn) const
{
  index_type n = 0;
  cn.assign_value(no_such_index, length());
  for (index_type k = 0; k < length(); k++) {
    index_type i = canonical(k);
    assert(i < length());
    if (cn[i] == no_such_index) {
      cn[i] = n++;
    }
    cn[k] = cn[i];
  }
  return n;
}

index_type equivalence::class_number_and_size
(index_vec& cn, index_vec& cs) const
{
  index_type n = 0;
  cn.assign_value(no_such_index, length());
  for (index_type k = 0; k < length(); k++) {
    index_type i = canonical(k);
    assert(i < length());
    if (cn[i] == no_such_index) {
      cn[i] = n++;
    }
    cn[k] = cn[i];
    //cs[cn[k]] += 1;
  }
  cs.assign_value(0, n);
  for (index_type k = 0; k < length(); k++) {
    assert(cn[k] < n);
    cs[cn[k]] += 1;
  }
  return n;
}

index_type equivalence::n_class_elements(index_type c) const
{
  index_type n = 0;
  for (index_type k = 0; k < length(); k++)
    if (canonical(k) == c) n += 1;
  return n;
}

void equivalence::make_map(index_vec& map) const
{
  index_set ce;
  canonical_elements(ce);
  map.assign_value(no_such_index, length());
  for (index_type k = 0; k < length(); k++) {
    index_type c = canonical(k);
    map[k] = ce.first(c);
  }
}

void set_hash_function::init(index_type n)
{
  resize(n);
  if (n == 0) return;
  (*this)[0] = 1;
  for (index_type k = 1; k < size(); k++) {
    (*this)[k] = ((2 * (*this)[k - 1]) % LARGE_PRIME);
  }
}

index_type set_hash_function::operator()(index_type& i, index_type v) const
{
  assert(i < size());
  return ((v + (*this)[i]) % LARGE_PRIME);
}

index_type set_hash_function::operator()
(const index_set& set, index_type v) const
{
  return ((v + (*this)(set)) % LARGE_PRIME);
}

index_type set_hash_function::operator()
(const bool_vec& set, index_type v) const
{
  return ((v + (*this)(set)) % LARGE_PRIME);
}

index_type set_hash_function::operator()(index_type& i) const
{
  assert(i < size());
  return (*this)[i];
}

index_type set_hash_function::operator()(const index_set& set) const
{
  index_type h = 0;
  for (index_type k = 0; k < set.length(); k++) {
    assert(set[k] < size());
    h = ((h + (*this)[set[k]]) % LARGE_PRIME);
  }
  return h;
}

index_type set_hash_function::operator()(const bool_vec& set) const
{
  assert(set.size() <= size());
  index_type h = 0;
  for (index_type k = 0; k < set.length(); k++) {
    if (set[k]) h = ((h + (*this)[k]) % LARGE_PRIME);
  }
  return h;
}

index_type set_hash_function::operator()(const bool* set, index_type n) const
{
  assert(n <= size());
  index_type h = 0;
  for (index_type k = 0; k < n; k++) {
    if (set[k]) h = ((h + (*this)[k]) % LARGE_PRIME);
  }
  return h;
}

END_HSPS_NAMESPACE
