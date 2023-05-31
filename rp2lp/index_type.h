#ifndef VECTOR_H
#define VECTOR_H

#include "config.h"
#include <limits.h>
#include <vector>
#include <utility>
#include <iostream>

BEGIN_HSPS_NAMESPACE

typedef unsigned int index_type;
const index_type index_type_max = (UINT_MAX - 1);
#define INDEX_TYPE_BITS 32
const index_type LARGE_PRIME = 2147483629U;
const index_type no_such_index = UINT_MAX;

typedef unsigned long count_type;
const count_type count_type_max = ULONG_MAX;

template<class T> class swapable_pair : public std::pair<T, T>
{
 public:
  swapable_pair()
    : std::pair<T, T>() { };
  swapable_pair(const T& v1, const T& v2)
    : std::pair<T, T>(v1, v2) { };
  swapable_pair(const T& v)
    : std::pair<T, T>(v, v) { };
  swapable_pair(const swapable_pair& p)
    : std::pair<T, T>(p) { };

  void swap();
};

template<class T> class comparable_pair : public swapable_pair<T>
{
 public:
  comparable_pair()
    : swapable_pair<T>() { };
  comparable_pair(const T& v1, const T& v2)
    : swapable_pair<T>(v1, v2) { };
  comparable_pair(const T& v)
    : swapable_pair<T>(v) { };
  comparable_pair(const comparable_pair& p)
    : swapable_pair<T>(p) { };

  void sort_ascending();
  void sort_descending();
};

typedef comparable_pair<index_type> index_pair;

template<class T> class zero_init_pair : public comparable_pair<T>
{
 public:
  zero_init_pair()
    : comparable_pair<T>(0) { };
  zero_init_pair(const T& v1, const T& v2)
    : comparable_pair<T>(v1, v2) { };
  zero_init_pair(const T& v)
    : comparable_pair<T>(v) { };
  zero_init_pair(const zero_init_pair& p)
    : comparable_pair<T>(p) { };
};


// forward declarations of the two representations of "set of indices"
// - they are used as arguments to some lvector methods
class index_set;
class bool_vec;

template<class T> class lvector : public std::vector<T>
{
 public:
  lvector() : std::vector<T>() { };
  lvector(const T& v, index_type l) : std::vector<T>(l, v) { };
  lvector(const lvector<T>& vec) : std::vector<T>(vec) { };
  // lvector(const T* arr, index_type n) : std::vector<T>(n) {
  //   for (index_type k = 0; k < n; k++) (*this)[k] = arr[k];
  // };

  // sub-vector constructors
  lvector(const lvector<T>& vec, const index_set& sel);
  lvector(const lvector<T>& vec, const bool_vec& sel);

  class element_reference {
    lvector*   _vec;
    index_type _pos;
  public:
    element_reference() : _vec(0), _pos(no_such_index) { };
    element_reference(lvector& v, index_type p) : _vec(&v), _pos(p) { };
    const T* operator->() const {
      assert(_vec != 0);
      return &((*_vec)[_pos]);
    };
    T* operator->() {
      assert(_vec != 0);
      return &((*_vec)[_pos]);
    };
    operator bool() const {
      if (_vec == 0) return false;
      if (_pos >= _vec->size()) return false;
      return true;
    };
  };

  class order {
   public:
    virtual bool operator()(const T& v0, const T& v1) const = 0;
  };

  class increasing_by_value {
    const lvector<T>& val;
   public:
    increasing_by_value(const lvector<T>& v) : val(v) { };
    virtual bool operator()(const index_type& i0, const index_type& i1) const {
      return (val[i0] < val[i1]);
    };
  };

  class decreasing_by_value {
    const lvector<T>& val;
   public:
    decreasing_by_value(const lvector<T>& v) : val(v) { };
    virtual bool operator()(const index_type& i0, const index_type& i1) const {
      return (val[i0] > val[i1]);
    };
  };

  index_type length() const;

#ifdef CHECK_VECTOR_INDEX
  typename std::vector<T>::reference operator[](typename std::vector<T>::size_type k) {
    assert(k < std::vector<T>::size());
    return std::vector<T>::operator[](k);
  };

  typename std::vector<T>::const_reference operator[](typename std::vector<T>::size_type k) const {
    assert(k < std::vector<T>::size());
    return std::vector<T>::operator[](k);
  };
#endif

  bool contains(const T& v) const;
  index_type first(const T& v) const;
  index_type last(const T& v) const;
  index_type next(const T& v, index_type i) const;
  index_type find(const T& v, bool_vec& s) const;
  index_type count(const T& v) const;
  index_type arg_max() const;
  index_type arg_max(const index_set& s) const;
  index_type arg_max(const bool_vec& s) const;
  index_type arg_min() const;
  index_type arg_min(const index_set& s) const;
  index_type arg_min(const bool_vec& s) const;
  index_type arg_first(const order& o) const;
  index_type arg_last(const order& o) const;

  index_pair first_common(const lvector<T>& vec) const;
  index_pair next_common(const lvector<T>& vec, index_pair p) const;
  void difference(const lvector& v1, lvector& d0, lvector& d1);

  bool operator==(const lvector& _vec) const;
  bool operator!=(const lvector& _vec) const;
  bool operator<(const lvector& vec) const;
  bool operator>(const lvector& vec) const;
  bool operator<=(const lvector& vec) const;
  bool operator>=(const lvector& vec) const;

  void assign_copy(const lvector& _vec);
  void assign_copy(const T* _arr, index_type n);
  void assign_value(const T& val);
  void assign_value(const T& val, index_type l);

  // note: map remaps vector indices
  void assign_remap(const lvector& vec, const lvector<index_type>& map);
  void remap(const lvector<index_type>& map);

  void assign_select(const lvector& _vec, const index_set& s);
  void assign_select(const lvector& _vec, const bool_vec& s);

  const lvector& operator=(const lvector& _vec);
  // const T& operator=(const T& _val);

  void set_length(index_type l);
  void set_length(index_type l, const T& v);
  void inc_length_to(index_type l);
  void inc_length_to(index_type l, const T& v);
  index_type inc_length() { return inc_length(1); };
  index_type inc_length(index_type d);
  index_type inc_length(index_type d, const T& v);
  index_type dec_length() { return dec_length(1); };
  index_type dec_length(index_type d);
  void clear();

  void append(const T& v);
  void append(const lvector& v);
  T&   append();
  void insert(const T& v, index_type p);
  index_type insert_ordered(const T& v, const order& o, index_type f = 0);
  index_type insert_ordered(const lvector& vec, const order& o);

  void remove(index_type p);
  void remove(index_type p0, index_type p1);
  void remove(const index_set& s);
  void remove(const index_set& s, lvector<index_type>& map);
  void remove(const bool_vec& s);
  void remove(const bool_vec& s, lvector<index_type>& map);
  void remove_duplicate_elements();

  void swap(index_type i, index_type j);
};

template<class T> class auto_expanding_vector : public lvector<T>
{
  T _default;
 public:
  auto_expanding_vector() : lvector<T>() { };
  auto_expanding_vector(const T& v) : lvector<T>(), _default(v) { };
  auto_expanding_vector(const T& v, index_type l)
    : lvector<T>(v, l), _default(v) { };
  auto_expanding_vector(const lvector<T>& vec)
    : lvector<T>(vec) { };
  auto_expanding_vector(const auto_expanding_vector<T>& vec)
    : lvector<T>(vec), _default(vec._default) { };

  typename std::vector<T>::reference
  operator[](typename std::vector<T>::size_type k)
  {
    this->inc_length_to(k + 1, _default);
    return lvector<T>::operator[](k);
  };

  typename std::vector<T>::const_reference
  operator[](typename std::vector<T>::size_type k) const
  {
    if (k >= std::vector<T>::size())
      return _default;
    else
      return lvector<T>::operator[](k);
  };

  void assign_value(const T& val)
  {
    _default = val;
    lvector<T>::assign_value(val);
  };

  void assign_value(const T& val, index_type l)
  {
    _default = val;
    lvector<T>::assign_value(val, l);
  };
};

template<class T> class read_only_vector_with_default
{
  const lvector<T>* _vec;
  T _default;
 public:
  read_only_vector_with_default(const lvector<T>* v, const T& d)
    : _vec(v), _default(d) { };

  typename std::vector<T>::const_reference
  operator[](typename std::vector<T>::size_type k) const
  {
    if (_vec == 0)
      return _default;
    if (k >= _vec->size())
      return _default;
    else
      return (*_vec)[k];
  };
};

typedef lvector<index_type> index_vec;
typedef lvector<index_pair> pair_vec;

class index_vec_util : public index_vec
{
 public:

  class decreasing_index_order : public index_vec::order {
  public:
    virtual bool operator()
      (const index_type& v0, const index_type& v1) const
      { return (v0 > v1); };
  };

  class increasing_index_order : public index_vec::order {
  public:
    virtual bool operator()
      (const index_type& v0, const index_type& v1) const
      { return (v0 < v1); };
  };

  class increasing_value_order : public index_vec::order {
    const index_vec& value;
  public:
    increasing_value_order(const index_vec& v) : value(v) { };
    virtual bool operator()
      (const index_type& v0, const index_type& v1) const
      {
	assert(v0 < value.length());
	assert(v1 < value.length());
	return (value[v0] < value[v1]);
      };
  };

  static class decreasing_index_order decreasing;
  static class increasing_index_order increasing;

  static void       fill(index_vec& vec, index_type max);

  static index_type min(const index_vec& vec, index_type def = no_such_index);
  static index_type max(const index_vec& vec, index_type def = no_such_index);
  static index_type min(const index_vec& vec, const index_set& sel,
			index_type def = no_such_index);
  static index_type max(const index_vec& vec, const index_set& sel,
			index_type def = no_such_index);
  static index_type min(const index_vec& vec, const bool_vec& sel,
			index_type def = no_such_index);
  static index_type max(const index_vec& vec, const bool_vec& sel,
			index_type def = no_such_index);
  static index_type sum(const index_vec& vec);
  static int        compare(const index_vec& v0, const index_vec& v1);
  static index_type hash(const index_vec& vec);

  void fill(index_type max);
  int  compare(const index_vec& v1) const;
  index_type hash() const;
};

class pair_vec_util : public pair_vec
{
 public:

  class decreasing_first_order : public pair_vec::order {
  public:
    virtual bool operator()
      (const index_pair& v0, const index_pair& v1) const
      { return (v0.first > v1.first); };
  };

  class decreasing_second_order : public pair_vec::order {
  public:
    virtual bool operator()
      (const index_pair& v0, const index_pair& v1) const
      { return (v0.second > v1.second); };
  };

  class increasing_first_order : public pair_vec::order {
  public:
    virtual bool operator()
      (const index_pair& v0, const index_pair& v1) const
      { return (v0.first < v1.first); };
  };

  class increasing_second_order : public pair_vec::order {
  public:
    virtual bool operator()
      (const index_pair& v0, const index_pair& v1) const
      { return (v0.second < v1.second); };
  };

  static class decreasing_first_order  decreasing_on_first;
  static class decreasing_second_order decreasing_on_second;
  static class increasing_first_order  increasing_on_first;
  static class increasing_second_order increasing_on_second;
};

template<class T> class svector : public lvector<T>
{
 public:
  svector() : lvector<T>() { };
  svector(const svector<T>& _svec) : lvector<T>(_svec) { };
  svector(const lvector<T>& _lvec) : lvector<T>() {
    for (index_type k = 0; k < _lvec.size(); k++) insert(_lvec[k]);
  };
  // svector(const T* _arr, index_type n) : lvector<T>() {
  //   for (index_type k = 0; k < n; k++) insert(_arr[k]);
  // };

  bool contains(const T& v) const;
  bool contains2(const T& v) const;
  bool contains(const svector& vec) const;
  bool subset(const svector& vec) const;

  index_pair first_common(const lvector<T>& vec) const;
  index_pair next_common(const lvector<T>& vec, index_pair p) const;
  index_pair first_common(const svector<T>& vec) const;
  index_pair next_common(const svector<T>& vec, index_pair p) const;
  index_type count_common(const svector& vec) const;

  void assign_singleton(const T& _val);
  void assign_values(const lvector<T>& vec);
  void insert(const T& v);
  void insert(const svector<T>& vec);
  void insert(const lvector<T>& vec);
  void intersect(const svector& vec);
  void difference(const svector& vec);
  void subtract(const svector& vec);
  void subtract(const T& v);

  void remove_greater_than(const T& _val);
  void remove_less_than(const T& _val);
};

class index_set : public svector<index_type>
{
 public:
  // standard svector constructors
  index_set()
    : svector<index_type>() { };
  index_set(const index_set& _svec)
    : svector<index_type>(_svec) { };
  index_set(const lvector<index_type>& _lvec)
    : svector<index_type>(_lvec) { };
  // index_set(const index_type* _arr, index_type n)
  //   : svector<index_type>(_arr, n) { };

  // conversion from array of bool and bool_vec
  index_set(const bool* _arr, index_type n);
  index_set(const bool_vec& _vec);

  // subset constructor
  index_set(const index_set& s0, const index_set& s);
  index_set(const index_set& s0, const bool_vec& s);

  // remap constructor
  index_set(const index_set& s0, const index_vec& map);

  // complement constructor: set becomes fill(n) - s
  index_set(index_type n, const index_set& s);

  // singleton set constructor
  index_set(index_type e);

  // two-element set constructor
  index_set(index_type e1, index_type e2);

  index_type first_common_element(const index_set& set) const;
  index_type first_common_element(const index_vec& vec) const;
  index_type first_common_element(const bool_vec& vec) const;
  index_type first_common_element(const bool* vec, index_type n) const;

  index_type count_common(const index_set& set) const;
  index_type count_common(const bool_vec& set) const;

  bool have_common_element(const index_set& set) const;
  bool have_common_element(const bool_vec& set) const;

  void insert(const index_type& v);
  void insert(const index_set& set);
  void insert(const index_vec& vec);
  void insert(const bool_vec& set);
  void intersect(const index_set& vec);
  void intersect(const bool_vec& set);
  void subtract(const index_vec& vec);
  void subtract(const bool_vec& set);
  void subtract(const index_type& v);

  bool contains(const index_type& v) const
    { return svector<index_type>::contains(v); };
  bool contains(const index_set& set) const
    { return svector<index_type>::contains(set); };
  bool contains(const bool_vec& set) const;
  bool contained_in(const lvector<bool>& set) const;

  bool* copy_to(bool* s, index_type n) const;

  void  fill(index_type to);
  // note: map here remaps elements
  void  assign_remap(const index_set& set, const index_vec& map);
  void  remap(const index_vec& map);
};

extern const index_set EMPTYSET;
typedef svector<index_pair> pair_set;

class bool_vec : public lvector<bool>
{
 public:
  // standard constructors
  bool_vec() : lvector<bool>() { };
  bool_vec(bool _val, index_type l) : lvector<bool>(_val, l) { };
  bool_vec(const lvector<bool>& _vec) : lvector<bool>(_vec) { };

  // construct from array + given size
  bool_vec(const bool* _arr, index_type n) : lvector<bool>(false, n) {
    for (index_type k = 0; k < n; k++) {
      if (_arr[k])
	(*this)[k] = true;
      else
	(*this)[k] = false;
    }
  };

  // conversion from index_set
  bool_vec(const index_set& set, index_type l);

  void complement();
  void insert(const lvector<bool>& vec);
  void insert(const index_set& set);
  void intersect(const lvector<bool>& vec);
  void intersect(const index_set& set);
  void subtract(const lvector<bool>& vec);
  void subtract(const index_set& set);
  bool subset(const lvector<bool>& vec) const;
  bool strict_subset(const lvector<bool>& vec) const;
  bool superset(const lvector<bool>& vec) const;
  bool strict_superset(const lvector<bool>& vec) const;
  bool contains(const bool& v) const;
  bool contains(const lvector<bool>& set) const;
  bool contains(const index_set& set) const;
  bool contains_any(const index_set& set) const;
  index_type first_common_element(const index_set& vec) const;
  index_type first_common_element(const lvector<bool>& vec) const;
  index_type count_common(const lvector<bool>& vec) const;
  index_type count_common(const index_set& set) const;
  index_set& copy_to(index_set& set) const;
  index_set& insert_into(index_set& set) const;
  index_set& subtract_from(index_set& set) const;
  bool* copy_to(bool* s, index_type n) const;
  int compare(const lvector<bool>& vec) const;
  index_type hash() const;
};

class index_set_vec : public lvector<index_set>
{
 public:
  index_set_vec()
    : lvector<index_set>() { };
  index_set_vec(const index_set& set, index_type l)
    : lvector<index_set>(set, l) { };
  index_set_vec(index_type l)
    : lvector<index_set>(EMPTYSET, l) { };
  index_set_vec(const index_set_vec& vec)
    : lvector<index_set>(vec) { };

  class decreasing_cardinality_order : public index_set_vec::order {
  public:
    virtual bool operator()
      (const index_set& v0, const index_set& v1) const
      { return (v0.size() > v1.size()); };
  };

  decreasing_cardinality_order decreasing_cardinality;

  index_type minimum_cardinality() const;
  index_type maximum_cardinality() const;
  index_type selected_minimum_cardinality(const index_set& sel) const;
  index_type selected_maximum_cardinality(const index_set& sel) const;
  index_type first_minimum_cardinality_set() const;
  index_type first_maximum_cardinality_set() const;

  index_type first_superset(const index_set& set) const;
  index_type first_superset(const bool_vec& set) const;
  index_type first_strict_superset(const index_set& set) const;
  index_type first_contains(index_type e) const;
  index_type next_contains(index_type e, index_type p) const;
  index_type first_subset(const index_set& set) const;
  index_type first_subset(const bool_vec& set) const;
  index_type first_strict_subset(const index_set& set) const;

  index_set& union_set(index_set& set) const;
  index_set& selected_union_set(const index_set& sel, index_set& set) const;
  index_set& intersection_set(index_set& set) const;

  void maximal(bool_vec& is_max);
  void minimal(bool_vec& is_min);

  void insert_maximal(const index_set& set);
  void insert_minimal(const index_set& set);
  void reduce_to_maximal();
  void reduce_to_minimal();

  void append_if_not_subset(const index_set& set);
  void append_if_not_superset(const index_set& set);
  void append_if_new(const index_set& set);
  void remove_sets_size_le(index_type l);
  void remove_empty_sets();

  void insert_in_all(index_type i);
  void insert_in_all(const index_set& set);
  void subtract_from_all(index_type i);
  void subtract_from_all(const index_set& set);

  void combinations_by_union(const index_set_vec& sv);
  void combinations_by_union(const index_set_vec& sv1,
			     const index_set_vec& sv2);
};

template<class T> class matrix : public lvector< lvector<T> >
{
 public:
  typedef lvector<T> row_type;

  matrix()
    : lvector<row_type>() { };
  matrix(const T& _val, index_type r, index_type c)
    : lvector<row_type>(row_type(_val, c), r) { };
  matrix(const matrix& _mat)
    : lvector<row_type>(_mat) { };

  // submatrix constructors
  matrix(const matrix& _mat, const bool_vec& rs, const bool_vec& cs);
  matrix(const matrix& _mat, const index_set& rs, const index_set& cs);

  index_type rows() const
  {
    return lvector<row_type>::length();
  };

  index_type columns() const
  {
    if (lvector<row_type>::length() == 0) return 0;
    else return (*this)[0].length();
  };

  void set_size(index_type r, index_type c);
  void assign_value(const T& _val);
  void assign_value(const T& _val, index_type r, index_type c);
  void append_row(const row_type& _row, const T& _vfill);
  void append_column(const row_type& _col, const T& _vfill);
  void delete_row(index_type i);
  void delete_column(index_type i);
  void delete_rows(const bool_vec& s);
  void delete_columns(const bool_vec& s);

  void column(index_type i, lvector<T>& c) const;
};

class bool_matrix : public matrix<bool> {
 public:
  bool_matrix() : matrix<bool>() { };
  bool_matrix(const bool& _val, index_type r, index_type c)
    : matrix<bool>(_val, r, c) { };
  bool_matrix(const bool_matrix& _mat) : matrix<bool>(_mat) { };
  bool_matrix(const bool_matrix& _mat, const bool_vec& rs, const bool_vec& cs)
    : matrix<bool>(_mat, rs, cs) { };
  bool_matrix(const bool_matrix& _mat, const index_set& rs, const index_set& cs)
    : matrix<bool>(_mat, rs, cs) { };


  void complement();
  // "insert" means union (can't use name "union" since it's a keyword,
  // and insert is consistent with bool_vec/index_set terminology)
  void insert(const bool_matrix& m);
  void intersect(const bool_matrix& m);
  void subtract(const bool_matrix& m);
  // this = m0 x m1
  void multiply(const bool_matrix& m0, const bool_matrix& m1);
  void transitive_closure();

  void and_of_rows(lvector<bool>& v) const;
  void and_of_rows(const bool_vec& s, lvector<bool>& v) const;
  void or_of_rows(lvector<bool>& v) const;
  void or_of_rows(const bool_vec& s, lvector<bool>& v) const;
  void and_of_columns(lvector<bool>& v) const;
  void and_of_columns(const bool_vec& s, lvector<bool>& v) const;
  void or_of_columns(lvector<bool>& v) const;
  void or_of_columns(const bool_vec& s, lvector<bool>& v) const;
};

typedef matrix<index_type> index_matrix;

class mapping : public index_vec
{
 public:
  // construct an identity mapping of size n
  static void identity_map
    (index_type n, index_vec& map)
    { index_vec_util::fill(map, n); };

  // extend domain of map to n (if map domain already goes to n or beyond,
  // this has no effect); new entries are identity.
  static void extend_domain(index_type n, index_vec& map);

  // construct inverse mapping, if possible
  static bool invert_map
    (const index_vec& map, index_vec& inv, index_type m = 0);

  // construct a map that deletes i:th of n elements
  // (maps 0..i-1 to identity, i to nil, and i+1..n to identity - 1)
  static void delete_index_map
    (index_type n, index_type i, index_vec& map);

  // construct a subset map (maps s[i] to i, j to nil for all j not in s)
  static void subset_map
    (index_type n, const index_set s, index_vec& map);
  
  // construct an all-to-i map
  static void all_to_one_map
    (index_type n, index_type i, index_vec& map);

  // compose two mappings: cm[x] = m1[m0[x]]; clearly, this assumes
  // the domain of m1 contains the range of m0; using the same object
  // for the first input (m0) and the result (cm) is safe
  static void compose
    (const index_vec& m0, const index_vec& m1, index_vec& cm);

  // remap a mapping: this is the result of composing m1^-1 * m0 * m1,
  // i.e., rm[x] = m1[m0[m1^-1[x]]]; the domain and range of rm is the
  // range of m1; the domain of m0 must contain the domain of m1, and
  // the domain of m1 must contain the range of m0, and m1 must be
  // invertible. The same objects can NOT be used for inputs (m0 or m1)
  // and the result (rm).
  static bool remap
    (const index_vec& m0, const index_vec& m1, index_vec& rm);

  // compute image under mapping of a vector of indices
  static void map_image
    (const index_vec& map, const index_vec& vec, index_vec& img);

  // compute image under mapping of a set of indices
  static void map_image
    (const index_vec& map, const index_set& s, index_set& img);

  // compute inverse image under mapping of a single index (the
  // inverse image may be a set, since mappings can be many to one)
  static void inverse_map_image
    (const index_vec& map, index_type x, index_set& img);

  // compute inverse image under mapping of a set of indices
  static void inverse_map_image
    (const index_vec& map, const index_set& x, index_set& img);

  // compute the range, i.e. the highest index in the map image + 1, or
  // equivalently, the number of distinct indices that can appear (but
  // not necessarily do appear) in the map image, corresponding to
  // domain 0..d - 1; map must be defined up to d - 1
  // (i.e. map.length() >= d).
  static index_type range(const index_vec& map, index_type d);

  // compute the set of values that actually appear in the range
  static void range_set(const index_vec& map, index_set& s);

  // construct a "compacted" map: this takes a map, renames the N
  // indices that appear in the range set to 0..N-1, and creates a
  // new mapping to 0..N-1
  static void compact(const index_vec& m0, index_vec& cm);

  mapping()
    : index_vec() { };
  mapping(index_type n) // identity map constructor
    : index_vec() { identity_map(n, *this); };
  // constructor for two kinds of map: if out == true, construct a
  // deleted index ("all but i") map; else construct an "all to i" map
  mapping(index_type n, index_type i, bool out) : index_vec() {
    if (out) delete_index_map(n, i, *this); else assign_value(i, n);
  };
  mapping(const index_vec& map)
    : index_vec(map) { };

  // assign identity mapping of size n to this
  void assign_identity(index_type n)
    { identity_map(n, *this); };

  // extend domain of this mapping to n; new entries are identity-mapped.
  void extend_domain_to(index_type n)
    { extend_domain(n, *this); };

  // apply mapping to a single index
  index_type operator()(index_type x) const
    { assert(x < size()); return (*this)[x]; };

  // apply mapping to a vector of indices
  index_vec operator()(const index_vec& vec) const
    { index_vec res; map_image(*this, vec, res); return res; };

  // apply inverse mapping to a single index
  index_vec& inverse(index_type x, index_set& res) const
    { inverse_map_image(*this, x, res); return res; };

  // apply inverse mapping to a set of indices
  index_vec& inverse(const index_set& x, index_set& res) const
    { inverse_map_image(*this, x, res); return res; };

  // assign inverse of this mapping to rmap
  bool invert(index_vec& rmap) const
    { return invert_map(*this, rmap); };

  // assign inverse of this mapping to this
  bool invert()
    { index_vec tmp(*this); return invert_map(tmp, *this); };

  index_type range() const
    { return range(*this, length()); };
};

class sparse_mapping : public pair_vec
{
 public:
  static void dense_to_sparse(const index_vec& dm, pair_vec sm);
  static void sparse_to_dense(const pair_vec& sm, index_vec dm);

  static index_type map_image
    (const pair_vec& map, index_type x);
  static void map_image
    (const pair_vec& map, const index_vec& vec, index_vec& img);
  static void inverse_map_image
    (const pair_vec& map, index_type x, index_set& img);
  static void inverse_map_image
    (const pair_vec& map, const index_set& x, index_set& img);

  sparse_mapping()
    : pair_vec() { };
  sparse_mapping(const pair_vec& m)
    : pair_vec(m) { };
  sparse_mapping(const index_vec& m)
    : pair_vec() { dense_to_sparse(m, *this); };

  index_type operator()(index_type x) const
    { return map_image(*this, x); };
  index_vec operator()(const index_vec& vec) const
    { index_vec res; map_image(*this, vec, res); return res; };
  index_set& inverse(index_type x, index_set& res) const
    { inverse_map_image(*this, x, res); return res; };
  index_set& inverse(const index_set& x, index_set& res) const
    { inverse_map_image(*this, x, res); return res; };
};


class equivalence : public index_vec
{
 public:
  equivalence()
    : index_vec() { };
  equivalence(index_type n)
    : index_vec(no_such_index, n) { index_vec_util::fill(*this, n); };
  equivalence(const equivalence& eq)
    : index_vec(eq) { };

  bool operator()(index_type a, index_type b) const;
  index_type canonical(index_type a) const;

  void extend(index_type a);
  void merge(index_type a, index_type b);
  void merge(const equivalence& eq);
  void merge(const index_set& set);
  void merge(const index_set& sa, const index_set& sb);
  void reset();
  void reset(index_type n);

  void canonical_set(index_set& set) const;
  void canonical_elements(index_set& set) const;
  void class_elements(index_type rep, index_set& set) const;
  index_type n_class_elements(index_type rep) const;
  void classes(index_set_vec& sets) const;
  index_type class_number(index_vec& cn) const;
  index_type class_number_and_size(index_vec& cn, index_vec& cs) const;
  void make_map(index_vec& map) const;
  index_type n_classes() const;
  index_type n_squeezed() const;
};


class set_hash_function : index_vec
{
 public:
  set_hash_function() : index_vec() { };
  set_hash_function(index_type n) : index_vec() { init(n); };

  void init(index_type n);

  index_type operator()(index_type& i) const;
  index_type operator()(const index_set& set) const;
  index_type operator()(const bool_vec& set) const;
  index_type operator()(const bool* set, index_type n) const;

  // incremental hashing
  index_type operator()(index_type& i, index_type v) const;
  index_type operator()(const index_set& set, index_type v) const;
  index_type operator()(const bool_vec& set, index_type v) const;
};


class s2index : public index_vec {
 public:
  s2index() : index_vec(0, 0) { };
  s2index(index_type n) : index_vec(0, n) { init(n); };

  void init(index_type n) {
    index_vec::assign_value(0, n);
    for (index_type i = 1; i < n; i++)
      (*this)[i] = ((*this)[i - 1] + i);
  };

  void extend_to(index_type n) {
    if (n <= index_vec::length()) return;
    index_type n0 = index_vec::length();
    index_vec::set_length(n);
    if (n0 == 0) {
      (*this)[0] = 0;
      n0 = 1;
    }
    for (index_type i = n0; i < n; i++)
      (*this)[i] = ((*this)[i - 1] + i);
  };

  index_type operator()(index_type i, index_type j) const {
    assert((i < index_vec::length()) && (j < index_vec::length()));
    if (i > j) {
      return ((*this)[i] + j);
    }
    else {
      return ((*this)[j] + i);
    }
  };

  index_pair inverse(index_type k) const {
    assert(k < n_pairs());
    for (index_type i = 0; (i + 1) < index_vec::length(); i++)
      if (((*this)[i] <= k) && (k < (*this)[i+1])) {
	return index_pair(i, k - (*this)[i]);
      }
    return index_pair(index_vec::length() - 1,
		      k - (*this)[index_vec::length() - 1]);
  };

  bool defined(index_type i) const {
    return (i < index_vec::length());
  };

  index_type n_pairs() const {
    if (index_vec::length() > 0)
      return (*this)[index_vec::length() - 1] + index_vec::length();
    else
      return 0;
  }
};

// Symmetric Binary Relation
class SBR : bool_vec {
  s2index ix;
 public:
  SBR() : bool_vec(false, 0), ix(0) { };
  SBR(index_type n) : ix(n) {
    bool_vec::assign_value(false, ix.n_pairs());
  };

  void extend_to(index_type n) {
    ix.extend_to(n);
    index_type n0 = bool_vec::length();
    bool_vec::set_length(ix.n_pairs());
    for (index_type i = n0; i < bool_vec::length(); i++)
      (*this)[i] = false;
  };

  bool operator()(index_type i, index_type j) const {
    return (*this)[ix(i, j)];
  };

  void set(index_type i, index_type j) {
    (*this)[ix(i, j)] = true;
  };

  void clear(index_type i, index_type j) {
    (*this)[ix(i, j)] = false;
  };

  void clear() {
    bool_vec::assign_value(false);
  };
};

///
/// old implement of s2index
///
// class s2index : public index_vec {
//  public:
//   s2index() : index_vec(0, 0) { };
//   s2index(index_type n) : index_vec(0, n) { init(n); };
// 
//   void init(index_type n) {
//     index_vec::assign_value(0, n);
//     for (HSPS::index_type i = 1; i < n; i++)
//       (*this)[i] = ((*this)[i - 1] + (n - i) + 1);
//   };
// 
//   index_type operator()(index_type i, index_type j) const {
//     assert((i < index_vec::length()) && (j < index_vec::length()));
//     if (i < j) {
//       return ((*this)[i] + (j - i));
//     }
//     else {
//       return ((*this)[j] + (i - j));
//     }
//   };
// 
//   index_pair inverse(index_type k) const {
//     assert(k < n_pairs());
//     for (index_type i = 0; (i + 1) < index_vec::length(); i++)
//       if (((*this)[i] <= k) && (k < (*this)[i+1])) {
// 	return index_pair(i, (k - (*this)[i]) + i);
//       }
//     assert(k == (*this)[index_vec::length() - 1]);
//     return index_pair(index_vec::length() - 1, index_vec::length() - 1);
//   };
// 
//   index_type n_pairs() const {
//     if (index_vec::length() > 0)
//       return (*this)[index_vec::length() - 1] + 1;
//     else
//       return 0;
//   }
// };

template<class T, class N> struct weighted
{
  T value;
  N weight;

  weighted() : weight(0) { };
  weighted(const T& v) : value(v), weight(0) { };
  weighted(const T& v, const N& w) : value(v), weight(w) { };
  weighted(const weighted& w) : value(w.value), weight(w.weight) { };
  ~weighted() { };

  weighted& operator=(const T& v)
  {
    value = v;
    weight = 0;
    return *this;
  };

  weighted& operator=(const weighted& w)
  {
    value = w.value;
    weight = w.weight;
    return *this;
  };

  bool operator==(const weighted& w) const
  {
    return (value == w.value);
  };

  bool operator!=(const weighted& w) const
  {
    return (value != w.value);
  };

  bool operator<(const weighted& w) const
  {
    return (value < w.value);
  };

  bool operator<=(const weighted& w) const
  {
    return (value <= w.value);
  };

  bool operator>(const weighted& w) const
  {
    return (value > w.value);
  };

  bool operator>=(const weighted& w) const
  {
    return (value >= w.value);
  };
};

template<class T, class N> class weighted_vec
: public lvector< weighted<T, N> >
{
 public:

  class decreasing_weight_order : public lvector< weighted<T,N> >::order {
  public:
    virtual bool operator()
      (const weighted<T,N>& v0, const weighted<T,N>& v1) const
      { return (v0.weight > v1.weight); };
  };

  class increasing_weight_order : public lvector< weighted<T,N> >::order {
  public:
    virtual bool operator()
      (const weighted<T,N>& v0, const weighted<T,N>& v1) const
      { return (v0.weight < v1.weight); };
  };

  static class decreasing_weight_order decreasing;
  static class increasing_weight_order increasing;

  void insert_increasing(const weighted<T,N>& v);
  void insert_decreasing(const weighted<T,N>& v);
  void insert_increasing(const T& v, const N& w);
  void insert_decreasing(const T& v, const N& w);
};

template<class T, class N> class weighted_set
: public svector< weighted<T,N> >
{
 public:
  void insert(const T& v, const N& w);
  void insert(const T& v);

  index_type arg_max();
  index_type arg_min();
};


// inlines

template<class T>
lvector<T>::lvector(const lvector<T>& vec, const index_set& sel)
  : std::vector<T>()
{
  for (index_type k = 0; k < sel.size(); k++) {
    assert(sel[k] < vec.size());
    append(vec[sel[k]]);
  }
}

template<class T>
lvector<T>::lvector(const lvector<T>& vec, const bool_vec& sel)
  : std::vector<T>()
{
  assert(sel.size() == vec.size());
  for (index_type k = 0; k < vec.size(); k++)
    if (sel[k])
      append(vec[k]);
}

template<class T>
bool lvector<T>::operator==(const lvector& _vec) const
{
  if (lvector<T>::size() != _vec.size()) return false;
  for (index_type k = 0; k < lvector<T>::size(); k++)
    if (!((*this)[k] == _vec[k])) return false;
  return true;
}

template<class T>
bool lvector<T>::operator!=(const lvector& _vec) const
{
  if (*this == _vec) return false;
  else return true;
}

template<class T>
bool lvector<T>::operator<(const lvector& vec) const
{
  if (lvector<T>::size() < vec.size()) return true;
  else if (lvector<T>::size() > vec.size()) return false;
  else {
    for (index_type k = 0; k < lvector<T>::size(); k++) {
      if ((*this)[k] < vec[k]) return true;
      else if ((*this)[k] > vec[k]) return false;
    }
    return false;
  }
}

template<class T>
bool lvector<T>::operator<=(const lvector& vec) const
{
  if (lvector<T>::size() < vec.size()) return true;
  else if (lvector<T>::size() > vec.size()) return false;
  else {
    for (index_type k = 0; k < lvector<T>::size(); k++) {
      if ((*this)[k] < vec[k]) return true;
      else if ((*this)[k] > vec[k]) return false;
    }
    return true;
  }
}

template<class T>
bool lvector<T>::operator>(const lvector& vec) const
{
  if (lvector<T>::size() < vec.size()) return false;
  else if (lvector<T>::size() > vec.size()) return true;
  else {
    for (index_type k = 0; k < lvector<T>::size(); k++) {
      if ((*this)[k] < vec[k]) return false;
      else if ((*this)[k] > vec[k]) return true;
    }
    return false;
  }
}

template<class T>
bool lvector<T>::operator>=(const lvector& vec) const
{
  if (lvector<T>::size() < vec.size()) return false;
  else if (lvector<T>::size() > vec.size()) return true;
  else {
    for (index_type k = 0; k < lvector<T>::size(); k++) {
      if ((*this)[k] < vec[k]) return false;
      else if ((*this)[k] > vec[k]) return true;
    }
    return true;
  }
}

template<class T>
bool lvector<T>::contains(const T& v) const
{
  for (index_type k = 0; k < lvector<T>::size(); k++)
    if ((*this)[k] == v) return true;
  return false;
}

template<class T>
index_type lvector<T>::first(const T& v) const
{
  for (index_type k = 0; k < lvector<T>::size(); k++)
    if ((*this)[k] == v) return k;
  return no_such_index;
}

template<class T>
index_type lvector<T>::last(const T& v) const
{
  for (index_type k = lvector<T>::size(); k > 0; k--)
    if ((*this)[k - 1] == v) return k;
  return no_such_index;
}

template<class T>
index_type lvector<T>::next(const T& v, index_type p) const
{
  for (index_type k = p + 1; k < lvector<T>::size(); k++)
    if ((*this)[k] == v) return k;
  return no_such_index;
}

template<class T>
index_type lvector<T>::find(const T& v, bool_vec& s) const
{
  index_type n = 0;
  s.assign_value(false, lvector<T>::size());
  for (index_type k = 0; k < lvector<T>::size(); k++)
    if ((*this)[k] == v) {
      s[k] = true;
      n += 1;
    }
  return n;
}

template<class T>
index_type lvector<T>::count(const T& v) const
{
  index_type c = 0;
  for (index_type k = 0; k < lvector<T>::size(); k++)
    if ((*this)[k] == v) c += 1;
  return c;
}

template<class T>
index_type lvector<T>::length() const
{
  return std::vector<T>::size();
}

template<class T>
index_type lvector<T>::arg_max() const
{
  if (lvector<T>::empty()) return no_such_index;
  index_type m = 0;
  for (index_type k = 1; k < lvector<T>::size(); k++)
    if ((*this)[k] > (*this)[m]) m = k;
  return m;
}

template<class T>
index_type lvector<T>::arg_max(const index_set& s) const
{
  if (s.empty()) return no_such_index;
  if (s[0] >= lvector<T>::size()) return no_such_index;
  index_type m = s[0];
  for (index_type k = 1; k < s.size(); k++) {
    if (s[k] >= lvector<T>::size()) return m;
    if ((*this)[s[k]] > (*this)[m]) m = s[k];
  }
  return m;
}

template<class T>
index_type lvector<T>::arg_max(const bool_vec& s) const
{
  assert(s.size() >= lvector<T>::size());
  index_type m = no_such_index;
  for (index_type k = 0; k < lvector<T>::size(); k++)
    if (s[k]) {
      if (m == no_such_index) m = k;
      else if ((*this)[k] > (*this)[m]) m = k;
    }
  return m;
}

template<class T>
index_type lvector<T>::arg_min() const
{
  if (lvector<T>::empty()) return no_such_index;
  index_type m = 0;
  for (index_type k = 1; k < lvector<T>::size(); k++)
    if ((*this)[k] < (*this)[m]) m = k;
  return m;
}

template<class T>
index_type lvector<T>::arg_min(const index_set& s) const
{
  if (s.empty()) return no_such_index;
  if (s[0] >= lvector<T>::size()) return no_such_index;
  index_type m = s[0];
  for (index_type k = 1; k < s.size(); k++) {
    if (s[k] >= lvector<T>::size()) return m;
    if ((*this)[s[k]] < (*this)[m]) m = s[k];
  }
  return m;
}

template<class T>
index_type lvector<T>::arg_min(const bool_vec& s) const
{
  assert(s.size() >= lvector<T>::size());
  index_type m = no_such_index;
  for (index_type k = 0; k < lvector<T>::size(); k++)
    if (s[k]) {
      if (m == no_such_index) m = k;
      else if ((*this)[k] < (*this)[m]) m = k;
    }
  return m;
}

template<class T>
index_type lvector<T>::arg_first(const order& o) const
{
  if (lvector<T>::empty()) return no_such_index;
  index_type m = 0;
  for (index_type k = 1; k < lvector<T>::size(); k++)
    if (o((*this)[k], (*this)[m])) m = k;
  return m;
}

template<class T>
index_type lvector<T>::arg_last(const order& o) const
{
  if (lvector<T>::empty()) return no_such_index;
  index_type m = 0;
  for (index_type k = 1; k < lvector<T>::size(); k++)
    if (o((*this)[m], (*this)[k])) m = k;
  return m;
}

template<class T>
index_pair lvector<T>::first_common(const lvector<T>& vec) const
{
  for (index_type i = 0; i < lvector<T>::size(); i++) {
    for (index_type j = 0; j < vec.size(); j++)
      if ((*this)[i] == vec[j]) return index_pair(i, j);
  }
  return index_pair(no_such_index, no_such_index);
}

template<class T>
index_pair lvector<T>::next_common(const lvector<T>& vec, index_pair p) const
{
  index_type i = p.first;
  index_type j = p.second + 1;
  while (j < vec.size()) {
    if ((*this)[i] == vec[j])
      return index_pair(i, j);
    j += 1;
  }
  i += 1;
  while (i < lvector<T>::size()) {
    j = 0;
    while (j < vec.size()) {
      if ((*this)[i] == vec[j])
	return index_pair(i, j);
      j += 1;
    }
    i += 1;
  }
  return index_pair(no_such_index, no_such_index);
}

template<class T>
void lvector<T>::difference
(const lvector& v1, lvector& d0, lvector& d1)
{
  d0.assign_copy(*this);
  d1.assign_copy(v1);
  index_type i0 = 0;
  while (i0 < d0.size()) {
    index_type i1 = d1.first(d0[i0]);
    if (i1 != no_such_index) {
      d0.remove(i0);
      d1.remove(i1);
    }
    else {
      i0 += 1;
    }
  }
}

template<class T>
void lvector<T>::assign_copy(const lvector& _vec)
{
  std::vector<T>::resize(_vec.size());
  for (index_type k = 0; k < _vec.size(); k++)
    (*this)[k] = _vec[k];
}

template<class T>
void lvector<T>::assign_copy(const T* _arr, index_type n)
{
  std::vector<T>::resize(n);
  for (index_type k = 0; k < n; k++)
    (*this)[k] = _arr[k];
}

template<class T>
void lvector<T>::assign_value(const T& _val)
{
  for (index_type k = 0; k < lvector<T>::size(); k++)
    (*this)[k] = _val;
}

template<class T>
void lvector<T>::assign_value(const T& _val, index_type l)
{
  std::vector<T>::resize(l);
  for (index_type k = 0; k < lvector<T>::size(); k++)
    (*this)[k] = _val;
}

template<class T>
void lvector<T>::assign_remap(const lvector<T>& vec, const index_vec& map)
{
  assert(map.length() == vec.length());
  index_type m = 0;
  for (index_type k = 0; k < vec.length(); k++)
    if (map[k] != no_such_index)
      if (map[k] > m) m = map[k];
  set_length(m + 1);
  for (index_type k = 0; k < vec.length(); k++)
    if (map[k] != no_such_index)
      (*this)[map[k]] = vec[k];
}

template<class T>
void lvector<T>::remap(const index_vec& map)
{
  lvector v0(*this);
  assign_remap(v0, map);
}

template<class T>
void lvector<T>::assign_select(const lvector& _vec, const index_set& s)
{
  set_length(s.length());
  for (index_type k = 0; k < s.length(); k++)
    (*this)[k] = _vec[s[k]];
}

template<class T>
void lvector<T>::assign_select(const lvector& _vec, const bool_vec& s)
{
  clear();
  for (index_type k = 0; k < _vec.length(); k++)
    if (s[k]) append(_vec[k]);
}

template<class T>
const lvector<T>& lvector<T>::operator=(const lvector<T>& _vec)
{
  assign_copy(_vec);
  return _vec;
}

// template<class T>
// const T& lvector<T>::operator=(const T& _val)
// {
//   assign_value(_val);
//   return _val;
// }

template<class T>
void lvector<T>::set_length(index_type l)
{
  std::vector<T>::resize(l);
}

template<class T>
void lvector<T>::set_length(index_type l, const T& v)
{
  std::vector<T>::resize(l, v);
}

template<class T>
void lvector<T>::inc_length_to(index_type l)
{
  if (std::vector<T>::size() < l)
    std::vector<T>::resize(l);
}

template<class T>
void lvector<T>::inc_length_to(index_type l, const T& v)
{
  if (std::vector<T>::size() < l)
    std::vector<T>::resize(l, v);
}

template<class T>
index_type lvector<T>::inc_length(index_type d)
{
  std::vector<T>::resize(std::vector<T>::size() + d);
  return std::vector<T>::size();
}

template<class T>
index_type lvector<T>::inc_length(index_type d, const T& v)
{
  std::vector<T>::resize(std::vector<T>::size() + d, v);
  return std::vector<T>::size();
}

template<class T>
index_type lvector<T>::dec_length(index_type d)
{
  assert(std::vector<T>::size() >= d);
  std::vector<T>::resize(std::vector<T>::size() - d);
  return std::vector<T>::size();
}

template<class T>
void lvector<T>::clear()
{
  std::vector<T>::clear();
}

template<class T>
void lvector<T>::append(const T& v)
{
  std::vector<T>::push_back(v);
}

template<class T>
void lvector<T>::append(const lvector<T>& v)
{
  for (index_type k = 0; k < v.size(); k++) append(v[k]);
}

template<class T>
T& lvector<T>::append()
{
  T v;
  std::vector<T>::push_back(v);
  // std::vector<T>::resize(std::vector<T>::size() + 1);
  return (*this)[std::vector<T>::size() - 1];
}

template<class T>
void lvector<T>::insert(const T& v, index_type p)
{
  if (p < lvector<T>::size()) {
    std::vector<T>::insert(std::vector<T>::begin() + p, v);
  }
  else {
    std::vector<T>::resize(p + 1);
    (*this)[p] = v;
  }
}

template<class T>
index_type lvector<T>::insert_ordered(const T& v, const order& o, index_type f)
{
  assert(f <= lvector<T>::size());
  for (index_type k = f; k < lvector<T>::size(); k++) {
    if (o(v, (*this)[k])) {
      insert(v, k);
      return k;
    }
  }
  append(v);
  return (lvector<T>::size() - 1);
}

template<class T>
index_type lvector<T>::insert_ordered(const lvector& vec, const order& o)
{
  if (vec.empty()) return no_such_index;
  index_type p0 = insert_ordered(vec[0], o);
  for (index_type k = 1; k < vec.size(); k++) {
    index_type p1 = insert_ordered(vec[k], o);
    if (p1 < p0) p0 = p1;
  }
  return p0;
}

template<class T>
void lvector<T>::remove(index_type p)
{
  if (p < lvector<T>::size())
    std::vector<T>::erase(std::vector<T>::begin() + p);
}

template<class T>
void lvector<T>::remove(index_type p0, index_type p1)
{
  assert(p0 < p1);
  if (p1 < lvector<T>::size())
    std::vector<T>::erase(std::vector<T>::begin() + p0,
			  std::vector<T>::begin() + p1);
  else
    std::vector<T>::erase(std::vector<T>::begin() + p0,
			  std::vector<T>::end());
}

template<class T>
void lvector<T>::remove(const bool_vec& s, index_vec& map)
{
  assert(s.size() >= lvector<T>::size());
  index_type scan_p = 0;
  index_type put_p = 0;
  index_vec rm_map(no_such_index, lvector<T>::size());
  while (scan_p < lvector<T>::size()) {
    if (!s[scan_p]) {
      if (put_p < scan_p) {
	(*this)[put_p] = (*this)[scan_p];
      }
      rm_map[scan_p] = put_p;
      put_p += 1;
    }
    else {
      rm_map[scan_p] = no_such_index;
    }
    scan_p += 1;
  }
  std::vector<T>::resize(put_p);

  for (index_type k = 0; k < map.size(); k++)
    if (map[k] != no_such_index) {
      assert(map[k] < rm_map.size());
      map[k] = rm_map[map[k]];
    }
}

template<class T>
void lvector<T>::remove(const bool_vec& s)
{
  assert(s.size() >= lvector<T>::size());
  index_type scan_p = 0;
  index_type put_p = 0;
  while (scan_p < lvector<T>::size()) {
    if (!s[scan_p]) {
      if (put_p < scan_p) {
	(*this)[put_p] = (*this)[scan_p];
      }
      put_p += 1;
    }
    scan_p += 1;
  }
  std::vector<T>::resize(put_p);
}

template<class T>
void lvector<T>::remove(const index_set& s)
{
  bool_vec s1(s, lvector<T>::size());
  remove(s1);
}

template<class T>
void lvector<T>::remove(const index_set& s, index_vec& map)
{
  bool_vec s1(s, std::vector<T>::size());
  lvector<T>::remove(s1, map);
}

template<class T>
void lvector<T>::remove_duplicate_elements()
{
  equivalence eq(lvector<T>::size());
  for (index_type i = 0; i < lvector<T>::size(); i++)
    for (index_type j = i+1; j < lvector<T>::size(); j++)
      if ((*this)[i] == (*this)[j])
	eq.merge(i, j);
  bool_vec s(false, lvector<T>::size());
  for (index_type i = 0; i < lvector<T>::size(); i++)
    if (eq.canonical(i) != i)
      s[i] = true;
  remove(s);
}

template<class T>
void lvector<T>::swap(index_type i, index_type j)
{
  T tmp = (*this)[i];
  (*this)[i] = (*this)[j];
  (*this)[j] = tmp;
}

template<class T>
void svector<T>::assign_singleton(const T& _val)
{
  lvector<T>::set_length(1);
  (*this)[0] = _val;
}

template<class T>
void svector<T>::assign_values(const lvector<T>& vec)
{
  lvector<T>::clear();
  for (index_type k = 0; k < vec.size(); k++)
    insert(vec[k]);
}

template<class T>
void svector<T>::insert(const T& v) {
  index_type i = 0;
  bool seeking = (i < std::vector<T>::size());
  while (seeking) {
    if ((*this)[i] < v) {
      i += 1;
      if (i >= std::vector<T>::size())
	seeking = false;
    }
    else {
      seeking = false;
    }
  }
  if (i < lvector<T>::size()) {
    if ((*this)[i] == v)
      return;
    else
      lvector<T>::insert(v, i);
  }
  else {
    lvector<T>::append(v);
  }
}

template<class T>
void svector<T>::insert(const svector<T>& vec)
{
  svector<T> tmp(*this);
  index_type r1 = 0;
  index_type r2 = 0;
  index_type w = 0;
  while ((r1 < tmp.size()) && (r2 < vec.size())) {
    if (tmp[r1] <= vec[r2]) {
      if (tmp[r1] == vec[r2]) r2 += 1;
      if (w < tmp.size())
	(*this)[w++] = tmp[r1++];
      else
	this->append(tmp[r1++]);
    }
    else {
      if (w < tmp.size())
	(*this)[w++] = vec[r2++];
      else
	this->append(vec[r2++]);
    }
  }
  if (r1 < tmp.size()) {
    assert(r2 >= vec.size());
    while (r1 < tmp.size()) {
      if (w < tmp.size())
	(*this)[w++] = tmp[r1++];
      else
	this->append(tmp[r1++]);
    }
  }
  else {
    while (r2 < vec.size()) {
      if (w < tmp.size())
	(*this)[w++] = vec[r2++];
      else
	this->append(vec[r2++]);
    }
  }
}

template<class T>
void svector<T>::insert(const lvector<T>& vec)
{
  for (index_type k = 0; k < vec.size(); k++) insert(vec[k]);
}

template<class T>
bool svector<T>::contains(const T& v) const
{
  index_type a = 0;
  index_type b = lvector<T>::size();
  while ((b - a) > 3) {
    index_type i = (a + b) / 2;
    if ((*this)[i] < v) a = i + 1;
    else if ((*this)[i] > v) b = i;
    else {
      // assert(contains2(v));
      return true;
    }
  }
  index_type i = a;
  while ((i < b) && ((*this)[i] < v)) i += 1;
  if (i < b)
    if ((*this)[i] == v) {
      // assert(contains2(v));
      return true;
    }
  // assert(!contains2(v));
  return false;
}

template<class T>
bool svector<T>::contains2(const T& v) const
{
  index_type i = 0;
  while ((i < lvector<T>::size()) &&
	 ((*this)[i] < v)) i += 1;
  if (i < lvector<T>::size())
    if ((*this)[i] == v) return true;
  return false;
}

template<class T>
bool svector<T>::contains(const svector& vec) const
{
  index_type v_i = 0;
  index_type i = 0;
  while (v_i < vec.size()) {
    if (i >= lvector<T>::size()) return false;
    if ((*this)[i] == vec[v_i]) {
      v_i += 1;
      i += 1;
    }
    else if ((*this)[i] > vec[v_i]) {
      return false;
    }
    else {
      while ((i < lvector<T>::size()) && ((*this)[i] < vec[v_i]))
	i += 1;
    }
  }
  return true;
}

template<class T>
bool svector<T>::subset(const svector& vec) const
{
  return vec.contains(*this);
}

template<class T>
void svector<T>::intersect(const svector& vec)
{
  index_type i = 0;
  while (i < lvector<T>::size()) {
    if (!vec.contains((*this)[i]))
      lvector<T>::remove(i);
    else
      i += 1;
  }
}

template<class T>
void svector<T>::difference(const svector& vec)
{
  svector d(vec);
  d.subtract(*this);
  subtract(vec);
  insert(d);
}

template<class T>
index_pair svector<T>::first_common(const lvector<T>& vec) const
{
  return lvector<T>::first_common(vec);
}

template<class T>
index_pair svector<T>::next_common(const lvector<T>& vec, index_pair p) const
{
  return lvector<T>::next_common(vec, p);
}

template<class T>
index_pair svector<T>::first_common(const svector<T>& vec) const
{
  index_type i = 0;
  index_type j = 0;
  while ((i < svector<T>::size()) && (j < vec.size())) {
    if ((*this)[i] == vec[j])
      return index_pair(i, j);
    else if
      ((*this)[i] < vec[j]) i += 1;
    else
      j += 1;
  }
  return index_pair(no_such_index, no_such_index);
}

template<class T>
index_pair svector<T>::next_common(const svector<T>& vec, index_pair p) const
{
  index_type i = p.first + 1;
  index_type j = p.second + 1;
  while ((i < svector<T>::size()) && (j < vec.size())) {
    if ((*this)[i] == vec[j])
      return index_pair(i, j);
    else if ((*this)[i] < vec[j])
      i += 1;
    else
      j += 1;
  }
  return index_pair(no_such_index, no_such_index);
}

template<class T>
index_type svector<T>::count_common(const svector& vec) const
{
  index_type i = 0;
  index_type j = 0;
  index_type c = 0;
  while ((i < svector<T>::size()) && (j < vec.size())) {
    if ((*this)[i] == vec[j]) {
      c += 1;
      i += 1;
      j += 1;
    }
    else if ((*this)[i] < vec[j])
      i += 1;
    else
      j += 1;
  }
  return c;
}

template<class T>
void svector<T>::subtract(const svector& vec)
{
  index_type i = 0;
  while (i < lvector<T>::size()) {
    if (vec.contains((*this)[i]))
      lvector<T>::remove(i);
    else
      i += 1;
  }
}

template<class T>
void svector<T>::subtract(const T& v)
{
  index_type i = 0;
  while (i < lvector<T>::size()) {
    if ((*this)[i] == v) {
      lvector<T>::remove(i);
      return;
    }
    else {
      i += 1;
    }
  }
}

template<class T>
void svector<T>::remove_greater_than(const T& _val)
{
  index_type i = lvector<T>::size();
  while (i > 0) {
    if ((*this)[i - 1] > _val)
      i -= 1;
    else {
      if (i < lvector<T>::size())
	lvector<T>::set_length(i);
      return;
    }
  }
  lvector<T>::set_length(0);
}

template<class T>
void svector<T>::remove_less_than(const T& _val)
{
  index_type i = 0;
  while (i < lvector<T>::size()) {
    if ((*this)[i] < _val)
      i += 1;
    else {
      if (i > 0) {
	for (index_type j = 0; (j + i) < lvector<T>::size(); j++)
	  (*this)[j] = (*this)[j + i];
	lvector<T>::set_length(lvector<T>::size() - i);
      }
      return;
    }
  }
  lvector<T>::set_length(0);
}

template<class T>
matrix<T>::matrix(const matrix& _mat, const bool_vec& rs, const bool_vec& cs)
{
  index_set rs1(rs);
  index_set cs1(cs);
  set_size(rs1.length(), cs1.length());
  for (index_type i = 0; i < rs1.length(); i++)
    for (index_type j = 0; j < cs1.length(); j++)
      (*this)[i][j] = _mat[rs1[i]][cs1[j]];
}

template<class T>
matrix<T>::matrix(const matrix& _mat, const index_set& rs, const index_set& cs)
{
  set_size(rs.length(), cs.length());
  for (index_type i = 0; i < rs.length(); i++)
    for (index_type j = 0; j < cs.length(); j++)
      (*this)[i][j] = _mat[rs[i]][cs[j]];
}

template<class T>
void matrix<T>::set_size(index_type r, index_type c)
{
  lvector<row_type>::set_length(r);
  for (index_type k = 0; k < lvector<row_type>::size(); k++)
    (*this)[k].set_length(c);
}

template<class T>
void matrix<T>::assign_value(const T& _val)
{
  for (index_type k = 0; k < lvector<row_type>::size(); k++)
    (*this)[k].assign_value(_val);
}

template<class T>
void matrix<T>::assign_value(const T& _val, index_type r, index_type c)
{
  lvector<row_type>::set_length(r);
  for (index_type k = 0; k < lvector<row_type>::size(); k++)
    (*this)[k].assign_value(_val, c);
}

template<class T>
void matrix<T>::append_row(const row_type& _row, const T& _vfill)
{
  if (lvector<row_type>::length() > 0) {
    index_type c = (*this)[0].length();
    if (c < _row.length()) {
      for (index_type k = 0; k < lvector<row_type>::length(); k++) {
	(*this)[k].set_length(_row.length());
	for (index_type i = c; i < _row.length(); i++)
	  (*this)[k][i] = _vfill;
      }
      append(_row);
    }
    else if (c > _row.length()) {
      index_type r = _row.length();
      append(_row);
      (*this)[lvector<row_type>::length() - 1].set_length(c);
      for (index_type i = r; i < c; i++)
	(*this)[lvector<row_type>::length() - 1][i] = _vfill;
    }
    else {
      append(_row);
    }
  }
  else {
    append(_row);
  }
}

template<class T>
void matrix<T>::append_column(const lvector<T>& _col, const T& _vfill)
{
  if (lvector<row_type>::length() > 0) {
    index_type c = (*this)[0].length();
    if (lvector<row_type>::length() < _col.length()) {
      index_type n = lvector<row_type>::length();
      set_length(_col.length());
      for (index_type k = n; k < lvector<row_type>::length(); k++)
	(*this)[k].assign_value(_vfill, c);
    }
    for (index_type k = 0; (k < lvector<row_type>::length()) &&
	   (k < _col.length()); k++)
      (*this)[k].append(_col[k]);
    for (index_type k = _col.length(); k < lvector<row_type>::length(); k++)
      (*this)[k].append(_vfill);
  }
  else {
    lvector<row_type>::set_length(_col.length());
    for (index_type k = 0; k < lvector<row_type>::length(); k++)
      (*this)[k].assign_value(_col[k], 1);
  }
}

template<class T>
void matrix<T>::delete_row(index_type i)
{
  lvector<row_type>::remove(i);
}

template<class T>
void matrix<T>::delete_column(index_type i)
{
  for (index_type k = 0; k < lvector<row_type>::length(); k++)
    (*this)[k].remove(i);
}

template<class T>
void matrix<T>::delete_rows(const bool_vec& s)
{
  lvector<row_type>::remove(s);
}

template<class T>
void matrix<T>::delete_columns(const bool_vec& s)
{
  for (index_type k = 0; k < lvector<row_type>::length(); k++)
    (*this)[k].remove(s);
}

template<class T>
void matrix<T>::column(index_type i, lvector<T>& c) const
{
  c.set_length(lvector<row_type>::length());
  for (index_type k = 0; k < lvector<row_type>::length(); k++)
    c[k] = (*this)[k][i];
}

template<class T>
void swapable_pair<T>::swap()
{
  T tmp = this->first;
  this->first = this->second;
  this->second = tmp;
}

template<class T>
void comparable_pair<T>::sort_ascending()
{
  if (this->first > this->second) swapable_pair<T>::swap();
}

template<class T>
void comparable_pair<T>::sort_descending()
{
  if (this->first < this->second) swapable_pair<T>::swap();
}

template<class T, class N>
class weighted_vec<T,N>::decreasing_weight_order
  weighted_vec<T,N>::decreasing;

template<class T, class N>
class weighted_vec<T,N>::increasing_weight_order
  weighted_vec<T,N>::increasing;

template<class T, class N>
void weighted_vec<T,N>::insert_increasing(const weighted<T,N>& v)
{
  this->insert_ordered(v, increasing);
}

template<class T, class N>
void weighted_vec<T,N>::insert_decreasing(const weighted<T,N>& v)
{
  this->insert_ordered(v, decreasing);
}

template<class T, class N>
void weighted_vec<T,N>::insert_increasing(const T& v, const N& w)
{
  this->insert_ordered(weighted<T,N>(v, w), increasing);
}

template<class T, class N>
void weighted_vec<T,N>::insert_decreasing(const T& v, const N& w)
{
  this->insert_ordered(weighted<T,N>(v, w), decreasing);
}

template<class T, class N>
void weighted_set<T,N>::insert(const T& v, const N& w)
{
  index_type p = svector< weighted<T,N> >::first(v);
  if (p == no_such_index) {
    svector< weighted<T,N> >::insert(weighted<T,N>(v, w));
  }
  else {
    (*this)[p].weight += w;
  }
}

template<class T, class N>
void weighted_set<T,N>::insert(const T& v)
{
  index_type p = svector< weighted<T,N> >::first(v);
  if (p == no_such_index) {
    svector< weighted<T,N> >::insert(weighted<T,N>(v, 1));
  }
  else {
    (*this)[p].weight += 1;
  }
}

template<class T, class N>
index_type weighted_set<T,N>::arg_max()
{
  if (weighted_set<T,N>::empty()) return no_such_index;
  index_type m = 0;
  for (index_type k = 1; k < weighted_set<T,N>::size(); k++)
    if ((*this)[k].weight > (*this)[m].weight) m = k;
  return m;
}

template<class T, class N>
index_type weighted_set<T,N>::arg_min()
{
  if (weighted_set<T,N>::empty()) return no_such_index;
  index_type m = 0;
  for (index_type k = 1; k < weighted_set<T,N>::size(); k++)
    if ((*this)[k].weight < (*this)[m].weight) m = k;
  return m;
}

template<class T>
inline std::ostream& operator<<(std::ostream& s, const swapable_pair<T>& p)
{
  return s << '(' << p.first << ',' << p.second << ')';
}

template<class T>
::std::ostream& operator<<(::std::ostream& s, const lvector<T>& _vec)
{
  s << '[';
  for (index_type k = 0; k < _vec.size(); k++) {
    if (k > 0) s << ',';
    s << _vec[k];
  }
  s << ']';
  return s;
}

inline std::ostream& operator<<(std::ostream& s, const mapping& m)
{
  s << '{';
  for (index_type k = 0; k < m.length(); k++) {
    if (k > 0) s << ',';
    s << k << '-' << '>';
    if (m[k] == no_such_index)
      s << '_';
    else
      s << m[k];
  }
  return s << '}';
}

inline std::ostream& operator<<(std::ostream& s, const equivalence& eq)
{
  s << '{';
  bool first = true;
  for (index_type k = 0; k < eq.length(); k++) {
    index_type c = eq.canonical(k);
    if (!first) {
      s << ',';
    }
    else {
      first = false;
    }
    s << k << '=' << c;
  }
  return s << '}';
}

template<class T, class N>
std::ostream& operator<<(::std::ostream& s, const weighted<T,N>& w)
{
  s << '<' << w.value << ':' << w.weight << '>';
}

END_HSPS_NAMESPACE

#endif
