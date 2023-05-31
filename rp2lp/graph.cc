
#include "graph.h"

BEGIN_HSPS_NAMESPACE

// construct empty graph
adjacency_list_graph::adjacency_list_graph()
{
}

// construct empty graph with s nodes
adjacency_list_graph::adjacency_list_graph(index_type s)
  : _size(s), _edges(EMPTYSET, s)
{
}

// copy constructor
adjacency_list_graph::adjacency_list_graph(const adjacency_list_graph& g)
  : _size(g._size), _edges(g._edges)
{
}

// subgraph constructor
adjacency_list_graph::adjacency_list_graph
(const adjacency_list_graph& g, const index_set& n)
  : _size(n.size()), _edges(EMPTYSET, n.size())
{
  assert(n.size() <= g._size);
  for (index_type i = 0; i < n.size(); i++) {
    index_type s = n[i];
    for (index_type j = 0; j < g._edges[s].size(); j++) {
      index_type d = n.first(g._edges[s][j]);
      if (d != no_such_index)
	_edges[i].append(d);
    }
  }
}

// copy-and-map constructor
adjacency_list_graph::adjacency_list_graph
(const adjacency_list_graph& g, const index_vec& m)
  : _size(0)
{
  assert(m.size() >= g._size);
  for (index_type i = 0; i < g._size; i++)
    if ((m[i] != no_such_index) && (m[i] > _size)) _size = m[i];
  _edges.assign_value(EMPTYSET, _size);
  for (index_type i = 0; i < g._size; i++)
    if (m[i] != no_such_index)
      _edges[m[i]].assign_remap(g._edges[i], m);
}

// quotient graph constructor
adjacency_list_graph::adjacency_list_graph
(const adjacency_list_graph& g, const equivalence& eq)
{
  assert(0);
}

// copy constructor from general graph
adjacency_list_graph::adjacency_list_graph(const graph& g)
  : _size(0), _edges(EMPTYSET, 0)
{
  copy(g);
}

// subgraph constructor from general graph
adjacency_list_graph::adjacency_list_graph
(const graph& g, const index_set& n)
{
  index_vec m;
  mapping::subset_map(g.size(), n, m);
  copy(g, m);
}

// return number of edges
index_type adjacency_list_graph::n_edges() const
{
  index_type n = 0;
  for (index_type i = 0; i < _size; i++)
    n += _edges[i].size();
  return n;
}

// return number of edges between nodes in two sets
index_type adjacency_list_graph::n_edges
(const index_set& from, const index_set& to) const
{
  index_type n = 0;
  for (index_type i = 0; i < from.size(); i++)
    n += _edges[from[i]].count_common(to);
  return n;
}

void adjacency_list_graph::distance(index_type s0, index_vec& d) const
{
  d.assign_value(no_such_index, size());
  d[s0] = 0;
  bool done = false;
  while (!done) {
    done = true;
    for (index_type i = 0; i < size(); i++)
      for (index_type j = 0; j < _edges[i].size(); j++)
	if (d[i] != no_such_index) {
	  if ((d[i] + 1) < d[_edges[i][j]]) {
	    d[_edges[i][j]] = d[i] + 1;
	    done = false;
	  }
	}
  }
}

index_type adjacency_list_graph::shortest_path
(index_type s, index_type t, index_vec& p) const
{
  index_vec d;
  distance(s, d);
  return extract_path(d, t, p);
}

index_type adjacency_list_graph::extract_path
(const index_vec& d, index_type n, index_vec& p) const
{
  index_type l = d[n];
  if (l == no_such_index) {
    p.clear();
    return l;
  }
  p.set_length(l + 1);
  p[l] = n;
  index_type r = l;
  while (r > 0) {
    index_type nn = no_such_index;
    for (index_type i = 0; (i < _size) && (nn == no_such_index); i++)
      if (d[i] == (r - 1))
	if (_edges[i].contains(n))
	  nn = i;
    assert(nn != no_such_index);
    n = nn;
    p[r - 1] = n;
    r = r - 1;
  }
  return l;
}

// assign empty graph with s nodes
void adjacency_list_graph::init(index_type s)
{
  _size = s;
  _edges.assign_value(EMPTYSET, _size);
}

// copy from general graph
void adjacency_list_graph::copy(const graph& g)
{
  _size = g.size();
  _edges.set_length(_size);
  for (index_type i = 0; i < _size; i++)
    _edges[i] = g.successors(i);
}

// mapped copy from general graph
void adjacency_list_graph::copy
(const graph& g, const index_vec& map)
{
  assert(map.size() >= g.size());
  _size = 0;
  for (index_type i = 0; i < g.size(); i++)
    if ((map[i] != no_such_index) && (map[i] > _size)) _size = map[i];
  _edges.assign_value(EMPTYSET, _size);
  for (index_type i = 0; i < g.size(); i++)
    if (map[i] != no_such_index)
      _edges[map[i]].assign_remap(g.successors(i), map);
}

// add an isolated node, return index of the new node (== new size - 1)
index_type adjacency_list_graph::add_node()
{
  index_type i = _size++;
  _edges.append(EMPTYSET);
  assert(_edges.size() == _size);
  return i;
}

// remove one or more nodes
void adjacency_list_graph::remove_node(index_type n)
{
  assert(n < _size);
  _edges.remove(n);
  for (index_type i = 0; i < _edges.size(); i++) {
    _edges[i].subtract(n);
    for (index_type j = 0; j < _edges[i].size(); j++)
      if (_edges[i][j] > n)
	_edges[i][j] = (_edges[i][j] - 1);
  }
  _size = _size - 1;
  assert(_size == _edges.size());
}

// add g as a subgraph; on return, m maps indices in g to indices in this
void adjacency_list_graph::add_graph
(const graph& g, mapping& m)
{
  m.set_length(g.size());
  _edges.set_length(_size + g.size());
  for (index_type i = 0; i < g.size(); i++) {
    _edges[_size + i] = g.successors(i);
    for (index_type j = 0; j < _edges[_size + i].size(); j++)
      _edges[_size + i][j] = (_edges[_size + i][j] + _size);
  }
  _size = (_size + g.size());
  assert(_size == _edges.size());
}

// add/remove edges:
void adjacency_list_graph::add_edge
(index_type src, index_type dst)
{
  assert(src < _size);
  assert(dst < _size);
  _edges[src].insert(dst);
}

void adjacency_list_graph::add_edge
(const index_set& srcs, index_type dst)
{
  assert(dst < _size);
  for (index_type i = 0; i < srcs.size(); i++) {
    assert(srcs[i] < _size);
    _edges[srcs[i]].insert(dst);
  }
}

void adjacency_list_graph::add_edge
(index_type src, const index_set& dsts)
{
  assert(src < _size);
  _edges[src].insert(dsts);
}

void adjacency_list_graph::add_edge_to_transitive_closure
(index_type src, index_type dst)
{
  assert(src < _size);
  assert(dst < _size);
  _edges[src].insert(dst);
  _edges[src].insert(_edges[dst]);
  for (index_type i = 0; i < _size; i++)
    if (_edges[i].contains(src)) {
      _edges[i].insert(dst);
      _edges[i].insert(_edges[dst]);
    }
}

void adjacency_list_graph::remove_edge
(index_type src, index_type dst)
{
  assert(src < _size);
  _edges[src].subtract(dst);
}

void adjacency_list_graph::remove_edges_from(index_type src)
{
  assert(src < _size);
  _edges[src].clear();
}

void adjacency_list_graph::remove_edges_to(index_type dst)
{
  for (index_type i = 0; i < _size; i++)
    _edges[i].subtract(dst);
}

void adjacency_list_graph::remove_edges_incident_on(index_type n)
{
  assert(n < _size);
  _edges[n].clear();
  for (index_type i = 0; i < _size; i++)
    _edges[i].subtract(n);
}

void adjacency_list_graph::clear_edges()
{
  for (index_type i = 0; i < _edges.size(); i++)
    _edges[i].clear();
}

void adjacency_list_graph::write_compact(::std::ostream& s) const
{
  s << "{";
  bool first = true;
  for (index_type i = 0; i < _edges.size(); i++)
    for (index_type j = 0; j < _edges[i].size(); j++) {
      if (!first) s << ",";
      s << i << "->" << _edges[i][j];
      first = false;
    }
  s << "}";
}

graph::graph()
  : _size(0), comp()
{
  // done
}

graph::graph(index_type s)
  : _size(0), comp()
{
  init(s);
}

graph::graph(const graph& g)
  : _size(0), comp()
{
  copy(g);
}

graph::graph(const graph& g, const index_set& n)
  : _size(0), comp()
{
  g.subgraph(*this, n);
}

graph::graph(const graph& g, const index_vec& m)
  : _size(0), comp()
{
  copy(g, m);
}

graph::graph(const graph& g, const equivalence& eq)
  : _size(0), comp()
{
  g.quotient(*this, eq);
}

graph::~graph()
{
  // done
}

bool graph::empty() const
{
  for (index_type i = 0; i < _size; i++)
    for (index_type j = 0; j < _size; j++)
      if (adj[i][j]) return false;
  return true;
}

bool graph::reachable(index_type n0, index_type n1) const
{
  assert(n0 < _size);
  assert(n1 < _size);
  bool_vec v(false, _size);
  reachable(n0, v);
  return v[n1];
}

index_type graph::count_reachable(index_type n0) const
{
  bool_vec v(false, _size);
  reachable(n0, v);
  return v.count(true);
}

void graph::reachable(bool_vec& v) const
{
  for (index_type k = 0; k < _size; k++)
    if (v[k])
      reachable(k, v);
}

void graph::descendants(index_type n0, bool_vec& s) const
{
  s.assign_value(false, _size);
  reachable(n0, s);
}

void graph::descendants(const index_set& s0, bool_vec& s) const
{
  s.assign_value(false, _size);
  for (index_type k = 0; k < s0.length(); k++)
    reachable(s0[k], s);
}

void graph::descendants(index_type n0, index_set& s) const
{
  bool_vec v(false, _size);
  reachable(n0, v);
  v.copy_to(s);
}

void graph::descendants(const index_set& s0, index_set& s) const
{
  bool_vec v(false, _size);
  for (index_type k = 0; k < s0.length(); k++)
    reachable(s0[k], v);
  v.copy_to(s);
}

void graph::nearest_common_descendants
(index_type n0, index_type n1, bool_vec& s) const
{
  descendants(n0, s);
  bool_vec s1;
  descendants(n1, s1);
  s.intersect(s1);
  s1.assign_copy(s);
  for (index_type i = 0; i < _size; i++)
    if (s[i]) {
      if (in[i].have_common_element(s1))
	s[i] = false;
    }
}

void graph::ancestors(index_type n0, bool_vec& s) const
{
  s.assign_value(false, _size);
  reverse_reachable(n0, s);
}

void graph::ancestors(const index_set& s0, bool_vec& s) const
{
  s.assign_value(false, _size);
  for (index_type k = 0; k < s0.length(); k++)
    reverse_reachable(s0[k], s);
}

void graph::ancestors(index_type n0, index_set& s) const
{
  bool_vec v(false, _size);
  reverse_reachable(n0, v);
  v.copy_to(s);
}

void graph::ancestors(const index_set& s0, index_set& s) const
{
  bool_vec v(false, _size);
  for (index_type k = 0; k < s0.length(); k++)
    reverse_reachable(s0[k], v);
  v.copy_to(s);
}

void graph::between
(index_type n0, index_type n1, bool_vec& s) const
{
  descendants(n0, s);
  bool_vec s1;
  ancestors(n1, s1);
  s.intersect(s1);
}

bool graph::acyclic() const
{
  bool_vec _reach(false, _size);
  for (index_type k = 0; k < _size; k++) {
    _reach.assign_value(false, _size);
    for (index_type i = 0; i < out[k].length(); i++)
      reachable(out[k][i], _reach);
    if (_reach[k]) return false;
  }
  return true;
}

bool graph::top_sort(index_vec& s) const
{
  s.clear();
  // std::cerr << "g = " << *this << std::endl;
  bool_vec rem(true, size());
  while (rem.count(true) > 0) {
    // std::cerr << "rem = " << rem << std::endl;
    index_type n = no_such_index;
    for (index_type i = 0; (i < size()) && (n == no_such_index); i++)
      if (rem[i] && (in[i].count_common(rem) == 0))
	n = i;
    // std::cerr << "n = " << n << std::endl;
    if (n == no_such_index)
      return false;
    s.append(n);
    rem[n] = false;
  }
  return true;
}

index_type graph::max_out_degree() const
{
  assert(_size > 0);
  index_type m = out[0].length();
  for (index_type k = 1; k < _size; k++)
    if (out[k].length() > m) m = out[k].length();
  return m;
}

index_type graph::max_in_degree() const
{
  assert(_size > 0);
  index_type m = in[0].length();
  for (index_type k = 1; k < _size; k++)
    if (in[k].length() > m) m = in[k].length();
  return m;
}

index_type graph::max_bi_degree() const
{
  assert(_size > 0);
  index_type m = bi[0].length();
  for (index_type k = 1; k < _size; k++)
    if (bi[k].length() > m) m = bi[k].length();
  return m;
}

index_type graph::max_bi_degree_node() const
{
  assert(_size > 0);
  index_type m = bi[0].length();
  index_type n = 0;
  for (index_type k = 1; k < _size; k++)
    if (bi[k].length() > m) {
      m = bi[k].length();
      n = k;
    }
  return n;
}

index_type graph::min_out_degree() const
{
  assert(_size > 0);
  index_type m = out[0].length();
  for (index_type k = 1; k < _size; k++)
    if (out[k].length() < m) m = out[k].length();
  return m;
}

index_type graph::min_in_degree() const
{
  assert(_size > 0);
  index_type m = in[0].length();
  for (index_type k = 1; k < _size; k++)
    if (in[k].length() < m) m = in[k].length();
  return m;
}

index_type graph::min_bi_degree() const
{
  assert(_size > 0);
  index_type m = bi[0].length();
  for (index_type k = 1; k < _size; k++)
    if (bi[k].length() < m) m = bi[k].length();
  return m;
}

index_type graph::min_bi_degree_node() const
{
  assert(_size > 0);
  index_type m = bi[0].length();
  index_type n = 0;
  for (index_type k = 1; k < _size; k++)
    if (bi[k].length() < m) {
      m = bi[k].length();
      n = k;
    }
  return n;
}

index_type graph::min_fillin_node() const
{
  assert(_size > 0);
  index_type m = fillin(bi[0]);
  index_type n = 0;
  for (index_type k = 1; k < _size; k++) {
    index_type f = fillin(bi[k]);
    if (f < m) {
      m = f;
      n = k;
    }
  }
  return n;
}

index_type graph::first_root() const
{
  for (index_type k = 0; k < size(); k++)
    if (in_degree(k) == 0) return k;
  return no_such_index;
}

index_type graph::first_leaf() const
{
  for (index_type k = 0; k < size(); k++)
    if (out_degree(k) == 0) return k;
  return no_such_index;
}

index_type graph::first_undirected_leaf() const
{
  for (index_type k = 0; k < size(); k++)
    if ((in[k].size() == 1) && (out[k].size() == 1))
      if (in[k][0] == out[k][0])
	return k;
  return no_such_index;
}

index_type graph::next_undirected_leaf(index_type l) const
{
  for (index_type k = l + 1; k < size(); k++)
    if ((in[k].size() == 1) && (out[k].size() == 1))
      if (in[0] == out[0])
	return k;
  return no_such_index;
}

index_type graph::first_undirected_leaf(const bool_vec& n) const
{
  for (index_type k = 0; k < size(); k++)
    if (n[k])
      if ((in[k].size() == 1) && (out[k].size() == 1))
	if (in[k][0] == out[k][0])
	  return k;
  return no_such_index;
}

index_type graph::first_undirected_leaf(const index_set& n) const
{
  for (index_type k = 0; k < n.size(); k++)
    if ((in[n[k]].size() == 1) && (out[n[k]].size() == 1))
      if (in[n[k]][0] == out[n[k]][0])
	return n[k];
  return no_such_index;
}

index_type graph::next_undirected_leaf(index_type l, const index_set& n) const
{
  index_type p = n.first(l);
  assert(p < n.size());
  for (index_type k = p + 1; k < n.size(); k++)
    if ((in[n[k]].size() == 1) && (out[n[k]].size() == 1))
      if (in[n[k]][0] == out[n[k]][0])
	return n[k];
  return no_such_index;
}

void graph::fringe(const index_set& n, index_set& fn) const
{
  fn.clear();
  for (index_type k = 0; k < n.length(); k++) {
    fn.insert(successors(n[k]));
  }
  fn.subtract(n);
}

void graph::bi_fringe(const index_set& n, index_set& fn) const
{
  fn.clear();
  for (index_type k = 0; k < n.length(); k++) {
    fn.insert(bidirectional(n[k]));
  }
  fn.subtract(n);
}

void graph::distance(index_type s0, index_vec& d) const
{
  index_set s0s;
  s0s.assign_singleton(s0);
  distance(s0s, d);
}

void graph::distance(const index_set& s0, index_vec& d) const
{
  d.assign_value(no_such_index, size());
  for (index_type k = 0; k < s0.length(); k++) {
    assert(s0[k] < size());
    d[s0[k]] = 0;
  }
  bool done = false;
  while (!done) {
    done = true;
    for (index_type i = 0; i < size(); i++)
      for (index_type j = 0; j < out[i].length(); j++)
	if (d[i] != no_such_index) {
	  if ((d[i] + 1) < d[out[i][j]]) {
	    d[out[i][j]] = d[i] + 1;
	    done = false;
	  }
	}
  }
}

index_type graph::distance(index_type s0, index_type s1) const
{
  index_vec d;
  distance(s0, d);
  return d[s1];
}

index_type graph::diameter() const
{
  index_type d_max = 0;
  index_vec d;
  for (index_type k = 0; k < size(); k++) {
    distance(k, d);
    index_type i = d.arg_max();
    if (d[i] > d_max)
      d_max = d[i];
  }
  return d_max;
}

index_type graph::extract_path
(const index_vec& d, index_type n, index_vec& p) const
{
  index_type l = d[n];
  if (l == no_such_index) {
    p.clear();
    return l;
  }
  p.set_length(l + 1);
  p[l] = n;
  index_type r = l;
  while (r > 0) {
    index_type nn = no_such_index;
    for (index_type i = 0; (i < in[n].size()) && (nn == no_such_index); i++)
      if (d[in[n][i]] == (r - 1))
	nn = in[n][i];
    assert(nn != no_such_index);
    n = nn;
    p[r - 1] = n;
    r = r - 1;
  }
  return l;
}

index_type graph::shortest_path
(index_type s, index_type t, index_vec& p) const
{
  index_vec d;
  distance(s, d);
  return extract_path(d, t, p);
}

index_type graph::shortest_path
(const index_set& s, const index_set& t, index_vec& p) const
{
  index_vec d;
  distance(s, d);
  index_type n = d.arg_min(t);
  assert(t.contains(n));
  return extract_path(d, n, p);
}

index_type graph::shortest_cycle(index_vec& p) const
{
  index_vec d;
  index_type l_min = no_such_index;
  for (index_type i = 0; i < _size; i++) {
    distance(out[i], d);
    if (d[i] != no_such_index)  {
      if ((l_min == no_such_index) || ((d[i] + 1) < l_min)) {
	extract_path(d, i, p);
	l_min = (d[i] + 1);
      }
    }
  }
  return l_min;
}

pair_set& graph::bidirectional_edges(pair_set& s) const
{
  s.clear();
  for (index_type i = 0; i < _size; i++)
    for (index_type j = i+1; j < _size; j++)
      if (adj[i][j] && adj[j][i]) s.insert(index_pair(i, j));
  return s;
}

bool graph::adjacent(index_type i, const index_set& n) const
{
  for (index_type j = 0; j < n.length(); j++)
    if (adjacent(i, n[j])) return true;
  return false;
}

bool graph::adjacent(const index_set& n, index_type i) const
{
  for (index_type j = 0; j < n.length(); j++)
    if (adjacent(n[j], i)) return true;
  return false;
}

bool graph::adjacent(const index_set& n0, const index_set& n1) const
{
  for (index_type i = 0; i < n0.length(); i++)
    for (index_type j = 0; j < n1.length(); j++)
      if (adjacent(n0[i], n1[j])) return true;
  return false;
}

bool graph::bi_adjacent(index_type i, const index_set& n) const
{
  for (index_type j = 0; j < n.length(); j++)
    if (bi_adjacent(i, n[j])) return true;
  return false;
}

index_type graph::n_edges() const
{
  index_type n = 0;
  for (index_type i = 0; i < _size; i++)
    for (index_type j = 0; j < _size; j++)
      if (adj[i][j]) n += 1;
  return n;
}

index_type graph::n_edges(const index_set& from, const index_set& to) const
{
  index_type n = 0;
  for (index_type i = 0; i < from.length(); i++)
    for (index_type j = 0; j < to.length(); j++)
      if (adj[from[i]][to[j]]) n += 1;
  return n;
}

pair_set& graph::edges(pair_set& s) const
{
  s.clear();
  for (index_type i = 0; i < _size; i++)
    for (index_type j = 0; j < _size; j++)
      if (adj[i][j]) s.insert(index_pair(i, j));
  return s;
}

index_type graph::n_induced_undirected_edges() const
{
  index_type n = 0;
  for (index_type i = 0; i < _size; i++)
    for (index_type j = i; j < _size; j++)
      if (adj[i][j] || adj[j][i]) n += 1;
  return n;
}

index_type graph::n_induced_undirected_edges
(const index_set& n0, const index_set& n1) const
{
  index_type n = 0;
  for (index_type i = 0; i < n0.length(); i++)
    for (index_type j = 0; j < n1.length(); j++)
      if (adj[n0[i]][n1[j]] || adj[n1[j]][n0[i]]) n += 1;
  return n;
}

index_type graph::n_bidirectional_edges() const
{
  index_type n = 0;
  for (index_type i = 0; i < _size; i++)
    for (index_type j = i; j < _size; j++)
      if (adj[i][j] && adj[j][i]) n += 1;
  return n;
}

index_type graph::n_bidirectional_edges
(const index_set& n0, const index_set& n1) const
{
  index_type n = 0;
  for (index_type i = 0; i < n0.length(); i++)
    for (index_type j = 0; j < n1.length(); j++)
      if (adj[n0[i]][n1[j]] && adj[n1[j]][n0[i]]) n += 1;
  return n;
}

bool graph::assign_node_level_top_down(index_vec& levels) const
{
  levels.assign_value(no_such_index, _size);
  for (index_type k = 0; k < _size; k++)
    if (in[k].empty())
      levels[k] = 0;
  bool done = false;
  while (!done) {
    done = true;
    for (index_type k = 0; k < _size; k++)
      if (levels[k] == no_such_index) {
	index_type l = 0;
	for (index_type i = 0; i < in[k].size(); i++)
	  if (levels[in[k][i]] == no_such_index)
	    l = no_such_index;
	  else if (levels[in[k][i]] > l)
	    l = levels[in[k][i]];
	if (l != no_such_index) {
	  levels[k] = l + 1;
	  done = false;
	}
      }
  }
  return (levels.first(no_such_index) == no_such_index);
}

index_type graph::component_node(index_type i) const
{
  for (index_type k = 0; k < _size; k++)
    if (comp[k] == i) return k;
  return no_such_index;
}

index_type graph::component_size(index_type i) const
{
  index_type n = 0;
  for (index_type k = 0; k < _size; k++)
    if (comp[k] == i) n += 1;
  return n;
}

void graph::component_node_set(index_type i, index_set& set) const
{
  set.clear();
  for (index_type k = 0; k < _size; k++)
    if (comp[k] == i) set.insert(k);
}

void graph::component_node_sets(index_set_vec& sets) const
{
  sets.assign_value(EMPTYSET, n_comp);
  for (index_type k = 0; k < _size; k++)
    sets[comp[k]].insert(k);
}

index_type graph::maximal_non_unit_component() const
{
  index_type k_max = no_such_index;
  index_type s_max = 1;
  for (index_type k = 0; k < n_components(); k++)
    if (component_size(k) > s_max) {
      k_max = k;
      s_max = component_size(k);
    }
  return k_max;
}

bool graph::equals(const graph& g) const
{
  if (g.size() != size()) return false;
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j) != g.adjacent(i, j)) return false;
  return true;
}

bool graph::equals(const graph& g, const index_vec& c) const
{
  if (g.size() != size()) return false;
  if (c.length() != g.size()) {
    ::std::cerr << "error: size of graph " << g
	      << " and permutation vector " << c
	      << " do not agree"
	      << ::std::endl;
    exit(255);
  }
  for (index_type i = 0; i < size(); i++) {
    assert(c[i] < g.size());
    for (index_type j = 0; j < size(); j++) {
      assert(c[j] < g.size());
      if (adjacent(i, j) != g.adjacent(c[i], c[j])) return false;
    }
  }
  return true;
}

void graph::difference
(const graph& g, const index_vec& c, pair_set& d0, pair_set& d1) const
{
  assert(g.size() == size());
  assert(c.length() == g.size());
  d0.clear();
  d1.clear();
  for (index_type i = 0; i < size(); i++) {
    assert(c[i] < g.size());
    for (index_type j = 0; j < size(); j++) {
      assert(c[j] < g.size());
      if (adjacent(i, j) && !g.adjacent(c[i], c[j]))
	d0.insert(index_pair(i, j));
      if (!adjacent(i, j) && g.adjacent(c[i], c[j]))
	d1.insert(index_pair(c[i], c[j]));
    }
  }
}

void graph::difference
(const graph& g, pair_set& d0, pair_set& d1) const
{
  assert(g.size() == size());
  d0.clear();
  d1.clear();
  for (index_type i = 0; i < size(); i++) {
    for (index_type j = 0; j < size(); j++) {
      if (adjacent(i, j) && !g.adjacent(i, j))
	d0.insert(index_pair(i, j));
      if (!adjacent(i, j) && g.adjacent(i, j))
	d1.insert(index_pair(i, j));
    }
  }
}

index_type graph::cardinality_of_difference(const graph& g) const
{
  assert(g.size() == size());
  index_type d = 0;
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adj[i][j] != g.adj[i][j]) d += 1;
  return d;
}

int graph::compare(const graph& g) const
{
  if (_size < g._size)
    return -1;
  else if (_size > g._size)
    return 1;
  else {
    for (index_type i = 0; i < _size; i++)
      for (index_type j = 0; j < _size; j++)
	if (!adj[i][j] && g.adj[i][j])
	  return -1;
	else if (!adj[i][j] && g.adj[i][j])
	  return 1;
    return 0;
  }
}

index_type graph::hash(const set_hash_function& f) const
{
  index_type h = 0;
  for (index_type i = 0; i < _size; i++)
    h = f(out[i], h);
  return h;
}

void graph::init(index_type size)
{
  _size = size;
  adj.assign_value(false, _size, _size);
  in.assign_value(EMPTYSET, _size);
  out.assign_value(EMPTYSET, _size);
  bi.assign_value(EMPTYSET, _size);
  comp.assign_value(0, _size);
  n_comp = 0;
}

void graph::copy(const graph& g)
{
  _size = g._size;
  adj.assign_copy(g.adj);
  in.assign_copy(g.in);
  out.assign_copy(g.out);
  bi.assign_copy(g.bi);
  comp.assign_copy(g.comp);
  n_comp = g.n_comp;
}

void graph::copy(const graph& g, const index_vec& map)
{
  assert(map.length() == g.size());
  index_type m = mapping::range(map, map.length());
  init(m);
  for (index_type i = 0; i < g.size(); i++)
    if (map[i] != no_such_index) {
      assert(map[i] < size());
      for (index_type j = 0; j < g.size(); j++)
	if (map[j] != no_such_index) {
	  assert(map[j] < _size);
	  if (g.adjacent(i, j))
	    add_edge(map[i], map[j]);
	}
    }
}

void graph::copy_and_rename(const graph& g, const index_vec& map)
{
  assert(map.length() == g.size());
  init(g.size());
  for (index_type i = 0; i < _size; i++) {
    assert(map[i] < _size);
    for (index_type j = 0; j < _size; j++) {
      assert(map[j] < _size);
      if (g.adjacent(i, j))
	add_edge(map[i], map[j]);
    }
  }
}

index_type graph::add_node()
{
  adj.set_size(_size + 1, _size + 1);
  for (index_type k = 0; k < _size; k++) {
    adj[k][_size] = false;
    adj[_size][k] = false;
  }
  adj[_size][_size] = false;

  in.set_length(_size + 1);
  in[_size].clear();
  out.set_length(_size + 1);
  out[_size].clear();
  bi.set_length(_size + 1);
  bi[_size].clear();

  comp.set_length(_size + 1);
  comp[_size] = 0;

  _size += 1;
  return _size - 1;
}

void graph::add_graph(const graph& g, mapping& m)
{
  index_type _new_size = _size + g.size();
  adj.set_size(_new_size, _new_size);
  in.set_length(_new_size);
  out.set_length(_new_size);
  bi.set_length(_new_size);
  comp.set_length(_new_size);

  m.assign_identity(g.size());
  for (index_type k = 0; k < g.size(); k++) {
    m[k] = _size + k;
    for (index_type i = 0; i < _new_size; i++) {
      adj[i][m[k]] = false;
      adj[m[k]][i] = false;
    }
    in[m[k]].clear();
    out[m[k]].clear();
    bi[m[k]].clear();
    comp[m[k]] = 0;
  }
  _size = _new_size;

  for (index_type i = 0; i < g.size(); i++)
    for (index_type j = 0; j < g.size(); j++)
      if (g.adjacent(i, j))
	add_edge(m[i], m[j]);
}

void graph::add_edge(index_type src, index_type dst)
{
  assert((src < _size) && (dst < _size));
  if (adj[src][dst]) return;
  adj[src][dst] = true;
  out[src].insert(dst);
  in[dst].insert(src);
  if (adj[dst][src]) {
    bi[src].insert(dst);
    bi[dst].insert(src);
  }
}

void graph::add_edge(const index_set& srcs, index_type dst)
{
  for (index_type k = 0; k < srcs.length(); k++)
    if (!adj[srcs[k]][dst])
      add_edge(srcs[k], dst);
}

void graph::add_edge(index_type src, const index_set& dsts)
{
  for (index_type k = 0; k < dsts.length(); k++)
    if (!adj[src][dsts[k]])
      add_edge(src, dsts[k]);
}

void graph::add_edge_to_transitive_closure
(index_type src, index_type dst)
{
  if (!adj[src][dst]) {
    add_edge(src, dst);
  }
  for (index_type i = 0; i < in[src].length(); i++) {
    if (!adj[in[src][i]][dst]) {
      add_edge(in[src][i], dst);
    }
    for (index_type j = 0; j < out[dst].length(); j++) {
      if (!adj[in[src][i]][out[dst][j]]) {
	add_edge(in[src][i], out[dst][j]);
      }
    }
  }
  for (index_type j = 0; j < out[dst].length(); j++) {
    if (!adj[src][out[dst][j]]) {
      add_edge(src, out[dst][j]);
    }
  }
}

void graph::add_edge_to_transitive_closure
(index_type src, index_type dst, pair_set& e)
{
  if (!adj[src][dst]) {
    add_edge(src, dst);
    e.insert(index_pair(src, dst));
  }
  for (index_type i = 0; i < in[src].length(); i++) {
    if (!adj[in[src][i]][dst]) {
      e.insert(index_pair(in[src][i], dst));
    }
    for (index_type j = 0; j < out[dst].length(); j++) {
      if (!adj[in[src][i]][out[dst][j]]) {
	e.insert(index_pair(in[src][i], out[dst][j]));
      }
    }
  }
  for (index_type j = 0; j < out[dst].length(); j++) {
    if (!adj[src][out[dst][j]]) {
      e.insert(index_pair(src, out[dst][j]));
    }
  }
  for (index_type k = 0; k < e.length(); k++)
    add_edge(e[k].first, e[k].second);
}

void graph::remove_edge(index_type src, index_type dst)
{
  assert((src < _size) && (dst < _size));
  if (adj[src][dst]) {
    adj[src][dst] = false;
    out[src].subtract(dst);
    in[dst].subtract(src);
    if (adj[dst][src]) {
      bi[src].subtract(dst);
      bi[dst].subtract(src);
    }
  }
}

void graph::remove_undirected_edge(index_type n0, index_type n1)
{
  assert((n0 < _size) && (n1 < _size));
  if (adj[n0][n1]) {
    adj[n0][n1] = false;
    out[n0].subtract(n1);
    in[n1].subtract(n0);
    bi[n0].subtract(n1);
  }
  if (adj[n1][n0]) {
    adj[n1][n0] = false;
    out[n1].subtract(n0);
    in[n0].subtract(n1);
    bi[n1].subtract(n0);
  }
}

void graph::remove_edges_from(index_type src)
{
  index_set ns(successors(src));
  for (index_type k = 0; k < ns.length(); k++)
    remove_edge(src, ns[k]);
}

void graph::remove_edges_to(index_type dst)
{
  index_set ns(predecessors(dst));
  for (index_type k = 0; k < ns.length(); k++)
    remove_edge(ns[k], dst);
}

void graph::remove_edges_incident_on(index_type n)
{
  remove_edges_from(n);
  remove_edges_to(n);
}

void graph::remove_edges(const pair_set& e)
{
  for (index_type k = 0; k < e.length(); k++)
    remove_edge(e[k].first, e[k].second);
}

void graph::remove_undirected_edges(const pair_set& e)
{
  for (index_type k = 0; k < e.length(); k++)
    remove_undirected_edge(e[k].first, e[k].second);
}

void graph::add_undirected_edge(index_type n0, index_type n1)
{
  add_edge(n0, n1);
  add_edge(n1, n0);
}

void graph::add_undirected_edges(const index_set& n0)
{
  for (index_type i = 0; i < n0.size(); i++)
    for (index_type j = i + 1; j < n0.size(); j++)
      add_undirected_edge(n0[i], n0[j]);
}

void graph::add_undirected_edges(const index_set& n0, const index_set& n1)
{
  for (index_type i = 0; i < n0.size(); i++)
    for (index_type j = 0; j < n1.size(); j++)
      if (n0[i] != n1[j])
	add_undirected_edge(n0[i], n1[j]);
}

void graph::add_undirected_edges(const pair_set& e)
{
  for (index_type k = 0; k < e.length(); k++)
    add_undirected_edge(e[k].first, e[k].second);
}

void graph::remove_node(index_type n)
{
  assert(n < _size);
  adj.delete_row(n);
  adj.delete_column(n);
  _size -= 1;
  recalculate();
//   graph g(*this);
//   init(g.size() - 1);
//   for (index_type i = 0; i < g.size(); i++)
//     for (index_type j = 0; j < g.size(); j++)
//       if (g.adjacent(i, j) && (i != n) && (j != n)) {
// 	index_type s = (i > n ? i - 1 : i);
// 	index_type d = (j > n ? j - 1 : j);
// 	add_edge(s, d);
//       }
}

void graph::remove_nodes(const index_set& ns)
{
  graph g(*this);
  // complement of ns
  index_set cns(size(), ns);
  g.subgraph(*this, cns);
}

void graph::clear_edges()
{
  adj.assign_value(false);
  in.assign_value(EMPTYSET);
  out.assign_value(EMPTYSET);
  bi.assign_value(EMPTYSET);
}

void graph::recalculate()
{
  for (index_type i = 0; i < _size; i++) {
    in[i].clear();
    out[i].clear();
  }
  for (index_type i = 0; i < _size; i++)
    for (index_type j = 0; j < _size; j++)
      if (adj[i][j]) {
	out[i].insert(j);
	in[j].insert(i);
      }
  for (index_type k = 0; k < _size; k++) {
    bi[k].assign_copy(out[k]);
    bi[k].intersect(in[k]);
  }
  comp.set_length(_size);
  comp.assign_value(0);
  n_comp = 0;
}

void graph::complement()
{
  for (index_type i = 0; i < _size; i++) {
    for (index_type j = 0; j < _size; j++)
      adj[i][j] = !adj[i][j];
    adj[i][i] = false;
  }
  recalculate();
}

void graph::complement_with_loops()
{
  for (index_type i = 0; i < _size; i++) {
    for (index_type j = 0; j < _size; j++)
      adj[i][j] = !adj[i][j];
  }
  recalculate();
}

void graph::remove_loops()
{
  for (index_type i = 0; i < _size; i++) if (adj[i][i]) {
    adj[i][i] = false;
    out[i].subtract(i);
    in[i].subtract(i);
    bi[i].subtract(i);
  }
}

void graph::make_undirected()
{
  for (index_type i = 0; i < _size; i++)
    for (index_type j = 0; j < _size; j++)
      if (adj[i][j] && !adj[j][i]) add_edge(j, i);
}

void graph::reverse()
{
  for (index_type i = 0; i < _size; i++)
    for (index_type j = i+1; j < _size; j++) {
      if (adj[i][j] && !adj[j][i]) {
	adj[i][j] = false;
	adj[j][i] = true;
      }
      else if (!adj[i][j] && adj[j][i]) {
	adj[i][j] = true;
	adj[j][i] = false;
      }
    }
  recalculate();
}

void graph::transitive_closure()
{
  for (index_type k = 0; k < _size; k++)
    for (index_type i = 0; i < _size; i++)
      if (adj[i][k])
	for (index_type j = 0; j < _size; j++)
	  if (adj[k][j]) adj[i][j] = true;
  recalculate();
}

void graph::missing_transitive_edges(pair_set& e) const
{
  e.clear();
  for (index_type i = 0; i < _size; i++) {
    bool_vec v(false, _size);
    reachable(i, v);
    for (index_type j = 0; j < _size; j++)
      if ((i != j) && v[j] && (!adj[i][j]))
	e.insert(index_pair(i, j));
  }
}

void graph::transitive_reduction()
{
  bool_matrix m2(adj);
  m2.transitive_closure();
  bool_matrix m3;
  m3.multiply(adj, m2);
  adj.subtract(m3);
  recalculate();
}

void graph::intersect(const graph& g)
{
  if (g._size != _size) {
    ::std::cerr << "error: can't intersect " << *this << " of size " << _size
	      << " with graph " << g << " of size" << g._size << ::std::endl;
    exit(255);
  }
  for (index_type i = 0; i < _size; i++)
    for (index_type j = 0; j < _size; j++)
      if (!g.adj[i][j]) adj[i][j] = false;
  recalculate();
}

graph& graph::subgraph(graph& sg, const index_set& nodes) const
{
  sg.init(nodes.length());
  for (index_type i = 0; i < nodes.length(); i++)
    for (index_type j = 0; j < nodes.length(); j++)
      if (adj[nodes[i]][nodes[j]]) sg.add_edge(i, j);
  return sg;
}

graph& graph::edge_subgraph(graph& sg, const index_set& nodes) const
{
  sg.init(size());
  for (index_type i = 0; i < nodes.length(); i++)
    for (index_type j = 0; j < nodes.length(); j++)
      if (adj[nodes[i]][nodes[j]]) sg.add_edge(nodes[i], nodes[j]);
  return sg;
}

graph& graph::undirected_edge_graph(graph& g, pair_set& es) const
{
  es.clear();
  for (index_type i = 0; i < size(); i++)
    for (index_type j = i + 1; j < size(); j++)
      if (bi_adjacent(i, j))
	es.insert(index_pair(i, j));
  g.init(es.size());
  for (index_type i = 0; i < es.size(); i++)
    for (index_type j = i + 1; j < es.size(); j++)
      if ((es[i].first == es[j].first) ||
	  (es[i].first == es[j].second) ||
	  (es[i].second == es[j].first) ||
	  (es[i].second == es[j].second))
	g.add_undirected_edge(i, j);
}

graph& graph::component_tree(graph& cg) const
{
  cg.init(n_comp);
  for (index_type i = 0; i < _size; i++)
    for (index_type j = 0; j < _size; j++)
      if (adj[i][j] && (comp[i] != comp[j]))
	cg.add_edge(comp[i], comp[j]);
  return cg;
}

equivalence& graph::component_partitioning(equivalence& eq) const
{
  eq.reset(_size);
  for (index_type i = 0; i < _size; i++)
    for (index_type j = i + 1; j < _size; j++)
      if ((comp[i] == comp[j])) {
	eq.merge(i, j);
      }
  return eq;
}

void graph::induced_partitioning(equivalence& eq) const
{
  eq.clear();
  eq.extend(_size);
  for (index_type i = 0; i < _size; i++)
    for (index_type j = 0; j < _size; j++)
      if (adj[i][j]) eq.merge(i, j);
}

void graph::induced_partitioning(index_set_vec& cs) const
{
  equivalence p;
  induced_partitioning(p);
  p.classes(cs);
}

graph& graph::induced_undirected_graph(graph& g) const
{
  g.copy(*this);
  g.make_undirected();
  return g;
}

graph& graph::minimal_equivalent_digraph(graph& g) const
{
  graph cg;
  component_tree(cg);
  cg.transitive_reduction();
  g.init(size());
  for (index_type k = 0; k < n_components(); k++) {
    index_set c;
    component_node_set(k, c);
    if (c.length() > 1) {
      for (index_type i = 0; i < (c.length() - 1); i++)
	g.add_edge(c[i], c[i+1]);
      g.add_edge(c[c.length() - 1], c[0]);
    }
  }
  for (index_type i = 0; i < cg.size(); i++)
    for (index_type j = 0; j < cg.size(); j++)
      if (cg.adjacent(i, j)) {
	index_type ci = component_node(i);
	index_type cj = component_node(j);
	assert(ci != no_such_index);
	assert(cj != no_such_index);
	g.add_edge(ci, cj);
      }
  return g;
}

graph& graph::minimal_distance_graph(graph& g, const index_set& s0) const
{
  index_vec d;
  distance(s0, d);
  g.init(size());
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j) && (d[i] != no_such_index) && (d[j] == (d[i] + 1)))
	g.add_edge(i, j);
  return g;
}

graph& graph::quotient(graph& g, const equivalence& eq) const
{
  g.init(eq.n_classes());
  index_vec m;
  eq.make_map(m);
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j) && !eq(i, j)) {
	assert(m[i] < eq.n_classes());
	assert(m[j] < eq.n_classes());
	g.add_edge(m[i], m[j]);
      }
}

graph& graph::complete_npart(graph& g, const index_set_vec& parts) const
{
  g.init(parts.length());
  for (index_type p = 0; p < parts.length(); p++)
    for (index_type q = p + 1; q < parts.length(); q++) {
      bool is_complete = true;
      for (index_type ip = 0; (ip < parts[p].length()) && is_complete; ip++)
	for (index_type iq = 0; (iq < parts[q].length()) && is_complete; iq++)
	  if (parts[p][ip] != parts[q][iq])
	    if (!bi_adjacent(parts[p][ip], parts[q][iq]))
	      is_complete = false;
      if (is_complete)
	g.add_undirected_edge(p, q);
    }
}

bool graph::is_clique(const index_set& nodes) const
{
  for (index_type i = 0; i < nodes.length(); i++)
    for (index_type j = i + 1; j < nodes.length(); j++)
      if (!adj[nodes[i]][nodes[j]] || !adj[nodes[j]][nodes[i]])
	return false;
  return true;
}

index_type graph::fillin(const index_set& nodes) const
{
  index_type f = 0;
  for (index_type i = 0; i < nodes.length(); i++)
    for (index_type j = i + 1; j < nodes.length(); j++)
      if (!bi_adjacent(nodes[i], nodes[j])) f += 1;
  return f;
}

bool graph::is_independent(const index_set& nodes) const
{
  for (index_type i = 0; i < nodes.length(); i++)
    for (index_type j = i + 1; j < nodes.length(); j++)
      if (adj[nodes[i]][nodes[j]] || adj[nodes[j]][nodes[i]])
	return false;
  return true;
}

bool graph::is_independent_range(index_type l, index_type u) const
{
  for (index_type i = l; i <= u; i++)
    for (index_type j = i + 1; j <= u; j++)
      if (adjacent(i, j) || adjacent(j, i)) return false;
  return true;
}

void graph::write_node_set(::std::ostream& s) const
{
  s << '{';
  for (index_type i = 0; i < _size; i++) {
    if (i > 0) s << ',';
    s << i << "[" << comp[i] << "]";
  }
  s << '}';
}

void graph::write_edge_set(::std::ostream& s) const
{
  s << '{';
  bool first = true;
  for (index_type i = 0; i < _size; i++)
    for (index_type j = 0; j < _size; j++)
      if (adj[i][j]) {
	if (!first) s << ',';
	first = false;
	s << i << "->" << j;
      }
  s << '}';
}

void graph::write_compact(::std::ostream& s) const
{
  s << '(';
  write_node_set(s);
  s << ',';
  write_edge_set(s);
  s << '}';
}

void graph::write_undirected_edge_set(::std::ostream& s) const
{
  s << '{';
  bool first = true;
  for (index_type i = 0; i < _size; i++)
    for (index_type j = i+1; j < _size; j++)
      if (adj[i][j]) {
	if (!first) s << ',';
	first = false;
	s << i << "-" << j;
      }
  s << '}';
}

void graph::write_adjacency_lists(::std::ostream& s) const
{
  for (index_type k = 0; k < _size; k++) {
    s << "node " << k << " (component " << comp[k] << "):" << ::std::endl;
    s << " in (" << in[k].length() << "):";
    for (index_type i = 0; i < _size; i++)
      if (adj[i][k]) s << " " << i;
    s << ::std::endl << " out (" << out[k].length() << "):"; 
    for (index_type i = 0; i < _size; i++)
      if (adj[k][i]) s << " " << i;
    s << ::std::endl;
  }
}

void graph::dump(::std::ostream& s) const
{
  s << "size = " << _size << std::endl;
  s << "adj = " << adj << std::endl;
  s << "in = " << in << std::endl;
  s << "out = " << out << std::endl;
  s << "bi = " << bi << std::endl;
  s << "comp = " << comp << std::endl;
}

void graph::randomize(count_type n, RNG& rnd)
{
  assert(_size > 1);
  for (index_type s = 0; s < n; s++) {
    index_type i = rnd.random_in_range(_size);
    index_type j = rnd.random_in_range(_size);
    if (i != j) {
      if (adj[i][j]) {
	adj[i][j] = false;
      }
      else {
	adj[i][j] = true;
      }
    }
  }
  recalculate();
}

void graph::randomize_undirected(count_type n, RNG& rnd)
{
  assert(_size > 1);
  for (index_type s = 0; s < n; s++) {
    index_type i = rnd.random_in_range(_size);
    index_type j = rnd.random_in_range(_size);
    if (i != j) {
      if (adj[i][j]) {
	adj[i][j] = false;
	adj[j][i] = false;
      }
      else {
	adj[i][j] = true;
	adj[j][i] = true;
      }
    }
  }
  recalculate();
}

void graph::randomize_connected(count_type n, RNG& rnd)
{
  assert(_size > 1);
  for (index_type s = 0; s < n; s++) {
    index_type i = rnd.random_in_range(_size);
    index_type j = rnd.random_in_range(_size, i);
    if (adj[i][j]) {
      adj[i][j] = false;
      if (!connected()) adj[i][j] = true;
    }
    else {
      adj[i][j] = true;
    }
  }
  recalculate();
}

void graph::randomize_undirected_connected(count_type n, RNG& rnd)
{
  assert(_size > 1);
  for (index_type s = 0; s < n; s++) {
    index_type i = rnd.random_in_range(_size);
    index_type j = rnd.random_in_range(_size, i);
    if (adj[i][j]) {
      adj[i][j] = false;
      adj[j][i] = false;
      if (!connected()) {
	adj[i][j] = true;
	adj[j][i] = true;
      }
    }
    else {
      adj[i][j] = true;
      adj[j][i] = true;
    }
  }
  recalculate();
}

void graph::randomize_strongly_connected(count_type n, RNG& rnd)
{
  assert(_size > 1);
  strongly_connected_components();
  index_type c = n_components();
  for (index_type s = 0; s < n; s++) {
    index_type i = rnd.random_in_range(_size);
    index_type j = rnd.random_in_range(_size, i);
    if (adj[i][j]) {
      adj[i][j] = false;
      strongly_connected_components();
      if (n_components() > c)
	adj[i][j] = true;
    }
    else {
      adj[i][j] = true;
    }
  }
  recalculate();
}

void graph::random_digraph(count_type n, RNG& rnd)
{
  clear_edges();
  randomize(n, rnd);
}

void graph::random_connected_digraph(count_type n, RNG& rnd)
{
  assert(_size > 1);
  clear_edges();
  for (index_type i = 0; i < _size - 1; i++) {
    adj[i][i+1] = true;
  }
  randomize_connected(n, rnd);
}

void graph::random_strongly_connected_digraph(count_type n, RNG& rnd)
{
  assert(_size > 1);
  clear_edges();
  for (index_type i = 0; i < _size - 1; i++) {
    adj[i][i+1] = true;
  }
  adj[_size - 1][0] = true;
  randomize_strongly_connected(n, rnd);
}

void graph::random_digraph_with_density(rational density, RNG& rnd)
{
  clear_edges();
  // assert(_size > 1);
  index_set e;
  index_type n = (_size * (_size - 1));
  index_type m = rational::floor_to(density*n, 1).numerator();
  rnd.select_fixed_set(e, m, n);
  for (index_type k = 0; k < e.length(); k++) {
    index_type s = e[k] / (_size - 1);
    index_type t = e[k] % (_size - 1);
    if (t >= s) t += 1;
    add_edge(s, t);
  }
}

void graph::random_tree(RNG& rnd)
{
  clear_edges();
  for (index_type k = 1; k < size(); k++) {
    index_type i = rnd.random_in_range(k);
    add_edge(i, k);
  }
}

void graph::random_tree(index_type b, index_type d, RNG& rnd)
{
  clear_edges();
  for (index_type k = 1; k < size(); k++) {
    // std::cerr << "k = " << k << ", g = " << *this << std::endl;
    index_set c;
    for (index_type i = 0; i < k; i++) {
      // std::cerr << i << ": "
      // << out_degree(i) << ", " << distance(0, i)
      // << std::endl;
      if ((out_degree(i) < b) && (distance(0, i) < d))
	c.insert(i);
    }
    // std::cerr << "c = " << c << std::endl;
    assert(!c.empty());
    index_type j = rnd.select_one_of(c);
    add_edge(j, k);
  }
}

void graph::max_clique
(index_set& sel, index_type next, index_set& clique) const
{
  if (next >= _size) {
    if (sel.length() > clique.length()) clique = sel;
  }
  else if (sel.contains(next)) {
    max_clique(sel, next + 1, clique);
  }
  else {
    if (bi[next].contains(sel)) {
      sel.insert(next);
      max_clique(sel, next + 1, clique);
      sel.subtract(next);
    }
    max_clique(sel, next + 1, clique);
  }
}

void graph::maximal_clique(index_set& clique) const
{
  assert(_size > 0);
  clique.clear();
  index_set sel;
  max_clique(sel, 0, clique);
}

void graph::maximal_clique_including
(index_type node, index_set& clique) const
{
  assert((0 <= node) && (node < _size) && (_size > 0));
  clique.clear();
  index_set sel;
  sel.assign_singleton(node);
  max_clique(sel, 0, clique);
}

void graph::maximal_clique_cover(index_set_vec& sets) const
{
  sets.clear();

  index_set clique;
  // max_clique(clique);
  // assert(!clique.empty());
  // sets.append(clique);

  index_set uncovered;
  uncovered.fill(_size);
  // uncovered.subtract(clique);
  while (!uncovered.empty()) {
    maximal_clique_including(uncovered[0], clique);
    assert(!clique.empty());
    sets.append(clique);
    uncovered.subtract(clique);
  }
}

void graph::all_max_cliques
(index_set& sel, index_type next, index_set_vec& cliques) const
{
  if (next >= _size) {
    if (sel.length() > cliques[0].length()) {
      cliques[0] = sel;
      cliques.set_length(1);
    }
    else if (sel.length() == cliques[0].length()) {
      cliques.append(sel);
    }
  }
  else if (sel.contains(next)) {
    all_max_cliques(sel, next + 1, cliques);
  }
  else {
    if (bi[next].contains(sel)) {
      sel.insert(next);
      all_max_cliques(sel, next + 1, cliques);
      sel.subtract(next);
    }
    all_max_cliques(sel, next + 1, cliques);
  }
}

void graph::all_maximal_cliques(index_set_vec& cliques) const
{
  assert(_size > 0);
  cliques.assign_value(EMPTYSET, 1);
  index_set sel;
  all_max_cliques(sel, 0, cliques);
}

void graph::all_maximal_cliques_including
(index_type node, index_set_vec& cliques) const
{
  assert(_size > 0);
  cliques.assign_value(EMPTYSET, 1);
  index_set sel;
  sel.assign_singleton(node);
  all_max_cliques(sel, 0, cliques);
}

void graph::scc_first_dfs
(index_type n, bool_vec& visited, index_vec& num) const
{
  visited[n] = true;
  for (index_type k = 0; k < _size; k++)
    if (adj[n][k] && !visited[k])
      scc_first_dfs(k, visited, num);
  num.append(n);
}

void graph::scc_second_dfs
(index_type n, bool_vec& visited, index_type c_id)
{
  visited[n] = true;
  comp[n] = c_id;
  for (index_type k = 0; k < _size; k++)
    if (adj[k][n] && !visited[k])
      scc_second_dfs(k, visited, c_id);
}

void graph::undirected_dfs
(index_type n, bool_vec& visited) const
{
  visited[n] = true;
  for (index_type k = 0; k < _size; k++)
    if ((adj[n][k] || adj[k][n]) && !visited[k])
      undirected_dfs(k, visited);
}

void graph::reachable(index_type n, bool_vec& visited) const
{
  index_vec q(n, _size);
  bool_vec inq(false, _size);
  index_type head = 0;
  index_type tail = 1;
  inq[n] = true;
  while (head < tail) {
    if (!visited[q[head]]) {
      visited[q[head]] = true;
      for (index_type k = 0; k < out[q[head]].size(); k++)
	if (!inq[out[q[head]][k]]) {
	  q[tail++] = out[q[head]][k];
	  inq[out[q[head]][k]] = true;
	}
    }
    head += 1; 
 }
  // visited[n] = true;
  // for (index_type k = 0; k < _size; k++)
  //   if (adj[n][k] && !visited[k])
  //     reachable(k, visited);
}

void graph::reverse_reachable(index_type n, bool_vec& visited) const
{
  visited[n] = true;
  for (index_type k = 0; k < _size; k++)
    if (adj[k][n] && !visited[k])
      reverse_reachable(k, visited);
}

void graph::strongly_connected_components()
{
  index_vec num(no_such_index, 0);
  bool_vec visited(false, _size);
  visited.set_length(_size);

  for (index_type k = 0; k < _size; k++) {
    if (!visited[k]) {
      scc_first_dfs(k, visited, num);
    }
  }
  assert(num.length() == _size);

  visited.assign_value(false);
  n_comp = 0;
  for (index_type k = num.length(); k > 0; k--) if (!visited[num[k - 1]]) {
    scc_second_dfs(num[k - 1], visited, n_comp);
    n_comp += 1;
  }
}

void graph::connected_components(equivalence& eq)
{
  eq.reset(_size);
  for (index_type i = 0; i < _size; i++)
    for (index_type j = i + 1; j < _size; j++)
      if (adj[i][j] || adj[j][i])
	eq.merge(i, j);
  eq.make_map(comp);
  n_comp = eq.n_classes();
}

void graph::connected_components()
{
  equivalence eq;
  connected_components(eq);
}

void graph::ramsey(const index_set& nodes, index_set& I, index_set& C) const
{
  I.clear();
  C.clear();

  if (nodes.empty()) return;

  index_type v = nodes[0];
  index_set n_v(nodes);
  n_v.intersect(bi[v]);
  n_v.subtract(v);
  ramsey(n_v, I, C);
  C.insert(v);

  index_set I2;
  index_set C2;
  n_v.assign_copy(nodes);
  n_v.subtract(out[v]);
  n_v.subtract(in[v]);
  n_v.subtract(v);
  ramsey(n_v, I2, C2);
  I2.insert(v);

  if (C2.length() > C.length()) C.assign_copy(C2);
  if (I2.length() > I.length()) I.assign_copy(I2);
}

void graph::apx_independent_set(const index_set& nodes, index_set& set) const
{
  set.clear();

  index_set n(nodes);
  index_set nextI;
  index_set nextC;

  while (!n.empty()) {
    ramsey(n, nextI, nextC);
    if (nextI.length() > set.length()) set.assign_copy(nextI);
    n.subtract(nextC);
  }
}

void graph::apx_independent_set(index_set& set) const
{
  index_set nodes;
  nodes.fill(_size);
  apx_independent_set(nodes, set);
}

void graph::apx_independent_set_including
(index_type node, index_set& set) const
{
  assert((0 <= node) && (node < _size) && (_size > 0));
  set.clear();

  index_set nodes;
  nodes.fill(_size);
  nodes.subtract(in[node]);
  nodes.subtract(out[node]);
  nodes.subtract(node);
  index_set nextI;
  index_set nextC;

  while (!nodes.empty()) {
    ramsey(nodes, nextI, nextC);
    if (nextI.length() > set.length()) set.assign_copy(nextI);
    nodes.subtract(nextC);
  }

  set.insert(node);
}

void graph::apx_independent_set_cover(index_set_vec& sets) const
{
  sets.clear();

  index_set I;
  apx_independent_set(I);
  assert(!I.empty());
  sets.append(I);

  index_set uncovered;
  uncovered.fill(_size);
  uncovered.subtract(I);
  while (!uncovered.empty()) {
    apx_independent_set_including(uncovered[0], I);
    assert(!I.empty());
    sets.append(I);
    uncovered.subtract(I);
  }
}

void graph::apx_independent_set_disjoint_cover(index_set_vec& sets) const
{
  sets.clear();
  index_set uncovered;
  uncovered.fill(_size);
  index_set I;
  while (!uncovered.empty()) {
    apx_independent_set(uncovered, I);
    assert(!I.empty());
    sets.append(I);
    uncovered.subtract(I);
  }
}

void graph::all_nondominated_cliques(index_set_vec &cliques) const
{
  /*
    Find all maximal (wrt set inclusion) cliques of the graph.
    The cliques are reported in no particular order.

    Implements the algorithm by Tomita, Tanaka and Takahashi
    ("The Worst-Case Time Complexity for Generating All Maximal
    Cliques"), which is optimal in the sense that it runs in
    O(3^(n/3)), which is about O(1.44^n), for an n-vertex graph, and
    there are graphs which have that many maximal cliques.

    Give or take a (low-order?) polynomial factor for inefficiencies
    in basic operations such as computing the set of neighbours
    from a certain set.
  */

  assert(cliques.empty());
  if (_size > 0) {
    for (index_type i = 0; i < _size; i++)
      assert(!adjacent(i, i)); // Self-loops confuse the algorithm.
    index_set current_clique;
    index_set candidates;
    candidates.fill(_size);
    all_nondominated_cliques_aux(cliques, current_clique, candidates, 1, false);
  }
}

void graph::all_cliques_geq
(index_type k, index_set_vec& cliques) const
{
  cliques.clear();
  if (_size > 0) {
    for (index_type i = 0; i < _size; i++)
      assert(!adjacent(i, i)); // Self-loops confuse the algorithm.
    index_set current_clique;
    index_set candidates;
    candidates.fill(_size);
    all_nondominated_cliques_aux(cliques, current_clique, candidates, k, false);
  }
}

void graph::one_maximal_clique(index_vec& clique) const
{
  index_set_vec cliques(EMPTYSET, 1);
  if (_size > 0) {
    cliques[0].assign_singleton(0);
    for (index_type i = 0; i < _size; i++)
      assert(!adjacent(i, i)); // Self-loops confuse the algorithm.
    index_set current_clique;
    index_set candidates;
    candidates.fill(_size);
    all_nondominated_cliques_aux(cliques, current_clique, candidates, 1, true);
  }
  clique = cliques[0];
}

void graph::all_nondominated_cliques_aux
(index_set_vec& cliques,
 index_set& current_clique,
 const index_set& candidates,
 index_type min,
 bool one_maximal_only) const
{
  if (candidates.empty()) {
    if (one_maximal_only) {
      if (current_clique.length() > cliques[0].length())
	cliques[0] = current_clique;
    }
    else {
      cliques.append_if_new(current_clique);
      // cliques.append(current_clique);
    }
  }
  else {
    // Find vertex with maximal number of successors in candidates.
    index_type best_vertex = no_such_index;
    int max_degree = 0;
    for (index_type i = 0; i < candidates.length(); i++) {
      index_type vertex = candidates[i];
      index_set neighbors = successors(vertex);
      neighbors.intersect(candidates);
      if (i == 0 || neighbors.length() > max_degree) {
        max_degree = neighbors.length();
        best_vertex = vertex;
      }
    }
    assert(best_vertex != no_such_index);
    // Iterate over candidates that are *not* adjacent to best_vertex.
    index_set chosen_set = candidates;
    chosen_set.subtract(successors(best_vertex));
    for (index_type i = 0; i < chosen_set.length(); i++) {
      index_type chosen = chosen_set[i];
      // Call recursively with chosen added to the current clique
      // and candidates reduced to chosen's neighbours.
      index_set new_candidates = candidates;
      new_candidates.intersect(successors(chosen));
      if ((current_clique.length() + new_candidates.length() + 1) >= min) {
	current_clique.insert(chosen);
	all_nondominated_cliques_aux(cliques,
				     current_clique,
				     new_candidates,
				     min, one_maximal_only);
	current_clique.subtract(chosen);
	if (one_maximal_only)
	  min = cliques[0].length();
      }
    }
  }
}

index_type graph::tree_decomposition_min_fillin(index_set_graph& td) const
{
  assert(size() > 0);
  if (size() == 1) {
    td.init(1);
    td.node_label(0).assign_singleton(0);
    return 1;
  }
  else {
    index_type n = min_fillin_node();
    mapping m(size(), n, true);
    graph g1(*this, m);
    assert(g1.size() == (size() - 1));
    for (index_type i = 0; i < bi[n].size(); i++)
      for (index_type j = i + 1; j < bi[n].size(); j++)
	if (!g1.bi_adjacent(m(bi[n][i]), m(bi[n][j])))
	  g1.add_undirected_edge(m(bi[n][i]), m(bi[n][j]));
    index_type w = g1.tree_decomposition_min_fillin(td);
    m.invert();
    td.remap_node_labels(m);
    index_type i = td.find_node_label_contains(bi[n]);
    assert(i != no_such_index);
    // short-cut: if the label of the existing node equals the
    // neighbourhood of n, the label of the new node will be a
    // strict superset of the existing one; thus, we merge the
    // two nodes into one (with the greater label)
    if (td.node_label(i) == bi[n]) {
      td.node_label(i).insert(n);
    }
    // the non-short-cut case:
    else {
      index_type np = td.add_node(bi[n]);
      td.node_label(np).insert(n);
      td.add_undirected_edge(np, i);
    }
    if ((bi[n].size() + 1) > w)
      return bi[n].size() + 1;
    else
      return w;
  }
}

bool graph::connected() const
{
  if (_size == 0) return true;
  bool_vec visited(false, _size);
  undirected_dfs(0, visited);
  return (visited.count(true) == _size);
}

bool graph::strongly_connected() const
{
  if (_size == 0) return true;
  bool_vec visited(false, _size);
  index_vec num(no_such_index, 0);
  scc_first_dfs(0, visited, num);
  return (visited.count(true) == _size);
}

void graph::write_digraph
(::std::ostream& s, bool with_node_indices, const char* name) const
{
  s << "digraph \"" << name << "\"" << "{" << ::std::endl;
  if (with_node_indices) {
    s << "node [shape=circle,width=0.5,height=0.5];" << ::std::endl;
    for (index_type k = 0; k < size(); k++)
      s << "\t" << k << " [label=\"" << k << "\"];" << std::endl;
  }
  else {
    s << "node [shape=point];" << ::std::endl;
  }
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j))
	s << "\t" << i << " -> " << j << ";" << ::std::endl;
  s << "}" << ::std::endl;
}

void graph::write_undirected_graph
(::std::ostream& s, bool with_node_indices, const char* name) const
{
  s << "graph \"" << name << "\"" << "{" << ::std::endl;
  if (with_node_indices) {
    s << "node [shape=circle,width=0.5,height=0.5];" << ::std::endl;
    for (index_type k = 0; k < size(); k++)
      s << "\t" << k << " [label=\"" << k << "\"];" << std::endl;
  }
  else {
    s << "node [shape=point];" << ::std::endl;
  }
  for (index_type i = 0; i < size(); i++)
    for (index_type j = i + 1; j < size(); j++)
      if (bi_adjacent(i, j))
	s << "\t" << i << " -- " << j << ";" << ::std::endl;
  s << "}" << ::std::endl;
}

void graph::write_component_labeled_digraph
(::std::ostream& s, const char* name) const
{
  write_labeled_digraph<index_vec>(s, *this, comp, false, name, no_such_index);
}

void graph::write_graph_correspondance
(::std::ostream& s,
 const graph& g,
 const index_vec& c,
 const char* name) const
{
  assert(g.size() == size());
  assert(c.length() == g.size());

  s << "digraph \"" << name << "\"" << ::std::endl << "{" << ::std::endl;
  s << "node [width=0,height=0];" << ::std::endl;

  s << "subgraph cluster0 {" << ::std::endl;
  for (index_type k = 0; k < size(); k++) {
    s << "\tG0_" << k << " [label=\"" << k << "\"]"
      << ::std::endl;
  }
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j)) {
	if (!g.adjacent(c[i], c[j])) {
	  s << "\tG0_" << i << " -> G0_" << j << " [style=bold];" << ::std::endl;
	}
	else {
	  s << "\tG0_" << i << " -> G0_" << j << ";" << ::std::endl;
	}
      }
  s << "}" << ::std::endl;

  s << "subgraph cluster1 {" << ::std::endl;
  for (index_type k = 0; k < g.size(); k++) {
    s << "\tG1_" << k << " [label=\"" << k << "\"]"
      << ::std::endl;
  }
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (g.adjacent(c[i], c[j])) {
	if (!adjacent(i, j)) {
	  s << "\tG1_" << i << " -> G1_" << j << " [style=bold];" << ::std::endl;
	}
	else {
	  s << "\tG1_" << i << " -> G1_" << j << ";" << ::std::endl;
	}
      }
  s << "}" << ::std::endl;

  for (index_type k = 0; k < size(); k++) {
    s << "G0_" << k << " -> G1_" << c[k] << " [style=dashed,dir=none];"
      << ::std::endl;
  }

  s << "}" << ::std::endl;
}

void index_graph::reverse()
{
  for (index_type i = 0; i < size(); i++)
    for (index_type j = i+1; j < size(); j++) {
      if (adjacent(i, j) && !adjacent(j, i)) {
	index_type l = (edge_has_label(i, j) ? edge_label(i, j) : 0);
	remove_edge(i, j);
	if ((l & EDGE_DIR) == ED_FORWARD)
	  l += (ED_BACK - ED_FORWARD);
	else if ((l & EDGE_DIR) == ED_BACK)
	  l -= (ED_BACK - ED_FORWARD);
	add_edge(j, i, l);
      }
      else if (!adjacent(i, j) && adjacent(j, i)) {
	index_type l = (edge_has_label(j, i) ? edge_label(j, i) : 0);
	remove_edge(j, i);
	if ((l & EDGE_DIR) == ED_FORWARD)
	  l += (ED_BACK - ED_FORWARD);
	else if ((l & EDGE_DIR) == ED_BACK)
	  l -= (ED_BACK - ED_FORWARD);
	add_edge(i, j, l);
      }
    }
}

void index_graph::reflect()
{
  for (index_type i = 0; i < size(); i++)
    for (index_type j = i + 1; j < size(); j++)
      if (adjacent(i, j)) {
	assert(!adjacent(j, i));
	add_edge(j, i);
	index_type l = (edge_has_label(i, j) ? edge_label(i, j) : 0);
	edge_label(j, i) = l;
	if ((l & EDGE_DIR) == ED_FORWARD)
	  edge_label(j, i) += (ED_BACK - ED_FORWARD);
	else if ((l & EDGE_DIR) == ED_BACK)
	  edge_label(j, i) -= (ED_BACK - ED_FORWARD);
      }
      else if (adjacent(j, i)) {
	add_edge(i, j);
	index_type l = (edge_has_label(j, i) ? edge_label(j, i) : 0);
	edge_label(i, j) = l;
	if ((l & EDGE_DIR) == ED_FORWARD)
	  edge_label(i, j) += (ED_BACK - ED_FORWARD);
	else if ((l & EDGE_DIR) == ED_BACK)
	  edge_label(i, j) -= (ED_BACK - ED_FORWARD);
      }
}

void index_graph::write_node_style
(std::ostream& s, index_type l)
{
  s << "shape=";
  if ((l & NODE_SHAPE) == NS_ELLIPSE)
    s << "ellipse";
  else if ((l & NODE_SHAPE) == NS_BOX)
    s << "box";
  else if ((l & NODE_SHAPE) == NS_POINT)
    s << "point";
  else if ((l & NODE_SHAPE) == NS_DIAMOND)
    s << "diamond";
  else if ((l & NODE_SHAPE) == NS_HEXAGON)
    s << "hexagon";
  else if ((l & NODE_SHAPE) == NS_OCTAGON)
    s << "octagon";
  else if ((l & NODE_SHAPE) == NS_PLAINTEXT)
    s << "plaintext";
  else
    s << "circle";
  if ((l & NODE_STYLE) == NS_FILLED)
    s << ",style=filled";
  else if ((l & NODE_STYLE) == NS_BOLD)
    s << ",style=bold";
  else if ((l & NODE_STYLE) == NS_DASHED)
    s << ",style=dashed";
  else if ((l & NODE_STYLE) == NS_DOTTED)
    s << ",style=dotted";
  if ((l & NS_DOUBLE) == NS_DOUBLE)
    s << ",peripheries=2";
}

void index_graph::write_edge_style
(std::ostream& s, index_type l)
{
  if ((l & EDGE_DIR) == ED_FORWARD)
    s << "dir=forward";
  else if ((l & EDGE_DIR) == ED_BACK)
    s << "dir=back";
  else if ((l & EDGE_DIR) == ED_BOTH)
    s << "dir=both";
  else
    s << "dir=none";
  if ((l & EDGE_STYLE) == ES_BOLD)
    s << ",style=bold";
  else if ((l & EDGE_STYLE) == ES_DASHED)
    s << ",style=dashed";
  else if ((l & EDGE_STYLE) == ES_DOTTED)
    s << ",style=dotted";
}

void index_graph::write_styled_digraph
(std::ostream& s,
 bool with_node_indices,
 const char* name,
 index_type c_id) const
{
  if (c_id != no_such_index) {
    s << "subgraph cluster" << c_id << " {" << std::endl;
    s << "node [width=0.5,height=0.5];" << ::std::endl;
  }
  else if (name) {
    s << "digraph \"" << name << "\" {" << std::endl;
    s << "node [width=0.5,height=0.5];" << ::std::endl;
  }
  for (index_type k = 0; k < size(); k++) {
    s << "\t" << k + (c_id != no_such_index ? c_id : 0) << " [";
    write_node_style(s, node_has_label(k) ? node_label(k) : 0);
    if (with_node_indices)
      s << ",label=\"" << k << "\"];" << std::endl;
    else
      s << ",label=\"\"];" << std::endl;
  }
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j)) {
	s << "\t" << i + (c_id != no_such_index ? c_id : 0)
	  << " -> " << j + (c_id != no_such_index ? c_id : 0) << " [";
	write_edge_style(s, edge_has_label(i, j) ? edge_label(i, j) : 0);
	s << "];" << std::endl;
      }
  if ((c_id != no_such_index) || (name != 0)) {
    s << "}" << ::std::endl;
  }
}

void index_graph::write_matrix
(std::ostream& s) const
{
  for (index_type i = 0; i < size(); i++) {
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j))
	s << ' ' << (edge_has_label(i, j) ? edge_label(i, j) : 0) + 1;
      else
	s << ' ' << 0;
    s << ' ' << (node_has_label(i) ? node_label(i) : 0) << std::endl;
  }
}

void index_graph::write_MATLAB
(std::ostream& s, const char* n, const char* t) const
{
  s << "defgraph('";
  if (n)
    s << n;
  else
    s << "NONAME";
  s << "', '";
  if (t)
    s << t;
  else
    s << "NOTYPE";
  s << "',[";
  for (index_type i = 0; i < size(); i++) {
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j))
	s << ' ' << (edge_has_label(i, j) ? edge_label(i, j) : 0) + 1;
      else
	s << ' ' << 0;
    if (i + 1 < size()) s << ";";
  }
  s << "], [";
  for (index_type i = 0; i < size(); i++) {
    s << ' ' << (node_has_label(i) ? node_label(i) : 0);
    if (i + 1 < size()) s << ";";
  }
  s << "]);" << std::endl;
}


// weighted_graph methods

weighted_graph& weighted_graph::quotient
(weighted_graph& g, const equivalence& eq) const
{
  g.init(eq.n_classes());
  index_vec m;
  eq.make_map(m);
  for (index_type i = 0; i < size(); i++) {
    g.increment_node_weight(i, weight(i));
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j) && !eq(i, j))
	g.increment_edge_weight(m[i], m[j], weight(i, j));
  }
}

void weighted_graph::copy(const weighted_graph& g, const index_vec& map)
{
  assert(map.length() == g.size());
  index_type m = mapping::range(map, map.length());
  init(m);
  for (index_type i = 0; i < g.size(); i++)
    if (map[i] != no_such_index) {
      assert(map[i] < size());
      increment_node_weight(map[i], g.weight(i));
      for (index_type j = 0; j < g.size(); j++)
	if (map[j] != no_such_index) {
	  assert(map[j] < size());
	  if (g.adjacent(i, j))
	    increment_edge_weight(map[i], map[j], g.weight(i, j));
	}
    }
}

NTYPE weighted_graph::apx_weighted_independent_set_1(index_set& set) const
{
  set.clear();
  NTYPE w_best = NEG_INF;

  NTYPE w_max = max_node_weight();
  index_type m = ilog(size());
  for (index_type k = 1; k <= m; k++) {
    NTYPE w_lb = (w_max / (1 << k));
    assert(w_lb > 0);
    NTYPE w_ub = (w_max / (1 << (k - 1)));

    index_set nodes;
    for (index_type i = 0; i < size(); i++)
      if ((w_lb < weight(i)) && (weight(i) <= w_ub))
	nodes.insert(i);
    index_set maxI;
    apx_independent_set(nodes, maxI);
    NTYPE w = 0;
    for (index_type i = 0; i < maxI.length(); i++)
      w += weight(maxI[i]);
    if (w > w_best) {
      set.assign_copy(maxI);
      w_best = w;
    }
  }
  return w_best;
}

NTYPE weighted_graph::apx_weighted_independent_set_2(index_set& set) const
{
  // sort nodes in order of increasing weighted degree
  weighted_vec<index_type,NTYPE> sorted;
  for (index_type k = 0; k < size(); k++) {
#ifdef NTYPE_RATIONAL
    NTYPE wd = safemul(weight(bidirectional(k)), weight(k).invert()).round();
#else
    NTYPE wd = (weight(bidirectional(k)) / weight(k));
#endif
    sorted.insert_increasing(k, wd);
  }

  set.clear();
  bool_vec rem(true, size());
  NTYPE w = 0;
  for (index_type k = 0; k < size(); k++) {
    index_type i = sorted[k].value;
    if (rem[i] && (weight(i) > 0)) {
      set.insert(i);
      rem.subtract(bidirectional(i));
      w += weight(i);
    }
  }
  return w;
}

NTYPE weighted_graph::apx_weighted_independent_set(index_set& set) const
{
  NTYPE v1 = apx_weighted_independent_set_1(set);
  NTYPE v2 = apx_weighted_independent_set_2(set);
  return MAX(v1, v2);
}

void weighted_graph::add_edge(index_type src, index_type dst)
{
  graph::add_edge(src, dst);
}

void weighted_graph::add_edge(index_type src, index_type dst, NTYPE w)
{
  graph::add_edge(src, dst);
  set_weight(src, dst, w);
}

void weighted_graph::increment_node_weight(index_type n, NTYPE w)
{
  if (!IS_ZERO(w))
    set_weight(n, weight(n) + w);
}

void weighted_graph::increment_edge_weight
(index_type src, index_type dst, NTYPE w)
{
  if (!adjacent(src, dst)) graph::add_edge(src, dst);
  if (!IS_ZERO(w))
    set_weight(src, dst, weight(src, dst) + w);
}

void weighted_graph::add_undirected_edge(index_type n0, index_type n1)
{
  graph::add_undirected_edge(n0, n1);
}

void weighted_graph::add_undirected_edge
(index_type n0, index_type n1, NTYPE w)
{
  graph::add_undirected_edge(n0, n1);
  set_weight(n0, n1, w);
  set_weight(n1, n0, w);
}

NTYPE weighted_graph::weight(index_type n) const
{
  if (node_has_label(n)) {
    return node_label(n);
  }
  else {
    return 0;
  }
}

NTYPE weighted_graph::weight(const index_set& ns) const
{
  NTYPE sw = 0;
  for (index_type k = 0; k < ns.length(); k++) {
    assert(ns[k] < size());
    sw += weight(ns[k]);
  }
  return sw;
}

NTYPE weighted_graph::bi_fringe_weight(const index_set& ns) const
{
  NTYPE sw = 0;
  for (index_type k = 0; k < size(); k++)
    if (!ns.contains(k) && bidirectional(k).have_common_element(ns))
      sw += weight(k);
  return sw;
}

NTYPE weighted_graph::weight(const bool_vec& ns) const
{
  NTYPE sw = 0;
  for (index_type k = 0; k < size(); k++)
    if (ns[k])
      sw += weight(k);
  return sw;
}

NTYPE weighted_graph::weight(index_type n0, index_type n1) const
{
  if (edge_has_label(n0, n1)) {
    return edge_label(n0, n1);
  }
  else {
    return 0;
  }
}

NTYPE weighted_graph::max_node_weight() const
{
  NTYPE w_max = NEG_INF;
  for (index_type k = 0; k < size(); k++)
    if (weight(k) > w_max) w_max = weight(k);
  return w_max;
}

void weighted_graph::set_weight(index_type n, NTYPE w)
{
  node_label(n) = w;
}

void weighted_graph::set_weight(index_type n0, index_type n1, NTYPE w)
{
  edge_label(n0, n1) = w;
}

void weighted_graph::set_node_weight(NTYPE w)
{
  for (index_type k = 0; k < size(); k++)
    node_label(k) = w;
}

void weighted_graph::set_edge_weight(NTYPE w)
{
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j))
	edge_label(i, j) = w;
}

void weighted_graph::transitive_closure()
{
  for (index_type k = 0; k < size(); k++)
    for (index_type i = 0; i < size(); i++)
      for (index_type j = 0; j < size(); j++)
	if (adjacent(i, k) && adjacent(k, j)) {
	  if (!adjacent(i, j)) {
	    add_edge(i, j, weight(i, k) + weight(k, j));
	  }
	  else if ((weight(i, k) + weight(k, j)) < weight(i, j)) {
	    set_weight(i, j, weight(i, k) + weight(k, j));
	  }
	}
  recalculate();
}

void weighted_graph::edge_matrix(cost_matrix& mat) const
{
  mat.assign_value(POS_INF, size(), size());
  for (index_type i = 0; i < size(); i++)
    mat[i][i] = 0;
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j))
	mat[i][j] = weight(i, j);
}

NTYPE weighted_graph::critical_path(cost_vec& s)
{
  cost_vec e(POS_INF, size());
  bool done = false;
  while (!done) {
    done = true;
    for (index_type k = 0; k < size(); k++) {
      NTYPE p = 0;
      for (index_type i = 0; i < predecessors(k).length(); i++)
	p = MAX(p, e[predecessors(k)[i]] + weight(predecessors(k)[i], k));
      if ((p + weight(k)) < e[k]) {
	e[k] = p + weight(k);
	done = false;
      }
    }
  }
  s.set_length(size());
  for (index_type k = 0; k < size(); k++)
    s[k] = (e[k] - weight(k));
  return cost_vec_util::max(e);
}


void weighted_graph::minimum_spanning_tree(graph& mst)
{
  // construct a list of edges, sorted in order of increasing weight
  weighted_vec<index_pair, NTYPE> e;
  for (index_type i = 0; i < size(); i++)
    for (index_type j = i + 1; j < size(); j++)
      if (bi_adjacent(i, j))
	e.insert_increasing(index_pair(i, j), weight(i, j));
  // mst is initialised to an empty graph (of same size)
  mst.init(size());
  // cc is used to track connected components in mst
  equivalence cc(size());
  for (index_type k = 0; k < e.size(); k++)
    if (!cc(e[k].value.first, e[k].value.second)) {
      mst.add_undirected_edge(e[k].value.first, e[k].value.second);
      cc.merge(e[k].value.first, e[k].value.second);
    }
}

void weighted_graph::minimum_spanning_tree(weighted_graph& mst)
{
  minimum_spanning_tree(mst);
  for (index_type i = 0; i < size(); i++)
    for (index_type j = i + 1; j < size(); j++)
      if (mst.adjacent(i, j))
	mst.set_weight(i, j, weight(i, j));
}

void weighted_graph::recursive_tree_decomposition
(NTYPE alpha, NTYPE delta, index_set_graph& t) const
{
  if (size() <= 2) {
    t.init(1);
    t.node_label(0).fill(size());
    return;
  }
  index_set s;
  index_set_vec p;
  min_vs(s, p, alpha, delta);
  t.init(1);
  t.node_label(0) = s;
  index_set todo;
  for (index_type k = 0; k < p.size(); k++) {
    index_type n = t.add_node(p[k]);
    t.add_undirected_edge(n, 0);
    todo.insert(n);
  }
  // std::cerr << "to-do = " << todo << std::endl;
  index_type ln = t.max_cardinality_undirected_leaf(todo);
  while (ln != no_such_index) {
    // std::cerr << "ln = " << ln
    //	      << ", |label(ln)| = " << t.cardinality(ln)
    //	      << std::endl;
    // ln is the unprocessed leaf with max card labelling set
    assert(t.bidirectional(ln).size() == 1);
    index_type pn = t.bidirectional(ln)[0]; // pn is the parent of ln
    // construct a mapping from graph nodes as follows:
    //  m[x] = 0  if x is in label(ln) and has neighbour in label(pn)
    //  m[x] = 1..u  for x in label(ln) with no neighbour in label(pn)
    //  m[x] = no_such_index  for all x not in label(ln)
    mapping m(size(), no_such_index, false);
    index_type u = 1;
    for (index_type k = 0; k < t.node_label(ln).size(); k++)
      if (bidirectional(t.node_label(ln)[k]).
	  have_common_element(t.node_label(pn)))
	m[t.node_label(ln)[k]] = 0;
      else
	m[t.node_label(ln)[k]] = u++;
    // if u <= 2, the graph constructed below will have too few nodes (<= 2)
    if (u > 2) {
      // construct new graph from this + mapping
      weighted_graph gl(*this, m);
      gl.remove_loops();
      // std::cerr << "gl = ";
      // gl.write_undirected_edge_set(std::cerr);
      // std::cerr << std::endl;
      // search for separator in this graph
      gl.min_vs(s, p, alpha, delta);
      // std::cerr << "|s| = " << s.size() << ", #p = " << p.size() << std::endl;
      // if successfull...
      if (p.size() >= 2) {
	// find where the "interface" (node 0) ended up
	if (s.contains(0)) {
	  assert(p.first_contains(0) == no_such_index);
	  // interface is the separator: set node_label(l) = s and make
	  // each partition a child of i
	  m.inverse(s, t.node_label(ln));
	  for (index_type k = 0; k < p.size(); k++) {
	    index_type nn = t.add_node(EMPTYSET);
	    m.inverse(p[k], t.node_label(nn));
	    t.add_undirected_edge(nn, ln);
	    todo.insert(nn);
	  }
	}
	else {
	  index_type i = p.first_contains(0);
	  assert(i != no_such_index);
	  assert(p.next_contains(0, i) == no_such_index);
	  // interface is in p[i]: set node_label(ln) = p[i], add a child
	  // ns of ln with label = s, then make each partition j != i a
	  // child of ns (note that ns is not put in to-do set, because
	  // it's not a leaf node)
	  m.inverse(p[i], t.node_label(ln));
	  index_type ns = t.add_node(EMPTYSET);
	  m.inverse(s, t.node_label(ns));
	  t.add_undirected_edge(ns, ln);
	  for (index_type k = 0; k < p.size(); k++)
	    if (k != i) {
	      index_type nn = t.add_node(EMPTYSET);
	      m.inverse(p[k], t.node_label(nn));
	      t.add_undirected_edge(nn, ns);
	      todo.insert(nn);
	    }
	}
      }
    }
    todo.subtract(ln);
    // std::cerr << "to-do = " << todo << std::endl;
    ln = t.max_cardinality_undirected_leaf(todo);
  }
}

NTYPE weighted_graph::min_vs
(index_set& s, index_set_vec& p, NTYPE alpha, NTYPE delta) const
{
  index_type i = 0;
  bool ok = false;
  while (!ok && (i < size())) {
    ok = vs_create_initial(s, i++);
  }
  if (!ok) {
    s.fill(size());
    p.clear();
    return weight(s);
  }
  // std::cerr << "initial s = " << s << std::endl;
  vs_compute_partitions(s, p);
  // std::cerr << "initial p = " << p << std::endl;
  NTYPE v = vs_improve(s, p, alpha, delta);

  ok = false;
  index_set s1;
  index_type j = size() - 1;
  while (!ok && (j > i)) {
    ok = vs_create_initial(s1, j--);
  }
  if (ok) {
    index_set_vec p1;
    vs_compute_partitions(s1, p1);
    NTYPE v1 = vs_improve(s1, p1, alpha, delta);
    if (v1 < v) {
      s = s1;
      p = p1;
      v = v1;
    }
  }

  return v;
}

bool weighted_graph::vs_create_initial(index_set& s, index_type i) const
{
  assert(size() > 2);
  assert(i < size());
  index_vec d;
  distance(i, d);
  index_type j = d.arg_max();
  // std::cerr << "i = " << i << ", d = " << d << ", j = " << j
  //	    << ", d[j] = " << d[j] << std::endl;
  if (d[j] < 2) return false;
  assert(!adjacent(i, j));
  bool_vec pi(false, size());
  pi[i] = true;
  bool_vec pj(false, size());
  pj[j] = true;
  bool done = false;
  while (!done) {
    done = true;
    bool_vec new_pi(false, size());
    bool_vec new_pj(false, size());
    for (index_type k = 0; k < size(); k++)
      if (!pi[k] && !pj[k]) {
	if (bidirectional(k).have_common_element(pi) &&
	    !bidirectional(k).have_common_element(pj) &&
	    !bidirectional(k).have_common_element(new_pj)) {
	  new_pi[k] = true;
	  done = false;
	}
	else if (bidirectional(k).have_common_element(pj) &&
		 !bidirectional(k).have_common_element(pi) &&
		 !bidirectional(k).have_common_element(new_pi)) {
	  new_pj[k] = true;
	  done = false;
	}
      }
    pi.insert(new_pi);
    pj.insert(new_pj);
  }
  s.clear();
  for (index_type k = 0; k < size(); k++)
    if (!pi[k] && !pj[k])
      s.insert(k);
  return true;
}

void weighted_graph::vs_compute_partitions
(const index_set& s, index_set_vec& p) const
{
  equivalence ceq(size());
  for (index_type i = 0; i < size(); i++)
    for (index_type j = i + 1; j < size(); j++)
      if (bi_adjacent(i, j) && (s.contains(i) == s.contains(j)))
	ceq.merge(i, j);
  // std::cerr << "ceq = " << ceq << std::endl;
  ceq.classes(p);
  bool_vec rmp(false, p.size());
  for (index_type k = 0; k < p.size(); k++)
    if (p[k].have_common_element(s)) {
      assert(s.contains(p[k]));
      rmp[k] = true;
    }
  p.remove(rmp);
}

NTYPE weighted_graph::vs_value
(const index_set& s, const index_set_vec& p, NTYPE alpha, NTYPE delta) const
{
  NTYPE wmin = POS_INF;
  NTYPE wmax = NEG_INF;
  for (index_type k = 0; k < p.size(); k++) {
    NTYPE w = weight(p[k]);
    wmin = MIN(wmin, w);
    wmax = MAX(wmax, w);
  }
  return (weight(s) +
	  (weight(s) * alpha * wmax/wmin) +
	  (delta * bi_fringe_weight(s)));
}

// NTYPE weighted_graph::vs_eval_move
// (const index_set& s, const index_set_vec& p,
//  const index_set& z, index_type d, NTYPE alpha, NTYPE delta) const
// {
//   // move z, which is a subset of s, to p[d]
//   assert(s.contains(z));
//   assert(d < p.size());
//   // set of nodes that will enter s: neighbours of z that are not
//   // already in s or in p[d]
//   index_set z_fringe;
//   bi_fringe(z, z_fringe);
//   index_set move_into_s(z_fringe);
//   move_into_s.subtract(s);
//   move_into_s.subtract(p[d]);
//   // compute wmin/wmax after move
//   NTYPE wmin = POS_INF;
//   NTYPE wmax = NEG_INF;
//   for (index_type k = 0; k < p.size(); k++) {
//     NTYPE w = weight(p[k]);
//     if (k == d) {
//       // p[d] will have weight(z) added
//       w += weight(z);
//     }
//     else {
//       // p[k] for k != d will have (move_into_s intersect p[k]) subtracted
//       index_set move_out_pk(move_into_s);
//       move_out_pk.intersect(p[k]);
//       w -= weight(move_out_pk);
//     }
//     wmin = MIN(wmin, w);
//     wmax = MAX(wmax, w);
//   }
//   index_set s_fringe;
//   bi_fringe(move_into_s, s_fringe);
//   s_fringe.subtract(s);
//   z_fringe.subtract(p[d]);
//   // weight(s) after move = (weight(s) - weight(z) + weight(move_into_s)):
//   NTYPE new_weight_s = weight(s) - weight(z) + weight(move_into_s);
//   // fringe(s) after move: -(fringe(z) - p[d]) + (fringe(move_into_s) - s)
//   return (new_weight_s +
// 	  (new_weight_s * alpha * wmax/wmin) +
// 	  (delta * (bi_fringe_weight(s) +
// 		    weight(s_fringe) -
// 		    weight(z_fringe))));
// }

NTYPE weighted_graph::vs_eval_move
(const index_set& s, const index_set_vec& p,
 const index_set& z, index_type d, NTYPE alpha, NTYPE delta) const
{
  // move z, which is a subset of s, to p[d]
  assert(s.contains(z));
  assert(d < p.size());
  index_set new_s(s);
  index_set_vec new_p(p);

  index_set into_s;
  bi_fringe(z, into_s);
  into_s.subtract(s);
  into_s.subtract(p[d]);
  new_s.subtract(z);
  new_s.insert(into_s);
  new_p[d].insert(z);
  bool bad_move = false;
  for (index_type k = 0; k < p.size(); k++)
    if (k != d) {
      new_p[k].subtract(into_s);
      if (new_p[k].empty())
	bad_move = true;
    }

  if (bad_move)
    return POS_INF;
  else
    return vs_value(new_s, new_p, alpha, delta);
}

NTYPE weighted_graph::vs_best_move
(const index_set& s, const index_set_vec& p,
 index_type first_i, index_set& z, index_type& d,
 NTYPE alpha, NTYPE delta) const
{
  assert(first_i < s.size());
  z.assign_singleton(s[first_i]);
  d = 0;
  NTYPE v = vs_eval_move(s, p, z, d, alpha, delta);
  for (index_type i = 1; i < p.size(); i++) {
    NTYPE new_v = vs_eval_move(s, p, z, i, alpha, delta);
    if (new_v < v) {
      v = new_v;
      d = i;
    }
  }
  index_type first_skip = no_such_index;
  for (index_type i = first_i + 1; i < s.size(); i++) {
    z.insert(s[i]);
    index_type min_new_d = 0;
    NTYPE min_new_v = vs_eval_move(s, p, z, min_new_d, alpha, delta);
    for (index_type j = 0; j < p.size(); j++) {
      NTYPE new_v = vs_eval_move(s, p, z, j, alpha, delta);
      if (new_v < min_new_v) {
	min_new_v = new_v;
	min_new_d = j;
      }
    }
    if (min_new_v < v) {
      v = min_new_v;
      d = min_new_d;
    }
    else {
      z.subtract(s[i]);
      if (first_skip == no_such_index)
	first_skip = i;
    }
  }
  if (first_skip != no_such_index) {
    index_set alt_z;
    index_type alt_d;
    NTYPE alt_v = vs_best_move(s, p, first_skip, alt_z, alt_d, alpha, delta);
    if (alt_v < v) {
      z.assign_copy(alt_z);
      d = alt_d;
      v = alt_v;
    }
  }
  return v;
}

NTYPE weighted_graph::vs_improve
(index_set& s, index_set_vec& p, NTYPE alpha, NTYPE delta) const
{
  assert(p.size() > 1);
  NTYPE best_val = POS_INF;
  NTYPE new_val = vs_value(s, p, alpha, delta);
  while (new_val < best_val) {
    best_val = new_val;
    // std::cerr << "current value = " << best_val << std::endl;
    index_set z;
    index_type d;
    NTYPE move_val = vs_best_move(s, p, 0, z, d, alpha, delta);
    if (move_val < best_val) {
      // do the move...
      index_set into_s;
      assert(d < p.size());
      bi_fringe(z, into_s);
      into_s.subtract(s);
      into_s.subtract(p[d]);
      s.subtract(z);
      s.insert(into_s);
      p[d].insert(z);
      for (index_type k = 0; k < p.size(); k++)
	if (k != d)
	  p[k].subtract(into_s);
      new_val = vs_value(s, p, alpha, delta);
      assert(new_val == move_val);
      assert(new_val < best_val);
    }
  }
  return best_val;
}

NTYPE weighted_graph::max_flow(index_type s, index_type t)
{
  cost_matrix f(0, size(), size());
  return max_flow(s, t, f);
}

NTYPE weighted_graph::max_flow(index_type s, index_type t, weighted_graph& rg)
{
  cost_matrix f(0, size(), size());
  NTYPE m = max_flow(s, t, f);
  rg.init(size());
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j) && ((weight(i, j) - f[i][j]) > 0))
	rg.add_edge(i, j, weight(i, j) - f[i][j]);
  return m;
}

NTYPE weighted_graph::max_flow(index_type s, index_type t, cost_matrix& f)
{
  // std::cerr << "max_flow(" << s << ", " << t << "):" << std::endl;
  f.assign_value(0, size(), size());
  pair_vec p;
  NTYPE m = 0;
  NTYPE a = augmenting_path(s, t, f, p);
  while (a > 0) {
    // std::cerr << "a = " << a << ", p = " << p << std::endl;
    for (index_type k = 0; k < p.length(); k++) {
      f[p[k].first][p[k].second] = 
	(f[p[k].first][p[k].second] + a);
      f[p[k].second][p[k].first] = 
	(f[p[k].second][p[k].first] - a);
    }
    m += a;
    a = augmenting_path(s, t, f, p);
    // assert(m < 3);
  }
  return m;
}

NTYPE weighted_graph::min_cut(index_type s, index_type t, bool_vec& s_set)
{
  weighted_graph rg;
  NTYPE c = max_flow(s, t, rg);
  rg.descendants(s, s_set);
}

NTYPE weighted_graph::min_cut(index_type s, index_type t, pair_set& e_set)
{
  weighted_graph rg;
  NTYPE c = max_flow(s, t, rg);
  bool_vec a_set;
  descendants(s, a_set);
  bool_vec s_set;
  rg.descendants(s, s_set);
  for (index_type i = 0; i < size(); i++) if (a_set[i])
    for (index_type j = 0; j < size(); j++) if (a_set[j])
      if (adjacent(i, j) && s_set[i] && !s_set[j])
	e_set.insert(index_pair(i, j));
  return c;
}

NTYPE weighted_graph::augmenting_path
(index_type s, index_type t, const cost_matrix& f, pair_vec& p)
{
  // std::cerr << "augmenting path " << s << "->" << t << "..." << std::endl;
  cost_vec m(NEG_INF, size());
  m[t] = POS_INF;
  index_vec d(no_such_index, size());
  bool done = false;
  while (!done) {
    // std::cerr << "m = " << m << std::endl;
    done = true;
    for (index_type i = 0; i < size(); i++) if (m[i] > 0) {
      for (index_type j = 0; j < predecessors(i).length(); j++) {
	NTYPE fji =
	  MIN(weight(predecessors(i)[j], i) - f[predecessors(i)[j]][i], m[i]);
	if (fji > m[predecessors(i)[j]]) {
	  m[predecessors(i)[j]] = fji;
	  d[predecessors(i)[j]] = i;
	  // std::cerr << "update: m[" << predecessors(i)[j] << "] = " << fji
	  //	    << ", d[" << predecessors(i)[j] << "] = " << i
	  //	    << std::endl;
	  done = false;
	}
      }
    }
  }
  p.set_length(0);
  if (m[s] > 0) {
    index_type k = s;
    while (k != t) {
      assert(m[k] >= m[s]);
      assert(d[k] != no_such_index);
      p.append(index_pair(k,d[k]));
      k = d[k];
    }
  }
  return m[s];
}

index_pair weighted_graph::max_weight_edge() const
{
  NTYPE w_max = NEG_INF;
  index_type n0 = no_such_index;
  index_type n1 = no_such_index;
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i,j) && ((n0 == no_such_index) || (weight(i, j) > w_max))) {
	n0 = i;
	n1 = j;
	w_max = weight(i, j);
      }
  return index_pair(n0, n1);
}

index_pair weighted_graph::min_weight_edge() const
{
  NTYPE w_min = POS_INF;
  index_type n0 = no_such_index;
  index_type n1 = no_such_index;
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i,j) && ((n0 == no_such_index) || (weight(i, j) < w_min))) {
	n0 = i;
	n1 = j;
	w_min = weight(i, j);
      }
  return index_pair(n0, n1);
}

void weighted_graph::min_and_max_edges
(const index_set& nodes,
 pair_set& e_min, NTYPE& w_min,
 pair_set& e_max, NTYPE& w_max) const
{
  e_min.clear();
  w_min = POS_INF;
  e_max.clear();
  w_max = NEG_INF;
  for (index_type i = 0; i < nodes.length(); i++)
    for (index_type j = 0; j < nodes.length(); j++)
      if ((i != j) && adjacent(nodes[i], nodes[j])) {
	if (weight(nodes[i], nodes[j]) < w_min) {
	  w_min = weight(nodes[i], nodes[j]);
	  e_min.clear();
	  e_min.insert(index_pair(nodes[i], nodes[j]));
	}
	else if (weight(nodes[i], nodes[j]) == w_min) {
	  e_min.insert(index_pair(nodes[i], nodes[j]));
	}
	if (weight(nodes[i], nodes[j]) > w_max) {
	  w_max = weight(nodes[i], nodes[j]);
	  e_max.clear();
	  e_max.insert(index_pair(nodes[i], nodes[j]));
	}
	else if (weight(nodes[i], nodes[j]) == w_max) {
	  e_max.insert(index_pair(nodes[i], nodes[j]));
	}
      }
}

void weighted_graph::write_node_set(::std::ostream& s) const
{
  s << '{';
  for (index_type i = 0; i < size(); i++) {
    if (i > 0) s << ',';
    s << i << "[" << weight(i) << "]";
  }
  s << '}';
}

void weighted_graph::write_edge_set(::std::ostream& s) const
{
  s << '{';
  bool first = true;
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j)) {
	if (!first) s << ',';
	first = false;
	if (edge_has_label(i, j)) {
	  s << i << "-[" << weight(i, j) << "]->" << j;
	}
	else {
	  s << i << "->" << j;
	}
      }
  s << '}';
}

void weighted_graph::write_compact(::std::ostream& s) const
{
  s << '(';
  write_node_set(s);
  s << ',';
  write_edge_set(s);
  s << '}';
}

void weighted_graph::write_matrix(::std::ostream& s) const
{
  s << '[';
  for (index_type i = 0; i < size(); i++) {
    if (i > 0) s << ' ';
    s << '[';
    for (index_type j = 0; j < size(); j++) {
      if (j > 0) s << ',';
      if (adjacent(i, j)) {
	if (edge_has_label(i, j)) {
	  s << PRINT_NTYPE(edge_label(i, j));
	}
	else {
	  s << "?";
	}
      }
      else {
	s << "INF";
      }
    }
    s << ']';
    if (i + 1 < size()) {
      s << ',' << '\n';
    }
    else {
      s << ']' << '\n';
    }
  }
}

NTYPE weighted_graph::maximal_matching(weighted_graph& matching)
{
  index_type ni = no_such_index;
  index_type nj = no_such_index;
  for (index_type i = 0; (i < size()) && (ni == no_such_index); i++)
    for (index_type j = i + 1; (j < size()) && (nj == no_such_index); j++)
      if (adjacent(i, j) && adjacent(j, i)) {
	ni = i;
	nj = j;
      }
  // if there exists an edge...
  if ((ni != no_such_index) && (nj != no_such_index)) {
    // try with edge ni->nj in (recurse on subgraph without nodes ni, nj)...
    index_set s;
    s.fill(size());
    s.subtract(ni);
    s.subtract(nj);
    weighted_graph* g1 = new weighted_graph(*this, s);
    weighted_graph* m1 = new weighted_graph();
    NTYPE v1 = weight(ni, nj) + g1->maximal_matching(*m1);
    delete g1;

    // try with edge ni->nj out...
    weighted_graph* g2 = new weighted_graph(*this);
    g2->remove_edge(ni, nj);
    NTYPE v2 = g2->maximal_matching(matching);
    delete g2;

    if (v1 > v2) {
      // construct the matching with ni->nj (if v2 >= v1, just return
      // the matching returned by second recursive call)
      matching.init(size());
      for (index_type i = 0; i < m1->size(); i++)
	for (index_type j = i + 1; j < m1->size(); j++)
	  if (m1->adjacent(i, j) && m1->adjacent(j, i))
	    matching.add_undirected_edge(s[i], s[j], weight(s[i], s[j]));
      matching.add_undirected_edge(ni, nj, weight(ni, nj));
    }
    delete m1;
    return MAX(v1, v2);
  }
  else {
    matching.init(size());
    return 0;
  }
}

NTYPE weighted_graph::apx_matching(bool_vec& matched)
{
  bool_vec rem(true, size());
  matched.assign_value(false, size());
  NTYPE val[2] = {0,0};
  index_type i = 0;
  index_type v = rem.first(true);
  while (v != no_such_index) {
    rem[v] = false;
    bool done = false;
    while (!done) {
      NTYPE w_max = NEG_INF;
      index_type v_next = no_such_index;
      for (index_type k = 0; k < bidirectional(v).length(); k++)
	if (rem[bidirectional(v)[k]] &&
	    (weight(v, bidirectional(v)[k]) > w_max)) {
	  w_max = weight(v, bidirectional(v)[k]);
	  v_next = bidirectional(v)[k];
	}
      if (v_next != no_such_index) {
	val[i] += w_max;
	i = ((i + 1) % 2);
	matched[v] = true;
	matched[v_next] = true;
	rem[v_next] = false;
	v = v_next;
      }
      else {
	done = true;
      }
    }
    v = rem.first(true);
  }
  return MAX(val[0], val[1]);
}


// index_set_graph methods

index_set_graph::index_set_graph
(const graph& g, const equivalence& eq)
  : labeled_graph<index_set,index_set>(eq.n_classes())
{
  index_vec m;
  eq.make_map(m);
  assert(m.length() == g.size());
  for (index_type i = 0; i < g.size(); i++)
    for (index_type j = 0; j < g.size(); j++)
      if (g.adjacent(i, j) && !eq(i, j))
	add_edge(m[i], m[j]);
  index_set ce;
  eq.canonical_elements(ce);
  assert(ce.length() == size());
  for (index_type i = 0; i < size(); i++)
    eq.class_elements(ce[i], node_label(i));
}

index_set_graph::index_set_graph
(const index_set_graph& g, const equivalence& eq)
{
  g.quotient(*this, eq);
}

void index_set_graph::add_edge
(index_type src, index_type dst)
{
  labeled_graph<index_set,index_set>::add_edge(src, dst);
}

void index_set_graph::add_edge
(index_type src, index_type dst, const index_set& lbl)
{
  labeled_graph<index_set,index_set>::add_edge(src, dst, lbl);
}

void index_set_graph::add_edge
(const index_set& srcs, index_type dst)
{
  labeled_graph<index_set,index_set>::add_edge(srcs, dst);
}

void index_set_graph::add_edge
(const index_set& srcs, index_type dst, const index_set& lbl)
{
  labeled_graph<index_set,index_set>::add_edge(srcs, dst, lbl);
}

void index_set_graph::add_edge
(index_type src, const index_set& dsts)
{
  labeled_graph<index_set,index_set>::add_edge(src, dsts);
}

void index_set_graph::add_edge
(index_type src, const index_set& dsts, const index_set& lbl)
{
  labeled_graph<index_set,index_set>::add_edge(src, dsts, lbl);
}

// add new label to edge if it already exists, else add edge with
// singleton set label
void index_set_graph::add_edge
(index_type src, index_type dst, index_type newlbl)
{
  graph::add_edge(src, dst);
  if (edge_has_label(src, dst))
    edge_label(src, dst).insert(newlbl);
  else {
    edge_label(src, dst).assign_singleton(newlbl);
  }
}

void index_set_graph::add_edge
(const index_set& srcs, index_type dst, index_type newlbl)
{
  for (index_type k = 0; k < srcs.length(); k++) {
    graph::add_edge(srcs[k], dst);
    if (edge_has_label(srcs[k], dst))
      edge_label(srcs[k], dst).insert(newlbl);
    else {
      edge_label(srcs[k], dst).assign_singleton(newlbl);
    }
  }
}

void index_set_graph::add_edge
(index_type src, const index_set& dsts, index_type newlbl)
{
  for (index_type k = 0; k < dsts.length(); k++) {
    graph::add_edge(src, dsts[k]);
    if (edge_has_label(src, dsts[k]))
      edge_label(src, dsts[k]).insert(newlbl);
    else {
      edge_label(src, dsts[k]).assign_singleton(newlbl);
    }
  }
}

void index_set_graph::remap_node_labels(const index_vec& map)
{
  for (index_type k = 0; k < size(); k++)
    if (node_has_label(k))
      node_label(k).remap(map);
}

void index_set_graph::remap_edge_labels(const index_vec& map)
{
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (edge_has_label(i, j))
	edge_label(i, j).remap(map);
}

index_type index_set_graph::find_node_label_contains(const index_set& s) const
{
  if (size() == 0) return no_such_index;
  if (s.empty()) return 0;
  for (index_type k = 0; k < size(); k++)
    if (node_has_label(k))
      if (node_label(k).contains(s))
	return k;
  return no_such_index;
}

index_type index_set_graph::cardinality(index_type n) const
{
  if (node_has_label(n))
    return node_label(n).size();
  else
    return 0;
}

index_type index_set_graph::max_node_cardinality() const
{
  index_type m = 0;
  for (index_type k = 0; k < size(); k++)
    if (cardinality(k) > m)
      m = cardinality(k);
  return m;
}

index_type index_set_graph::max_cardinality_undirected_leaf
(const index_set& n) const
{
  index_type i = first_undirected_leaf(n);
  if (i == no_such_index) return no_such_index; // graph has no leaf in n
  index_type i_max = i;
  index_type c_max = cardinality(i);
  i = next_undirected_leaf(i, n);
  while (i != no_such_index) {
    index_type c = cardinality(i);
    if (c > c_max) {
      i_max = i;
      c_max = c;
    }
    i = next_undirected_leaf(i, n);
  }
  return i_max;
}

void index_set_graph::contract(index_type n0, index_type n1)
{
  // merge n1's label into n0's label
  node_label(n0).insert(node_label(n1));
  // for each incoming edge to n1:
  const index_set& preds = predecessors(n1);
  for (index_type i = 0; i < preds.size(); i++) {
    // if the edge already exists to n0, merge edge labels:
    if (adjacent(preds[i], n0)) {
      if (edge_has_label(preds[i], n1))
	edge_label(preds[i], n0).insert(edge_label(preds[i], n1));
    }
    // otherwise, add new edge with n1's edge label:
    else {
      add_edge(preds[i], n0, edge_label(preds[i], n1));
    }
  }
  // for each outgoing edge from n1:
  const index_set& succs = successors(n1);
  for (index_type i = 0; i < succs.size(); i++) {
    // if the edge already exists from n0, merge edge labels:
    if (adjacent(n0, succs[i])) {
      if (edge_has_label(n1, succs[i]))
	edge_label(n1, succs[i]).insert(edge_label(n1, succs[i]));
    }
    // otherwise, add new edge with n1's edge label:
    else {
      add_edge(n0, succs[i], edge_label(n1, succs[i]));
    }
  }
  // finally, remove n1:
  remove_node(n1);
}

index_set_graph& index_set_graph::quotient
(index_set_graph& g, const equivalence& eq) const
{
  graph::quotient(g, eq);
  g.clear_node_labels();
  g.clear_edge_labels();
  index_vec m;
  eq.make_map(m);
  for (index_type i = 0; i < g.size(); i++)
    g.node_label(i) = EMPTYSET;
  for (index_type i = 0; i < size(); i++)
    g.node_label(m[i]).insert(node_label(i));
}

index_set_graph& index_set_graph::edge_label_quotient
(index_set_graph& g, const equivalence& eq) const
{
  graph::quotient(g, eq);
  g.clear_node_labels();
  g.clear_edge_labels();
  index_vec m;
  eq.make_map(m);
  // assign empty sets to all nodes/edges in g
  for (index_type i = 0; i < g.size(); i++) {
    g.node_label(i) = EMPTYSET;
    for (index_type j = 0; j < g.size(); j++)
      if (g.adjacent(i, j))
	g.edge_label(i, j) = EMPTYSET;
  }
  // construct label sets
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j)) {
	if (m[i] == m[j]) {
	  if (edge_has_label(i, j))
	    g.node_label(m[i]).insert(edge_label(i, j));
	}
	else {
	  assert(g.adjacent(i, j));
	  if (edge_has_label(i, j))
	    g.edge_label(m[i], m[j]).insert(edge_label(i, j));
	}
      }
}

index_set_graph& index_set_graph::union_reachable(index_set_graph& g) const
{
  g.copy(*this);
  for (index_type k = 0; k < size(); k++) {
    index_set n;
    g.descendants(k, n);
    for (index_type i = 0; i < n.length(); i++)
      g.node_label(k).insert(node_label(n[i]));
  }
}

index_type index_set_graph::union_of_edges_on_path
(index_type s, index_type t, index_set& u) const
{
  u.clear();
  index_vec p;
  index_type l = shortest_path(s, t, p);
  if (l == no_such_index) return l;
  // double check...
  assert(p.length() > 0);
  assert(p[0] == s);
  assert(p[p.length() - 1] == t);
  //std::cerr << "path from " << s << " to " << t << ": " << p << std::endl;
  for (index_type i = 1; i < p.length(); i++) {
    if (edge_has_label(p[i - 1], p[i]))
      u.insert(edge_label(p[i - 1], p[i]));
  }
  return l;
}

bool index_set_graph::union_of_edges_between
(index_type s, index_type t, index_set& u) const
{
  u.clear();
  bool_vec b;
  between(s, t, b);
  if (b.count(true) == 0) return false;
  // double check...
  assert(b[s] && b[t]);
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (b[i] && b[j] && adjacent(i, j)) {
	if (edge_has_label(i, j))
	  u.insert(edge_label(i, j));
      }
  return true;
}

bool index_set_graph::union_of_edges_strictly_between
(index_type s, index_type t, index_set& u) const
{
  u.clear();
  bool_vec b;
  between(s, t, b);
  if (b.count(true) == 0) return false;
  // double check...
  assert(b[s] && b[t]);
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (b[i] && b[j] && adjacent(i, j) && ((i != s) || (j != t))) {
	if (edge_has_label(i, j))
	  u.insert(edge_label(i, j));
      }
  return true;
}

void index_set_graph::edge_label_preserving_transitive_reduction()
{
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j)) {
	index_set u;
	bool conn = union_of_edges_strictly_between(i, j, u);
	assert(conn);
	edge_label(i, j).subtract(u);
	if (edge_label(i, j).empty())
	  remove_edge(i, j);
      }
}

void index_set_graph::merge_labels(const index_set& ns)
{
  if (ns.empty()) return;
  assert(ns[0] < size());
  for (index_type i = 1; i < ns.length(); i++) {
    assert(ns[i] < size());
    node_label(ns[0]).insert(node_label(ns[i]));
  }
  for (index_type i = 1; i < ns.length(); i++) {
    node_label(ns[i]).assign_copy(node_label(ns[0]));
  }
}

void index_set_graph::merge_labels_upwards()
{
  strongly_connected_components();
  equivalence eq;
  component_partitioning(eq);
  // cdag: DAG of strongly connected components of *this
  index_set_graph cdag(*((graph*)this), eq);
  // 1. union node labels among nodes in same component
  for (index_type k = 0; k < cdag.size(); k++) {
    assert(!cdag.node_label(k).empty());
    merge_labels(cdag.node_label(k));
  }
  // 2. propagate upwards...
  while (!cdag.empty()) {
    // cdag.write_digraph(std::cerr, "CDAG");
    index_type l = cdag.first_leaf();
    // std::cerr << "l = " << l << std::endl;
    assert(l != no_such_index);
    for (index_type i = 0; i < cdag.predecessors(l).length(); i++) {
      index_type n_i = cdag.predecessors(l)[i];
      assert(n_i < size());
      for (index_type j = 0; j < cdag.node_label(n_i).length(); j++) {
	index_type n_j = cdag.node_label(n_i)[j];
	assert(n_j < size());
	assert(!cdag.node_label(l).empty());
	assert(cdag.node_label(l)[0] < size());
	node_label(n_j).insert(node_label(cdag.node_label(l)[0]));
      }
    }
    cdag.remove_node(l);
  }
}

void index_set_graph::merge_labels_downwards()
{
  strongly_connected_components();
  equivalence eq;
  component_partitioning(eq);
  // cdag: DAG of strongly connected components of *this
  index_set_graph cdag(*((graph*)this), eq);
  // 1. union node labels among nodes in same component
  for (index_type k = 0; k < cdag.size(); k++) {
    assert(!cdag.node_label(k).empty());
    merge_labels(cdag.node_label(k));
  }
  // 2. propagate downwards...
  while (!cdag.empty()) {
    index_type l = cdag.first_root();
    assert(l != no_such_index);
    for (index_type i = 0; i < cdag.successors(l).length(); i++) {
      index_type n_i = cdag.successors(l)[i];
      for (index_type j = 0; j < cdag.node_label(n_i).length(); j++) {
	index_type n_j = cdag.node_label(n_i)[j];
	node_label(n_j).insert(node_label(cdag.node_label(l)[0]));
      }
    }
    cdag.remove_node(l);
  }
}

index_set_graph& index_set_graph::subgraph_set_size_gt
(index_set_graph& g, index_type l)
{
  index_set r;
  for (index_type k = 0; k < size(); k++)
    if (node_has_label(k)) {
      if (node_label(k).length() > l) r.insert(k);
    }
  subgraph(g, r);
  return g;
}

void index_set_graph::write_edge_set(::std::ostream& s) const
{
  s << '{';
  bool first = true;
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j)) {
	if (!first) s << ',';
	first = false;
	s << i << ":";
	if (node_has_label(i))
	  s << node_label(i);
	else
	  s << "[]";
	if (edge_has_label(i, j))
	  s << "-" << edge_label(i,j) << "->" << j << ":";
	else
	  s << "-[]->" << j << ":";
	if (node_has_label(j))
	  s << node_label(j);
	else
	  s << "[]";
      }
  s << '}';
}

void index_set_graph::write_undirected_edge_set(::std::ostream& s) const
{
  s << '{';
  bool first = true;
  for (index_type i = 0; i < size(); i++)
    for (index_type j = i + 1; j < size(); j++)
      if (bi_adjacent(i, j)) {
	if (!first) s << ',';
	first = false;
	s << i << ":";
	if (node_has_label(i))
	  s << node_label(i);
	else
	  s << "[]";
	s << "-" << j << ":";
	if (node_has_label(j))
	  s << node_label(j);
	else
	  s << "[]";
      }
  s << '}';
}

void index_set_graph::write_DOT
(std::ostream& s, bool undirected, const char* name) const
{
  if (undirected)
    s << "graph";
  else if (strncmp(name, "cluster", 7) == 0)
    s << "subgraph";
  else
    s << "digraph";
  s << " \"" << name << "\"" << ::std::endl << "{" << ::std::endl;
  for (index_type k = 0; k < size(); k++) {
    s << "\t" << k << " [label=\"" << k << ":{";
    if (node_has_label(k)) {
      for (index_type i = 0; i < node_label(k).length(); i++) {
	if (i > 0) s << ",";
	s << node_label(k)[i];
      }
    }
    s << "}\"];" << ::std::endl;
  }
  if (undirected) {
    for (index_type i = 0; i < size(); i++)
      for (index_type j = i + 1; j < size(); j++)
	if (bi_adjacent(i, j)) {
	  s << "\t" << i << " -- " << j;
	  if (edge_has_label(i, j))
	    s << " [label=\"" << edge_label(i, j) << "\"]";
	  s << ";" << ::std::endl;
	}
  }
  else {
    for (index_type i = 0; i < size(); i++)
      for (index_type j = 0; j < size(); j++)
	if (adjacent(i, j)) {
	  s << "\t" << i << " -> " << j;
	  if (edge_has_label(i, j))
	    s << " [label=\"" << edge_label(i, j) << "\"]";
	  s << ";" << ::std::endl;
	}
  }
  s << "}" << ::std::endl;
}

END_HSPS_NAMESPACE
