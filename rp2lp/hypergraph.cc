
#include "hypergraph.h"

BEGIN_HSPS_NAMESPACE

hypergraph::hypergraph()
  : _size(0)
{
  // done
}

hypergraph::hypergraph(index_type n)
  : _size(n), _in(EMPTYSET, n)
{
  // done
}

hypergraph::hypergraph(index_type n, const index_set_vec e)
  : _size(n), _edges(e), _in(EMPTYSET, n)
{
  for (index_type k = 0; k < _edges.length(); k++) {
    for (index_type i = 0; k < _edges[k].length(); i++) {
      assert(_edges[k][i] < _size);
      _in[_edges[k][i]].insert(k);
    }
  }
}

hypergraph::hypergraph(const hypergraph& h)
  : _size(h._size), _edges(h._edges), _in(h._in)
{
  // done
}

hypergraph::~hypergraph()
{
  // done
}

index_type hypergraph::size() const
{
  return _size;
}

index_type hypergraph::n_edges() const
{
  return _edges.length();
}

const index_set& hypergraph::edge(index_type k) const
{
  assert(k < _edges.length());
  return _edges[k];
}

const index_set_vec& hypergraph::edges() const
{
  return _edges;
}

const index_set& hypergraph::adjacent_edges(index_type k) const
{
  assert(k < _size);
  return _in[k];
}

bool hypergraph::is_minimal_edge(index_type k) const
{
  assert(k < _edges.size());
  index_type i = _edges.first_strict_subset(_edges[k]);
  return (i == no_such_index);
}

bool hypergraph::is_maximal_edge(index_type k) const
{
  assert(k < _edges.size());
  index_type i = _edges.first_strict_superset(_edges[k]);
  return (i == no_such_index);
}

void hypergraph::add_edge(const index_set& e)
{
  _edges.append(e);
  for (index_type k = 0; k < e.length(); k++) {
    assert(e[k] < _size);
    _in[e[k]].insert(_edges.length() - 1);
  }
}

void hypergraph::add_singleton_edge(index_type n)
{
  index_set e;
  e.assign_singleton(n);
  add_edge(e);
}

void hypergraph::add_binary_edge(index_type n1, index_type n2)
{
  index_set e;
  e.insert(n1);
  e.insert(n2);
  add_edge(e);
}

void hypergraph::reduce_to_minimal()
{
  _edges.reduce_to_minimal();
  for (index_type k = 0; k < _size; k++)
    _in[k].clear();
  for (index_type k = 0; k < _edges.size(); k++)
    for (index_type i = 0; i < _edges[k].size(); i++) {
      assert(_edges[k][i] < _size);
      _in[_edges[k][i]].insert(k);
    }
}

bool hypergraph::is_independent(const index_set& ns) const
{
  bool_vec t(false, _edges.length());
  for (index_type k = 0; k < ns.length(); k++) {
    assert(ns[k] < _size);
    for (index_type i = 0; i < _in[ns[k]].length(); i++)
      if (!t[_in[ns[k]][i]]) {
	if (ns.contains(_edges[_in[ns[k]][i]]))
	  return false;
	t[_in[ns[k]][i]] = true;
      }
  }
  return true;
}

bool hypergraph::is_independent(const bool_vec& ns) const
{
  bool_vec t(false, _edges.length());
  for (index_type k = 0; (k < ns.length()) && (k < _size); k++)
    if (ns[k])
      for (index_type i = 0; i < _in[k].length(); i++)
	if (!t[_in[k][i]]) {
	  if (ns.contains(_edges[_in[k][i]]))
	    return false;
	  t[_in[k][i]] = true;
	}
  return true;
}

void hypergraph::covering_sets(index_set_vec& cs) const
{
  cs.clear();
  cs.append(EMPTYSET);
  bool_vec o(false, _edges.length());
  for (index_type k = 0; k < _size; k++)
    for (index_type i = 0; i < _in[k].length(); i++) {
      index_type e = _in[k][i];
      if (!o[e]) {
	// std::cerr << "uncovered edge: " << _edges[e] << std::endl;
	index_type m = cs.length();
	for (index_type j = 0; j < m; j++) {
	  // std::cerr << " set " << cs[j] << "...";
	  if (cs[j].first_common_element(_edges[e]) == no_such_index) {
	    // std::cerr << " does not cover this edge" << std::endl;
	    for (index_type l = 1; l < _edges[e].length(); l++) {
	      index_set x(cs[j]);
	      x.insert(_edges[e][l]);
	      cs.append(x);
	    }
	    cs[j].insert(_edges[e][0]);
	  }
	  // else {
	  //  std::cerr << " already covers this edge" << std::endl;
	  // }
	}
	cs.reduce_to_minimal();
	o[_in[k][i]] = true;
      }
    }
  // std::cerr << "edges: " << _edges << std::endl;
  // std::cerr << "minimal covering sets: " << cs << std::endl;
}

void hypergraph::covering_sets_2(index_set_vec& cs) const
{
  cs.clear();
  cs.append(EMPTYSET);
  // csc[i] = indices of sets in cs that contain element i
  index_set_vec csc(EMPTYSET, _size);
  for (index_type k = 0; k < _edges.length(); k++) {
    index_type m = cs.length();
    for (index_type j = 0; j < m; j++) {
      if (cs[j].first_common_element(_edges[k]) == no_such_index) {
	// set cs[j] must be extended to cover edge[k]
	for (index_type i = 1; i < _edges[k].length(); i++) {
	  index_type e = _edges[k][i];
	  index_set x(cs[j]);
	  x.insert(e);
	  // only add new set x to collection if it's not dominated
	  // by (i.e., contains) any set already in cs
	  bool non_dom = true;
	  for (index_type p = 0; (p < csc[e].size()) && non_dom; p++)
	    if (x.contains(cs[csc[e][p]]))
	      non_dom = false;
	  if (non_dom) {
	    cs.append(x);
	    // update csc
	    for (index_type p = 0; p < x.size(); p++)
	      csc[x[p]].insert(cs.size() - 1);
	  }
	}
	// dominance check not done for this element/set
	cs[j].insert(_edges[k][0]);
	// update csc
	csc[_edges[k][0]].insert(j);
      }
    }
  }
  // cs may now contain some non-minimal sets, so must be reduced
  cs.reduce_to_minimal();
}

void hypergraph::independent_sets(index_set_vec& is) const
{
  index_set_vec cs;
  covering_sets(cs);
  is.assign_value(EMPTYSET, cs.length());
  for (index_type k = 0; k < is.length(); k++) {
    is[k].fill(_size);
    is[k].subtract(cs[k]);
  }
  // std::cerr << "maximal independent sets: " << is << std::endl;
}

void hypergraph::write_edge_set(std::ostream& s) const
{
  s << '{';
  for (index_type k = 0; k < _edges.size(); k++) {
    if (k > 0) s << ',';
    s << _edges[k];
  }
  s << '}';
}


END_HSPS_NAMESPACE
