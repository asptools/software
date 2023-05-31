#ifndef HYPERGRAPH_H
#define HYPERGRAPH_H

#include "index_type.h"

BEGIN_HSPS_NAMESPACE

class hypergraph {
  index_type    _size;
  index_set_vec _edges;
  index_set_vec _in;

 public:
  hypergraph();
  hypergraph(index_type n);
  hypergraph(index_type n, const index_set_vec e);
  hypergraph(const hypergraph& h);
  ~hypergraph();

  index_type size() const;
  index_type n_edges() const;
  const index_set& edge(index_type k) const;
  const index_set_vec& edges() const;
  const index_set& adjacent_edges(index_type k) const;
  bool is_minimal_edge(index_type k) const;
  bool is_maximal_edge(index_type k) const;

  void add_edge(const index_set& e);
  void add_singleton_edge(index_type n);
  void add_binary_edge(index_type n1, index_type n2);
  void reduce_to_minimal();

  bool is_independent(const index_set& ns) const;
  bool is_independent(const bool_vec& ns) const;
  void covering_sets(index_set_vec& cs) const;
  void covering_sets_2(index_set_vec& cs) const;
  void independent_sets(index_set_vec& is) const;

  void write_edge_set(std::ostream& s) const;
};

inline ::std::ostream& operator<<(::std::ostream& s, const hypergraph& g)
{
  g.write_edge_set(s);
  return s;
};

END_HSPS_NAMESPACE

#endif
