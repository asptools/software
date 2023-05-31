#ifndef GRAPH_H
#define GRAPH_H

#include "config.h"
#include "index_type.h"
#include "numeric_type.h"
#ifndef NTYPE_RATIONAL
#include "rational.h"
#endif
#include "rng.h"

#include <map>
#include <iostream>

BEGIN_HSPS_NAMESPACE

// forward declarations
class graph;
class index_set_graph;
class weighted_graph;

class adjacency_list_graph {
 private:
  index_type    _size;
  index_set_vec _edges;

 public:
  // empty graph
  adjacency_list_graph();
  // empty graph with s nodes
  adjacency_list_graph(index_type s);
  // copy constructor
  adjacency_list_graph(const adjacency_list_graph& g);
  // subgraph constructor
  adjacency_list_graph(const adjacency_list_graph& g, const index_set& n);
  // copy-and-map constructor
  adjacency_list_graph(const adjacency_list_graph& g, const index_vec& m);
  // quotient graph constructor
  adjacency_list_graph(const adjacency_list_graph& g, const equivalence& eq);

  // copy constructor from general graph
  adjacency_list_graph(const graph& g);
  // subgraph constructor from general graph
  adjacency_list_graph(const graph& g, const index_set& n);

  ~adjacency_list_graph() { };

  // return size (== number of nodes)
  index_type size() const { return _size; };

  // query existence of directed edge (i,j)
  bool adjacent(index_type i, index_type j) const {
    return _edges[i].contains(j);
  };

  // query existence of UN-directed edge {i,j}
  bool bi_adjacent(index_type i, index_type j) const {
    return (adjacent(i, j) && adjacent(j, i));
  };

  // return number of edges
  index_type n_edges() const;

  // return number of edges between nodes in two sets
  index_type n_edges(const index_set& from, const index_set& to) const;

  // return set/number of successor nodes to node i.
  const index_set& successors(index_type i) const { return _edges[i]; };
  index_type out_degree(index_type i) const { return _edges[i].size(); };

  void distance(index_type s0, index_vec& d) const;
  index_type shortest_path(index_type s, index_type t, index_vec& p) const;

  // assign empty graph with s nodes
  void init(index_type s);
  // copy from general graph
  void copy(const graph& g);
  // mapped copy from general graph
  void copy(const graph& g, const index_vec& map);
  // add an isolated node, return index of the new node (== new size - 1)
  index_type add_node();
  // remove a node
  void remove_node(index_type n);

  // add g as a subgraph; on return, m maps indices in g to indices in this
  void add_graph(const graph& g, mapping& m);

  // add/remove edges:
  void add_edge(index_type src, index_type dst);
  void add_edge(const index_set& srcs, index_type dst);
  void add_edge(index_type src, const index_set& dsts);

  void add_edge_to_transitive_closure(index_type src, index_type dst);

  void remove_edge(index_type src, index_type dst);
  void remove_edges_from(index_type src);
  void remove_edges_to(index_type dst);
  void remove_edges_incident_on(index_type n);

  void clear_edges();

  void write_compact(::std::ostream& s) const;

 private:
  index_type extract_path(const index_vec& d, index_type n, index_vec& p) const;
};

class graph {
 private:
  index_type    _size;
  bool_matrix   adj;
  index_set_vec in;
  index_set_vec out;
  index_set_vec bi;
  index_vec     comp;
  index_type    n_comp;

 public:
  // empty graph
  graph();
  // empty graph with s nodes
  graph(index_type s);
  // copy constructor
  graph(const graph& g);
  // subgraph constructor
  graph(const graph& g, const index_set& n);
  // copy-and-map constructor
  graph(const graph& g, const index_vec& m);
  // quotient graph constructor
  graph(const graph& g, const equivalence& eq);
  ~graph();

  // retrieve basic graph information
  index_type size() const { return _size; };
  bool adjacent(index_type i, index_type j) const { return adj[i][j]; };
  bool adjacent(index_type i, const index_set& n) const;
  bool adjacent(const index_set& n, index_type i) const;
  bool adjacent(const index_set& n0, const index_set& n1) const;
  bool bi_adjacent(index_type i, index_type j) const
    { return (adj[i][j] && adj[j][i]); };
  bool bi_adjacent(index_type i, const index_set& n) const;
  index_type n_edges() const;
  index_type n_edges(const index_set& from, const index_set& to) const;
  pair_set& edges(pair_set& s) const;
  index_type n_induced_undirected_edges() const;
  index_type n_induced_undirected_edges(const index_set& n0,
					const index_set& n1) const;
  index_type n_bidirectional_edges() const;
  index_type n_bidirectional_edges(const index_set& n0,
				   const index_set& n1) const;
  const index_set& successors(index_type i) const { return out[i]; };
  index_type out_degree(index_type i) const { return out[i].length(); };
  const index_set& predecessors(index_type i) const { return in[i]; };
  index_type in_degree(index_type i) const { return in[i].length(); };
  const index_set& bidirectional(index_type i) const { return bi[i]; };
  index_type bi_degree(index_type i) const { return bi[i].length(); };
  pair_set& bidirectional_edges(pair_set& s) const;

  void descendants(index_type n0, bool_vec& s) const;
  void descendants(const index_set& s0, bool_vec& s) const;
  void descendants(index_type n0, index_set& s) const;
  void descendants(const index_set& s0, index_set& s) const;
  void ancestors(index_type n0, bool_vec& s) const;
  void ancestors(const index_set& s0, bool_vec& s) const;
  void ancestors(index_type n0, index_set& s) const;
  void ancestors(const index_set& s0, index_set& s) const;

  void nearest_common_descendants
    (index_type n0, index_type n1, bool_vec& s) const;

  // between(n0, n1) = descendants(n0) intersect ancestors(n1)
  void between(index_type n0, index_type n1, bool_vec& s) const;

  index_type max_out_degree() const;
  index_type max_in_degree() const;
  index_type max_bi_degree() const;
  index_type min_out_degree() const;
  index_type min_in_degree() const;
  index_type min_bi_degree() const;

  index_type max_bi_degree_node() const;
  index_type min_bi_degree_node() const;
  index_type min_fillin_node() const;

  bool empty() const;
  bool connected() const;
  bool strongly_connected() const;
  bool reachable(index_type n0, index_type n1) const;
  void reachable(bool_vec& v) const;
  index_type count_reachable(index_type n0) const;
  bool acyclic() const;
  bool top_sort(index_vec& s) const;

  index_type first_root() const;
  index_type first_leaf() const;
  index_type first_undirected_leaf() const;
  index_type next_undirected_leaf(index_type l) const;
  index_type first_undirected_leaf(const bool_vec& n) const;
  index_type first_undirected_leaf(const index_set& n) const;
  index_type next_undirected_leaf(index_type l, const index_set& n) const;
  void fringe(const index_set& n, index_set& fn) const;
  void bi_fringe(const index_set& n, index_set& fn) const;

  void distance(index_type s0, index_vec& d) const;
  void distance(const index_set& s0, index_vec& d) const;
  index_type distance(index_type s0, index_type s1) const;
  index_type diameter() const;

  // find a shortest path s->t; returns path length (#edges), no_such_index
  // if no path exists; vector p contains path nodes, including s and t.
  index_type shortest_path(index_type s, index_type t, index_vec& p) const;
  // shortest path from any node in set s to any node in set t
  index_type shortest_path(const index_set& s,
			   const index_set& t,
			   index_vec& p) const;

  // find a shortest cycle in the graph; if a cycle is found, p contains
  // the cycle nodes, in sequence and without repetition, and the function
  // returns the cycle length; otherwise returns no_such_index.
  index_type shortest_cycle(index_vec& p) const;

  bool assign_node_level_top_down(index_vec& levels) const;

  // retrieve information about graph components
  // note: strongly_connected_components()/connected_components() must be
  //  called first to compute the components; methods are non-constant,
  //  because components are cached; thus need to be called again if the
  //  graph changes
  void strongly_connected_components();
  // use first method to keep the "belong to same component" equivalence
  // relation, second to ignore it
  void connected_components(equivalence& eq);
  void connected_components();

  index_type component(index_type i) const { return comp[i]; };
  const index_vec& component_map() const { return comp; };
  index_type n_components() const { return n_comp; };
  index_type component_node(index_type i) const;
  index_type component_size(index_type i) const;
  void component_node_set(index_type i, index_set& set) const;
  void component_node_sets(index_set_vec& sets) const;
  graph& component_tree(graph& cg) const;
  equivalence& component_partitioning(equivalence& eq) const;
  index_type maximal_non_unit_component() const;

  // extract subgraphs/induced graphs
  // note: methods below assign to first argument and returns same.
  // note: induced_partitioning is the same as the component partitioning
  //       on the induced undirected graph.
  graph& subgraph(graph& g, const index_set& n) const;
  graph& edge_subgraph(graph& g, const index_set& nodes) const;
  void induced_partitioning(equivalence& eq) const;
  void induced_partitioning(index_set_vec& cs) const;
  graph& induced_undirected_graph(graph& g) const;
  graph& minimal_equivalent_digraph(graph& g) const;
  graph& minimal_distance_graph(graph& g, const index_set& s0) const;
  graph& quotient(graph& g, const equivalence& eq) const;
  graph& complete_npart(graph& g, const index_set_vec& parts) const;
  graph& undirected_edge_graph(graph& g, pair_set& es) const;

  // compare graphs
  // note: the difference and distance methods assert that the graphs are
  //       of equal size.
  bool equals(const graph& g) const;
  bool equals(const graph& g, const index_vec& c) const;
  void difference(const graph& g,
		  const index_vec& c,
		  pair_set& d0,
		  pair_set& d1) const;
  void difference(const graph& g,
		  pair_set& d0,
		  pair_set& d1) const;
  index_type cardinality_of_difference(const graph& g) const;

  // standard order compare on graphs
  int compare(const graph& g) const;

  // hash the graph using a set hash function
  index_type hash(const set_hash_function& f) const;

  // basic graph modification
  void init(index_type s); // assign empty graph with s nodes
  void copy(const graph& g);
  // mapped copy
  void copy(const graph& g, const index_vec& map);
  // renamed copy: map must be a permutation
  void copy_and_rename(const graph& g, const index_vec& map);
  // add an isolated node, return index of the new node (== new size - 1)
  index_type add_node();
  void remove_node(index_type n);
  void remove_nodes(const index_set& ns);

  // add g as a subgraph; on return, m maps indices in g to indices in this
  void add_graph(const graph& g, mapping& m);

  void add_edge(index_type src, index_type dst);
  void add_edge(const index_set& srcs, index_type dst);
  void add_edge(index_type src, const index_set& dsts);
  void add_edge_to_transitive_closure(index_type src,
				      index_type dst);
  void add_edge_to_transitive_closure(index_type src,
				      index_type dst,
				      pair_set& e);

  void remove_edge(index_type src, index_type dst);
  void remove_edges_from(index_type src);
  void remove_edges_to(index_type dst);
  void remove_edges_incident_on(index_type n);
  void remove_edges(const pair_set& e);

  void add_undirected_edge(index_type n0, index_type n1);
  // add undirected edge between every pair of nodes in n0
  void add_undirected_edges(const index_set& n0);
  // add undirected edge between every node in n0 and every node in n1
  void add_undirected_edges(const index_set& n0, const index_set& n1);
  // add a set of undirected edges
  void add_undirected_edges(const pair_set& e);
  void remove_undirected_edge(index_type n0, index_type n1);
  void remove_undirected_edges(const pair_set& e);

  void clear_edges();
  void recalculate();

  // more graph modifications
  void complement();
  void complement_with_loops();
  void remove_loops();
  void make_undirected();
  void reverse();
  void transitive_closure();
  void missing_transitive_edges(pair_set& e) const;
  void transitive_reduction();
  void intersect(const graph& g);

  void randomize(count_type n, RNG& rnd);
  void randomize_undirected(count_type n, RNG& rnd);
  void randomize_connected(count_type n, RNG& rnd);
  void randomize_undirected_connected(count_type n, RNG& rnd);
  void randomize_strongly_connected(count_type n, RNG& rnd);

  void random_digraph(count_type n, RNG& rnd);
  void random_connected_digraph(count_type n, RNG& rnd);
  void random_strongly_connected_digraph(count_type n, RNG& rnd);
  void random_digraph_with_density(rational density, RNG& rnd);
  void random_tree(RNG& rnd);
  void random_tree(index_type b, index_type d, RNG& rnd);

  // clique/independent set finding
  bool is_clique(const index_set& nodes) const;
  index_type fillin(const index_set& nodes) const;
  bool is_independent(const index_set& nodes) const;
  bool is_independent_range(index_type l, index_type u) const;
  void maximal_clique(index_set& clique) const;
  void maximal_clique_including(index_type node, index_set& clique) const;
  void maximal_clique_cover(index_set_vec& sets) const;
  void all_maximal_cliques(index_set_vec& cliques) const;
  void all_maximal_cliques_including(index_type node, index_set_vec& cliques)
    const;
  void apx_independent_set(index_set& set) const;
  void apx_independent_set(const index_set& nodes, index_set& set) const;
  void apx_independent_set_including(index_type node, index_set& set) const;
  void apx_independent_set_cover(index_set_vec& sets) const;
  void apx_independent_set_disjoint_cover(index_set_vec& sets) const;

  void all_nondominated_cliques(index_set_vec &cliques) const;
    // I would have called this "all_maximal_cliques", but that name
    // already exists and does something else. :-(
  void all_cliques_geq(index_type k, index_set_vec& cliques) const;
  void one_maximal_clique(index_vec& clique) const;

  // tree decomposition methods
  index_type tree_decomposition_min_fillin(index_set_graph& td) const;

  // write methods
  void write_node_set(::std::ostream& s) const;
  void write_edge_set(::std::ostream& s) const;
  void write_compact(::std::ostream& s) const;
  void write_undirected_edge_set(::std::ostream& s) const;
  void write_adjacency_lists(::std::ostream& s) const;

  void dump(::std::ostream& s) const;

  // DOT write methods (for fancy formatting, use template method below)
  void write_digraph(::std::ostream& s,
		     bool with_node_indices,
		     const char* name) const;
  void write_undirected_graph(::std::ostream& s,
			      bool with_node_indices,
			      const char* name) const;
  void write_component_labeled_digraph(::std::ostream& s,
				       const char* name) const;
  void write_graph_correspondance(::std::ostream& s,
				  const graph& g,
				  const index_vec& c,
				  const char* name) const;


  // internal submethods; methods that are constant are declared public for
  // convenience, but not part of graph's "official" interface.
 public:
  void max_clique(index_set& sel,
		  index_type next,
		  index_set& clique) const;
  void all_max_cliques(index_set& sel,
		       index_type next,
		       index_set_vec& cliques) const;
  void all_nondominated_cliques_aux(index_set_vec &cliques,
                                    index_set &current_clique,
                                    const index_set &candidates,
				    index_type min,
				    bool one_maximal_only) const;
  void ramsey(const index_set& nodes, index_set& I, index_set& C) const;

  void undirected_dfs(index_type n, bool_vec& visited) const;

  void reachable(index_type n, bool_vec& v) const;
  void reverse_reachable(index_type n, bool_vec& v) const;

 private:
  void scc_first_dfs(index_type n, bool_vec& visited, index_vec& num) const;
  void scc_second_dfs(index_type n, bool_vec& visited, index_type c_id);
  index_type extract_path(const index_vec& d, index_type n, index_vec& p) const;
};

// write graph in DOT format with specified node labels; the type of the
// node label vector is a template argument: it must support subscripting
// (by index_type) and its elements must be printable (with <<).
// if c_id != no_such_index, the graph is printed as a cluster subgraph,
// and c_id is added to all node indices;
// else, if name != 0, the graph is printed as a digraph;
// else (i.e., c_id == no_such_index AND name == 0), only nodes and edges
// are printed (without enclosing graph declaration).
template<class LS>
void write_labeled_digraph
(std::ostream& s,
 const graph& g,
 const LS& ls,
 bool with_node_indices = false,
 const char* name = 0,
 index_type  c_id = no_such_index)
{
  if (c_id != no_such_index) {
    s << "subgraph cluster" << c_id << "{" << std::endl;
  }
  else if (name) {
    s << "digraph \"" << name << "\" {" << ::std::endl;
  }
  s << "node [width=0,height=0];" << ::std::endl;
  for (index_type k = 0; k < g.size(); k++) {
    if (with_node_indices) {
      s << "\t" << k + (c_id == no_such_index ? 0 : c_id)
	<< " [label=\"(" << k << ") " << ls[k] << "\"];"
	<< ::std::endl;
    }
    else {
      s << "\t" << k + (c_id == no_such_index ? 0 : c_id)
	<< " [label=\"" << ls[k] << "\"];"
	<< ::std::endl;
    }
  }
  for (index_type i = 0; i < g.size(); i++)
    for (index_type j = 0; j < g.size(); j++)
      if (g.adjacent(i, j))
	s << "\t" << i + (c_id == no_such_index ? 0 : c_id)
	  << " -> " << j + (c_id == no_such_index ? 0 : c_id)
	  << ";" << ::std::endl;
  if ((c_id != no_such_index) || (name != 0)) {
    s << "}" << ::std::endl;
  }
}

template<class LS>
void write_labeled_undirected_graph
(std::ostream& s,
 const graph& g,
 const LS& ls,
 bool with_node_indices = false,
 const char* name = 0)
{
  if (name) {
    s << "graph \"" << name << "\" {" << ::std::endl;
  }
  s << "node [width=0,height=0];" << ::std::endl;
  for (index_type k = 0; k < g.size(); k++) {
    if (with_node_indices) {
      s << "\t" << k
	<< " [label=\"(" << k << ") " << ls[k] << "\"];"
	<< ::std::endl;
    }
    else {
      s << "\t" << k
	<< " [label=\"" << ls[k] << "\"];"
	<< ::std::endl;
    }
  }
  for (index_type i = 0; i < g.size(); i++)
    for (index_type j = i + 1; j < g.size(); j++)
      if (g.bi_adjacent(i, j))
	s << "\t" << i
	  << " -- " << j
	  << ";" << ::std::endl;
  if (name != 0) {
    s << "}" << ::std::endl;
  }
}

typedef lvector<graph> graph_vec;

template<class N, class E> class labeled_graph : public graph {
 public:
  typedef std::map<index_type, N> node_label_map;
  typedef std::map<index_pair, E> edge_label_map;

  typedef lvector<N> node_label_vec;
  typedef lvector<E> edge_label_vec;

 protected:
  node_label_map _node_label;
  edge_label_map _edge_label;

 public:
  labeled_graph();
  labeled_graph(index_type size);
  labeled_graph(const graph& g);
  labeled_graph(const labeled_graph& g);
  labeled_graph(const graph& g, const index_set& nodes);
  labeled_graph(const labeled_graph& g, const index_set& nodes);
  labeled_graph(const graph& g, const equivalence& eq);
  ~labeled_graph();

  N& node_label(index_type n);
  E& edge_label(index_type i, index_type j);
  const N& node_label(index_type n) const;
  const E& edge_label(index_type i, index_type j) const;
  bool node_has_label(index_type n) const;
  bool edge_has_label(index_type i, index_type j) const;

  labeled_graph& subgraph(labeled_graph& g, const index_set& n) const;

  index_type node_with_label(const N& l) const;
  index_pair edge_with_label(const E& l) const;

  void init(index_type size);
  void init(index_type size, const N& n, const E& e);
  void clear_node_labels();
  void clear_edge_labels();
  void clear_labels();
  void copy(const graph& g);
  void copy(const labeled_graph& g);

  void add_graph(const graph& g, mapping& m);
  void add_graph(const labeled_graph& g, mapping& m);

  void add_edge(index_type src, index_type dst);
  void add_edge(index_type src, index_type dst, const E& lbl);
  void add_edge(const index_set& srcs, index_type dst);
  void add_edge(const index_set& srcs, index_type dst, const E& lbl);
  void add_edge(index_type src, const index_set& dsts);
  void add_edge(index_type src, const index_set& dsts, const E& lbl);

  index_type add_node();
  index_type add_node(const N& l);

  void remove_node(index_type n);
  void remove_edge(index_type src, index_type dst);
  void remove_edges_from(index_type src);
  void remove_edges_to(index_type dst);
  void remove_edges_incident_on(index_type n);
  void remove_undirected_edge(index_type n0, index_type n1);
  void remove_edges(const pair_set& e);
  void remove_undirected_edges(const pair_set& e);
  void clear_edges();

  void remove_edges_with_label(const E& l);

  void write_digraph(::std::ostream& s,
		     bool with_node_indices,
		     bool with_node_labels,
		     bool with_edge_labels,
		     bool compact_edges,
		     const char* name) const;
  void write_matrix(::std::ostream& s,
		    const char* unlabeled_edge,
		    const char* missing_edge) const;
};


// inlines, labeled_graph

template<class N, class E>
labeled_graph<N,E>::labeled_graph()
{
  // done
}

template<class N, class E>
labeled_graph<N,E>::labeled_graph(index_type s)
  : graph(s)
{
  // done
}

template<class N, class E>
labeled_graph<N,E>::labeled_graph(const graph& g)
  : graph(g)
{
  // done
}

template<class N, class E>
labeled_graph<N,E>::labeled_graph(const labeled_graph& g)
  : graph(g), _node_label(g._node_label), _edge_label(g._edge_label)
{
  // done
}

template<class N, class E>
labeled_graph<N,E>::labeled_graph(const graph& g, const index_set& n)
{
  g.subgraph(*this, n);
}

template<class N, class E>
labeled_graph<N,E>::labeled_graph
(const labeled_graph& g, const index_set& n)
{
  g.subgraph(*this, n);
}

template<class N, class E>
labeled_graph<N,E>::labeled_graph
(const graph& g, const equivalence& eq)
{
  g.quotient(*this, eq);
}

template<class N, class E>
labeled_graph<N,E>::~labeled_graph()
{
  // done
}

template<class N, class E>
N& labeled_graph<N,E>::node_label(index_type n)
{
  assert(n < size());
  return _node_label[n];
}

template<class N, class E>
E& labeled_graph<N,E>::edge_label(index_type i, index_type j)
{
  assert((i < size()) && (j < size()));
  return _edge_label[index_pair(i, j)];
}

template<class N, class E>
const N& labeled_graph<N,E>::node_label(index_type n) const
{
  assert(node_has_label(n));
  return (_node_label.find(n)->second);
}

template<class N, class E>
const E& labeled_graph<N,E>::edge_label(index_type i, index_type j) const
{
  assert(edge_has_label(i, j));
  return (_edge_label.find(index_pair(i, j))->second);
}

template<class N, class E>
bool labeled_graph<N,E>::node_has_label(index_type n) const
{
  assert(n < size());
  return (_node_label.find(n) != _node_label.end());
}

template<class N, class E>
bool labeled_graph<N,E>::edge_has_label(index_type i, index_type j) const
{
  assert((i < size()) && (j < size()));
  return (_edge_label.find(index_pair(i, j)) != _edge_label.end());
}

template<class N, class E>
void labeled_graph<N,E>::init(index_type s)
{
  graph::init(s);
  _node_label.clear();
  _edge_label.clear();
}

template<class N, class E>
void labeled_graph<N,E>::clear_node_labels()
{
  _node_label.clear();
}

template<class N, class E>
void labeled_graph<N,E>::clear_edge_labels()
{
  _edge_label.clear();
}

template<class N, class E>
void labeled_graph<N,E>::clear_labels()
{
  _node_label.clear();
  _edge_label.clear();
}

template<class N, class E>
void labeled_graph<N,E>::copy(const labeled_graph& g)
{
  graph::copy(g);
  _node_label = g._node_label;
  _edge_label = g._edge_label;
}

template<class N, class E>
labeled_graph<N,E>& labeled_graph<N,E>::subgraph
(labeled_graph& g, const index_set& n) const
{
  g.init(n.length());
  for (index_type k = 0; k < n.length(); k++) {
    assert(n[k] < size());
    if (node_has_label(n[k])) {
      g.node_label(k) = node_label(n[k]);
    }
  }
  for (index_type i = 0; i < n.length(); i++)
    for (index_type j = 0; j < n.length(); j++)
      if (adjacent(n[i], n[j])) {
	g.add_edge(i, j);
	if (edge_has_label(n[i], n[j])) {
	  g.edge_label(i, j) = edge_label(n[i], n[j]);
	}
      }
}

template<class N, class E>
index_type labeled_graph<N,E>::node_with_label(const N& l) const
{
  for (index_type k = 0; k < size(); k++)
    if (node_has_label(k))
      if (node_label(k) == l)
	return k;
  return no_such_index;
}

template<class N, class E>
index_pair labeled_graph<N,E>::edge_with_label(const E& l) const
{
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j) && edge_has_label(i, j))
	if (edge_label(i, j) == l)
	  return index_pair(i, j);
  return no_such_index;
}

template<class N, class E>
void labeled_graph<N,E>::add_graph
(const graph& g, mapping& m)
{
  graph::add_graph(g, m);
}

template<class N, class E>
void labeled_graph<N,E>::add_graph
(const labeled_graph& g, mapping& m)
{
  graph::add_graph(g, m);
  for (index_type i = 0; i < g.size(); i++) {
    if (g.node_has_label(i))
      node_label(m[i]) = g.node_label(i);
    for (index_type j = 0; j < g.size(); j++)
      if (g.adjacent(i, j) && g.edge_has_label(i, j))
	edge_label(m[i], m[j]) = g.edge_label(i, j);
  }
}

template<class N, class E>
void labeled_graph<N,E>::add_edge
(index_type src, index_type dst)
{
  graph::add_edge(src, dst);
}

template<class N, class E>
void labeled_graph<N,E>::add_edge
(index_type src, index_type dst, const E& lbl)
{
  graph::add_edge(src, dst);
  edge_label(src, dst) = lbl;
}

template<class N, class E>
void labeled_graph<N,E>::add_edge
(const index_set& srcs, index_type dst)
{
  graph::add_edge(srcs, dst);
}

template<class N, class E>
void labeled_graph<N,E>::add_edge
(const index_set& srcs, index_type dst, const E& lbl)
{
  graph::add_edge(srcs, dst);
  for (index_type k = 0; k < srcs.length(); k++)
    edge_label(srcs[k], dst) = lbl;
}

template<class N, class E>
void labeled_graph<N,E>::add_edge
(index_type src, const index_set& dsts)
{
  graph::add_edge(src, dsts);
}

template<class N, class E>
void labeled_graph<N,E>::add_edge
(index_type src, const index_set& dsts, const E& lbl)
{
  graph::add_edge(src, dsts);
  for (index_type k = 0; k < dsts.length(); k++)
    edge_label(src, dsts[k]) = lbl;
}

template<class N, class E>
index_type labeled_graph<N,E>::add_node()
{
  return graph::add_node();
}

template<class N, class E>
index_type labeled_graph<N,E>::add_node(const N& l)
{
  index_type n = graph::add_node();
  node_label(n) = l;
  return n;
}

template<class N, class E>
void labeled_graph<N,E>::remove_node(index_type n)
{
  assert(n < size());
  labeled_graph g(*this);
  index_set ns;
  ns.fill(size());
  ns.subtract(n);
  g.subgraph(*this, ns);
//  // move node labels
//  for (index_type k = n + 1; k < size(); k++)
//    if (node_has_label(k))
//      node_label(k - 1) = node_label(k);
//  // erase node label now out of range
//  if (node_has_label(size() - 1))
//    _node_label.erase(size() - 1);
//  for (index_type i = 0; i < size(); i++)
//    for (index_type j = 0; j < size(); j++)
//      if (adjacent(i, j) && edge_has_label(i, j)) {
//	if ((i == n) || (j == n)) {
//	  _edge_label.erase(index_pair(i, j));
//	}
//	else {
//	  index_type s = (i > n ? i - 1 : i);
//	  index_type d = (j > n ? j - 1 : j);
//	  edge_label(s, d) = edge_label(i, j);
//	}
//      }
//  graph::remove_node(n);
}

template<class N, class E>
void labeled_graph<N,E>::remove_edge(index_type src, index_type dst)
{
  _edge_label.erase(index_pair(src, dst));
  graph::remove_edge(src, dst);
}

template<class N, class E>
void labeled_graph<N,E>::remove_undirected_edge(index_type n0, index_type n1)
{
  _edge_label.erase(index_pair(n0, n1));
  _edge_label.erase(index_pair(n1, n0));
  graph::remove_undirected_edge(n0, n1);
}

template<class N, class E>
void labeled_graph<N,E>::remove_edges_from(index_type src)
{
  index_set ns(successors(src));
  for (index_type k = 0; k < ns.length(); k++)
    remove_edge(src, ns[k]);
}

template<class N, class E>
void labeled_graph<N,E>::remove_edges_to(index_type dst)
{
  index_set ns(predecessors(dst));
  for (index_type k = 0; k < ns.length(); k++)
    remove_edge(ns[k], dst);
}

template<class N, class E>
void labeled_graph<N,E>::remove_edges_incident_on(index_type n)
{
  remove_edges_from(n);
  remove_edges_to(n);
}

template<class N, class E>
void labeled_graph<N,E>::remove_edges(const pair_set& e)
{
  for (index_type k = 0; k < e.length(); k++) {
    remove_edge(e[k].first, e[k].second);
  }
}

template<class N, class E>
void labeled_graph<N,E>::remove_undirected_edges(const pair_set& e)
{
  for (index_type k = 0; k < e.length(); k++) {
    remove_undirected_edge(e[k].first, e[k].second);
  }
}

template<class N, class E>
void labeled_graph<N,E>::clear_edges()
{
  _edge_label.clear();
  graph::clear_edges();
}

template<class N, class E>
void labeled_graph<N,E>::remove_edges_with_label(const E& l)
{
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j) && (edge_label(i, j) == l))
	remove_edge(i, j);
}

template<class N, class E>
void labeled_graph<N,E>::write_digraph
(::std::ostream& s,
 bool with_node_indices,
 bool with_node_labels,
 bool with_edge_labels,
 bool compact_edges,
 const char* name) const
{
  if (strncmp(name, "cluster", 7) == 0)
    s << "subgraph";
  else
    s << "digraph";
  s << " \"" << name << "\"" << ::std::endl << "{" << ::std::endl;
  if (with_node_indices || with_node_labels) {
    s << "\tnode [shape=ellipse];" << ::std::endl;
  }
  else {
    s << "\tnode [shape=point];" << ::std::endl;
  }
  for (index_type i = 0; i < size(); i++) {
    s << "\t" << i;
    if (with_node_indices || with_node_labels) {
      s << " [label=\"";
      if (with_node_indices) {
	if (with_node_labels) {
	  s << i << ": ";
	  if (node_has_label(i)) {
	    s << node_label(i);
	  }
	}
	else {
	  s << i;
	}
      }
      else {
	if (node_has_label(i)) {
	  s << node_label(i);
	}
      }
      s << "\"]";
    }
    s << ";" << std::endl;
  }
  for (index_type i = 0; i < size(); i++)
    for (index_type j = 0; j < size(); j++)
      if (adjacent(i, j)) {
	if (adjacent(j, i) && compact_edges) {
	  if (i < j) {
	    s << "\t" << i << " -> " << j;
	    s << " [dir=both";
	    if (with_edge_labels &&
		(edge_has_label(i, j) || edge_has_label(j, i))) {
	      if (!edge_has_label(i, j)) {
		s << ",label=\"" << edge_label(j, i) << "\"";
	      }
	      else if (!edge_has_label(j, i)) {
		s << ",label=\"" << edge_label(i, j) << "\"";
	      }
	      else if (edge_label(i, j) == edge_label(j, i)) {
		s << ",label=\"" << edge_label(i, j) << "\"";
	      }
	      else {
		s << ",label=\"" << edge_label(i, j) << ", "
		  << edge_label(j, i) << "\"";
	      }
	    }
	    s << "]" << ::std::endl;
	  }
	}
	else {
	  s << "\t" << i << " -> " << j;
	  if (with_edge_labels && edge_has_label(i, j)) {
	    s << " [label=\"" << edge_label(i, j) << "\"]";
	  }
	  s << ";" << ::std::endl;
	}
      }
  s << "}" << ::std::endl;
}

template<class N, class E>
void labeled_graph<N,E>::write_matrix
(::std::ostream& s, const char* unlabeled_edge, const char* missing_edge) const
{
  s << '[';
  for (index_type i = 0; i < size(); i++) {
    if (i > 0) s << ' ';
    s << '[';
    for (index_type j = 0; j < size(); j++) {
      if (j > 0) s << ',';
      if (adjacent(i, j)) {
	if (edge_has_label(i, j)) {
	  s << edge_label(i, j);
	}
	else {
	  s << unlabeled_edge;
	}
      }
      else {
	s << missing_edge;
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

template<class N, class E, class F>
void apply_to_nodes(labeled_graph<N,E>& g, F& f)
{
  for (index_type k = 0; k < g.size(); k++)
    if (g.node_has_label(k))
      g.node_label(k) = f(g.node_label(k));
}

// class index_graph

class index_graph : public labeled_graph<index_type,index_type> {
  static const index_type NODE_SHAPE = 2 + 4 + 8 + 16;
  static const index_type NODE_STYLE = 64 + 128 + 256;
  static const index_type EDGE_STYLE = 64 + 128;
  static const index_type EDGE_DIR = 512 + 1024;

 public:
  static const index_type NS_CIRCLE    = 0;
  static const index_type NS_ELLIPSE   = 2;
  static const index_type NS_BOX       = 4;
  static const index_type NS_POINT     = 6;
  static const index_type NS_DIAMOND   = 8;
  static const index_type NS_HEXAGON   = 10;
  static const index_type NS_OCTAGON   = 12;
  static const index_type NS_PLAINTEXT = 14;

  static const index_type NS_NORMAL    = 0;
  static const index_type NS_DOUBLE    = 32;
  static const index_type NS_BOLD      = 64;
  static const index_type NS_DASHED    = 128;
  static const index_type NS_DOTTED    = 192;
  static const index_type NS_FILLED    = 256;

  static const index_type ED_NONE      = 0;
  static const index_type ED_FORWARD   = 512;
  static const index_type ED_BACK      = 1024;
  static const index_type ED_BOTH      = ED_FORWARD + ED_BACK;

  static const index_type ES_NORMAL    = NS_NORMAL;
  static const index_type ES_BOLD      = NS_BOLD;
  static const index_type ES_DASHED    = NS_DASHED;
  static const index_type ES_DOTTED    = NS_DOTTED;

  static const index_type STYLE_MAX = 2048;

  index_graph()
    : labeled_graph<index_type, index_type>() { };
  index_graph(index_type size)
    : labeled_graph<index_type, index_type>(size) { };
  index_graph(const graph& g)
    : labeled_graph<index_type, index_type>(g) { };
  index_graph(const index_graph& g)
    : labeled_graph<index_type, index_type>(g) { };
  index_graph(const graph& g, const index_set& nodes)
    : labeled_graph<index_type, index_type>(g, nodes) { };
  index_graph(const index_graph& g, const index_set& nodes)
    : labeled_graph<index_type, index_type>(g, nodes) { };
  index_graph(const graph& g, const equivalence& eq)
    : labeled_graph<index_type, index_type>(g, eq) { };
  ~index_graph() { };

  // reverse edges but simulatneously reverse ED label
  void reverse();

  // "reflect around main diagonal": for every edge (i, j) add edge (j, i)
  // with opposite ED label; asserts that (j, i) does not already exist
  void reflect();

  static void write_node_style(std::ostream& s, index_type l);
  static void write_edge_style(std::ostream& s, index_type l);

  void write_styled_digraph(std::ostream& s,
			    bool with_node_indices = false,
			    const char* name = 0,
			    index_type c_id = no_such_index) const;
  void write_matrix(std::ostream& s) const;
  void write_MATLAB(std::ostream& s,
		    const char* n,
		    const char* t) const;
};


template<class LS>
void write_styled_digraph
(std::ostream& s,
 const index_graph& g,
 const LS& ls,
 bool with_node_indices = false,
 const char* name = 0,
 index_type  c_id = no_such_index)
{
  if (c_id != no_such_index) {
    s << "subgraph cluster" << c_id << " {" << std::endl;
    s << "node [width=0.5,height=0.5];" << ::std::endl;
  }
  else if (name) {
    s << "digraph \"" << name << "\" {" << std::endl;
    s << "node [width=0.5,height=0.5];" << ::std::endl;
  }
  for (index_type k = 0; k < g.size(); k++) {
    s << "\t" << k + (c_id != no_such_index ? c_id : 0) << " [";
    index_graph::write_node_style(s, g.node_has_label(k) ? g.node_label(k) : 0);
    if (with_node_indices)
      s << ",label=\"(" << k << ") " << ls[k] << "\"];" << std::endl;
    else
      s << ",label=\"" << ls[k] << "\"];" << std::endl;
  }
  for (index_type i = 0; i < g.size(); i++)
    for (index_type j = 0; j < g.size(); j++)
      if (g.adjacent(i, j)) {
	s << "\t" << i + (c_id != no_such_index ? c_id : 0)
	  << " -> " << j + (c_id != no_such_index ? c_id : 0) << " [";
	index_graph::write_edge_style(s, g.edge_has_label(i, j) ? g.edge_label(i, j) : 0);
	s << "];" << std::endl;
      }
  if ((c_id != no_such_index) || (name != 0)) {
    s << "}" << ::std::endl;
  }
}


// class weighted_graph

class weighted_graph : public labeled_graph<NTYPE,NTYPE> {
 public:
  weighted_graph() { };
  weighted_graph(index_type s)
    : labeled_graph<NTYPE,NTYPE>(s) { };
  weighted_graph(const graph& g)
    : labeled_graph<NTYPE,NTYPE>(g) { };
  weighted_graph(const graph& g, NTYPE w)
    : labeled_graph<NTYPE,NTYPE>(g) { set_node_weight(w); };
  weighted_graph(const weighted_graph& g)
    : labeled_graph<NTYPE,NTYPE>(g) { };
  weighted_graph(const graph& g, const index_set& n)
    : labeled_graph<NTYPE,NTYPE>(g, n) { };
  weighted_graph(const weighted_graph& g, const index_set& n)
    : labeled_graph<NTYPE,NTYPE>(g, n) { };
  weighted_graph(const graph& g, const equivalence& eq)
    : labeled_graph<NTYPE,NTYPE>(g, eq) { };
  weighted_graph(const weighted_graph& g, const equivalence& eq);
  weighted_graph(const weighted_graph& g, const index_vec& m)
    : labeled_graph<NTYPE,NTYPE>() { copy(g, m); };
  ~weighted_graph() { };

  // weighted quotient construction: sums weights of merged nodes and edges
  weighted_graph& quotient(weighted_graph& g, const equivalence& eq) const;

  void copy(const graph& g) { labeled_graph<NTYPE,NTYPE>::copy(g); };
  void copy(const weighted_graph& g) { labeled_graph<NTYPE,NTYPE>::copy(g); };
  void copy(const graph& g, const index_vec& map) {
    clear_labels();
    graph::copy(g, map);
  };
  void copy(const weighted_graph& g, const index_vec& map);

  // approximate maximal weighted independent set (considering node weights),
  // two different algorithms; last method returns max of both.
  NTYPE apx_weighted_independent_set_1(index_set& set) const;
  NTYPE apx_weighted_independent_set_2(index_set& set) const;
  NTYPE apx_weighted_independent_set(index_set& set) const;

  void add_edge(index_type src, index_type dst);
  void add_edge(index_type src, index_type dst, NTYPE w);
  void add_undirected_edge(index_type n0, index_type n1);
  void add_undirected_edge(index_type n0, index_type n1, NTYPE w);

  // get/set node/edge weight; weight defaults to 0 if not set
  NTYPE weight(index_type n) const;
  NTYPE weight(index_type n0, index_type n1) const;
  // weight of a set of nodes = sum of weights of nodes in set
  NTYPE weight(const index_set& ns) const;
  NTYPE weight(const bool_vec& ns) const;
  // weight of the bidirectional fringe around ns:
  NTYPE bi_fringe_weight(const index_set& ns) const;
  void set_weight(index_type n, NTYPE w);
  void set_weight(index_type n0, index_type n1, NTYPE w);
  void increment_node_weight(index_type n, NTYPE w);
  void increment_edge_weight(index_type src, index_type dst, NTYPE w);
  // set weight of every node/edge
  void set_node_weight(NTYPE w);
  void set_edge_weight(NTYPE w);
  NTYPE max_node_weight() const;

  // get the edge weight matrix (+INF where no edge)
  void edge_matrix(cost_matrix& mat) const;

  // transitive closure adjusts edge weights to sum along least-weight path
  void transitive_closure();

  // simple CPM scheduler, considers both node and edge weights as times;
  // return value is the makespan
  NTYPE critical_path(cost_vec& s);

  // MST implements Kruskal's algorithm; assumes that the graph is
  // undirected, with symmetric edge weights; second form copies weights
  // of edges in the tree
  void minimum_spanning_tree(graph& mst);
  void minimum_spanning_tree(weighted_graph& mst);

  // recursive tree decomposition (splitting on min weight separating
  // vertex sets)
  void recursive_tree_decomposition
    (NTYPE alpha, NTYPE delta, index_set_graph& t) const;

  // (approx) min weight separating vertex set, with soft balance
  // constraint (alpha = 0: no balance constraint; alpha = POS_INF:
  // perfect balance required) and counting separator weight as
  // weight(s) + delta * weight(neighbourhood(s))
  NTYPE min_vs(index_set& s, index_set_vec& p, NTYPE alpha, NTYPE delta) const;

  // max s->t flow, interpreting edge weights as capacity (node weights
  // are ignored); second and third method return the residual graph and
  // the flow assignment, respectively.
  NTYPE max_flow(index_type s, index_type t);
  NTYPE max_flow(index_type s, index_type t, weighted_graph& rg);
  NTYPE max_flow(index_type s, index_type t, cost_matrix& f);

  // min s-t cut (based on max flow); first method outputs the
  // s-partition, second outputs edges in the cut
  NTYPE min_cut(index_type s, index_type t, bool_vec& s_set);
  NTYPE min_cut(index_type s, index_type t, pair_set& e_set);

  index_pair max_weight_edge() const;
  index_pair min_weight_edge() const;
  void min_and_max_edges(const index_set& nodes,
			 pair_set& e_min, NTYPE& w_min,
			 pair_set& e_max, NTYPE& w_max) const;

  // note: implementation of max match is STUPID and EXPONENTIAL!
  NTYPE maximal_matching(weighted_graph& matching);
  NTYPE apx_matching(bool_vec& nodes);

  void write_node_set(::std::ostream& s) const;
  void write_edge_set(::std::ostream& s) const;
  void write_compact(::std::ostream& s) const;

  void write_matrix(::std::ostream& s) const;

 private:
  NTYPE augmenting_path(index_type s, index_type t,
			const cost_matrix& f, pair_vec& p);
  bool vs_create_initial(index_set& s, index_type i) const;
  void vs_compute_partitions(const index_set& s, index_set_vec& p) const;
  NTYPE vs_value(const index_set& s, const index_set_vec& p,
		 NTYPE alpha, NTYPE delta) const;
  NTYPE vs_eval_move(const index_set& s, const index_set_vec& p,
		     const index_set& z, index_type d,
		     NTYPE alpha, NTYPE delta) const;
  NTYPE vs_best_move(const index_set& s, const index_set_vec& p,
		     index_type first_i, index_set& z, index_type& d,
		     NTYPE alpha, NTYPE delta) const;
  NTYPE vs_improve(index_set& s, index_set_vec& p,
		   NTYPE alpha, NTYPE delta) const;
};

class index_set_graph : public labeled_graph<index_set,index_set> {
 public:
  index_set_graph() { };
  index_set_graph(index_type s)
    : labeled_graph<index_set,index_set>(s) { };
  index_set_graph(const graph& g)
    : labeled_graph<index_set,index_set>(g) { };
  index_set_graph(const index_set_graph& g)
    : labeled_graph<index_set,index_set>(g) { };
  index_set_graph(const graph& g, const index_set& n)
    : labeled_graph<index_set,index_set>(g, n) { };
  index_set_graph(const index_set_graph& g, const index_set& n)
    : labeled_graph<index_set,index_set>(g, n) { };
  // note: quotient construction from standard graph is specialised
  // (assigns sets of merged nodes to node labels)
  index_set_graph(const graph& g, const equivalence& eq);
  // quotient construction from index_set_graph assigns to each node
  // label union of node labels in the same equivalence class (calls
  // quotient method below)
  index_set_graph(const index_set_graph& g, const equivalence& eq);
  ~index_set_graph() { };

  void add_edge(index_type src, index_type dst);
  void add_edge(index_type src, index_type dst, const index_set& lbl);
  void add_edge(const index_set& srcs, index_type dst);
  void add_edge(const index_set& srcs, index_type dst, const index_set& lbl);
  void add_edge(index_type src, const index_set& dsts);
  void add_edge(index_type src, const index_set& dsts, const index_set& lbl);

  // add new label to edge if it already exists, else add edge with
  // singleton set label
  void add_edge(index_type src, index_type dst, index_type newlbl);
  void add_edge(const index_set& srcs, index_type dst, index_type newlbl);
  void add_edge(index_type src, const index_set& dsts, index_type newlbl);

  // apply a mapping to node/edge labels
  void remap_node_labels(const index_vec& map);
  void remap_edge_labels(const index_vec& map);

  // cardinality of set labelling node n; 0 if node has no label
  index_type cardinality(index_type n) const;
  index_type max_node_cardinality() const;

  index_type max_cardinality_undirected_leaf() const;
  index_type max_cardinality_undirected_leaf(const index_set& n) const;

  // return index of (first) node whose label contains s (or no_such_index)
  index_type find_node_label_contains(const index_set& s) const;

  // union_reachable: assigns to every nodes label the union of labels
  // of nodes reachable from that node (done in place).
  void union_reachable();

  // compute the union of edge labels on a shortest path from s to t;
  // returns the length of the path (== no_such_index iff t not
  // reachable from s) and the union in 'u'.
  index_type union_of_edges_on_path
    (index_type s, index_type t, index_set& u) const;

  // compute the union of edge labels on all paths from s to t;
  // returns true iff t reachable from s and the union in 'u'.
  // the strict version excludes from u the label on the edge
  // from s to t (if any).
  bool union_of_edges_between
    (index_type s, index_type t, index_set& u) const;
  bool union_of_edges_strictly_between
    (index_type s, index_type t, index_set& u) const;

  void edge_label_preserving_transitive_reduction();

  void merge_labels(const index_set& ns);
  void merge_labels_upwards();
  void merge_labels_downwards();

  // contract two nodes: all labels are assigned union.
  void contract(index_type n0, index_type n1);

  // assign g the quotient graph of this under equivalence eq. node
  // labels in g are assigned the union of node labels in the
  // corresponding equivalence class. edge labels in the original graph
  // (this) are ignored.
  index_set_graph& quotient(index_set_graph& g, const equivalence& eq) const;
  // edge_label_quotient also assigns g the quotient graph of this under
  // equivalence relation eq, but computes node and edge labels differently:
  // node labels are assigned the union of all edge labels between pairs of
  // nodes in the corresponding equivalence class; edge labels are assigned
  // the union of labels on edges between nodes in corresponding classes.
  // node labels in original graph (this) are ignored.
  index_set_graph& edge_label_quotient(index_set_graph& g,
				       const equivalence& eq) const;
  index_set_graph& union_reachable(index_set_graph& g) const;
  index_set_graph& subgraph_set_size_gt(index_set_graph& g, index_type l);

  void write_edge_set(::std::ostream& s) const;
  void write_undirected_edge_set(::std::ostream& s) const;
  void write_DOT(::std::ostream& s, bool undirected, const char* name) const;

  void write_digraph(::std::ostream& s, const char* name) const
  { write_DOT(s, false, name); }
  void write_undirected_graph(::std::ostream& s, const char* name) const
  { write_DOT(s, true, name); }
};


// stream operators

inline ::std::ostream& operator<<(::std::ostream& s, const graph& g)
{
  g.write_compact(s);
  return s;
};

inline ::std::ostream& operator<<(::std::ostream& s, const adjacency_list_graph& g)
{
  g.write_compact(s);
  return s;
};

inline ::std::ostream& operator<<(::std::ostream& s, const weighted_graph& g)
{
  g.write_compact(s);
  return s;
};

inline ::std::ostream& operator<<(::std::ostream& s, const index_set_graph& g)
{
  g.write_edge_set(s);
  return s;
};


END_HSPS_NAMESPACE

#endif
