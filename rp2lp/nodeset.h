#ifndef NODESET_H
#define NODESET_H

#include "config.h"
#include "index_type.h"
#include "search.h"
// #include "stats.h"
#include "name.h"
#include <map>
#include "function.h"

BEGIN_HSPS_NAMESPACE

class Node;
typedef Node* nodep;
typedef lvector<nodep> node_vec;

struct Link {
  Node* node;
  NTYPE delta;

  Link() : node(0), delta(0) { };
  Link(Node* s, NTYPE d) : node(s), delta(d) { };
};

typedef lvector<Link> link_vec;
typedef std::map<Transition*, Link, compare_as_less_ptr<Transition*> > link_map;

class Node {
 public:
  index_type id;
  State*     state;
  link_map   succ;
  bool       closed;

  NTYPE      acc;
  NTYPE      est;
  NTYPE      val;
  NTYPE      opt;
  index_type pos;
  count_type exp;

  Node*       bp_pre;
  Transition* bp_trans;
  NTYPE       bp_delta;
  link_vec*   all_pre;

  Node() : id(0), state(0), closed(false), acc(0), est(0), val(0),
    opt(POS_INF), pos(no_such_index), exp(0), bp_pre(0), bp_delta(0),
    bp_trans(0), all_pre(0) { };
  ~Node();

  NTYPE min_delta_to(Node* n) const;
  bool has_successor(Node* n) const;
  bool has_successor(index_type id) const;

  // backup performs a bellman backup of the est value based on
  // the est values of successors and link deltas; if the node's
  // est value increases, it recursively calls backup on its
  // parent (bp_pre), if any. This could be used as a hook to
  // trigger some kind of heuristic learning. Future work...
  void backup();

  bool solved() const { return FINITE(opt); };
  void cache_optimal_path(NTYPE c);
  bool back_path_contains(Node* n) const;
  void add_predecessor(Node* n, NTYPE d);

  void write(std::ostream& s, const Name* p = 0);
  void write_short(std::ostream& s, const Name* p = 0);

  void write_back_path(std::ostream& s);

  void write_graph_node(std::ostream& s);
  void write_graph_edges(std::ostream& s);
};

// default node order is the standard A* order: increasing on f-value,
// tie-breaking on decreasing h-value.
class DefaultNodeOrder : public node_vec::order {
 public:
  DefaultNodeOrder() { };
  virtual bool operator()(const nodep& v0, const nodep& v1) const;
};

class NodeQueue : public node_vec {
  static DefaultNodeOrder default_node_order;

  const node_vec::order& before;
 public:
  NodeQueue(const node_vec::order& b = default_node_order);
  ~NodeQueue();

  void  check_queue();

  Node* peek();
  void  enqueue(Node*);
  Node* dequeue();
  void  clear();

  // these should be private, but BFS uses them when updating costs/paths
  void shift_up(index_type i);
  void shift_down(index_type i);
};

class NodeSet {
 protected:
  static index_type next_id;
  node_vec roots;

 public:
  static bool write_state_in_graph_node;

  virtual ~NodeSet();

  virtual Node* insert_node(State& s) = 0;
  virtual Node* find_node(State& s) = 0;
  virtual void  clear() = 0;
  virtual void  collect_nodes(node_vec& ns) = 0;

  Node*  insert_root_node(State& s);
  void   make_root(Node* n);
  node_vec& root_nodes();

  void   compute_reverse_links(node_vec& v);
  void   compute_reverse_links();
  void   mark_solved_apsp();
  void   mark_solved();

  void   set_back_path_solution_cost(Node* n, NTYPE c_sol);
  void   cache_pg(NTYPE c_sol);

  void   back_path_to_sequence(Node* n, node_vec& ns);

  void  write_short(std::ostream& s, const Name* p = 0);
  void  write_graph(std::ostream& s);
  void  write_graph_compact(std::ostream& s);
  void  write_graph_rainbow(std::ostream& s);
};

class TreeNodeSet : public Node, public NodeSet {
  TreeNodeSet* left;
  TreeNodeSet* right;
 public:
  TreeNodeSet();
  virtual ~TreeNodeSet();

  virtual Node* insert_node(State& s);
  virtual Node* find_node(State& s);
  virtual void  clear();

  virtual void collect_nodes(node_vec& ns);
};

class HashNodeSet : public NodeSet {
  index_type    size;
  TreeNodeSet** tab;
 public:
  HashNodeSet(index_type s = 31337);
  virtual ~HashNodeSet();

  virtual Node* insert_node(State& s);
  virtual Node* find_node(State& s);
  virtual void  clear();

  virtual void collect_nodes(node_vec& ns);
};

END_HSPS_NAMESPACE

#endif
