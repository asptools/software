#ifndef BFS_H
#define BFS_H

#include "config.h"
#include "search_base.h"
#include "nodeset.h"
#include "problem.h"

BEGIN_HSPS_NAMESPACE

class BFS : public SingleSearchAlgorithm {

 protected:
  static const count_type TRACE_LEVEL_2_NOTIFY = 10000;

  HashNodeSet graph;
  NodeQueue   queue;
  Node* current_node;
  Node* current_sol;
  NTYPE best_node_cost;

  count_type acc_succ;
  count_type new_succ;

  void update_path(Node* n, Node* p, Transition* t, NTYPE d);

  NTYPE weight;

 public:
  static Instance* trace_print_instance;
  bool greedy;

  NodeSet& state_space() { return graph; };
  NodeQueue& open_queue() { return queue; };

  BFS(Statistics& s, SearchResult& r);
  BFS(Statistics& s, SearchResult& r, index_type nt_size);
  virtual ~BFS();

  void set_weight(NTYPE w) { weight = w; };

  virtual NTYPE start(State& s, NTYPE b);
  virtual NTYPE start(State& s);
  virtual NTYPE resume();

  virtual NTYPE main();
  virtual NTYPE new_state(State& s, NTYPE bound);

  virtual NTYPE cost() const;
  virtual bool done() const;
};

class BFS_PX : public BFS {
  NTYPE threshold;
 public:
  BFS_PX(Statistics& s, SearchResult& r, NTYPE thresh)
    : BFS(s, r), threshold(thresh) { };
  BFS_PX(Statistics& s, SearchResult& r, index_type nt_size, NTYPE thresh)
    : BFS(s, r, nt_size), threshold(thresh) { };
  virtual ~BFS_PX() { };

  virtual NTYPE main();
};

END_HSPS_NAMESPACE

#endif
