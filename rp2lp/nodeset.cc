
#include "nodeset.h"

#include <cmath>

BEGIN_HSPS_NAMESPACE

Node::~Node()
{
  if (state) delete state;
  if (bp_trans) delete bp_trans;
}

NTYPE Node::min_delta_to(Node* n) const
{
  NTYPE d = POS_INF;
  for (link_map::const_iterator p = succ.begin(); p != succ.end(); p++)
    if (p->second.node == n)
      d = MIN(d, p->second.delta);
  return d;
}

bool Node::has_successor(Node* n) const
{
  for (link_map::const_iterator p = succ.begin(); p != succ.end(); p++)
    if (p->second.node == n)
      return true;
  return false;
}

bool Node::has_successor(index_type id) const
{
  for (link_map::const_iterator p = succ.begin(); p != succ.end(); p++)
    if (p->second.node->id == id)
      return true;
  return false;
}

void Node::backup()
{
  // std::cerr << "backup ";
  // write(std::cerr);
  NTYPE c_min = POS_INF;
  for (link_map::const_iterator p = succ.begin(); p != succ.end(); p++) {
    NTYPE c_p = p->second.delta + p->second.node->est;
    // std::cerr << " " << p->second.node->id << " : "
    // 	      << p->second.delta << " + " << p->second.node->est
    // 	      << std::endl;
    c_min = MIN(c_min, c_p);
  }
  if (c_min > est) {
    est = c_min;
    // std::cerr << "update: new est = " << est
    // 	      << ", h = " << state->est_cost()
    // 	      << std::endl;
    if (bp_pre) {
      bp_pre->backup();
    }
  }
}

void Node::cache_optimal_path(NTYPE c)
{
  if (c < opt) {
    opt = c;
    est = opt;
  }
  if (bp_pre) {
    bp_pre->cache_optimal_path(c + bp_delta);
  }
}

bool Node::back_path_contains(Node* n) const
{
  if (n == this) return true;
  if (bp_pre) {
    return bp_pre->back_path_contains(n);
  }
  return false;
}

void Node::add_predecessor(Node* n, NTYPE d)
{
  if (all_pre == 0) all_pre = new link_vec;
  for (index_type k = 0; k < all_pre->length(); k++)
    if ((*all_pre)[k].node == n) {
      (*all_pre)[k].delta = MIN((*all_pre)[k].delta, d);
      return;
    }
  all_pre->append(Link(n, d));
}

void Node::write(std::ostream& s, const Name* p)
{
  assert(state);
  if (p) s << p << " ";
  s << "NODE: " << id
    << ' ' << state->is_max()
    << ' ' << state->is_final()
    << ' ' << acc
    << ' ' << opt
    << ' ' << est
    << ' ' << val
    << ' ' << exp
    << ' ' << succ.size();
  state->write_eval(s, " H:", false);
  s << " STATE: " << *state
    << std::endl;
}

void Node::write_short(std::ostream& s, const Name* p)
{
  assert(state);
  if (p) s << p << " ";
  s << "NODE: " << id
    << ' ' << state->is_max()
    << ' ' << state->is_final()
    << ' ' << acc
    << ' ' << opt
    << ' ' << est
    << ' ' << val
    << ' ' << exp
    << ' ' << succ.size();
  state->write_eval(s, " H:", true);
}

void Node::write_back_path(std::ostream& s)
{
  if (bp_pre) {
    s << "N" << id << " -[";
    if (bp_trans) {
      bp_trans->write(s);
    }
    else {
      s << "?";
    }
    s << " : " << bp_delta << "]-> N" << bp_pre->id;
    if (bp_trans) {
      if (bp_trans->predecessor() != bp_pre->bp_trans) {
	s << "  BROKEN CHAIN!" << std::endl;
      }
    }
    s << std::endl;
    bp_pre->write_back_path(s);
  }
}

void Node::write_graph_node(std::ostream& s)
{
  assert(state);
  s << "N" << id << " [label=\"" << id << ": ";
  if (NodeSet::write_state_in_graph_node) {
    s << *state << " / ";
  }
  s << acc << " + " << est << " / " << opt << "\"";
  if (state->is_max()) {
    s << ",shape=box";
  }
  else {
    s << ",shape=ellipse";
  }
  if (state->is_final()) {
    s << ",peripheries=2";
  }
  else if (exp == 0) {
    s << ",style=dashed";
  }
  s << "];" << std::endl;
}

void Node::write_graph_edges(std::ostream& s)
{
  for (link_map::const_iterator i = succ.begin(); i != succ.end(); i++) {
    s << "N" << id << " -> " << " N" << i->second.node->id << " [label=\"";
    if (i->first != 0)
      i->first->write(s);
    else
      s << "?";
    s << ":" << i->second.delta << "\"";
    if (FINITE(opt) && ((i->second.node->opt + i->second.delta) == opt)) {
      s << ",style=\"bold\"";
    }
    s << "];" << std::endl;
  }
}

bool DefaultNodeOrder::operator()(const nodep& v0, const nodep& v1) const
{
  if (v0->val < v1->val) return true;
  if ((v0->val == v1->val) && (v0->est < v1->est)) return true;
  return false;
}

DefaultNodeOrder NodeQueue::default_node_order;

NodeQueue::NodeQueue(const node_vec::order& b)
  : node_vec(0, 0), before(b)
{
  // done
}

NodeQueue::~NodeQueue()
{
  // done
}

#define ROOT(i)        (i == 0)
#define PARENT(i)      (((i+1) >> 1) - 1)
#define LEFT_CHILD(i)  (((i+1)*2) - 1)
#define RIGHT_CHILD(i) ((i+1)*2)

void NodeQueue::check_queue()
{
  for (index_type i = 0; i < (length() / 2); i++) {
    if (LEFT_CHILD(i) < length()) {
      if (before((*this)[LEFT_CHILD(i)], (*this)[i])) {
	std::cerr << "error in queue:";
	for (index_type k = 0; k < length(); k++)
	  std::cerr << " #" << (*this)[k]->id << ":" << (*this)[k]->val;
	std::cerr << " at position " << i << " / " << LEFT_CHILD(i)
		  << std::endl;
	exit(255);
      }
    }
    if (RIGHT_CHILD(i) < length()) {
      if (before((*this)[RIGHT_CHILD(i)], (*this)[i])) {
	std::cerr << "error in queue:";
	for (index_type k = 0; k < length(); k++)
	  std::cerr << " #" << (*this)[k]->id << ":" << (*this)[k]->val;
	std::cerr << " at position " << i << " / " << RIGHT_CHILD(i)
		  << std::endl;
	exit(255);
      }
    }
  }
}

void NodeQueue::shift_up(index_type i) {
  while (!ROOT(i)) {
    if (before((*this)[i], (*this)[PARENT(i)])) {
      swap(i, PARENT(i));
      (*this)[i]->pos = i;
      (*this)[PARENT(i)]->pos = PARENT(i);
      i = PARENT(i);
    }
    else return;
  }
}

void NodeQueue::shift_down(index_type i) {
  index_type j;
  while (i < (length() / 2)) {
    assert(LEFT_CHILD(i) < length());
    if (LEFT_CHILD(i) == (length() - 1)) {
      j = LEFT_CHILD(i);
    }
    else {
      assert(RIGHT_CHILD(i) < length());
      if (before((*this)[LEFT_CHILD(i)], (*this)[RIGHT_CHILD(i)]))
	j = LEFT_CHILD(i);
      else
	j = RIGHT_CHILD(i);
    }
    if (before((*this)[j], (*this)[i])) {
      swap(i, j);
      (*this)[i]->pos = i;
      (*this)[j]->pos = j;
      i = j;
    }
    else return;
  }
}

void NodeQueue::enqueue(Node* n)
{
  append(n);
  n->pos = length() - 1;
  shift_up(length() - 1);
#ifdef CHECK_QUEUE_OPERATIONS
  std::cerr << "checking queue after ENQUEUE..." << std::endl;
  check_queue();
#endif
}

Node* NodeQueue::dequeue()
{
  if (length() == 0) {
    return 0;
  }
  Node* n = (*this)[0];
  n->pos = no_such_index;
  if (length() > 1) {
    (*this)[0] = (*this)[length() - 1];
    (*this)[0]->pos = 0;
    dec_length();
    shift_down(0);
  }
  else {
    clear();
  }
#ifdef CHECK_QUEUE_OPERATIONS
  std::cerr << "checking queue after DEQUEUE..." << std::endl;
  check_queue();
#endif
  return n;
}

Node* NodeQueue::peek()
{
#ifdef CHECK_QUEUE_OPERATIONS
  std::cerr << "checking queue before PEEK..." << std::endl;
  check_queue();
#endif
  if (length() == 0) return 0;
  return (*this)[0];
}

void NodeQueue::clear()
{
  for (index_type i = 0; i < size(); i++)
    (*this)[i]->pos = no_such_index;
  node_vec::clear();
}


index_type NodeSet::next_id = 1;
bool NodeSet::write_state_in_graph_node = true;

NodeSet::~NodeSet()
{
  // can't call clear because no implement?
}

void NodeSet::back_path_to_sequence(Node* n, node_vec& ns)
{
  if (n->bp_pre) {
    back_path_to_sequence(n->bp_pre, ns);
  }
  ns.append(n);
}

void NodeSet::set_back_path_solution_cost(Node* n, NTYPE c_sol)
{
  n->opt = c_sol;
  if (n->bp_pre) {
    set_back_path_solution_cost(n->bp_pre, c_sol + n->bp_delta);
  }
}

Node* NodeSet::insert_root_node(State& s)
{
  Node* n = insert_node(s);
  make_root(n);
  return n;
}

void NodeSet::make_root(Node* n)
{
  if (roots.first(n) == no_such_index)
    roots.append(n);
}

node_vec& NodeSet::root_nodes()
{
  return roots;
}

void NodeSet::cache_pg(NTYPE c_sol)
{
  node_vec v(0, 0);
  collect_nodes(v);

  for (index_type k = 0; k < v.length(); k++) if (v[k]->closed) {
    if (FINITE(v[k]->opt)) {
      v[k]->state->store(v[k]->opt, true);
    }
    else {
      v[k]->state->store(c_sol - v[k]->acc, false);
    }
  }
}

void NodeSet::mark_solved_apsp()
{
  node_vec v(0, 0);
  collect_nodes(v);
  index_type n = v.length();
  NTYPE** d = new NTYPE*[n];
  for (index_type i = 0; i < n; i++)
    d[i] = new NTYPE[n];

  std::cerr << "initializing..." << std::endl;
  for (index_type i = 0; i < n; i++)
    for (index_type j = 0; j < n; i++)
      d[i][j] = POS_INF;
  std::cerr << "adding links..." << std::endl;
  for (index_type k = 0; k < v.length(); k++)
    for (link_map::const_iterator i = v[k]->succ.begin();
	 i != v[k]->succ.end(); i++) {
      index_type p = v.first(i->second.node);
      assert(p != no_such_index);
      d[k][p] = MIN(d[k][p], i->second.delta);
    }
  std::cerr << "computing shortest paths..." << std::endl;
  for (index_type k = 0; k < n; k++)
    for (index_type i = 0; i < n; i++)
      for (index_type j = 0; j < n; i++)
	d[i][j] = MIN(d[i][j], d[i][k] + d[k][j]);
      
  std::cerr << "initializing..." << std::endl;
  for (index_type k = 0; k < v.length(); k++)
    v[k]->opt = POS_INF;
  std::cerr << "computing min distances..." << std::endl;
  for (index_type k = 0; k < v.length(); k++)
    if (v[k]->state->is_final()) {
      v[k]->opt = 0;
      for (index_type i = 0; i < v.length(); i++)
	v[i]->opt = MIN(v[i]->opt, d[i][k]);
    }
}

class NodeOrderByOpt : public node_vec::order {
public:
  NodeOrderByOpt() { };
  virtual bool operator()(const nodep& v0, const nodep& v1) const;
};

bool NodeOrderByOpt::operator()(const nodep& v0, const nodep& v1) const
{
  if (v0->opt < v1->opt) return true;
  return false;
}

void NodeSet::compute_reverse_links(node_vec& v)
{
  for (index_type k = 0; k < v.length(); k++)
    for (link_map::iterator i = v[k]->succ.begin(); i != v[k]->succ.end(); i++)
      i->second.node->add_predecessor(v[k], i->second.delta);
}

void NodeSet::compute_reverse_links()
{
  node_vec v(0, 0);
  collect_nodes(v);
  compute_reverse_links(v);
}

void NodeSet::mark_solved()
{
  node_vec v(0, 0);
  // std::cerr << "collecting nodes..." << std::endl;
  collect_nodes(v);
  // std::cerr << "computing reverse links..." << std::endl;
  compute_reverse_links(v);

  NodeOrderByOpt opt_is_less;
  NodeQueue q(opt_is_less);

  // std::cerr << "initializing..." << std::endl;
  for (index_type k = 0; k < v.length(); k++) {
    if (v[k]->state->is_final()) {
      v[k]->opt = 0;
      q.enqueue(v[k]);
    }
    else {
      v[k]->opt = POS_INF;
    }
  }

  index_type nn = v.length();
  index_type i = 0;
  // index_type l = 0;

  // std::cerr << "marking..." << std::endl;
  while (!q.empty() > 0) {
    Node* n = q.dequeue();
    if (n->all_pre) {
      for (index_type k = 0; k < n->all_pre->length(); k++) {
	NTYPE d = (*(n->all_pre))[k].delta;
	if ((n->opt + d) < (*(n->all_pre))[k].node->opt) {
	  (*(n->all_pre))[k].node->opt = n->opt + d;
	  q.enqueue((*(n->all_pre))[k].node);
	}
      }
    }
    i += 1;
    // if ((i - l) > (nn / 20)) {
    // std::cerr << rational(i, nn).decimal()*100.0 << "%..." << std::endl;
    // l = i;
    // }
  }
}

void NodeSet::write_short(std::ostream& s, const Name* p)
{
  node_vec v(0, 0);
  collect_nodes(v);
  for (index_type k = 0; k < v.length(); k++) {
    v[k]->write_short(s, p);
  }
}

void NodeSet::write_graph(std::ostream& s)
{
  node_vec v(0, 0);
  collect_nodes(v);

  s << "digraph BFS_SEARCH_SPACE {" << std::endl;
  s << "node [width=0,height=0,shape=box];" << std::endl;
  for (index_type k = 0; k < v.length(); k++) {
    v[k]->write_graph_node(s);
  }
  for (index_type k = 0; k < v.length(); k++) {
    v[k]->write_graph_edges(s);
  }
  s << "}" << std::endl;
}

void NodeSet::write_graph_compact(std::ostream& s)
{
  node_vec v(0, 0);
  collect_nodes(v);

  s << "digraph BFS_SEARCH_SPACE {" << std::endl;
  s << "node [width=0.2,height=0.2,shape=circle];" << std::endl;
  for (index_type k = 0; k < v.length(); k++) {
    s << "N" << v[k]->id << " [label=\"\"";
    if (roots.contains(v[k])) {
      s << ",style=filled,fillcolor=black";
    }
    if (v[k]->state->is_final()) {
      s << ",peripheries=2";
    }
    else if (INFINITE(v[k]->opt)) {
      s << ",style=dashed";
    }
    s << "];" << std::endl;
  }
  for (index_type k = 0; k < v.length(); k++)
    for (link_map::const_iterator i = v[k]->succ.begin();
	 i != v[k]->succ.end(); i++) {
      if (i->second.node->has_successor(v[k])) {
	if (v[k]->id < i->second.node->id)
	  s << "N" << v[k]->id << " -> N" << i->second.node->id
	    << " [dir=both];" << std::endl;
      }
      else {
	s << "N" << v[k]->id << " -> N" << i->second.node->id
	  << ";" << std::endl;
      }
    }
  s << "}" << std::endl;
}

void NodeSet::write_graph_rainbow(std::ostream& s)
{
  node_vec v(0, 0);
  collect_nodes(v);

  NTYPE e_min = POS_INF;
  NTYPE e_max = NEG_INF;
  for (index_type k = 0; k < v.length(); k++)
    if (FINITE(v[k]->est)) {
      e_min = MIN(e_min, v[k]->est);
      e_max = MAX(e_max, v[k]->est);
    }
  double e_range = N_TO_D(e_max - e_min);

  std::cerr << "e_max = " << e_max
	    << ", e_min = " << e_min
	    << ", e_range = " << e_range
	    << std::endl;

  s << "digraph BFS_SEARCH_SPACE {" << std::endl;
  s << "node [width=0.2,height=0.2];" << std::endl;
  for (index_type k = 0; k < v.length(); k++) {
    s << "N" << v[k]->id << " [label=\"\"";
    if (roots.contains(v[k])) {
      s << ",shape=box";
    }
    else if (INFINITE(v[k]->opt)) {
      s << ",shape=diamond";
    }
    else {
      s << ",shape=circle";
    }
    if (v[k]->state->is_final()) {
      s << ",peripheries=2";
    }
    if (INFINITE(v[k]->est)) {
      s << ",style=filled,fillcolor=black];" << std::endl;
    }
    else {
      double e_val = N_TO_D(v[k]->est - e_min);
      //std::cerr << "e_val = " << e_val << std::endl;
      double h = (300 - ((e_val/e_range)*300)) + 180;
      if (h >= 360) h = (h - 360);
      //std::cerr << "h = " << h << std::endl;
      double hp = h/60;
      double tmp = hp - floor(hp / 2);
      if (tmp > 2) tmp = 2;
      //std::cerr << "hp = " << hp << ", tmp = " << tmp << std::endl;
      double x = 1 - ((tmp - 1) < 0 ? -1 * (tmp - 1) : (tmp - 1));
      //std::cerr << "x = " << x << std::endl;
      double r, g, b;
      if (hp < 1) {
	r = 1; g = x; b = 0;
      }
      else if (hp < 2) {
	r = x; g = 1; b = 0;
      }
      else if (hp < 3) {
	r = 0; g = 1; b = x;
      }
      else if (hp < 4) {
	r = 0; g = 1; b = x;
      }
      else if (hp < 5) {
	r = x; g = 0; b = 1;
      }
      else {
	r = 1; g = 0; b = x;
      }
      s << ",style=filled,fillcolor=\"" << h/360 << ", 1, 1\"];" << std::endl;
    }
  }
  for (index_type k = 0; k < v.length(); k++)
    for (link_map::const_iterator i = v[k]->succ.begin();
	 i != v[k]->succ.end(); i++) {
      if (i->second.node->has_successor(v[k])) {
	if ((v[k]->acc < i->second.node->acc) ||
	    ((v[k]->acc == i->second.node->acc) &&
	     (v[k]->id < i->second.node->id)))
	  s << "N" << v[k]->id << " -> N" << i->second.node->id
	    << " [dir=both];" << std::endl;
      }
      else {
	s << "N" << v[k]->id << " -> N" << i->second.node->id
	  << ";" << std::endl;
      }
    }
  s << "}" << std::endl;
}

TreeNodeSet::TreeNodeSet()
  : left(0), right(0)
{
  id = next_id++;
}

TreeNodeSet::~TreeNodeSet()
{
  clear();
}

Node* TreeNodeSet::insert_node(State& s)
{
  if (!state) return this;
  int d = state->compare(s);
  if (d < 0) {
    if (!right) right = new TreeNodeSet();
    return right->insert_node(s);
  }
  else if (d > 0) {
    if (!left) left = new TreeNodeSet();
    return left->insert_node(s);
  }
  else return this;
}

Node* TreeNodeSet::find_node(State& s)
{
  if (!state) return 0;
  int d = state->compare(s);
  if (d < 0) {
    if (!right) return 0;
    else return right->find_node(s);
  }
  else if (d > 0) {
    if (!left) return 0;
    else return left->find_node(s);
  }
  else return this;
}

void TreeNodeSet::clear()
{
  if (state) {
    delete state;
    state = NULL;
  }
  if (bp_trans) {
    delete bp_trans;
    bp_trans = NULL;
  }
  if (left) {
    delete left;
    left = 0;
  }
  if (right) {
    delete right;
    right = 0;
  }
}

void TreeNodeSet::collect_nodes(node_vec& ns)
{
  ns.append(this);
  if (left) {
    left->collect_nodes(ns);
  }
  if (right) {
    right->collect_nodes(ns);
  }
}

HashNodeSet::HashNodeSet(index_type s)
  : size(s), tab(0)
{
  tab = new TreeNodeSet*[size];
  for (index_type k = 0; k < size; k++) tab[k] = 0;
}

HashNodeSet::~HashNodeSet() {
  clear();
  delete tab;
}

Node* HashNodeSet::insert_node(State& s)
{
  index_type i = (s.hash() % size);
  if (tab[i]) {
    return tab[i]->insert_node(s);
  }
  else {
    tab[i] = new TreeNodeSet();
    return tab[i];
  }
}

Node* HashNodeSet::find_node(State& s)
{
  index_type i = (s.hash() % size);
  if (tab[i]) {
    return tab[i]->find_node(s);
  }
  else {
    return 0;
  }
}

void HashNodeSet::clear()
{
  for (index_type k = 0; k < size; k++) if (tab[k]) {
    delete tab[k];
    tab[k] = 0;
  }
}

void HashNodeSet::collect_nodes(node_vec& ns)
{
  for (index_type k = 0; k < size; k++) {
    if (tab[k]) tab[k]->collect_nodes(ns);
  }
}

END_HSPS_NAMESPACE
