
#include "bss.h"

BEGIN_HSPS_NAMESPACE

//#define TRACE_PRINT_LOTS

BeamStackSearch::BeamStackSearch
(Statistics& s, SearchResult& r, index_type w)
  : SearchAlgorithm(s, r), beam_width(w)
{
  // done
}

BeamStackSearch::~BeamStackSearch()
{
}

NTYPE BeamStackSearch::start(State& s, NTYPE b)
{
  set_cost_limit(b);
  start(s);
}

NTYPE BeamStackSearch::start(State& s)
{
  stats.start();
  NTYPE val = main(s);
  stats.stop();
  return val;
}

void BeamStackSearch::make_new_layer(index_type l)
{
  layers.inc_length_to(l + 1);
  assert(layers[l].nodes == NULL);
  layers[l].nodes = new HashNodeSet((3 * beam_width) - 1);
  layers[l].open.set_length(0);
  layers[l].next = 0;
  layers[l].L = signature(ZERO);
  layers[l].U = signature(get_cost_limit());
}

void BeamStackSearch::reset_layer(index_type l)
{
  assert(layers.size() >= l + 1);
  assert(layers[l].nodes == NULL);
  layers[l].nodes = new HashNodeSet((3 * beam_width) - 1);
  layers[l].open.set_length(0);
  layers[l].next = 0;
}

void BeamStackSearch::delete_layer(index_type l)
{
  assert(l < layers.size());
  assert(layers[l].nodes != NULL);
  delete layers[l].nodes;
  layers[l].nodes = NULL;
  layers[l].open.clear();
  layers[l].next = 0;
}

void BeamStackSearch::delete_all_layers()
{
  for (index_type i = 0; i < layers.size(); i++)
    if (layers[i].nodes != NULL)
      delete_layer(i);
}

NTYPE BeamStackSearch::main(State& s)
{
  // initialise the beam stack
  make_new_layer(0);
  Node* root = layers[0].nodes->insert_node(s);
  root->state = s.copy();
  root->acc = ZERO;
  root->est = s.est_cost();
  root->val = root->est;
  root->bp_pre = NULL;
  root->bp_trans = NULL;
  layers[0].next = 0;
  layers[0].open.append(root);
  current = 0;
  make_new_layer(1);
  // main loop
  while (true) { // exit is through return in backtracking loop
    // expand nodes in current layers
    std::cerr << "expanding layer " << current
	      << " : " << layers[current].next << " / "
	      << layers[current].open.size()
	      << ", " << layers[current].min_open()
	      << "--" << layers[current].max_open()
	      << " into layer " << current + 1
	      << " " << layers[current + 1].L << "--" << layers[current + 1].U
	      << std::endl;
    while (!layers[current].empty()) {
      // "pop" next open node in current layer
      current_node = layers[current].open[layers[current].next];
      layers[current].next += 1;
      // if this is a goal state, it must be a new best plan:
      if (current_node->state->is_final()) {
	// ignore goal state if cost >= limit
	if (current_node->acc < get_cost_limit()) {
	  std::cerr << "found plan with cost "
		    << current_node->acc
		    << std::endl;
	  set_solved(true, false); // not optimal!
	  set_cost_limit(current_node->acc);
	  result.solution(*(current_node->state), current_node->acc);
	  if (done()) {
	    std::cerr << "exiting because done() == true" << std::endl;
	    return current_node->acc;
	  }
	}
      }
      // if not a goal state, generate the admitted successors
      // (all pruning is done in new_state):
      else {
#ifdef TRACE_PRINT_LOTS
	std::cerr << "expanding " << *(current_node->state)
	 	  << " (" << (void*)current_node << ")" << std::endl;
#endif
	stats.expand_node(*(current_node->state));
	current_node->state->expand(*this, layers[current + 1].bound());
      }
      if (stats.break_signal_raised()) return POS_INF;
    }
    std::cerr << "complete: layer " << current + 1
	      << " : " << layers[current + 1].next << " / "
	      << layers[current + 1].open.size()
	      << ", " << layers[current + 1].min_open()
	      << "--" << layers[current + 1].max_open()
	      << ", " << layers[current + 1].L
	      << "--" << layers[current + 1].U
	      << std::endl;
    std::cerr << stats << std::endl;
    if (!layers[current + 1].empty()) {
      current += 1;
      make_new_layer(current + 1);
    }
    else {
      std::cerr << "backtracking..." << std::endl;
      // current layer expanded but generated no accepted successors
      // in the [L,U] range of next layer.
      layers[current + 1].L = layers[current + 1].U.clone();
      layers[current + 1].U = signature(get_cost_limit());
      delete_layer(current + 1);
      std::cerr << "updated layer " << current + 1 << " to [ "
		<< layers[current + 1].L << " , " << layers[current + 1].U
		<< " ]" << std::endl;
      // next layer is exhausted: we must backtrack
      while (!(layers[current + 1].L < layers[current + 1].U)) {
	// if current layer is zero, we cannot backtrack any more
	if (current == 0) {
	  result.no_more_solutions(get_cost_limit());
	  return POS_INF;
	}
	current -= 1;
	layers[current + 1].L = layers[current + 1].U.clone();
	layers[current + 1].U = signature(get_cost_limit());
	delete_layer(current + 1);
	std::cerr << "updated layer " << current + 1 << " to [ "
		  << layers[current + 1].L << " , " << layers[current + 1].U
		  << " ]" << std::endl;
      }
      // done backtracking: the next layer now has updated bounds, so
      // we need to re-expand the open list in the (now) current layer.
      layers[current].next = 0;
      // std::cerr << "is " << layers[current + 1].L
      // 		<< " < " << layers[current + 1].U
      // 		<< " ? "
      // 		<< (layers[current + 1].L < layers[current + 1].U)
      // 		<< std::endl;
      assert(layers[current + 1].L < layers[current + 1].U);
      reset_layer(current + 1);
    }
  }
}

NTYPE BeamStackSearch::main2(State& s)
{
  NTYPE cnew = main(s);
  std::cerr << "cnew = " << cnew << ", done = " << done() << std::endl;
  while (!done() && (cnew < POS_INF)) {
    delete_all_layers();
    set_cost_limit(cnew);
    cnew = main(s);
    std::cerr << "cnew = " << cnew << ", done = " << done() << std::endl;
  }
}

index_type BeamStackSearch::insert_pos
(const signature& nss, node_vec& open, index_type last)
{
  if (last == no_such_index) {
    if (open.size() == 0) return 0;
    last = open.size();
  }
  assert(last <= open.size());
  while (last > 0) {
    signature slast(*(open[last - 1]));
    if (slast < nss)
      return last;
    last -= 1;
  }
  return 0;
}

NTYPE BeamStackSearch::new_state(State& s, NTYPE bound)
{
#ifdef TRACE_PRINT_LOTS
  std::cerr << "generated " << s << std::endl;
#endif
  NTYPE s_acc = s.acc_cost();
  NTYPE s_est = s.est_cost();
  NTYPE s_val = s_acc + s_est;
  if (s_val >= get_cost_limit()) {
    return POS_INF;
  }
  // is the state already in the stack with a cheaper or equal path?
  for (index_type i = 0; i <= current + 1; i++) {
    Node* n = layers[i].nodes->find_node(s);
    if (n != NULL)
      if ((n->acc <= s_acc) && (n->pos != no_such_index)) {
#ifdef TRACE_PRINT_LOTS
	std::cerr << "- discarded by loop check" << std::endl;
#endif
	return POS_INF;
      }
  }
  Transition* s_trans = s.transition();
  signature nss(s_val, s_est, current_node, s_trans);
#ifdef TRACE_PRINT_LOTS
  std::cerr << "- signature: " << nss << std::endl;
#endif
  // is the state below the next layer lower limit?
  if (nss < layers[current + 1].L) {
    delete s_trans;
    return POS_INF;
  }
  // is the state above the next layer upper limit?
  if (layers[current + 1].U < nss) {
#ifdef TRACE_PRINT_LOTS
    std::cerr << "- discarded as > U = " << layers[current + 1].U << std::endl;
#endif
    delete s_trans;
    return POS_INF;
  }
  // otherwise, insert new node in next layer open:
  Node* new_node = layers[current + 1].nodes->insert_node(s);
  // is this a cheaper path to a state already in the next layer?
  if (new_node->state != NULL) {
    if (s_acc >= new_node->acc) {
      assert(new_node->pos == no_such_index);
      return POS_INF;
    }
#ifdef TRACE_PRINT_LOTS
    if (s_acc > new_node->acc) {
      std::cerr << "identical state: " << s << " with acc = " << s_acc
		<< " reached by path:" << std::endl;
      Transition* p = s.transition_path();
      p->write_path(std::cerr);
      std::cerr << "previous state: " << *(new_node->state)
		<< " with acc = " << new_node->state->acc_cost()
		<< " / " << new_node->acc << " reached by path:" << std::endl;
      Transition* q = new_node->state->transition_path();
      q->write_path(std::cerr);
      std::cerr << "pos = " << new_node->pos << std::endl;
      std::cerr << "pos in open = " << layers[current + 1].open.first(new_node)
		<< std::endl;
    }
#endif
    delete new_node->state;
    if (new_node->bp_trans != NULL)
      delete new_node->bp_trans;
  }
  stats.create_node(s);
  new_node->state = s.copy();
  new_node->acc = s_acc;
  new_node->est = s_est;
  new_node->val = s_val;
  new_node->bp_pre = current_node;
  new_node->bp_trans = s_trans;
  index_type pos_in_open = layers[current + 1].open.first(new_node);
  // is the "new" node already in open?
  if (pos_in_open != no_such_index) {
#ifdef TRACE_PRINT_LOTS
    std::cerr << "- found at " << pos_in_open << std::endl;
#endif
    index_type new_pos = insert_pos(nss, layers[current + 1].open, pos_in_open);
    assert(new_pos <= pos_in_open);
#ifdef TRACE_PRINT_LOTS
    std::cerr << "- move to " << new_pos << std::endl;
#endif
    for (index_type k = pos_in_open; k > new_pos; k--)
      layers[current + 1].open[k] = layers[current + 1].open[k-1];
    layers[current + 1].open[new_pos] = new_node;
    new_node->pos = 1;
  }
  else {
    index_type new_pos =
      insert_pos(nss, layers[current + 1].open, no_such_index);
    assert(0 <= new_pos && new_pos < layers[current + 1].open.size() + 1);
#ifdef TRACE_PRINT_LOTS
    std::cerr << "- insert at " << new_pos << std::endl;
#endif
    layers[current + 1].open.insert(new_node, new_pos);
    new_node->pos = 1;
  }
  // for (index_type k = 1; k < layers[current + 1].open.size(); k++) {
  //   signature tmp1(*(layers[current + 1].open[k - 1]));
  //   signature tmp2(*(layers[current + 1].open[k]));
  //   if (!(tmp1 < tmp2)) {
  //     for (index_type i = 0; i < layers[current + 1].open.size(); i++) {
  // 	signature tmp(*(layers[current + 1].open[i]));
  // 	std::cerr << " " << i << " : " << tmp << std::endl;
  //     }
  //   }
  //   assert(tmp1 < tmp2);
  // }
  if (layers[current + 1].open.size() > beam_width) {
#ifdef TRACE_PRINT_LOTS
    std::cerr << "pruning: |open| = " << layers[current + 1].open.size()
	      << std::endl;
#endif
    for (index_type k = beam_width; k < layers[current + 1].open.size(); k++) {
      layers[current + 1].open[k]->pos = no_such_index;
    }
    signature newU(*(layers[current + 1].open[beam_width]));
#ifdef TRACE_PRINT_LOTS
    std::cerr << " newU = " << newU << std::endl;
#endif
    assert(newU < layers[current + 1].U);
    layers[current + 1].U = newU;
    layers[current + 1].open.set_length(beam_width);
  }
  return POS_INF;
}

NTYPE BeamStackSearch::cost() const
{
  if (solved())
    return get_cost_limit();
  else
    return POS_INF;
}

END_HSPS_NAMESPACE
