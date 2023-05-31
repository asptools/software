
#include "bfs.h"
#include "plans.h"

BEGIN_HSPS_NAMESPACE

Instance* BFS::trace_print_instance = 0;

void BFS::update_path(Node* n, Node* p, Transition* t, NTYPE d) {
  NTYPE c_new = p->acc + d;
  if (c_new < n->acc) {
    n->acc = c_new;
    // if (n->state->is_final()) {
    //   if (current_sol) {
    // 	if (n->acc < current_sol->acc) {
    // 	  std::cerr << "best solution: " << PRINT_NTYPE(n->acc)
    // 		    << " (" << stats << ")" << std::endl;
    // 	  current_sol = n;
    // 	}
    //   }
    //   else {
    // 	std::cerr << "best solution: " << PRINT_NTYPE(n->acc)
    // 		  << " (" << stats << ")" << std::endl;
    // 	current_sol = n;
    //   }
    // }
    n->val = c_new + (n->est * weight);
    n->bp_pre = p;
    n->bp_trans = t;
    assert(t);
    t->set_predecessor(p->bp_trans);
    n->bp_delta = d;
    if (n->pos != no_such_index) {
      queue.shift_up(n->pos);
    }
    for (link_map::iterator i = n->succ.begin(); i != n->succ.end(); i++)
      update_path(i->second.node, n, i->first, i->second.delta);
  }
}

NTYPE BFS::new_state(State& s, NTYPE bound) {
  NTYPE s_est = s.est_cost();
  if (INFINITE(s_est)) return POS_INF;

  acc_succ += 1;

  NTYPE s_delta = s.delta_cost();
  Transition* s_trans = s.transition();

  Node* n = graph.insert_node(s);
  if (n->state) {
    if (greedy) return POS_INF;
    assert(current_node);
    assert(s_trans);
    link_map::iterator i = current_node->succ.find(s_trans);

    // if we have found a cheaper way to the new state..
    if (n->acc > (current_node->acc + s_delta)) {
      // the successor links cannot already exist
      if (i != current_node->succ.end()) {
	std::cerr << "assert fail:" << std::endl;
	std::cerr << i->first << " / ";
	if (i->first) i->first->write(std::cerr);
	std::cerr << std::endl;
	std::cerr << s_trans << " / ";
	s_trans->write(std::cerr);
	std::cerr << std::endl;
	int c1 = s_trans->compare(*(i->first));
	int c2 = i->first->compare(*s_trans);
	std::cerr << c1 << " / "  << c2 << std::endl;
      }
      assert(i == current_node->succ.end());
      if (trace_level > 3) {
	std::cerr << "update #" << n->id << " " << s
		  << ": " << n->acc << "/" << n->est << "/" << n->val
		  << " -> "
		  << current_node->acc + s.delta_cost() << "/"
		  << n->est << "/"
		  << current_node->acc + s.delta_cost() + n->est
		  << ", #" << current_node->id << std::endl;
	Node* p = n->bp_pre;
	while (p) {
	  std::cerr << " #" << p->id
		    << ": " << p->acc << "/" << p->est << "/" << p->val
		    << " <- ";
	  p = p->bp_pre;
	}
	std::cerr << "0" << std::endl;
	p = current_node;
	while (p) {
	  std::cerr << " #" << p->id
		    << ": " << p->acc << "/" << p->est << "/" << p->val
		    << " <- ";
	  p = p->bp_pre;
	}
	std::cerr << "0" << std::endl;
      }
      //delete n->state;
      //n->state = s.copy();
      //n->state->State::set_predecessor(0);
      // update costs and back-pointers
      update_path(n, current_node, s_trans, s_delta);
    }
    // else, if the successor link does not exist in map, we have found
    // a new, but more (or equally) expensive path to the node
    else if (i == current_node->succ.end()) {
      // record new link in current_node's successors
      current_node->succ.insert(link_map::value_type(s_trans, Link(n, s_delta)));
      if (trace_level > 3) {
	std::cerr << "add link #" << current_node->id
		  << " to " << n->id << " ";
	s_trans->write(std::cerr);
	std::cerr << std::endl;
      }
    }
    else {
      delete s_trans;
    }
  }

  else {
    new_succ += 1;
    stats.create_node(s);
    n->state = s.copy();
    n->state->set_predecessor(0);
    if (current_node) {
      assert(s_trans);
      current_node->succ.insert(link_map::value_type(s_trans, Link(n, s_delta)));
      n->bp_pre = current_node;
      n->bp_trans = s_trans;
      s_trans->set_predecessor(current_node->bp_trans);
      n->bp_delta = s_delta;
      n->acc = current_node->acc + s_delta;
    }
    else {
      graph.make_root(n);
      n->bp_pre = 0;
      n->bp_trans = 0;
      n->bp_delta = s_delta;
      if (n->bp_delta > 0) {
	std::cerr << "warning: root node with non-zero delta cost"
		  << std::endl;
      }
      n->acc = s_delta;
    }
    n->est = s_est;
    if (greedy)
      n->val = n->est;
    else
      n->val = n->acc + (n->est * weight);
    n->exp = 0;
    // if (s.is_final()) {
    //   if (current_sol) {
    // 	if (n->acc < current_sol->acc) {
    // 	  std::cerr << "best solution: " << PRINT_NTYPE(n->acc)
    // 		    << " (" << stats << ")" << std::endl;
    // 	  current_sol = n;
    // 	}
    //   }
    //   else {
    // 	std::cerr << "best solution: " << PRINT_NTYPE(n->acc)
    // 		  << " (" << stats << ")" << std::endl;
    // 	current_sol = n;
    //   }
    // }
    queue.enqueue(n);

    if (trace_level > 3) {
      std::cerr << "new #" << n->id << " " << s
		<< ": " << n->acc << "/" << n->est << "/" << n->val;
      if (trace_level > 4) {
	std::cerr << " >> ";
	s.write_path(std::cerr);
      }
      std::cerr << std::endl;
    }
  }

  return POS_INF;
}

BFS::BFS(Statistics& s, SearchResult& r)
  : SingleSearchAlgorithm(s, r),
    graph(31337),
    current_node(0),
    current_sol(0),
    best_node_cost(0),
    weight(1),
    greedy(false)
{
  // done
}

BFS::BFS(Statistics& s, SearchResult& r, index_type nt_size)
  : SingleSearchAlgorithm(s, r),
    graph(nt_size),
    current_node(0),
    current_sol(0),
    best_node_cost(0),
    weight(1),
    greedy(false)
{
  // done
}

BFS::~BFS() {
  /* delete graph; */
}

NTYPE BFS::main() {
  NTYPE last_h = 0;
  bool first_iteration = true;
  while ((!solved() || result.more()) && !queue.empty()) {
    if (stats.break_signal_raised()) return best_node_cost;
    current_node = queue.peek();
    if ((current_node->val > best_node_cost) || first_iteration) {
      if (!solved()) stats.current_lower_bound(current_node->val);
      result.no_more_solutions(best_node_cost);
      if (!solved()) result.lower_bound(current_node->val);
      if (solved() && !result.more()) return best_node_cost;
      best_node_cost = current_node->val;
      if (trace_level > 0) {
	std::cerr << "f = " << PRINT_NTYPE(current_node->val)
		  << ", q = " << queue.length()
		  << ", " << stats << std::endl;
      }
    }
    else if (greedy && (trace_level > 0 ) && (current_node->est != last_h)) {
      std::cerr << "h = " << PRINT_NTYPE(current_node->est)
		<< ", q = " << queue.length()
		<< ", " << stats << std::endl;
      last_h = current_node->est;
    }
    if (best_node_cost > get_cost_limit()) {
      return best_node_cost;
    }
    // if (current_sol) {
    //   if ((current_sol->acc <= current_node->val) || greedy) {
    // 	if (trace_level > 0) {
    // 	  std::cerr << "solution (cost = " << PRINT_NTYPE(current_node->acc)
    // 		    << " (" << current_node->state->acc_cost()
    // 		    << "), depth = " << current_node->state->depth() << ")"
    // 		    << std::endl;
    // 	}
    // 	set_solved(true, !greedy);
    // 	result.solution(*(current_sol->state),
    // 			current_sol->bp_trans,
    // 			current_sol->acc);
    //   }
    // }
    current_node = queue.dequeue();
    assert(!current_node->state->is_max());
    if (current_node->state->is_final()) {
      if (trace_level > 0) {
	std::cerr << "solution (cost = " << PRINT_NTYPE(current_node->acc)
		  << " (" << current_node->state->acc_cost()
		  << "), depth = " << current_node->state->depth() << ")"
		  << std::endl;
      }
      set_solved(true, !greedy);
      //graph.set_back_path_solution_cost(current_node, current_node->acc);
      //State* f_state = graph.build_path(current_node);
      result.solution(*(current_node->state),
		      current_node->bp_trans,
		      current_node->acc);
    }
    if (!done()) {
      assert(current_node->exp == 0); // node should only ever be expanded once
      stats.expand_node(*(current_node->state));
      if (trace_level > 2) {
	std::cerr << "expanding #" << current_node->id
		  << " " << *(current_node->state)
		  << ": " << current_node->acc
		  << "/" << current_node->est
		  << "/" << current_node->val;
	if (current_node->bp_pre) {
	  std::cerr << ", pre = #" << current_node->bp_pre->id << " by ";
	  if (trace_print_instance) {
	    PrintActions tprinter(*trace_print_instance, std::cerr);
	    current_node->bp_trans->insert(tprinter);
	  }
	  else {
	    current_node->bp_trans->write(std::cerr);
	  }
	}
	std::cerr << std::endl;
      }
      else if (trace_level > 1) {
	if ((stats.nodes() % TRACE_LEVEL_2_NOTIFY) == 0)
	  std::cerr << stats << std::endl;
      }
      acc_succ = 0;
      new_succ = 0;
      index_type l = queue.length();
      current_node->state->expand(*this, POS_INF);
      current_node->exp += 1;
      if (!done()) current_node->closed = true;
      assert(new_succ == (queue.length() - l));
      if (trace_level > 3) {
	std::cerr << "done expanding #" << current_node->id
		  << " (" << acc_succ << " acc. / " << new_succ
		  << " new succs.), closed = " << current_node->closed
		  << std::endl;
      }
      // After a node has been fully expanded, we can back up new
      // estimated costs from its successors (possibly updating the
      // node's estimated cost, and therefore also that of its parent).
      // At present, this is completely redundant because the updated
      // cost estimates are not used for anything. This is a
      // place-holder for possible future development.
      //current_node->backup();
    }
    first_iteration = false;
    current_node = 0;
  }
  if (queue.empty()) {
    result.no_more_solutions(POS_INF);
    if (!solved()) best_node_cost = POS_INF;
  }
  return best_node_cost;
}

NTYPE BFS::start(State& s, NTYPE b)
{
  return start(s);
}

NTYPE BFS::start(State& s)
{
  reset();
  graph.clear();
  queue.set_length(0);
  best_node_cost = 0;

  stats.start();
  current_node = 0;
  current_sol = 0;
  new_state(s, POS_INF);

  NTYPE val = main();
  stats.stop();
  return val;
}

NTYPE BFS::resume() {
  stats.start();
  NTYPE val = main();
  stats.stop();
  return val;
}

NTYPE BFS::cost() const {
  return best_node_cost;
}

bool BFS::done() const {
  return ((solved() && !result.more()) || stats.break_signal_raised());
}

NTYPE BFS_PX::main() {
  while ((!solved() || result.more()) && !queue.empty()) {
    if (stats.break_signal_raised()) return best_node_cost;
    current_node = queue.peek();
    if (current_node->val > best_node_cost) {
      if (!solved()) stats.current_lower_bound(current_node->val);
      result.no_more_solutions(best_node_cost);
      if (!solved()) result.lower_bound(current_node->val);
      if (solved() && !result.more()) return best_node_cost;
      best_node_cost = current_node->val;
      if (trace_level > 0) {
	std::cerr << "f = " << current_node->val
		  << ", q = " << queue.length()
		  << ", " << stats << std::endl;
      }
    }
    if (best_node_cost > get_cost_limit()) {
      return best_node_cost;
    }
    current_node = queue.dequeue();
    assert(!current_node->state->is_max());
    if (current_node->state->is_final()) {
      set_solved(true, !greedy);
      result.solution(*(current_node->state),
		      current_node->bp_trans,
		      current_node->acc);
    }
    else {
      stats.expand_node(*(current_node->state));
      if (trace_level > 2) {
	if (current_node->exp == 0)
	  std::cerr << "expanding ";
	else
	  std::cerr << "re-expanding ";
	std::cerr << *(current_node->state) << ": "
		  << current_node->acc << "/"
		  << current_node->est << "/"
		  << current_node->val
		  << " with bound " << current_node->val + threshold
		  << std::endl;
	acc_succ = 0;
	new_succ = 0;
      }
      else if (trace_level > 1) {
	if ((stats.nodes() % TRACE_LEVEL_2_NOTIFY) == 0)
	  std::cerr << stats << std::endl;
      }
      NTYPE new_val =
	current_node->state->expand(*this, current_node->val + threshold);
      if (trace_level > 2) {
	std::cerr << "new value is " << new_val << ", " << acc_succ
		  << " successors accepted, " << new_succ << " are new"
		  << std::endl;
      }
      if (FINITE(new_val)) {
	current_node->val = new_val;
	queue.enqueue(current_node);
      }
      current_node->exp += 1;
    }
    current_node = 0;
  }
  if (queue.empty()) {
    result.no_more_solutions(POS_INF);
    if (!solved()) best_node_cost = POS_INF;
  }
  return best_node_cost;
}

END_HSPS_NAMESPACE
