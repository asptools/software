
#include "bfhs.h"

BEGIN_HSPS_NAMESPACE

NTYPE BFHS::new_state(State& s, NTYPE bound) {
  if (trace_level > 3) {
    std::cerr << "new " << s << ": " << s.delta_cost() << "/" << s.est_cost();
    if (trace_level > 4) {
      std::cerr << " >> ";
      s.write_path(std::cerr);
    }
    NTYPE e = (current_layer_cost + s.delta_cost() + s.est_cost()) - get_cost_limit();
    std::cerr << "; limit = " << get_cost_limit() << ", e = " << e;
    std::cerr << std::endl;
  }

  NTYPE s_est = s.est_cost();
  if (INFINITE(s_est)) return POS_INF;

  NTYPE e = (current_layer_cost + s.delta_cost() + s_est) - get_cost_limit();
  if (e > 0) {
    least_over_bound = MIN(least_over_bound, e);
    return POS_INF;
  }

  Node* n = 0;
  if (previous_layer) {
    n = previous_layer->find_node(s);
    if (n) return POS_INF;
  }
  n = current_layer->find_node(s);
  if (n) return POS_INF;

  acc_succ += 1;

  n = next_layer->insert_node(s);
  if (!n->state) {
    new_succ += 1;
    stats.create_node(s);
    n->state = s.copy();
    n->state->set_predecessor(0);
    n->acc = current_layer_cost + s.delta_cost();
    n->est = s_est;
    n->val = n->acc + n->est;
    next_open.append(n);
  }

  return POS_INF;
}

BFHS::BFHS(Statistics& s, SearchResult& r, index_type lt_size)
  : SearchAlgorithm(s, r),
    layer_table_size(lt_size),
    previous_layer(0),
    current_layer(0),
    next_layer(0),
    current_open(0, 0),
    next_open(0, 0),
    current_layer_cost(0),
    least_over_bound(POS_INF)
{
  // done
}

BFHS::~BFHS()
{
  // nada
}

NTYPE BFHS::main(State& s) {
  current_layer_cost = 0;
  least_over_bound = POS_INF;
  previous_layer = 0;
  current_layer = new HashNodeSet(layer_table_size);
  current_open.clear();
  stats.create_node(s);
  Node* n = current_layer->insert_node(s);
  n->state = s.copy();
  n->state->set_predecessor(0);
  n->acc = 0;
  n->est = s.est_cost();
  n->val = n->acc + n->est;
  current_open.append(n);
  next_layer = new HashNodeSet(layer_table_size);
  next_open.clear();

  while (!solved() &&
	 (current_open.length() > 0) &&
	 !stats.break_signal_raised()) {
    if (trace_level > 1) {
      std::cerr << "current layer = " << current_layer_cost
		<< " (open = " << current_open.length()
		<< ", " << stats
		<< ")" << std::endl;
    }
    for (index_type k = 0;
	 (k < current_open.length()) &&
	   !solved() &&
	   !stats.break_signal_raised();
	 k++) {
      if (trace_level > 2) {
	std::cerr << "k = " << k
		  << ", #open = " << current_open.length()
		  << ", open[k].state = " << current_open[k]->state
		  << ", open[k] = " << current_open[k]->acc
		  << "/" << current_open[k]->est
		  << "/" << current_open[k]->val
		  << std::endl;
	std::cerr << "expanding " << *(current_open[k]->state)
		  << " (final = " << current_open[k]->state->is_final()
		  << ", val = " << current_open[k]->val
		  << ")" << std::endl;
      }
      assert(current_open[k]->state);
      assert(!current_open[k]->state->is_max());
      if (current_open[k]->state->is_final()) {
	set_solved(true, false);
      }
      else {
	stats.expand_node(*(current_open[k]->state));
	current_open[k]->state->expand(*this, POS_INF);
      }
    }

    if (!solved() && !stats.break_signal_raised()) {
      if (previous_layer) delete previous_layer;
      previous_layer = current_layer;
      current_layer = next_layer;
      current_layer_cost += 1;
      next_layer = new HashNodeSet(layer_table_size);
      current_open.clear();
      for (index_type k = 0; k < next_open.length(); k++)
	current_open.append(next_open[k]);
      // current_open.assign_copy(next_open);
      next_open.clear();
    }
  }

  if (previous_layer) delete previous_layer;
  delete current_layer;
  delete next_layer;
  return current_layer_cost;
}

NTYPE BFHS::start(State& s, NTYPE b)
{
  reset();
  set_cost_limit(b);
  stats.start();
  NTYPE val = main(s);
  stats.stop();
  return val;
}

NTYPE BFHS::start(State& s)
{
  return start(s, get_cost_limit());
}

NTYPE BFHS::cost() const
{
  return current_layer_cost;
}

bool BFHS::done() const
{
  return (solved() || stats.break_signal_raised());
}

BFIDA::BFIDA(Statistics& s, SearchResult& r, index_type lt_size)
  : BFHS(s, r, lt_size)
{
}

BFIDA::~BFIDA()
{
  // done
}

NTYPE BFIDA::start(State& s, NTYPE b)
{
  reset();
  set_cost_limit(b);
  stats.start();

  while (!solved() && FINITE(get_cost_limit()) && !stats.break_signal_raised()) {
    if (trace_level > 0) {
      std::cerr << "cost bound = " << get_cost_limit()
		<< " (" << stats << ")"
		<< std::endl;
    }
    NTYPE c = main(s);
    if (!solved() && !stats.break_signal_raised()) {
      // std::cerr << "least over bound = " << least_over_bound << std::endl;
      stats.current_lower_bound(get_cost_limit() + 1);
      set_cost_limit(get_cost_limit() + 1);
    }
  }
  stats.stop();
  return get_cost_limit();
}

NTYPE BFIDA::start(State& s)
{
  return start(s, s.est_cost());
}

BFHS_XD::BFHS_XD(Statistics& s, SearchResult& r, index_type lt_size)
  : BFHS(s, r, lt_size)
{
}

BFHS_XD::~BFHS_XD()
{
  // done
}

NTYPE BFHS_XD::start(State& s, NTYPE b)
{
  reset();
  set_cost_limit(b);
  stats.start();

  NTYPE d = 1;

  while (!solved() && FINITE(get_cost_limit()) && !stats.break_signal_raised()) {
    if (trace_level > 0) {
      std::cerr << "cost bound = " << get_cost_limit()
		<< " (" << stats << ")"
		<< std::endl;
    }
    NTYPE c = main(s);
    if (!solved() && !stats.break_signal_raised()) {
      set_cost_limit(s.est_cost() + d);
      d = (d * 2);
    }
  }
  stats.stop();
  return get_cost_limit();
}

NTYPE BFHS_XD::start(State& s)
{
  return start(s, 1);
}

END_HSPS_NAMESPACE
