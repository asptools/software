
#include "mutex.h"

BEGIN_HSPS_NAMESPACE

void StaticMutex::compute(const index_set& init, const bool_vec* aa)
{
  read_only_vector_with_default<bool> aa_def(aa, true);

  pairix.init(instance.n_atoms());
  tab.assign_value(false, pairix.n_pairs());
  for (index_type i = 0; i < init.size(); i++)
    for (index_type j = i; j < init.size(); j++)
      tab[pairix(init[i], init[j])] = true;

  bool done = false;
  while (!done) {
#ifdef TRACE_PRINT_LOTS
    std::cerr << "begin iteration" << std::endl;
#endif
    done = true;
    for (index_type a = 0; a < instance.n_actions(); a++)
      if (aa_def[a]) {
	const Instance::Action& act = instance.actions[a];
	if (!mutex(act.pre)) {
#ifdef TRACE_PRINT_LOTS
	  std::cerr << "applying action " << a << "."
		    << instance.actions[a].name << std::endl;
#endif
	  for (index_type i = 0; i < instance.n_atoms(); i++) {
	    // if i in act.add...
	    if (act.add.contains(i)) {
	      // update i
	      if (!tab[pairix(i, i)]) {
#ifdef TRACE_PRINT_LOTS
		std::cerr << "updating " << i << "."
			  << instance.atoms[i].name << std::endl;
#endif
		tab[pairix(i, i)] = true;
		done = false;
	      }
	      for (index_type j = 0; j < act.add.length(); j++)
		// update {i,q} for every q in act.add
		if (act.add[j] > i) // q < i have already been updated (as i)
		  if (!tab[pairix(act.add[j], i)]) {
#ifdef TRACE_PRINT_LOTS
		    std::cerr << "updating " << i << "."
			      << instance.atoms[i].name
			      << ", " << act.add[j]
			      << instance.atoms[act.add[j]].name
			      << std::endl;
#endif
		    tab[pairix(act.add[j], i)] = true;
		    done = false;
		  }
	    }
	    // else if i not in act.del (i.e. case act + noop(i))...
	    else if (tab[pairix(i, i)] && !act.del.contains(i)) {
	      // compute H(act.pre U {i})
	      bool ok = true;
	      for (index_type k = 0; (k < act.pre.length()) && ok; k++)
		if (!tab[pairix(i, act.pre[k])])
		  ok = false;
	      // update {i,q} for every q in act.add
	      if (ok) {
		for (index_type j = 0; j < act.add.length(); j++)
		  if (!tab[pairix(i, act.add[j])]) {
#ifdef TRACE_PRINT_LOTS
		    std::cerr << "updating " << i << "."
			      << instance.atoms[i].name
			      << ", " << act.add[j]
			      << instance.atoms[act.add[j]].name
			      << std::endl;
#endif
		    tab[pairix(i, act.add[j])] = true;
		    done = false;
		  }
	      }
	    }
	  } // for each atom i
	} // if (eval(act.pre))
      } // for each action
#ifdef TRACE_PRINT_LOTS
    std::cerr << "end iteration, done = " << done << std::endl;
#endif
  }
}

StaticMutex::StaticMutex
(Instance& ins, bool compute_on_construction)
  : Heuristic(ins)
{
  if (compute_on_construction)
    compute(instance.init_atoms);
}

StaticMutex::StaticMutex
(Instance& ins, const index_set& init)
  : Heuristic(ins)
{
  compute(init);
}

StaticMutex::StaticMutex
(Instance& ins, const index_set& init, const bool_vec& aa)
  : Heuristic(ins)
{
  compute(init, &aa);
}

StaticMutex::StaticMutex
(Instance& ins, const index_set& init, const bool_vec* aa)
  : Heuristic(ins)
{
  compute(init, aa);
}

StaticMutex::StaticMutex(Instance& ins, const StaticMutex& mx)
  : Heuristic(ins)
{
  pairix.init(instance.n_atoms());
  tab.assign_value(false, pairix.n_pairs());
  assert(mx.tab.size() >= pairix.n_pairs());
  for (index_type i = 0; i < pairix.n_pairs(); i++)
    tab[i] = mx.tab[i];
}

void StaticMutex::recompute()
{
  compute(instance.init_atoms);
}

void StaticMutex::recompute(const index_set& init)
{
  compute(init);
}

void StaticMutex::recompute(const bool_vec& aa)
{
  compute(instance.init_atoms, &aa);
}

void StaticMutex::recompute(const index_set& init, const bool_vec& aa)
{
  compute(init, &aa);
}

bool StaticMutex::unreachable(index_type i) const
{
  if (!pairix.defined(i)) return false;
  return !tab[pairix(i, i)];
}

bool StaticMutex::mutex(index_type i, index_type j) const
{
  if (!pairix.defined(i)) return false;
  if (!pairix.defined(j)) return false;
  return !tab[pairix(i, j)];
}

bool StaticMutex::mutex(const index_set& s) const
{
#ifdef TRACE_PRINT_LOTS
  std::cerr << "eval " << s << std::endl;
#endif
  for (index_type i = 0; i < s.size(); i++)
    for (index_type j = i; j < s.size(); j++)
      if (pairix.defined(s[i]) && pairix.defined(s[j]))
	if (!tab[pairix(s[i], s[j])]) {
#ifdef TRACE_PRINT_LOTS
	  std::cerr << " - (" << s[i] << "," << s[j] << ") = pair #"
		    << pairix(s[i], s[j]) << " is false" << std::endl;
#endif
	  return true;
	}
#ifdef TRACE_PRINT_LOTS
  std::cerr << " - return true" << std::endl;
#endif
  return false;
}

bool StaticMutex::mutex(const index_set& s, index_type i) const
{
  if (!pairix.defined(i)) return false;
  if (!tab[pairix(i, i)]) return true;
  if (mutex(s)) return true;
  for (index_type k = 0; k < s.size(); k++)
    if (pairix.defined(s[k]))
      if (!tab[pairix(s[k], i)])
	return true;
  return false;
}

bool StaticMutex::mutex(const index_set& s1, const index_set& s2) const
{
  if (mutex(s1)) return true;
  if (mutex(s2)) return true;
  for (index_type k1 = 0; k1 < s1.size(); k1++)
    for (index_type k2 = 0; k2 < s2.size(); k2++)
      if (pairix.defined(s1[k1]) && pairix.defined(s2[k2]))
	if (!tab[pairix(s1[k1], s2[k2])])
	  return true;
  return false;
}

NTYPE StaticMutex::eval(const index_set& s)
{
  return (mutex(s) ? POS_INF : 0);
}

NTYPE StaticMutex::eval(const bool_vec& s)
{
  return (mutex(s) ? POS_INF : 0);
}

NTYPE StaticMutex::incremental_eval(const index_set& s, index_type i_new)
{
  return (mutex(s, i_new) ? POS_INF : 0);
}

NTYPE StaticMutex::incremental_eval(const bool_vec& s, index_type i_new)
{
  return (mutex(s, i_new) ? POS_INF : 0);
}

NTYPE StaticMutex::eval(index_type atom)
{
  return (mutex(atom, atom) ? POS_INF : 0);
}

NTYPE StaticMutex::eval(index_type atom1, index_type atom2)
{
  return (mutex(atom1, atom2) ? POS_INF : 0);
}

void StaticMutex::write(std::ostream& s) const
{
  s << '{';
  bool first = true;
  for (index_type i = 0; i < instance.n_atoms(); i++) {
    if (unreachable(i)) {
      if (!first)
	s << ',';
      else
	first = false;
      s << instance.atoms[i].name;
    }
    else {
      for (index_type j = i + 1; j < instance.n_atoms(); j++)
	if (!unreachable(j) && mutex(i, j)) {
	  if (!first)
	    s << ',';
	  else
	    first = false;
	  s << '(' << instance.atoms[i].name << ','
	    << instance.atoms[j].name << ')';
	}
    }
  }
  s << '}';
}

//#define TRACE_PRINT_LOTS

void Reachability::compute_main
(read_only_vector_with_default<bool>& aa)
{
#ifdef TRACE_PRINT_LOTS
  std::cerr << "initial reached = " << index_set(tab) << std::endl;
#endif
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (aa[k] && (rem_pre[k] == 0)) {
#ifdef TRACE_PRINT_LOTS
      std::cerr << "applying " << k << "..." << std::endl;
#endif
      //atab[k] = true;
      for (index_type i = 0; i < instance.actions[k].add.size(); i++)
	if (!tab[instance.actions[k].add[i]]) {
	  tab[instance.actions[k].add[i]] = true;
	  q[q_tail++] = instance.actions[k].add[i];
	}
    }
  while (q_head < q_tail) {
#ifdef TRACE_PRINT_LOTS
    std::cerr << "reached = " << index_set(tab) << std::endl;
    std::cerr << "q = " << q << " (" << q_head << ", " << q_tail << ")"
      	      << std::endl;
#endif
    assert(q_tail <= instance.n_atoms());
    index_type p = q[q_head++];
    for (index_type k = 0; k < instance.atoms[p].req_by.size(); k++) {
      index_type a = instance.atoms[p].req_by[k];
#ifdef TRACE_PRINT_LOTS
      std::cerr << " - req by " << a << std::endl;
#endif
      assert(rem_pre[a] > 0);
      rem_pre[a] -= 1;
      if (aa[a] && (rem_pre[a] == 0)) {
#ifdef TRACE_PRINT_LOTS
	std::cerr << "applying " << a << "..." << std::endl;
#endif
	//atab[a] = true;
	for (index_type i = 0; i < instance.actions[a].add.size(); i++)
	  if (!tab[instance.actions[a].add[i]]) {
	    tab[instance.actions[a].add[i]] = true;
	    q[q_tail++] = instance.actions[a].add[i];
	  }
      }
    }
  }
#ifdef TRACE_PRINT_LOTS
  std::cerr << "final reached = " << index_set(tab) << std::endl;
#endif
}

void Reachability::compute(const index_set& init, const bool_vec* aa)
{
  read_only_vector_with_default<bool> aa_def(aa, true);
  // index_type n = 0;
  // for (index_type i = 0; i < instance.n_actions(); i++)
  //   if (aa_def[i]) n += 1;
  // std::cerr << "Reachability::compute: " << n << "/" << instance.n_actions()
  // 	    << std::endl;
  tab.assign_value(false, instance.n_atoms());
  //atab.assign_value(false, instance.n_actions());
  for (index_type k = 0; k < instance.n_actions(); k++)
    rem_pre[k] = instance.actions[k].pre.size();
  q_head = 0;
  q_tail = 0;
  for (index_type i = 0; i < init.size(); i++) {
    tab[init[i]] = true;
    q[q_tail++] = init[i];
  }
  compute_main(aa_def);
}

void Reachability::compute(const bool_vec& init, const bool_vec* aa)
{
  assert(init.size() == instance.n_atoms());
  read_only_vector_with_default<bool> aa_def(aa, true);
  tab.assign_value(false, instance.n_atoms());
  //atab.assign_value(false, instance.n_actions());
  for (index_type k = 0; k < instance.n_actions(); k++)
    rem_pre[k] = instance.actions[k].pre.size();
  q_head = 0;
  q_tail = 0;
  for (index_type i = 0; i < instance.n_atoms(); i++)
    if (init[i]) {
      tab[i] = true;
      q[q_tail++] = i;
    }
  compute_main(aa_def);
}

void Reachability::init_structs()
{
  tab.assign_value(false, instance.n_atoms());
  //atab.assign_value(false, instance.n_actions());
  rem_pre.assign_value(0, instance.n_actions());
  q.assign_value(no_such_index, instance.n_atoms());
}

void Reachability::save_state(reachability_state& s) const
{
  s.tab.assign_copy(tab);
  //s.atab.assign_copy(atab);
  s.rem_pre.assign_copy(rem_pre);
  s.q_head = q_head;
}

void Reachability::restore_state(const reachability_state& s)
{
  tab.assign_copy(s.tab);
  //atab.assign_copy(s.atab);
  rem_pre.assign_copy(s.rem_pre);
  q_head = s.q_head;
  q_tail = s.q_head;
}

Reachability::Reachability
(Instance& ins, bool compute_on_construction)
  : Heuristic(ins)
{
  init_structs();
  if (compute_on_construction)
    compute(ins.init_atoms, 0);
}

Reachability::Reachability
(Instance& ins, const index_set& init)
  : Heuristic(ins)
{
  init_structs();
  compute(init, 0);
}

Reachability::Reachability
(Instance& ins, const index_set& init, const bool_vec& aa)
  : Heuristic(ins)
{
  init_structs();
  compute(init, &aa);
}

Reachability::Reachability
(Instance& ins, const index_set& init, const bool_vec* aa)
  : Heuristic(ins)
{
  init_structs();
  compute(init, aa);
}

void Reachability::recompute()
{
  compute(instance.init_atoms, 0);
}

void Reachability::recompute(const index_set& init)
{
  compute(init, 0);
}

void Reachability::recompute(const bool_vec& aa)
{
  compute(instance.init_atoms, &aa);
}

void Reachability::recompute(const index_set& init, const bool_vec& aa)
{
  compute(init, &aa);
}

void Reachability::recompute(const index_set& init, const index_set& aa)
{
  bool_vec bv_aa(aa, instance.n_actions());
  compute(init, &bv_aa);
}

void Reachability::recompute_bv_init(const bool_vec& init)
{
  compute(init, 0);
}

void Reachability::update(const bool_vec& aa)
{
  read_only_vector_with_default<bool> aa_def(&aa, true);
  compute_main(aa_def);
}

void Reachability::update(const bool_vec& aa, index_type new_aa)
{
  if (rem_pre[new_aa] == 0) {
    for (index_type i = 0; i < instance.actions[new_aa].add.size(); i++)
      if (!tab[instance.actions[new_aa].add[i]]) {
	tab[instance.actions[new_aa].add[i]] = true;
	q[q_tail++] = instance.actions[new_aa].add[i];
      }
  }
  while (q_head < q_tail) {
    assert(q_tail <= instance.n_atoms());
    index_type p = q[q_head++];
    for (index_type k = 0; k < instance.atoms[p].req_by.size(); k++) {
      index_type a = instance.atoms[p].req_by[k];
      assert(rem_pre[a] > 0);
      rem_pre[a] -= 1;
      if (aa[a] && (rem_pre[a] == 0)) {
	for (index_type i = 0; i < instance.actions[a].add.size(); i++)
	  if (!tab[instance.actions[a].add[i]]) {
	    tab[instance.actions[a].add[i]] = true;
	    q[q_tail++] = instance.actions[a].add[i];
	  }
      }
    }
  }
}

bool Reachability::order_relaxed_plan
(const index_set& init, const index_set& goal,
 const index_set& aa, index_vec& ao)
{
  bool_vec raa(aa, instance.n_actions());
  tab.assign_value(false, instance.n_atoms());
  for (index_type k = 0; k < instance.n_actions(); k++)
    rem_pre[k] = instance.actions[k].pre.size();
  q_head = 0;
  q_tail = 0;
  for (index_type i = 0; i < init.size(); i++) {
    tab[init[i]] = true;
    q[q_tail++] = init[i];
  }
  ao.clear();
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (raa[k] && (rem_pre[k] == 0)) {
      raa[k] = false;
      ao.append(k);
      for (index_type i = 0; i < instance.actions[k].add.size(); i++)
	if (!tab[instance.actions[k].add[i]]) {
	  tab[instance.actions[k].add[i]] = true;
	  q[q_tail++] = instance.actions[k].add[i];
	}
      if (tab.contains(goal)) return true;
    }
  while (q_head < q_tail) {
    assert(q_tail <= instance.n_atoms());
    index_type p = q[q_head++];
    for (index_type k = 0; k < instance.atoms[p].req_by.size(); k++) {
      index_type a = instance.atoms[p].req_by[k];
      assert(rem_pre[a] > 0);
      rem_pre[a] -= 1;
      if (raa[a] && (rem_pre[a] == 0)) {
	raa[a] = false;
	ao.append(a);
	for (index_type i = 0; i < instance.actions[a].add.size(); i++)
	  if (!tab[instance.actions[a].add[i]]) {
	    tab[instance.actions[a].add[i]] = true;
	    q[q_tail++] = instance.actions[a].add[i];
	  }
	if (tab.contains(goal)) return true;
      }
    }
  }
  return false;
}

bool Reachability::unreachable(index_type i) const
{
  // if (tab[i])
  //   return false;
  // else
  //   return true;
  return !(tab[i]);
}

bool Reachability::unreachable(const index_set& s) const
{
  for (index_type k = 0; k < s.size(); k++)
    if (!tab[s[k]]) return true;
  return false;
}

bool Reachability::unreachable(const bool_vec& s) const
{
  for (index_type k = 0; k < instance.n_atoms(); k++)
    if (s[k])
      if (!tab[k]) return true;
  return false;
}

void Reachability::unreachable_subset(const index_set& s, index_set& u) const
{
  u.clear();
  for (index_type k = 0; k < s.size(); k++)
    if (!tab[s[k]])
      u.insert(s[k]);
}

void Reachability::write(std::ostream& s) const
{
  bool_vec u(tab);
  u.complement();
  instance.write_atom_set(s, u);
}

NTYPE Reachability::eval(const index_set& s)
{
  return (unreachable(s) ? POS_INF : 0);
}

NTYPE Reachability::eval(const bool_vec& s)
{
  return (unreachable(s) ? POS_INF : 0);
}

NTYPE Reachability::incremental_eval(const index_set& s, index_type i_new)
{
  return ((unreachable(i_new) || unreachable(s)) ? POS_INF : 0);
}

NTYPE Reachability::incremental_eval(const bool_vec& s, index_type i_new)
{
  return ((unreachable(i_new) || unreachable(s)) ? POS_INF : 0);
}

NTYPE Reachability::eval(index_type atom)
{
  return (unreachable(atom) ? POS_INF : 0);
}

NTYPE Reachability::eval(index_type atom1, index_type atom2)
{
  return ((unreachable(atom1) || unreachable(atom2)) ? POS_INF : 0);
}


void landmark_graph_viaP2
(const bool_vec& atms, const bool_vec& acts,
 Instance& insP2, s2index& pair_map, index_vec& act_map,
 graph& g)
{
  g.init(atms.size());
  bool_vec aa1(true, insP2.n_actions());
  for (index_type k = 0; k < insP2.n_actions(); k++)
    if (!acts[act_map[k]])
      aa1[k] = false;
  Reachability inc2(insP2, false);
  for (index_type i = 0; i < atms.size(); i++)
    if (atms[i]) {
      bool_vec aa(aa1);
      for (index_type k = 0; k < insP2.atoms[pair_map(i, i)].req_by.size(); k++)
	aa[insP2.atoms[pair_map(i, i)].req_by[k]] = false;
      inc2.recompute(aa);
      for (index_type j = 0; j < atms.size(); j++)
	if ((i != j) && atms[j]) {
	  if (inc2.unreachable(pair_map(j, j)))
	    g.add_edge(i, j);
	}
    }
}

void landmark_graph_viaP2
(Instance& ins, const bool_vec& atms, const bool_vec& acts,
 Heuristic* inc, graph& g)
{
  Instance insP2(ins.name);
  s2index   pair_map;
  index_vec act_map;
  insP2.makeP2(ins, inc, pair_map, act_map);
  insP2.cross_reference();
  landmark_graph_viaP2(atms, acts, insP2, pair_map, act_map, g);
}

void landmark_graph_viaP2
(Instance& ins, const bool_vec& atms, const bool_vec& acts, graph& g)
{
  StaticMutex mx(ins, acts);
  landmark_graph_viaP2(ins, atms, acts, &mx, g);
}


void landmark_graph_triggered_edges
(Instance& ins, graph& lmg, set_edge_vec& trev)
{
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (!ins.atoms[i].init) {
      index_set cand;
      for (index_type k = 0; k < ins.atoms[i].add_by.size(); k++)
	cand.insert(ins.actions[ins.atoms[i].add_by[k]].pre);
      cand.subtract(lmg.predecessors(i));
      for (index_type j = 0; j < cand.size(); j++) {
	index_set act_to_rem;
	for (index_type k = 0; k < ins.atoms[i].add_by.size(); k++)
	  if (!ins.actions[ins.atoms[i].add_by[k]].pre.contains(cand[j]))
	    act_to_rem.insert(ins.atoms[i].add_by[k]);
	assert(act_to_rem.size() < ins.atoms[i].add_by.size());
	trev.append(std::pair<index_set, index_pair>
		    (act_to_rem, index_pair(cand[j], i)));
      }
    }
}

ForwardReachabilityCheck::ForwardReachabilityCheck
(Instance& i, const index_set& g)
  : Reachability(i, false), goals(g)
{
  // done
}

ForwardReachabilityCheck::~ForwardReachabilityCheck()
{
  // done
}

NTYPE ForwardReachabilityCheck::eval(const index_set& s)
{
  recompute(s);
  return Reachability::eval(goals);
}

NTYPE ForwardReachabilityCheck::eval(const bool_vec& s)
{
  recompute_bv_init(s);
  return Reachability::eval(goals);
}

END_HSPS_NAMESPACE
