
#include "additive.h"
#include "preprocess.h"
#include "enumerators.h"

BEGIN_HSPS_NAMESPACE

count_type AH::Hmax_wins = 0;
count_type AH::Hsum_wins = 0;
count_type AH::draws = 0;

bool AH::use_linear_scan_eval = false;

AH::AH(Instance& i, Stopwatch& s)
  : Heuristic(i), Hmax(0), H_vec(0, 0), stats(s)
{
  // done
}

AH::~AH()
{
  if (Hmax) delete Hmax;
  for (index_type k = 0; k < H_vec.length(); k++)
    if (H_vec[k]) delete H_vec[k];
}

void AH::compute_additive(const ACF& cost, const index_set_vec& p, bool useH2)
{
  stats.start();

  index_type i_first = H_vec.length();
  index_type n_useless = 0;

  for (index_type k = 0; k < p.length(); k++) if (!p[k].empty()) {
    if (stats.break_signal_raised()) return;
    // bool_vec dis(p[k], instance.n_actions());
    // dis.complement();
    if (trace_level > 0) {
      std::cerr << "counting " << p[k].length()
		<< " of " << instance.n_actions() << " actions..."
		<< std::endl;
      if (trace_level > 1) {
	instance.write_action_set(std::cerr, p[k]);
	std::cerr << std::endl;
      }
    }
    // DiscountACF d_cost(cost, dis);
    FracACF d_cost(cost, instance.n_actions(), 0);
    d_cost.set(p[k], 1);
    CostTable* Hd = (use_linear_scan_eval ?
		     new LinearScanH2Eval(instance, stats) :
		     new CostTable(instance, stats));
    if (useH2)
      Hd->compute_H2(d_cost);
    else
      Hd->compute_H1(d_cost);
    if (Hd->max_finite() == 0) {
      n_useless += 1;
    }
    else {
      H_vec.append(Hd);
    }
  }

  std::cerr << H_vec.length() << " additive h^" << (useH2 ? 2 : 1)
	    << " built" << std::endl;
  if (n_useless > 0) {
    std::cerr << n_useless << " useless action sets detected" << std::endl;
  }

  if (use_linear_scan_eval) {
    for (index_type i = i_first; i < H_vec.length(); i++) {
      if (stats.break_signal_raised()) return;
      ((LinearScanH2Eval*)H_vec[i])->compile_finite();
    }
  }

  stats.stop();
}

void AH::compute_fractional(const ACF& cost, index_set_vec& p, bool useH2)
{
  stats.start();

  index_vec occ(0, instance.n_actions());
  for (index_type k = 0; k < p.length(); k++)
    for (index_type i = 0; i < p[k].length(); i++)
      occ[p[k][i]] += 1;

  index_type i_first = H_vec.length();
  index_set  useless_p;

  for (index_type k = 0; k < p.length(); k++) if (!p[k].empty()) {
    if (stats.break_signal_raised()) return;
    if (trace_level > 0) {
      std::cerr << "counting " << p[k].length()
		<< " of " << instance.n_actions() << " actions..."
		<< std::endl;
      if (trace_level > 1) {
	instance.write_action_set(std::cerr, p[k]);
	std::cerr << std::endl;
      }
    }
    FracACF d_cost(cost, instance.n_actions(), 0);
    for (index_type i = 0; i < p[k].length(); i++)
      d_cost.set(p[k][i], R_TO_N(1, occ[p[k][i]]));
    CostTable* Hd = (use_linear_scan_eval ?
		     new LinearScanH2Eval(instance, stats) :
		     new CostTable(instance, stats));
    if (useH2)
      Hd->compute_H2(d_cost);
    else
      Hd->compute_H1(d_cost);
    H_vec.append(Hd);
    if (H_vec[H_vec.length() - 1]->max_finite() == 0)
      useless_p.insert(k);
  }

  std::cerr << H_vec.length() << " additive h^" << (useH2 ? 2 : 1)
	    << " built" << std::endl;

  if ((useless_p.length() > 0) && ((H_vec.length() - i_first) > 0)) {
    std::cerr << "useless action sets detected, recomputing..." << std::endl;
    for (index_type i = i_first; i < H_vec.length(); i++)
      delete H_vec[i];
    H_vec.set_length(i_first);
    for (index_type i = 0; i < useless_p.length(); i++)
      p[useless_p[i]].clear();
    compute_fractional(cost, p, useH2);
  }
  else if (use_linear_scan_eval) {
    for (index_type i = i_first; i < H_vec.length(); i++) {
      if (stats.break_signal_raised()) return;
      ((LinearScanH2Eval*)H_vec[i])->compile_finite();
    }
  }

  stats.stop();
}

void AH::compute_max(const ACF& cost, bool useH2)
{
  if (Hmax) delete Hmax;
  std::cerr << "computing h^2..." << std::endl;
  Hmax = new CostTable(instance, stats);
  Hmax->compute_H2(cost);  
}

void AH::disable_max()
{
  if (Hmax) delete Hmax;
  Hmax = 0;
}

class decreasing_cost : public index_vec::order {
  CostTable* h;
 public:
  decreasing_cost(CostTable* _h) : h(_h) { };
  virtual bool operator()(const index_type& v0, const index_type& v1) const;
};

bool decreasing_cost::operator()
(const index_type& v0, const index_type& v1) const
{
  return (h->eval(v0) >= h->eval(v1));
}

class increasing_cost : public index_vec::order {
  CostTable* h;
 public:
  increasing_cost(CostTable* _h) : h(_h) { };
  virtual bool operator()(const index_type& v0, const index_type& v1) const;
};

bool increasing_cost::operator()
(const index_type& v0, const index_type& v1) const
{
  return (h->eval(v0) <= h->eval(v1));
}

void AH::compute_with_relevance_partitioning
(const ACF& cost, const index_set& g)
{
  stats.start();
  Hmax = new CostTable(instance, stats);
  Hmax->compute_H2(cost);

  // note: instance is const, and preprocessor does not take a const
  // instance as argument, because some preprocessor methods modify
  // the instance, but we won't call any of those here...
  Preprocessor* prep = new Preprocessor((Instance&)instance, stats);
  bool_vec rem(true, instance.n_actions());
  index_type n_rem = instance.n_actions();

  index_vec g_sort;
  g_sort.insert_ordered(g, decreasing_cost(Hmax));

  for (index_type k = 0; (k < g_sort.length()) && (rem.count(true) > 0); k++) {
    bool_vec arg(false, instance.n_actions()); // Actions Relevant to Goal
    prep->strictly_relevant_actions(g_sort[k], Hmax->eval(g_sort[k]),
				    *Hmax, cost, arg);
    arg.intersect(rem);
    if (arg.count(true) > 0) {
      if (trace_level > 0) {
	std::cerr << "goal " << instance.atoms[g_sort[k]].name
		  << ": counting " << arg.count(true)
		  << " of " << instance.n_actions() << " actions..."
		  << std::endl;
	if (trace_level > 1) {
	  instance.write_action_set(std::cerr, arg);
	  std::cerr << std::endl;
	}
      }
      bool_vec d(arg);
      d.complement();
      DiscountACF d_cost(cost, d);
      CostTable* Hd = new CostTable(instance, stats);
      Hd->compute_H2(d_cost);
      H_vec.append(Hd);
      rem.subtract(arg);
    }
  }

  if (rem.count(true) > 0) {
    if (trace_level > 0) {
      std::cerr << "remaining: " << rem.count(true)
		<< " of " << instance.n_actions() << " actions..."
		<< std::endl;
      if (trace_level > 1) {
	instance.write_action_set(std::cerr, rem);
	std::cerr << std::endl;
      }
    }
    bool_vec d(rem);
    d.complement();
    DiscountACF d_cost(cost, d);
    CostTable* Hd = new CostTable(instance, stats);
    Hd->compute_H2(d_cost);
    H_vec.append(Hd);
  }

  if (trace_level > 0) {
    std::cerr << H_vec.length() << " additive H2 built" << std::endl;
  }

  delete prep;
  stats.stop();
}

// void AH::compute_with_random_relevance_partitioning
// (const ACF& cost, const index_set& g, RNG& rnd, bool useH2)
// {
//   stats.start();
// 
//   Hmax = new CostTable(instance, stats);
//   if (useH2) {
//     Hmax->compute_H2(cost);
//   }
//   else {
//     Hmax->compute_H1(cost);
//   }
// 
//   index_set_vec p; // action partitions
//   p.assign_value(EMPTYSET, g.length());
// 
//   Preprocessor* prep = new Preprocessor(instance, stats);
// 
//   index_set_vec r; // relevant actions
//   r.assign_value(EMPTYSET, g.length());
// 
//   NTYPE offset = 0;
//   bool_vec rem(true, instance.n_actions());
//   bool conv = false;
// 
//   while (!conv && (rem.count(true) > 0)) {
//     conv = true;
//     if (trace_level > 0) {
//       std::cerr << "checking relevance at H() + " << offset << "..."
// 		<< std::endl;
//     }
// 
//     for (index_type k = 0; k < g.length(); k++) {
//       bool_vec rel(false, instance.n_actions());
//       prep->strictly_relevant_actions(g[k], Hmax->eval(g[k]) + offset,
// 				      *Hmax, cost, rel);
//       rel.copy_to(r[k]);
//       if (trace_level > 0) {
// 	std::cerr << r[k].length() << " actions relevant to "
// 		  << instance.atoms[g[k]].name << " at H() + " << offset
// 		  << std::endl;
//       }
//     }
// 
//     for (index_type k = 0; k < instance.n_actions(); k++) if (rem[k]) {
//       index_set s;
//       for (index_type i = 0; i < g.length(); i++)
// 	if (r[i].contains(k)) s.insert(i);
//       if (s.length() == 1) {
// 	p[s[0]].insert(k);
// 	if (trace_level > 0) {
// 	  std::cerr << "action " << instance.actions[k].name
// 		    << " relevant to single goal: assigned to "
// 		    << instance.atoms[g[s[0]]].name << " (" << s[0] << ")"
// 		    << std::endl;
// 	}
// 	rem[k] = false;
// 	conv = false;
//       }
//       else if (s.length() > 1) {
// 	index_type l = (rnd.random() % s.length());
// 	p[s[l]].insert(k);
// 	if (trace_level > 0) {
// 	  std::cerr << "action " << instance.actions[k].name
// 		    << " relevant to " << s.length() << " goals: assigned to "
// 		    << instance.atoms[g[s[l]]].name << " (" << s[l] << ")"
// 		    << std::endl;
// 	}
// 	rem[k] = false;
// 	conv = false;
//       }
//     }
// 
//     offset += 1;
//     if (trace_level > 0) {
//       std::cerr << rem.count(true) << " actions remaining" << std::endl;
//     }
//   }
// 
// //   for (index_type k = 0; k < instance.n_actions(); k++) if (rem[k]) {
// //     index_type l = (rnd.random() % g.length());
// //     p[l].insert(k);
// //     std::cerr << "action " << instance.actions[k].name
// // 	      << " not relevant to any goal: assigned to "
// // 	      << instance.atoms[g[l]].name << " (" << l << ")"
// // 	      << std::endl;
// //   }
// 
//   compute_additive(cost, p, useH2);
// 
//   delete prep;
//   stats.stop();
// }

#define ROLLBACK

#ifdef ROLLBACK
void AH::compute_with_iterative_assignment
(const ACF& cost, const index_set& g, index_set_vec& app,
 bool useH2, bool fractional, const index_set& g_prio)
{
  stats.start();

  std::cerr << "computing h^1..." << std::endl;
  Hmax = new CostTable(instance, stats);
  Hmax->compute_H1(cost);
  NTYPE h_limit = Hmax->eval(g_prio);

  index_set_vec p; // action partition
  p.assign_value(EMPTYSET, g.length());
  bool_vec rem(true, instance.n_actions());
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (IS_ZERO(cost(k))) rem[k] = false;

  std::cerr << app.length() << " action pre-partitions..." << std::endl;
  for (index_type k = 0; k < app.length(); k++) if (!app[k].empty()) {
    if (!fractional) {
      if (app[k].count_common(rem) != app[k].length()) {
	std::cerr << "error: action pre-partitioning is not a partitioning!"
		  << std::endl;
	exit(255);
      }
    }
    if (trace_level > 0) {
      std::cerr << "classifying set of " << app[k].length() << " actions..."
		<< std::endl;
    }
    cost_vec l(0, g.length());
    for (index_type i = 0; i < g.length(); i++) {
      DiscountACF d_cost(cost, instance.n_actions());
      d_cost.discount(app[k]);
      CostTable* Hd = new CostTable(instance, stats);
      Hd->compute_H1(d_cost);
      l[i] = (Hmax->eval(g[i]) - Hd->eval(g[i]));
      if (trace_level > 0) {
	std::cerr << "loss for goal " << instance.atoms[g[i]].name
		  << " = " << l[i] << std::endl;
      }
      delete Hd;
    }
    index_type mi = l.arg_max();
    assert(mi != no_such_index);
    if (fractional) {
      for (index_type i = 0; i < g.length(); i++) if (l[i] == l[mi]) {
	if (trace_level > 0) {
	  std::cerr << "set assigned to goal " << instance.atoms[g[i]].name
		    << std::endl;
	}
	p[i].insert(app[k]);
      }
    }
    else {
      index_type mn = p[mi].length();
      NTYPE      mv = Hmax->eval(g[mi]);
      bool is_goal = g_prio.contains(g[mi]);
      for (index_type i = mi + 1; i < g.length(); i++)
	if (l[i] == l[mi]) {
	  NTYPE vi = Hmax->eval(g[i]);
	  bool gi_is_goal = (!is_goal ? g_prio.contains(g[i]) : false);
	  // alt. 1: break ties towards higher (but not too high) values:
	  //  if ((vi > mv) && (vi <= h_limit))
	  // alt. 2: break ties towards lower (but greater-than-zero) values:
	  if (vi < mv) {
	    mi = i;
	    mn = p[i].length();
	    mv = Hmax->eval(g[mi]);
	  }
	  else if ((vi == mv) && !is_goal && gi_is_goal) {
	    mi = i;
	    mn = p[i].length();
	    is_goal = true;
	  }
	  // alt. 1: prefer distribution:
	  else if ((vi == mv) && (p[i].length() < mn)) {
	  // alt. 2: prefer concentration:
	  // else if ((vi == mv) && (p[i].length() > mn)) {
	    mi = i;
	    mn = p[i].length();
	  }
	}
      assert(mi < g.length());
      if (trace_level > 0) {
	std::cerr << "set assigned to goal " << instance.atoms[g[mi]].name
		  << std::endl;
      }
      p[mi].insert(app[k]);
    }
    rem.subtract(app[k]);
  }

  bool done = (rem.count(true) == 0);
  if (!done) {
    index_set_vec p_new;
    p_new.assign_value(EMPTYSET, g.length());
    while (!done) {
      std::cerr << rem.count(true) << " actions remain..." << std::endl;
      done = true;
      for (index_type k = 0; k < instance.n_actions(); k++) if (rem[k]) {
	bool assigned = false;
	for (index_type i = 0; (i < g.length()) && !assigned; i++) {
	  DiscountACF d_cost(cost, instance.n_actions());
	  for (index_type j = 0; j < g.length(); j++)
	    if (j != i) d_cost.discount(p[j]);
	  d_cost.discount(k);
	  CostTable* Hd = new CostTable(instance, stats);
	  Hd->compute_H1(d_cost);
	  if (Hd->eval(g[i]) < Hmax->eval(g[i])) {
	    p_new[i].insert(k);
	    rem[k] = false;
	    if (!fractional)
	      assigned = true;
	    done = false;
	  }
	  delete Hd;
	}
      }
      for (index_type i = 0; i < g.length(); i++) {
	p[i].insert(p_new[i]);
	p_new[i].clear();
      }
    }
  }

  if (rem.count(true) > 0) {
    std::cerr << rem.count(true) << " actions still remain..." << std::endl;
    if (fractional) {
      for (index_type i = 0; i < g.length(); i++)
	if (!p[i].empty())
	  p[i].insert(rem);
    }
  }

  if (fractional)
    compute_fractional(cost, p, useH2);
  else
    compute_additive(cost, p, useH2);
  if (useH2)
    compute_max(cost, true);

  stats.stop();
}
#else
void AH::compute_with_iterative_assignment
(const ACF& cost, const index_set& g, index_set_vec& app,
 bool useH2, bool fractional, const index_set& g_prio)
{
  stats.start();
  assert(!fractional);

  std::cerr << "computing h^1..." << std::endl;
  Hmax = new CostTable(instance, stats);
  Hmax->compute_H1(cost);
  NTYPE h_limit = Hmax->eval(g_prio);

  // extend pre-partitioning to cover all actions
  bool_vec rem(true, instance.n_actions());
  for (index_type k = 0; k < app.length(); k++) {
    if (app[k].count_common(rem) != app[k].length()) {
      std::cerr << "error: action pre-partitioning is not a partitioning!"
		<< std::endl;
      exit(255);
    }
    rem.subtract(app[k]);
  }
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (rem[k] && (cost(k) > 0)) {
      app.append(EMPTYSET);
      app[app.length() - 1].assign_singleton(k);
    }

  // action partitioning
  index_set_vec p(EMPTYSET, g.length());

  // actions not assigned to any partition
  bool_vec unass(true, instance.n_actions());
  // action pre-part sets remaining
  rem.assign_value(true, app.length());
  for (index_type k = 0; k < app.length(); k++)
    if (app[k].empty()) rem[k] = false;
  // partitions already known to be useless
  bool_vec useless(false, g.length());
  bool done = false;
  while (!done) { //  && (rem.count(true) > 0)
    // set-up:
    done = true;
    // h^1 counting all and only unassigned actions
    CostTable* Hu = new CostTable(instance, stats);
    DiscountACF d_cost(cost, instance.n_actions());
    d_cost.count_only(unass);
    Hu->compute_H1(d_cost);
    bool hu_is_useless = (Hu->max_finite() <= 0);
    // h^1 counting actions assigned to p[k] + unassigned actions
    lvector<CostTable*> Hp(0, g.length());
    for (index_type i = 0; i < g.length(); i++) if (!useless[i]) {
      if (p[i].empty()) {
	if (hu_is_useless)
	  useless[i] = true;
	else
	  Hp[i] = Hu;
      }
      else {
	d_cost.count_only(unass);
	d_cost.count(p[i]);
	Hp[i] = new CostTable(instance, stats);
	Hp[i]->compute_H1(d_cost);
	if (Hp[i]->max_finite() <= 0) {
	  useless[i] = true;
	  delete Hp[i];
	  Hp[i] = 0;
	}
      }
    }
    std::cerr << rem.count(true) << " unassigned action pre-partitions, "
	      << useless.count(false) << " non-useless partitions"
	      << std::endl;
    // main testing loop:
    for (index_type k = 0; k < app.length(); k++) if (rem[k]) {
      cost_vec loss(0, g.length());
      for (index_type i = 0; i < g.length(); i++) if (!useless[i]) {
	d_cost.count_only(unass);
	d_cost.count(p[i]);
	d_cost.discount(app[k]);
	CostTable* Hd = new CostTable(instance, stats);
	Hd->compute_H1(d_cost);
	loss[i] = (Hp[i]->eval(g[i]) - Hd->eval(g[i]));
	delete Hd;
      }
      index_type imax = loss.arg_max();
      assert(imax < g.length());
      NTYPE lmax = loss[imax];
      bool is_goal = g_prio.contains(g[imax]);
      if (lmax > 0) {
	// tie-breaking
	for (index_type i = imax + 1; i < g.length(); i++)
	  if (!useless[i] && (loss[i] == lmax)) {
	    // bool gi_is_goal = (!is_goal ? g_prio.contains(g[i]) : false);
	    // if (!is_goal && gi_is_goal) {
	    //   imax = i;
	    // }
	    // else
	    if (Hmax->eval(g[i]) < Hmax->eval(g[imax])) {
	      imax = i;
	    }
	    else if ((Hmax->eval(g[i]) == Hmax->eval(g[imax])) &&
		     (p[i].length() < p[imax].length())) {
	      imax = i;
	    }
	  }
	assert(imax < g.length());
	p[imax].insert(app[k]);
	rem[k] = false;
	unass.subtract(app[k]);
	// std::cerr << "set " << k << " assigned to goal " << imax << std::endl;
	done = false;
      }
    }
    // clean-up:
    for (index_type i = 0; i < g.length(); i++)
      if ((Hp[i] != 0) && (Hp[i] != Hu))
	delete Hp[i];
    delete Hu;
    for (index_type i = 0; i < g.length(); i++)
      if (useless[i] && !p[i].empty()) {
	for (index_type k = 0; k < app.length(); k++)
	  if (!rem[k] && !app[k].empty())
	    if (p[i].contains(app[k])) {
	      rem[k] = true;
	      unass.insert(app[k]);
	    }
	p[i].clear();
	done = false;
      }

    if (done && (rem.count(true) > 1)) {
      index_type k0 = rem.first(true);
      assert(k0 < app.length());
      assert(rem[k0]);
      for (index_type k = k0 + 1; k < app.length(); k++) if (rem[k]) {
	app[k0].insert(app[k]);
	app[k].clear();
	rem[k] = false;
      }
      done = false;
    }
  }
  std::cerr << "construction finished: "
	    << useless.count(false) << " non-useless partitions, "
	    << unass.count(true) << " actions remain unassigned"
	    << std::endl;

  compute_additive(cost, p, useH2);
  if (useH2)
    compute_max(cost, true);

  stats.stop();
}
#endif

void AH::compute_with_iterative_assignment_1
(const ACF& cost, const index_set& g,
 bool useH2, bool fractional, bool optimal,
 const index_set& g_prio)
{
  stats.start();

  index_set_vec app; // action pre-partition

  index_set v;
  for (index_type k = 0; k < instance.n_invariants(); k++)
    if ((instance.invariants[k].lim == 1) && instance.invariants[k].exact)
      v.insert(k);
  if (v.empty()) {
    std::cerr << "warning: no invariants usable for action partitioning"
	      << std::endl;
  }
  else {
    std::cerr << v.length() << " variables in problem..." << std::endl;

    std::cerr << "associating actions to variables..." << std::endl;
    index_set_vec a(v.length());
    for (index_type k = 0; k < instance.n_actions(); k++)
      if (!IS_ZERO(cost(k))) {
	index_set e(instance.actions[k].add);
	e.insert(instance.actions[k].del);
	for (index_type i = 0; i < v.length(); i++)
	  if (instance.invariants[v[i]].set.count_common(e) > 0)
	    a[i].insert(k);
      }

    std::cerr << "building graph of conflicts..." << std::endl;
    graph c(v.length());
    for (index_type i = 0; i < v.length(); i++)
      for (index_type j = i+1; j < v.length(); j++)
	if (a[i].count_common(a[j]))
	  c.add_undirected_edge(i, j);

    index_set s;
    if (optimal) {
      std::cerr << "searching for maximal conflict-free set..." << std::endl;
      c.complement();
      c.maximal_clique(s);
    }
    else {
      std::cerr << "searching for conflict-free set..." << std::endl;
      c.apx_independent_set(s);
    }
    std::cerr << s.length() << " variables in maximal additive set"
	      << std::endl;

    app.assign_select(a, s);
  }

  compute_with_iterative_assignment(cost, g, app, useH2, fractional, g_prio);
  stats.stop();
}


void AH::compute_with_iterative_assignment_2
(const ACF& cost, const index_set& g,
 bool useH2, bool fractional, const index_set& g_prio)
{
  stats.start();

  index_set_vec app; // action pre-partition

#ifdef ROLLBACK
  graph cag;
  instance.coadd_graph(cag);
  index_set s;
  cag.apx_independent_set(s);
  app.assign_value(EMPTYSET, s.length());
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (!IS_ZERO(cost(k))) {
      index_pair ij = s.first_common(instance.actions[k].add);
      if (ij.first != no_such_index)
	app[ij.first].insert(k);
    }
#else
  graph cag(instance.n_atoms());
  bool_vec ra(false, instance.n_actions());
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (cost(k) > 0) {
      for (index_type i = 0; i < instance.actions[k].add.length(); i++)
	for (index_type j = i + 1; j < instance.actions[k].add.length(); j++)
	  cag.add_undirected_edge(instance.actions[k].add[i],
				  instance.actions[k].add[j]);
      ra[k] = true;
    }
  index_set rn;
  rn.fill(instance.n_atoms());

  while (!rn.empty()) {
    index_set s;
    cag.apx_independent_set(rn, s);
    index_set_vec new_app(EMPTYSET, s.length());
    for (index_type k = 0; k < instance.n_actions(); k++) if (ra[k]) {
      index_pair ij = s.first_common(instance.actions[k].add);
      if (ij.first != no_such_index) {
	new_app[ij.first].insert(k);
	ra[k] = false;
      }
    }
    for (index_type k = 0; k < new_app.length(); k++)
      if (!new_app[k].empty()) app.append(new_app[k]);
    rn.subtract(s);
    for (index_type i = 0; i < rn.length(); i++)
      for (index_type j = i + 1; j < rn.length(); j++)
	if (cag.adjacent(i, j)) {
	  index_pair p = instance.atoms[rn[i]].add_by.first_common(instance.atoms[rn[j]].add_by);
	  bool edge_rem = false;
	  while ((p.first != no_such_index) && edge_rem) {
	    if (ra[instance.atoms[rn[i]].add_by[p.first]])
	      edge_rem = true;
	    else
	      p = instance.atoms[rn[i]].add_by.next_common(instance.atoms[rn[j]].add_by, p);
	  }
	  if (!edge_rem)
	    cag.remove_undirected_edge(rn[i], rn[j]);
	}
  }

  for (index_type k = 0; k < instance.n_actions(); k++)
    if (ra[k] && !instance.actions[k].add.empty())
      std::cerr << "oops: action " << instance.actions[k].name
		<< " remains after loop"
		<< std::endl;

  for (index_type i = 0; i < app.length(); i++)
    for (index_type j = i + 1; j < app.length(); j++)
      if (app[i].have_common_element(app[j]))
	std::cerr << "oops: action partitioning is not a partitioning!"
		  << std::endl;
#endif
  compute_with_iterative_assignment(cost, g, app, useH2, fractional, g_prio);
  stats.stop();
}


void AH::compute_with_new_decomposition
(const ACF& cost, const index_set& g, bool useH2)
{
  stats.start();

  index_set_graph cg;
  bool_vec  a(true, instance.n_atoms());
  for (index_type k = 0; k < instance.n_invariants(); k++)
    if (instance.invariants[k].lim == 1) {
      cg.add_node();
      cg.node_label(cg.size() - 1).assign_copy(instance.invariants[k].set);
      a.subtract(instance.invariants[k].set);
    }
  for (index_type k = 0; k < instance.n_atoms(); k++)
    if (a[k]) {
      cg.add_node();
      cg.node_label(cg.size() - 1).assign_singleton(k);
    }
  std::cerr << cg.size() << " state variables" << std::endl;

  for (index_type k = 0; k < instance.n_actions(); k++) {
    index_set ps;
    for (index_type i = 0; i < cg.size(); i++)
      if (instance.actions[k].pre.have_common_element(cg.node_label(i)))
	ps.insert(i);
    for (index_type i = 0; i < cg.size(); i++)
      if (instance.actions[k].add.have_common_element(cg.node_label(i)))
	cg.add_edge(ps, i);
  }

  // cg.write_digraph(std::cerr, "CG");

  index_set_graph ccg(cg);
  ccg.merge_labels_downwards();

  // ccg.write_digraph(std::cerr, "CCG");

  index_set gv;
  for (index_type k = 0; k < cg.size(); k++)
    if (g.have_common_element(cg.node_label(k)))
      gv.insert(k);

  // std::cerr << "gv = " << gv << std::endl;

  index_set_vec p(EMPTYSET, gv.length());
  for (index_type k = 0; k < instance.n_actions(); k++) {
    //bool_vec aff(false, gv.length());
    //for (index_type i = 0; i < gv.length(); i++)
    //  if (instance.actions[k].add.have_common_element(cg.node_label(gv[i])))
    //aff[i] = true;
    //if (aff.count(true)) {
    //  for (index_type i = 0; i < gv.length(); i++)
    //    if (aff[i]) p[i].insert(k);
    //}
    //else {
    for (index_type i = 0; i < gv.length(); i++)
      if (instance.actions[k].add.have_common_element(ccg.node_label(gv[i])))
	p[i].insert(k);
    //}
  }

  if (trace_level > 0) {
    for (index_type k = 0; k < gv.length(); k++) {
      std::cerr << "goal variable #" << k + 1 << ": ";
      instance.write_atom_set(std::cerr, cg.node_label(gv[k]));
      std::cerr << std::endl << " actions: ";
      instance.write_action_set(std::cerr, p[k]);
      std::cerr << std::endl;
    }
  }

  compute_fractional(cost, p, useH2);
  stats.stop();
}


void AH::compute_bottom_up_2
(const ACF& cost, const index_set& g, index_type p_max, bool useH2)
{
  stats.start();
  compute_max(cost, true);
  std::cerr << "H2max computed in " << stats.time() << " seconds" << std::endl;

  PairIndexFunction pi(instance.n_atoms());
  graph prg(pi.max_out() + 1);
  for (index_type k = 0; k < instance.n_actions(); k++) {
    index_set pre_pairs;
    for (index_type i = 0; i < instance.actions[k].pre.length(); i++)
      for (index_type j = i; j < instance.actions[k].pre.length(); j++)
	pre_pairs.insert(pi(instance.actions[k].pre[i],
			    instance.actions[k].pre[j]));
    for (index_type i = 0; i < instance.n_atoms(); i++) {
      if (instance.actions[k].add.contains(i)) {
	prg.add_edge(pre_pairs, pi(i, i));
	for (index_type j = 0; j < instance.actions[k].add.length(); j++)
	  if (instance.actions[k].add[j] > i)
	    prg.add_edge(pre_pairs, pi(i, instance.actions[k].add[j]));
      }
      else if (!instance.actions[k].del.contains(i)) {
	if (FINITE(Hmax->incremental_eval(instance.actions[k].pre, i))) {
	  for (index_type j = 0; j < instance.actions[k].add.length(); j++) {
	    prg.add_edge(pre_pairs, pi(instance.actions[k].add[j], i));
	    if (!instance.actions[k].pre.contains(i))
	      for (index_type l = 0; l < instance.actions[k].pre.length(); l++)
		prg.add_edge(pi(instance.actions[k].pre[l], i),
			     pi(instance.actions[k].add[j], i));
	  }
	}
      }
    }
//     for (index_type i = 0; i < instance.actions[k].add.length(); i++) {
//       for (index_type j = i; j < instance.actions[k].add.length(); j++)
// 	prg.add_edge(pre_pairs, pi(instance.actions[k].add[i],
// 				   instance.actions[k].add[j]));
//       for (index_type j = 0; j < instance.n_atoms(); j++)
// 	if (!instance.actions[k].del.contains(j) &&
// 	    !instance.actions[k].add.contains(j) &&
// 	    FINITE(Hmax->incremental_eval(instance.actions[k].pre, j))) {
// 	  for (index_type l = 0; l < instance.actions[k].add.length(); l++) {
// 	    prg.add_edge(pre_pairs, pi(instance.actions[k].add[l], j));
// 	    for (index_type m = 0; m < instance.actions[k].pre.length(); m++)
// 	      prg.add_edge(pi(instance.actions[k].pre[m], j),
// 			   pi(instance.actions[k].add[l], j));
// 	  }
// 	}
//    }
  }
  std::cerr << "prg construction finished in " << stats.time()
	    << " seconds" << std::endl;

  bool_vec rp(false, prg.size());
  prg.ancestors(g, rp);

  std::cerr << "relevant pairs: " << rp.count(true)
	    << " of " << rp.length() << ", " << stats.time() << " seconds"
	    << std::endl;

  graph nmg(instance.n_atoms());
  for (index_type i = 0; i < instance.n_atoms(); i++)
    for (index_type j = i + 1; j < instance.n_atoms(); j++)
      if (rp[pi(i, j)])
	nmg.add_undirected_edge(i, j);

//   pair_vec q;
//   for (index_type i = 0; i < g.length(); i++)
//     for (index_type j = i + 1; j < g.length(); j++)
//       q.append(index_pair(g[i], g[j]));
//
//   while (q.length() > 0) {
//     index_pair p = q[0];
//     p.sort_ascending();
//     q.remove(0);
//     if (!nmg.adjacent(p.first, p.second) &&
// 	FINITE(Hmax->eval(p.first, p.second))) {
//       nmg.add_undirected_edge(p.first, p.second);
//       const index_set& adds1st = instance.atoms[p.first].add_by;
//       for (index_type k = 0; k < adds1st.length(); k++) {
// 	if (instance.actions[adds1st[k]].add.contains(p.second)) {
// 	  for (index_type i = 0; i < instance.actions[adds1st[k]].pre.length(); i++)
// 	    for (index_type j = i + 1; j < instance.actions[adds1st[k]].pre.length(); j++)
// 	      if (!nmg.adjacent(instance.actions[adds1st[k]].pre[i],
// 				instance.actions[adds1st[k]].pre[j]))
// 		q.append(index_pair(instance.actions[adds1st[k]].pre[i],
// 				    instance.actions[adds1st[k]].pre[j]));
// 	}
// 	else if (!instance.actions[adds1st[k]].del.contains(p.second)) {
// 	  for (index_type i = 0; i < instance.actions[adds1st[k]].pre.length(); i++) {
// 	    for (index_type j = i + 1; j < instance.actions[adds1st[k]].pre.length(); j++)
// 	      if (!nmg.adjacent(instance.actions[adds1st[k]].pre[i],
// 				instance.actions[adds1st[k]].pre[j]))
// 		q.append(index_pair(instance.actions[adds1st[k]].pre[i],
// 				    instance.actions[adds1st[k]].pre[j]));
// 	    if (!nmg.adjacent(instance.actions[adds1st[k]].pre[i], p.second))
// 	      q.append(index_pair(instance.actions[adds1st[k]].pre[i],
// 				  p.second));
// 	  }
// 	}
//       }
//       const index_set& adds2nd = instance.atoms[p.second].add_by;
//       for (index_type k = 0; k < adds2nd.length(); k++) {
// 	if (!instance.actions[adds2nd[k]].del.contains(p.first) &&
// 	    !instance.actions[adds2nd[k]].add.contains(p.first)) {
// 	  for (index_type i = 0; i < instance.actions[adds2nd[k]].pre.length(); i++) {
// 	    for (index_type j = i + 1; j < instance.actions[adds2nd[k]].pre.length(); j++)
// 	      if (!nmg.adjacent(instance.actions[adds2nd[k]].pre[i],
// 				instance.actions[adds2nd[k]].pre[j]))
// 		q.append(index_pair(instance.actions[adds2nd[k]].pre[i],
// 				    instance.actions[adds2nd[k]].pre[j]));
// 	    if (!nmg.adjacent(instance.actions[adds2nd[k]].pre[i], p.first))
// 	      q.append(index_pair(instance.actions[adds2nd[k]].pre[i],
// 				  p.first));
// 	  }
// 	}
//       }
//     }
//   }

  std::cerr << "computed graph: " << nmg.n_edges() / 2
	    << " of " << (nmg.size() * (nmg.size() - 1)) / 2
	    << " edges present, " << stats.time() << " seconds"
	    << std::endl;

  index_set_vec p_atms;
  nmg.apx_independent_set_disjoint_cover(p_atms);
  std::cerr << "cover by " << p_atms.length() << " independent sets"
	    << std::endl;

  if (p_atms.length() > p_max) {
    cost_matrix mc(0, p_atms.length(), p_atms.length());
    for (index_type i = 0; i < p_atms.length(); i++)
      for (index_type j = i + 1; j < p_atms.length(); j++)
	mc[i][j] = nmg.n_induced_undirected_edges(p_atms[i], p_atms[j]);

    while (p_atms.length() > p_max) {
      index_type last = p_atms.length() - 1;
      index_type i_min = 0;
      for (index_type i = 1; i < last; i++)
	if (mc[i][last] < mc[i_min][last])
	  i_min = i;
      p_atms[i_min].insert(p_atms[last]);
      p_atms.dec_length();
      for (index_type i = 0; i < p_atms.length(); i++)
	if (i < i_min)
	  mc[i][i_min] =
	    nmg.n_induced_undirected_edges(p_atms[i], p_atms[i_min]);
	else if (i > i_min)
	  mc[i_min][i] =
	    nmg.n_induced_undirected_edges(p_atms[i], p_atms[i_min]);
    }
  }

  index_set_vec p_acts(EMPTYSET, p_atms.length());

  for (index_type k = 0; k < p_atms.length(); k++)
    for (index_type i = 0; i < p_atms[k].length(); i++)
      p_acts[k].insert(instance.atoms[p_atms[k][i]].add_by);

  compute_fractional(cost, p_acts, useH2);
  stats.stop();
}

void AH::compute_bottom_up
(const ACF& cost, const index_set& g, bool do_glue, bool useH2)
{
  stats.start();

  index_set_graph cg;
  bool_vec  rem(true, instance.n_atoms());
  for (index_type k = 0; k < instance.n_invariants(); k++)
    if (instance.invariants[k].lim == 1) {
      cg.add_node();
      cg.node_label(cg.size() - 1).assign_copy(instance.invariants[k].set);
      rem.subtract(instance.invariants[k].set);
    }
  for (index_type k = 0; k < instance.n_atoms(); k++)
    if (rem[k]) {
      cg.add_node();
      cg.node_label(cg.size() - 1).assign_singleton(k);
    }
  std::cerr << cg.size() << " state variables" << std::endl;

  for (index_type k = 0; k < instance.n_actions(); k++) {
    index_set ps;
    for (index_type i = 0; i < cg.size(); i++)
      if (instance.actions[k].pre.have_common_element(cg.node_label(i)))
	ps.insert(i);
    for (index_type i = 0; i < cg.size(); i++)
      if (instance.actions[k].add.have_common_element(cg.node_label(i)))
	cg.add_edge(ps, i);
  }

  equivalence p_eq(instance.n_atoms());

  rem.assign_value(true, instance.n_atoms());
  for (index_type k = 0; k < cg.size(); k++)
    if (cg.out_degree(k) == 0) {
      p_eq.merge(cg.node_label(k));
      rem.subtract(cg.node_label(k));
    }

  std::cerr << p_eq.n_classes() << " initial partitions" << std::endl;

  compute_max(cost, true);
  graph hg0;
  Hmax->compute_heuristic_graph(cost, hg0);

  if (do_glue) {
    bool done = false;
    while (!done) {
      done = true;
      index_set_graph hg(hg0, p_eq);
      index_set gn;
      for (index_type i = 0; i < hg.size(); i++)
	if (hg.node_label(i).have_common_element(instance.goal_atoms))
	  gn.insert(i);
      bool_vec rn;
      hg.ancestors(gn, rn);
      for (index_type i = 0; i < hg.size(); i++) {
	if (rn[i]) {
	  index_type nc = hg.predecessors(i).count_common(rn);
	  if (nc == 1) {
	    index_type fc = hg.predecessors(i).first_common_element(rn);
	    p_eq.merge(hg.node_label(i), hg.node_label(fc));
	    done = false;
	  }
	}
	else {
	  if (hg.in_degree(i) == 1) {
	    index_type j = hg.predecessors(i)[0];
	    p_eq.merge(hg.node_label(i), hg.node_label(j));
	    done = false;
	  }
	}
      }
      std::cerr << p_eq.n_classes() << " partitions after glueing"
		<< std::endl;
    }
  }

  index_set_vec p_atms;
  p_eq.classes(p_atms);
  // bool_vec useless(false, p_atms.length());
  // for (index_type k = 0; k < p_atms.length(); k++)
  //   if (instance.init_atoms.contains(p_atms[k]))
  //     useless[k] = true;
  // p_atms.remove(useless);
  // std::cerr << p_atms.length() << " final partitions" << std::endl;

  index_set_vec p_acts(EMPTYSET, p_atms.length());
  for (index_type k = 0; k < p_atms.length(); k++)
    for (index_type i = 0; i < p_atms[k].length(); i++)
      p_acts[k].insert(instance.atoms[p_atms[k][i]].add_by);

  compute_fractional(cost, p_acts, useH2);
  stats.stop();
}

void AH::compute_with_k_cuts
(const ACF& cost, const index_set& g, index_type k, bool useH2)
{
  stats.start();
  weighted_graph lg(instance.n_actions());
  for (index_type i = 0; i < instance.n_actions(); i++)
    for (index_type j = 0; j < instance.n_actions(); j++)
      if ((i != j) &&
	  instance.actions[i].add.have_common_element(instance.actions[j].add))
	lg.add_undirected_edge(i, j, 1);

  equivalence a_eq;
  lg.induced_partitioning(a_eq);
  index_type b = a_eq.n_classes();
  std::cerr << b << " initial partitions" << std::endl;

  if (b < k) {
    weighted_vec<index_pair, NTYPE> cuts;
    for (index_type i = 0; i < instance.n_actions(); i++)
      for (index_type j = i+1; j < instance.n_actions(); j++)
	if (lg.adjacent(i, j)) {
	  NTYPE c = lg.max_flow(i, j); // == min cut cost
	  std::cerr << "cut " << i << "-" << j << ": " << c << std::endl;
	  cuts.insert_increasing(index_pair(i, j), c);
	}

    for (index_type i = 0; i < cuts.length(); i++)
      std::cerr << cuts[i].value << ", cost = " << cuts[i].weight
		<< std::endl;

    pair_set e_all;
    index_type i = 0;
    while (b < k) {
      assert(i < cuts.length());
      index_type s = cuts[i].value.first;
      index_type t = cuts[i].value.second;
      assert(a_eq.canonical(s) == a_eq.canonical(t));
      pair_set e_cut;
      lg.min_cut(s, t, e_cut);
      if (!e_all.contains(e_cut)) {
	b += 1;
	e_all.insert(e_cut);
	std::cerr << "selected cut " << e_cut
		  << " with cost " << cuts[i].weight << std::endl;
      }
      i += 1;
    }

    lg.remove_undirected_edges(e_all);
    lg.induced_partitioning(a_eq);
    index_type b = a_eq.n_classes();
    std::cerr << b << " final partitions" << std::endl;
  }

  index_set_vec p_acts;
  a_eq.classes(p_acts);
  compute_additive(cost, p_acts, useH2);
  stats.stop();
}

void AH::compute_with_layered_partitioning
(const ACF& cost, const index_set& g)
{
  stats.start();

  std::cerr << "computing H2..." << std::endl;
  Hmax = new CostTable(instance, stats);
  Hmax->compute_H2(cost);

  index_set_vec p; // action partition
  bool_vec rem(true, instance.n_actions());

  std::cerr << "computing layered action partition..." << std::endl;
  compute_layered_action_partition(2, g, p, no_such_index, rem);

  std::cerr << "computing additive H2 ("
	    << p.length() << " components)..."
	    << std::endl;
  compute_additive(cost, p, true);
  stats.stop();
}

void AH::max_cost_mset
(index_type m, const index_set& g, index_set& s)
{
  assert(Hmax); // Hmax must be computed/set

  NTYPE v_max = NEG_INF;
  mSubsetEnumerator se(g.length(), m);
  bool more = se.first();
  while (more) {
    index_set cs;
    se.current_set(g, cs);
    if (Hmax->eval(cs) > v_max) {
      v_max = Hmax->eval(cs);
      s.assign_copy(cs);
    }
    more = se.next();
  }
}

void AH::goal_relevant_remaining_actions
(const index_set& g, index_set& a, bool_vec& rem)
{
  a.clear();
  for (index_type k = 0; k < instance.n_actions(); k++)
    if (rem[k]) {
      const Instance::Action& act = instance.actions[k];
      if ((act.add.count_common(g) > 0) && (act.del.count_common(g) == 0))
	a.insert(k);
    }
}

void AH::compute_layered_action_partition
(index_type m, const index_set& g, index_set_vec& p, index_type pg,
 bool_vec& rem)
{
  // split g into size m partitions
  std::cerr << "splitting " << g << "..." << std::endl;
  index_set g0(g);
  index_set_vec sub_g;
  while (!g0.empty()) {
    index_set next_p;
    max_cost_mset(m, g0, next_p);
    sub_g.append(next_p);
    g0.subtract(next_p);
  }
  std::cerr << "result: " << sub_g << std::endl;
  // assign immediately relevant actions to each goal partition
  index_vec sub_p(no_such_index, sub_g.length());
  index_set_vec new_g(sub_g.length());
  for (index_type k = 0; k < sub_g.length(); k++) {
    index_set rra;
    goal_relevant_remaining_actions(sub_g[k], rra, rem);
    std::cerr << "goal subset " << k << " = " << sub_g[k]
	      << ", rra = " << rra << std::endl;
    if (!rra.empty()) {
      // if pg != no_such_index, assign rra to p[pg]
      if (pg != no_such_index) {
	p[pg].insert(rra);
	sub_p[k] = pg;
	pg = no_such_index;
      }
      // otherwise, create a new partition
      else {
	p.inc_length();
	p[p.length() - 1].insert(rra);
	sub_p[k] = p.length() - 1;
      }
      // remove rra's from rem, and add precs to new_g
      for (index_type i = 0; i < rra.length(); i++) {
	rem[rra[i]] = false;
	new_g[k].insert(instance.actions[rra[i]].pre);
      }
    }
  }
  std::cerr << "action partitions: " << p << std::endl;
  // recurse with collected preconditions of assigned actions
  std::cerr << "new_g = " << new_g << std::endl;
  for (index_type k = 0; k < sub_g.length(); k++) if (!new_g[k].empty()) {
    compute_layered_action_partition(m, new_g[k], p, sub_p[k], rem);
  }
}

NTYPE AH::eval(const index_set& s)
{
  NTYPE v_max = (Hmax ? Hmax->eval(s) : 0);
  if (INFINITE(v_max)) {
#ifdef AH_EXTRA_STATS
    // draws += 1;
#endif
    return v_max;
  }
  NTYPE v_sum = 0;
  if (trace_level > 2) {
    std::cerr << "AH: v_max = " << v_max << ", v_sum = ";
  }
  for (index_type k = 0; k < H_vec.length(); k++) {
    assert(H_vec[k]);
    v_sum += H_vec[k]->eval(s);
    if (trace_level > 2) {
      std::cerr << " + " << H_vec[k]->eval(s);
    }
  }
  if (trace_level > 2) {
    std::cerr << " = " << v_sum << std::endl;
  }
#ifdef PRINT_EXTRA_STATS
  std::cerr << "AH: " << v_max << " " << v_sum << std::endl;
#endif
#ifdef AH_EXTRA_STATS
  if (v_sum > v_max) {
    Hsum_wins += 1;
    return v_sum;
  }
  else if (v_max > v_sum) {
    Hmax_wins += 1;
    return v_max;
  }
  else {
    draws += 1;
    return v_max;
  }
#else
  return MAX(v_max, v_sum);
#endif
}

NTYPE AH::eval(const bool_vec& s)
{
  NTYPE v_max = (Hmax ? Hmax->eval(s) : 0);
  if (INFINITE(v_max)) {
#ifdef AH_EXTRA_STATS
    // draws += 1;
#endif
    return v_max;
  }
  NTYPE v_sum = 0;
  if (trace_level > 2) {
    std::cerr << "AH: v_max = " << v_max << ", v_sum = ";
  }
  for (index_type k = 0; k < H_vec.length(); k++) {
    assert(H_vec[k]);
    v_sum += H_vec[k]->eval(s);
    if (trace_level > 2) {
      std::cerr << " + " << H_vec[k]->eval(s);
    }
  }
  if (trace_level > 2) {
    std::cerr << " = " << v_sum << std::endl;
  }
#ifdef PRINT_EXTRA_STATS
  std::cerr << "AH: " << v_max << " " << v_sum << std::endl;
#endif
#ifdef AH_EXTRA_STATS
  if (v_sum > v_max) {
    Hsum_wins += 1;
    return v_sum;
  }
  else if (v_max > v_sum) {
    Hmax_wins += 1;
    return v_max;
  }
  else {
    draws += 1;
    return v_max;
  }
#else
  return MAX(v_max, v_sum);
#endif
}

void AH::write_eval(const index_set& s, std::ostream& st, char* p, bool e)
{
  NTYPE v_max = (Hmax ? Hmax->eval(s) : 0);
  NTYPE v_sum = 0;
  if (p) st << p << " (";
  for (index_type k = 0; k < H_vec.length(); k++) {
    assert(H_vec[k]);
    v_sum += H_vec[k]->eval(s);
    if (k > 0) st << " ";
    st << H_vec[k]->eval(s);
  }
  st << ") " << v_sum << " " << v_max;
  if (e) st << std::endl;
}

void AH::write_eval(const bool_vec& s, std::ostream& st, char* p, bool e)
{
  NTYPE v_max = (Hmax ? Hmax->eval(s) : 0);
  NTYPE v_sum = 0;
  if (p) st << p << " (";
  for (index_type k = 0; k < H_vec.length(); k++) {
    assert(H_vec[k]);
    v_sum += H_vec[k]->eval(s);
    if (k > 0) st << " ";
    st << H_vec[k]->eval(s);
  }
  st << ") " << v_sum << " " << v_max;
  if (e) st << std::endl;
}

NTYPE AH::incremental_eval(const index_set& s, index_type i_new)
{
  NTYPE v_max = (Hmax ? Hmax->incremental_eval(s, i_new) : 0);
  NTYPE v_sum = 0;
  for (index_type k = 0; k < H_vec.length(); k++) {
    assert(H_vec[k]);
    v_sum += H_vec[k]->incremental_eval(s, i_new);
  }
#ifdef AH_EXTRA_STATS
  if (v_sum > v_max) {
    Hsum_wins += 1;
    return v_sum;
  }
  else if (v_max > v_sum) {
    Hmax_wins += 1;
    return v_max;
  }
  else {
    draws += 1;
    return v_max;
  }
#else
  return MAX(v_max, v_sum);
#endif
}

NTYPE AH::incremental_eval(const bool_vec& s, index_type i_new)
{
  NTYPE v_max = (Hmax ? Hmax->incremental_eval(s, i_new) : 0);
  NTYPE v_sum = 0;
  for (index_type k = 0; k < H_vec.length(); k++) {
    assert(H_vec[k]);
    v_sum += H_vec[k]->incremental_eval(s, i_new);
  }
#ifdef AH_EXTRA_STATS
  if (v_sum > v_max) {
    Hsum_wins += 1;
    return v_sum;
  }
  else if (v_max > v_sum) {
    Hmax_wins += 1;
    return v_max;
  }
  else {
    draws += 1;
    return v_max;
  }
#else
  return MAX(v_max, v_sum);
#endif
}

NTYPE AH::eval_to_bound(const index_set& s, NTYPE bound)
{
  NTYPE v_max = (Hmax ? Hmax->eval_to_bound(s, bound) : 0);
  NTYPE v_sum = 0;
  for (index_type k = 0; (k < H_vec.length()) && (v_sum < bound); k++) {
    assert(H_vec[k]);
    v_sum += H_vec[k]->eval(s);
  }
#ifdef AH_EXTRA_STATS
  if (v_sum > v_max) {
    Hsum_wins += 1;
    return v_sum;
  }
  else if (v_max > v_sum) {
    Hmax_wins += 1;
    return v_max;
  }
  else {
    draws += 1;
    return v_max;
  }
#else
  return MAX(v_max, v_sum);
#endif
}

NTYPE AH::eval_to_bound(const bool_vec& s, NTYPE bound)
{
  NTYPE v_max = (Hmax ? Hmax->eval_to_bound(s, bound) : 0);
  NTYPE v_sum = 0;
  for (index_type k = 0; (k < H_vec.length()) && (v_sum < bound); k++) {
    assert(H_vec[k]);
    v_sum += H_vec[k]->eval(s);
  }
#ifdef AH_EXTRA_STATS
  if (v_sum > v_max) {
    Hsum_wins += 1;
    return v_sum;
  }
  else if (v_max > v_sum) {
    Hmax_wins += 1;
    return v_max;
  }
  else {
    draws += 1;
    return v_max;
  }
#else
  return MAX(v_max, v_sum);
#endif
}

void AH::store(const index_set& s, NTYPE v, bool opt)
{
  if (Hmax) Hmax->store(s, v, opt);
}

void AH::store(const bool_vec& s, NTYPE v, bool opt)
{
  if (Hmax) Hmax->store(s, v, opt);
}

END_HSPS_NAMESPACE
