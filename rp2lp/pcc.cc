
#include "pcc.h"

BEGIN_HSPS_NAMESPACE

void ppc::write(std::ostream& s, const Instance& ins) const
{
  s << "(preference ";
  if (name)
    s << name << " ";
  write_PDDL_constraint(s, ins);
  s << ") ;; weight = " << PRINT_NTYPE(weight) << std::endl;
}

void ppc::write_PDDL_constraint(std::ostream& s, const Instance& ins) const
{
  switch (pct) {
    // unary constraints
  case pc_at_end:
  case pc_always:
  case pc_sometime:
  case pc_at_most_once:
    switch (pct) {
    case pc_at_end:
      s << "(at end";
      break;
    case pc_always:
      s << "(always";
      break;
    case pc_sometime:
      s << "(sometime";
      break;
    case pc_at_most_once:
      s << "(at-most-once";
      break;
    }
    if (s_c.size() > 1) {
      if (c_dis)
	s << " (or";
      else 
	s << " (and";
    }
    for (index_type i = 0; i < s_c.size(); i++)
      s << " " << ins.atoms[s_c[i]].name;
    if (s_c.size() > 1) s << ")";
    s << ")";
    break;
    // binary constraints
  case pc_sometime_before:
  case pc_sometime_after:
    if (pct == pc_sometime_before)
      s << "(sometime-before";
    else
      s << "(sometime-after";
    if (s_t.size() > 1) {
      if (t_dis)
	s << " (or";
      else
	s << " (and";
    }
    for (index_type i = 0; i < s_t.size(); i++)
      s << " " << ins.atoms[s_t[i]].name;
    if (s_t.size() > 1) s << ")";
    if (s_c.size() > 1) {
      if (c_dis)
	s << " (or";
      else
	s << " (and";
    }
    for (index_type i = 0; i < s_c.size(); i++)
      s << " " << ins.atoms[s_c[i]].name;
    if (s_c.size() > 1) s << ")";
    s << ")";
    break;
  default:
    s << "<<unprintable plan constraint type " << pct << ">>";
  }
}

void ppc::print(std::ostream& s, const Instance& ins) const
{
  switch (pct) {
  case pc_at_end:
    assert(s_c.length() == 1);
    s << "AtEnd(";
    ins.atoms[s_c[0]].name->write(s, Name::NC_INSTANCE);
    s << ")";
    break;
  case pc_always:
    assert(s_c.length() == 1);
    s << "Always(";
    ins.atoms[s_c[0]].name->write(s, Name::NC_INSTANCE);
    s << ")";
    break;
  case pc_sometime:
    assert(s_c.length() == 1);
    s << "Sometime(";
    ins.atoms[s_c[0]].name->write(s, Name::NC_INSTANCE);
    s << ")";
    break;
  case pc_at_most_once:
    assert(s_c.length() == 1);
    s << "AtMostOnce(";
    ins.atoms[s_c[0]].name->write(s, Name::NC_INSTANCE);
    s << ")";
    break;
  case pc_sometime_before:
    assert(s_t.length() == 1);
    assert(s_c.length() == 1);
    s << "SometimeBefore(";
    ins.atoms[s_t[0]].name->write(s, Name::NC_INSTANCE);
    s << ",";
    ins.atoms[s_c[0]].name->write(s, Name::NC_INSTANCE);
    s << ")";
    break;
  case pc_sometime_after:
    assert(s_t.length() == 1);
    assert(s_c.length() == 1);
    s << "SometimeAfter(";
    ins.atoms[s_t[0]].name->write(s, Name::NC_INSTANCE);
    s << ",";
    ins.atoms[s_c[0]].name->write(s, Name::NC_INSTANCE);
    s << ")";
    break;
  default:
    s << "<<unprintable plan constraint type " << pct << ">>";
  }
}

void ppc::enforce(Instance& ins, index_vec& map) const
{
  if (pct == pc_at_end) {
    if (c_dis)
      ins.enforce_disjunctive_goal(s_c, name, map);
    else
      for (index_type i = 0; i < s_c.size(); i++)
	ins.atoms[s_c[i]].goal = true;
  }
  else if (pct == pc_always) {
    if (c_dis)
      ins.enforce_pc_always_disjunction(s_c, name, map);
    else
      ins.enforce_pc_always_conjunction(s_c, name, map);
  }
  else if (pct == pc_sometime) {
    if (c_dis)
      ins.enforce_pc_sometime_disjunction(s_c, name);
    else
      ins.enforce_pc_sometime_conjunction(s_c, name, map);
  }
  else if (pct == pc_at_most_once) {
    if (c_dis)
      ins.enforce_pc_at_most_once_disjunction(s_c, name, map);
    else
      ins.enforce_pc_at_most_once_conjunction(s_c, name, map);
  }
  else if ((pct == pc_sometime_before) && !c_dis && !t_dis) {
    ins.enforce_pc_sometime_before_cc(s_t, s_c, name, map);
  }
  else if ((pct == pc_sometime_after) &&
	   (s_t.size() == 1) &&
	   (s_c.size() == 1)) {
    index_type g = ins.compile_pc_sometime_after(s_t[0], s_c[0], name, &map);
    ins.atoms[g].goal = true;
  }
  else {
    std::cerr << "error: can't enforce plan constraint ";
    write_PDDL_constraint(std::cerr, ins);
    std::cerr << std::endl;
    exit(255);
  }
  assert(map.size() == ins.n_actions());
}

index_type ppc::compile(Instance& ins, index_vec& map) const
{
  if (pct == pc_always) {
    if (c_dis)
      return ins.compile_pc_always_disjunction(s_c, name, &map);
    else
      return ins.compile_pc_always_conjunction(s_c, name);
  }
  else if (pct == pc_sometime) {
    if (c_dis)
      return ins.compile_pc_sometime_disjunction(s_c, name);
    else
      return ins.compile_pc_sometime_conjunction(s_c, name, &map);
  }
  else if (pct == pc_at_most_once) {
    if (c_dis)
      return ins.compile_pc_at_most_once_disjunction(s_c, name);
    else
      return ins.compile_pc_at_most_once_conjunction(s_c, name);
  }
  else if ((pct == pc_sometime_before) && !c_dis && !t_dis) {
    return ins.compile_pc_sometime_before_cc(s_t, s_c, name);
  }
  else if ((pct == pc_sometime_after) &&
	   (s_t.size() == 1) &&
	   (s_c.size() == 1)) {
    return ins.compile_pc_sometime_after(s_t[0], s_c[0], name);
  }
  else {
    std::cerr << "error: can't compile plan constraint ";
    write_PDDL_constraint(std::cerr, ins);
    std::cerr << std::endl;
    exit(255);
  }
}

bool ppc::test(ExecTrace* t) const
{
  if (pct == pc_at_end) {
    if (c_dis)
      return t->final_state()->test_disjunction(s_c);
    else
      return t->final_state()->test_conjunction(s_c);
  }
  else if (pct == pc_always) {
    if (c_dis)
      return t->test_always_disjunction(s_c);
    else
      return t->test_always_conjunction(s_c);
  }
  else if (pct == pc_sometime) {
    if (c_dis)
      return t->test_sometime_disjunction(s_c);
    else
      return t->test_sometime_conjunction(s_c);
  }
  else if (pct == pc_at_most_once) {
    if (c_dis)
      return t->test_at_most_once_disjunction(s_c);
    else
      return t->test_at_most_once_conjunction(s_c);
  }
  else if ((pct == pc_sometime_before) &&
	   (s_c.length() == 1) && (s_t.length() == 1)) {
    return t->test_sometime_before(s_t[0], s_c[0]);
  }
  else if ((pct == pc_sometime_after) &&
	   (s_c.length() == 1) && (s_t.length() == 1)) {
    return t->test_sometime_after(s_t[0], s_c[0]);
  }
  else {
    std::cerr << "error: can't test ";
    write_PDDL_constraint(std::cerr, t->get_instance());
    std::cerr << std::endl;
    exit(255);
  }
}

// bool extract_ppc
// (PDDL_Base* b,
//  PDDL_Base::Preference* p,
//  Name* n,
//  PDDL_Base::Goal* g1,
//  PDDL_Base::Goal* g2,
//  Instance& ins,
//  plan_constraint_type t,
//  index_type k,
//  ppc_vec& ppcv)
// {
//   PDDL_Base::signed_atom_vec a;
//   bool g1_is_dis;
//   index_set s1;
//   if (b->goal_to_atom_vec(g1, a, g1_is_dis)) {
//     if (!g1_is_dis) {
//       Name* n = (p->name ?
// 		 (Name*)new StringName(p->name->print_name) :
// 		 (Name*)new EnumName("ppc", k));
//       NTYPE w = (b->metric ? p->value(b->metric_type, b->metric) : 1);
//       b->instantiate_atom_set(ins, a, s1);
//       for (index_type i = 0; i < s1.length(); i++)
// 	ins.atoms[s1[i]].goal = true;
//       if (g2) {
// 	PDDL_Base::signed_atom_vec a2;
// 	bool g2_is_dis;
// 	index_set s2;
// 	if (b->goal_to_atom_vec(g2, a2, g2_is_dis)) {
// 	  if (!g2_is_dis) {
// 	    b->instantiate_atom_set(ins, a2, s2);
// 	    for (index_type i = 0; i < s2.length(); i++)
// 	      ins.atoms[s2[i]].goal = true;
// 	    ppcv.append(ppc(n, w, t, s1, s2, p));
// 	  }
// 	  else {
// 	    std::cerr << "error: can't handle goal ";
// 	    p->print(std::cerr);
// 	    std::cerr << std::endl;
// 	    return false;
// 	  }
// 	}
// 	else {
// 	  std::cerr << "failed to extract atoms from goal ";
// 	  g2->print(std::cerr);
// 	  std::cerr << std::endl;
// 	  return false;
// 	}
//       }
//       else {
// 	ppcv.append(ppc(n, w, t, s1, EMPTYSET, p));
//       }
//       return true;
//     }
//     else {
//       std::cerr << "error: can't handle goal ";
//       p->print(std::cerr);
//       std::cerr << std::endl;
//       return false;
//     }
//   }
//   else {
//     std::cerr << "failed to extract atoms from goal ";
//     g1->print(std::cerr);
//     std::cerr << std::endl;
//     return false;
//   }
// }

// bool extract_ppcs
// (PDDL_Base* b,
//  Instance& ins,
//  ppc_vec& ppcv)
// {
//   bool ok = true;
//   for (index_type k = 0; k < b->dom_preferences.size(); k++) {
//     if (b->dom_preferences[k]->goal->g_class == PDDL_Base::goal_always) {
//       if (!extract_ppc(b, b->dom_preferences[k],
// 		       ((PDDL_Base::SimpleSequenceGoal*)
// 			b->dom_preferences[k]->goal)->constraint,
// 		       0, ins, pc_always, k, ppcv))
// 	ok = false;
//     }
//     else if (b->dom_preferences[k]->goal->g_class
// == PDDL_Base::goal_sometime) {
//       if (!extract_ppc(b, b->dom_preferences[k],
// 		       ((PDDL_Base::SimpleSequenceGoal*)
// 			b->dom_preferences[k]->goal)->constraint,
// 		       0, ins, pc_sometime, k, ppcv))
// 	ok = false;
//     }
//     else if (b->dom_preferences[k]->goal->g_class
// 	     == PDDL_Base::goal_at_most_once) {
//       if (!extract_ppc(b, b->dom_preferences[k],
// 		       ((PDDL_Base::SimpleSequenceGoal*)
// 			b->dom_preferences[k]->goal)->constraint,
// 		       0, ins, pc_at_most_once, k, ppcv))
// 	ok = false;
//     }
//     else if (b->dom_preferences[k]->goal->g_class
// 	     == PDDL_Base::goal_sometime_before) {
//       if (!extract_ppc(b, b->dom_preferences[k],
// 		       ((PDDL_Base::TriggeredSequenceGoal*)
// 			b->dom_preferences[k]->goal)->constraint,
// 		       ((PDDL_Base::TriggeredSequenceGoal*)
// 			b->dom_preferences[k]->goal)->trigger,
// 		       ins, pc_sometime_before, k, ppcv))
// 	ok = false;
//     }
//     else if (b->dom_preferences[k]->goal->g_class
// 	     == PDDL_Base::goal_sometime_after) {
//       if (!extract_ppc(b, b->dom_preferences[k],
// 		       ((PDDL_Base::TriggeredSequenceGoal*)
// 			b->dom_preferences[k]->goal)->constraint,
// 		       ((PDDL_Base::TriggeredSequenceGoal*)
// 			b->dom_preferences[k]->goal)->trigger,
// 		       ins, pc_sometime_after, k, ppcv))
// 	ok = false;
//     }
//     else {
//       if (!extract_ppc(b, b->dom_preferences[k],
// b->dom_preferences[k]->goal, 0, ins, pc_at_end, k, ppcv))
// 	ok = false;
//     }
//   }
//   return ok;
// }

bool extract_atoms
(PDDL_Base* b,
 PDDL_Base::Goal* g,
 Instance& ins,
 index_set& s,
 bool& is_dis)
{
  PDDL_Base::signed_atom_vec a;
  if (b->goal_to_atom_vec(g, a, is_dis)) {
    assert(a.size() > 0);
    b->instantiate_atom_set(ins, a, s);
    for (index_type i = 0; i < s.length(); i++)
      ins.atoms[s[i]].goal = true;
    if (s.size() <= 1) is_dis = false;
    return true;
  }
  else {
    std::cerr << "failed to extract atoms from goal ";
    g->print(std::cerr);
    std::cerr << std::endl;
    return false;
  }
}

bool extract_ppc
(PDDL_Base* b,
 PDDL_Base::Preference* p,
 const Name* n,
 Instance& ins,
 ppc_vec& ppcv)
{
  // std::cerr << "extracting from ";
  // p->print(std::cerr);
  // std::cerr << "..." << std::endl;
  if (p->goal->g_class == PDDL_Base::goal_always) {
    PDDL_Base::Goal* g1 =
      ((PDDL_Base::SimpleSequenceGoal*)p->goal)->constraint;
    bool g1_is_dis;
    index_set s1;
    if (extract_atoms(b, g1, ins, s1, g1_is_dis)) {
      NTYPE w = (b->metric ? p->value(b->metric_type, b->metric) : 1);
      ppcv.append(ppc(n, w, pc_always, s1, g1_is_dis, p));
      return true;
    }
  }
  if (p->goal->g_class == PDDL_Base::goal_sometime) {
    PDDL_Base::Goal* g1 =
      ((PDDL_Base::SimpleSequenceGoal*)p->goal)->constraint;
    bool g1_is_dis;
    index_set s1;
    if (extract_atoms(b, g1, ins, s1, g1_is_dis)) {
      NTYPE w = (b->metric ? p->value(b->metric_type, b->metric) : 1);
      ppcv.append(ppc(n, w, pc_sometime, s1, g1_is_dis, p));
      return true;
    }
  }
  if (p->goal->g_class == PDDL_Base::goal_at_most_once) {
    PDDL_Base::Goal* g1 =
      ((PDDL_Base::SimpleSequenceGoal*)p->goal)->constraint;
    bool g1_is_dis;
    index_set s1;
    if (extract_atoms(b, g1, ins, s1, g1_is_dis)) {
      NTYPE w = (b->metric ? p->value(b->metric_type, b->metric) : 1);
      ppcv.append(ppc(n, w, pc_at_most_once, s1, g1_is_dis, p));
      return true;
    }
  }
  else if (p->goal->g_class == PDDL_Base::goal_sometime_before) {
    PDDL_Base::Goal* g1 =
      ((PDDL_Base::TriggeredSequenceGoal*)p->goal)->constraint;
    PDDL_Base::Goal* g2 =
      ((PDDL_Base::TriggeredSequenceGoal*)p->goal)->trigger;
    bool g1_is_dis;
    index_set s1;
    bool g2_is_dis;
    index_set s2;
    if (extract_atoms(b, g1, ins, s1, g1_is_dis) &&
	extract_atoms(b, g2, ins, s2, g2_is_dis)) {
      NTYPE w = (b->metric ? p->value(b->metric_type, b->metric) : 1);
      ppcv.append(ppc(n, w, pc_sometime_before,
		      s1, g1_is_dis, s2, g2_is_dis, p));
      return true;
    }
  }
  else if (p->goal->g_class == PDDL_Base::goal_sometime_after) {
    PDDL_Base::Goal* g1 =
      ((PDDL_Base::TriggeredSequenceGoal*)p->goal)->constraint;
    PDDL_Base::Goal* g2 =
      ((PDDL_Base::TriggeredSequenceGoal*)p->goal)->trigger;
    bool g1_is_dis;
    index_set s1;
    bool g2_is_dis;
    index_set s2;
    if (extract_atoms(b, g1, ins, s1, g1_is_dis) &&
	extract_atoms(b, g2, ins, s2, g2_is_dis)) {
      NTYPE w = (b->metric ? p->value(b->metric_type, b->metric) : 1);
      ppcv.append(ppc(n, w, pc_sometime_after,
		      s1, g1_is_dis, s2, g2_is_dis, p));
      return true;
    }
  }
  else {
    bool g1_is_dis;
    index_set s1;
    if (extract_atoms(b, p->goal, ins, s1, g1_is_dis)) {
      NTYPE w = (b->metric ? p->value(b->metric_type, b->metric) : 1);
      ppcv.append(ppc(n, w, pc_at_end, s1, g1_is_dis, p));
      return true;
    }
  }
  std::cerr << "error: can't handle goal ";
  p->print(std::cerr);
  std::cerr << std::endl;
  return false;
}

bool extract_ppc
(PDDL_Base* b,
 PDDL_Base::Preference* p,
 index_type& k,
 Instance& ins,
 ppc_vec& ppcv)
{
  if (p->param.size() > 0) {
    assert(p->name != NULL);
    name_vec nv;
    PDDL_Base::preference_vec pv;
    p->unquantify(0, nv, pv);
    assert(nv.size() == pv.size());
    bool ok = true;
    for (index_type k = 0; k < pv.size(); k++)
      if (!extract_ppc(b, pv[k],
		       new NameWithContext(nv[k], Name::NC_INSTANCE, 0),
		       ins, ppcv))
	ok = false;
    return ok;
  }
  else {
    const Name* n = (p->name ?
		     (Name*)new StringName(p->name->print_name) :
		     (Name*)new EnumName("ppc", k++));
    return extract_ppc(b, p, n, ins, ppcv);
  }
}

bool extract_ppcs
(PDDL_Base* b,
 Instance& ins,
 ppc_vec& ppcv)
{
  bool ok = true;
  index_type n = 0;
  for (index_type k = 0; k < b->dom_preferences.size(); k++) {
    if (!extract_ppc(b, b->dom_preferences[k], n, ins, ppcv))
      ok = false;
  }
  return ok;
}

/// virtual Propagator implementation

InferenceTracker::InferenceTracker(Propagator& p)
  : propagator(p), verbose(false), debug_mode(false)
{
  // done?
}

InferenceTracker::~InferenceTracker()
{
  // done
}

void InferenceTracker::mark_steps(index_type i, bool_vec& m) const
{
  assert(i < steps.size());
  if (!m[i]) {
    m[i] = true;
    for (index_type k = 0; k < steps[i].pre.size(); k++)
      mark_steps(steps[i].pre[k], m);
  }
}

void InferenceTracker::print_assertion
(const assertion& a, std::ostream& to) const
{
  switch (a.pc) {
  case pc_always:
    to << "A";
    propagator.print_state_formula(a.a, to);
    break;
  case pc_sometime:
    to << "E";
    propagator.print_state_formula(a.a, to);
    break;
  case pc_at_most_once:
    to << "O";
    propagator.print_state_formula(a.a, to);
    break;
  case pc_sometime_before:
    to << "(";
    propagator.print_state_formula(a.b, to);
    to << " SB ";
    propagator.print_state_formula(a.a, to);
    to << ")";
    break;
  case pc_sometime_after:
    to << "(";
    propagator.print_state_formula(a.b, to);
    to << " SA ";
    propagator.print_state_formula(a.a, to);
    to << ")";
    break;
  case pc_never:
    to << "N";
    propagator.print_state_formula(a.a, to);
    break;
  case pc_disallowed:
    to << "D" << a.a << "." << propagator.ins.actions[a.a].name;
    break;
  case pc_contradiction:
    to << "CONTRADICTION";
    break;
  default:
    to << "ERROR: " << a.pc << " NOT PRINTABLE";
  }
}

void InferenceTracker::print_step(index_type i, std::ostream& to) const
{
  assert(i < steps.size());
  to << "(" << i << ") ";
  print_assertion(steps[i].ass, to);
  if (steps[i].rule)
    to << " " << steps[i].rule << " " << steps[i].pre;
  else
    to << " ??? " << steps[i].pre;
  if (!steps[i].valid) to << " - INVALID";
  to << std::endl;
}

void InferenceTracker::print_proof(std::ostream& to) const
{
  print_proof(assertion(pc_contradiction), to);
}

void InferenceTracker::print_proof(const assertion& a, std::ostream& to) const
{
  std::map<assertion, index_type>::const_iterator i = index.find(a);
  if (i == index.end()) {
    print_assertion(a, to);
    to << ": NO PROOF!" << std::endl;
  }
  else {
    print_assertion(a, to);
    to << ": BEGIN PROOF" << std::endl;
    bool_vec used(false, steps.size());
    mark_steps(i->second, used);
    for (index_type k = 0; k < steps.size(); k++)
      if (used[k])
	print_step(k, to);
    to << "== END PROOF" << std::endl;
  }
}

void InferenceTracker::print_all(std::ostream& to) const
{
  to << "== BEGIN PROOF ==" << std::endl;
  for (index_type k = 0; k < steps.size(); k++)
    print_step(k, to);
  to << "== END PROOF ==" << std::endl;  
}

void InferenceTracker::extract_proof
(const assertion& a, std::vector<step>& p) const
{
  p.clear();
  std::map<assertion, index_type>::const_iterator i = index.find(a);
  if (i != index.end()) {
    bool_vec used(false, steps.size());
    mark_steps(i->second, used);
    index_set s_used(used);
    index_vec s_map;
    bool ok = mapping::invert_map(s_used, s_map);
    assert(ok);
    for (index_type k = 0; k < s_used.size(); k++) {
      step s(steps[s_used[k]].ass, steps[s_used[k]].rule);
      for (index_type i = 0; i < steps[s_used[k]].pre.size(); i++) {
	assert(s_map[steps[s_used[k]].pre[i]] != no_such_index);
	assert(s_map[steps[s_used[k]].pre[i]] < k);
	s.pre.append(s_map[steps[s_used[k]].pre[i]]);
      }
      s.valid = steps[s_used[k]].valid;
      p.push_back(s);
    }
  }
}

bool InferenceTracker::assertion::match
(plan_constraint_type t, index_type f1, index_type f2) const
{
  if (pc != t) return false;
  switch (pc) {
  case pc_at_end:
  case pc_always:
  case pc_sometime:
  case pc_at_most_once:
  case pc_never:
  case pc_disallowed:
  case pc_contradiction:
    return (f1 == a);
  case pc_sometime_before:
  case pc_sometime_after:
    return ((f1 == a) && (f2 == b));
  }
  return false;
}

bool InferenceTracker::extract_proof_premises
(const assertion& a, index_set& p) const
{
  p.clear();
  std::map<assertion, index_type>::const_iterator i = index.find(a);
  if (i == index.end()) return false; // no proof of assertion
  bool_vec used(false, steps.size());
  mark_steps(i->second, used);
  bool ok = true;
  for (index_type k = 0; k < steps.size(); k++)
    if (used[k]) {
      // std::cerr << "used step: ";
      // print_step(k, std::cerr);
      if (strcmp(steps[k].rule, "given") == 0) {
	// std::cerr << " - is a premise" << std::endl;
	bool found = false;
	for (index_type j = 0; (j < propagator.pct.size()) && !found; j++)
	  if (steps[k].ass.match((propagator.pct[j] == pc_at_end ?
				  pc_sometime : propagator.pct[j]),
				 propagator.pcf[j].first,
				 propagator.pcf[j].second)) {
	    p.insert(j);
	    found = true;
	  }
	if (!found)
	  ok = false;
      }
    }
  return ok;
}

void InferenceTracker::validate() const
{
  std::cerr << "validate called" << std::endl;
  //std::cout << "sb = " << propagator.sometime_before << std::endl;
  // for (index_type i = 0; i < propagator.sometime_before.size(); i++) {
  //   assertion a(pc_sometime_before,
  // 		propagator.sometime_before[i].first,
  // 		propagator.sometime_before[i].second);
  //   std::cout << "checking ";
  //   print_assertion(a, std::cout);
  //   std::cout << "..." << std::endl;
  //   print_proof(a, std::cout);
  // }
  exit(0);
}

// premise: A alpha
void InferenceTracker::premise_A
(index_type alpha, const char* rule)
{
  step s(assertion(pc_always, alpha), rule);
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
}

// premise: E alpha
void InferenceTracker::premise_E
(index_type alpha, const char* rule)
{
  step s(assertion(pc_sometime, alpha), rule);
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
}

// premise: sometime-before(alpha, beta) (or: beta SB alpha)
void InferenceTracker::premise_SB
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_sometime_before, alpha, beta), rule);
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
}

// premise: sometime-after(alpha, beta) (or: beta SA alpha)
void InferenceTracker::premise_SA
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_sometime_after, alpha, beta), rule);
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
}

// premise: AMO alpha
void InferenceTracker::premise_AMO
(index_type alpha, const char* rule)
{
  step s(assertion(pc_at_most_once, alpha), rule);
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
}

void InferenceTracker::infer_E_from_A
(index_type alpha, const char* rule)
{
  step s(assertion(pc_sometime, alpha), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_always, alpha));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// A alpha & mutex(alpha, beta) => N beta
void InferenceTracker::infer_N_from_A_mutex
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_never, alpha), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_always, beta));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  if (!propagator.mutex(alpha, beta))
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// A alpha & (del(act) -> !alpha) => D act
void InferenceTracker::infer_D_from_A
(index_type act, index_type alpha, const char* rule)
{
  step s(assertion(pc_disallowed, act), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_always, alpha));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  if (!propagator.act_implies_not[alpha].contains(act))
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// N alpha & (pre(act)|add(act) -> alpha) => D act
void InferenceTracker::infer_D_from_N
(index_type act, index_type alpha, const char* rule)
{
  step s(assertion(pc_disallowed, act), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_never, alpha));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  if (!propagator.act_implies[alpha].contains(act))
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

void InferenceTracker::infer_D_from_AMO
(index_type act, index_type alpha, const char* rule)
{
  step s(assertion(pc_disallowed, act), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_at_most_once, alpha));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  if (!propagator.act_change[2*alpha].contains(act))
    s.valid = false;
  if (!propagator.sf_initial[alpha])
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// // sometime-before(alpha, gamma) | alpha -> gamma;
// // sometime-before(gamma, beta) | gamma -> beta
// // => sometime-before(alpha, beta)
// void InferenceTracker::infer_transitive_SB
// (index_type alpha, index_type beta, const char* rule)
// {
//   step s(assertion(pc_sometime_before, alpha, beta), rule);
//   index_vec p;
//   index_type l = propagator.prec.shortest_path(alpha, beta, p);
//   if (l != no_such_index) {
//     assert(p.size() >= 2);
//     for (index_type k = 0; (k + 1) < (p.size()); k++) {
//       std::map<assertion, index_type>::iterator i =
// 	index.find(assertion(pc_sometime_before, p[k], p[k + 1]));
//       if (i != index.end()) {
// 	s.pre.append(i->second);
//       }
//       else if (!propagator.sf_implies[p[k]].contains(p[k + 1])) {
// 	s.valid = false;
// 	if (debug_mode) {
// 	  std::cerr << "infer transitive SB: alpha = " << alpha
// 		    << ", beta = " << beta << ", path = " << p
// 		    << " - failed at step " << k
// 		    << std::endl;
// 	  std::cerr << "prec = " << propagator.prec << std::endl;
// 	}
//       }
//     }
//     if (s.pre.size() == 0)
//       s.valid = false;
//   }
//   else {
//     s.valid = false;
//     if (debug_mode) {
//       std::cerr << "infer transitive SB: alpha = " << alpha
// 		<< ", beta = " << beta << " - no path" << std::endl;
//     }
//   }
//   steps.push_back(s);
//   index[s.ass] = steps.size() - 1;
//   if (verbose) print_step(steps.size() - 1, std::cerr);
//   if (!s.valid && debug_mode) validate();
// }

// // (gamma SA alpha) | (alpha -> gamma) &
// // (beta SA gamma) | (gamma -> beta) => beta SA alpha
// void InferenceTracker::infer_transitive_SA
// (index_type alpha, index_type beta, const char* rule)
// {
//   step s(assertion(pc_sometime_after, alpha, beta), rule);
//   index_vec p;
//   index_type l = propagator.wprec.shortest_path(alpha, beta, p);
//   //std::cerr << "l = " << l << ", p = " << p << std::endl;
//   if (l != no_such_index) {
//     assert(p.size() >= 2);
//     for (index_type k = 0; (k + 1) < (p.size()); k++) {
//       std::map<assertion, index_type>::iterator i =
// 	index.find(assertion(pc_sometime_after, p[k], p[k + 1]));
//       if (i != index.end())
// 	s.pre.append(i->second);
//       else if (!propagator.sf_implies[p[k]].contains(p[k + 1]))
// 	s.valid = false;
//     }
//     if (s.pre.size() == 0)
//       s.valid = false;
//   }
//   else {
//     s.valid = false;
//   }
//   steps.push_back(s);
//   index[s.ass] = steps.size() - 1;
//   if (verbose) print_step(steps.size() - 1, std::cerr);
//   if (!s.valid && debug_mode) validate();
// }

// find a path from alpha to beta in union of g and sf_implies, such
// that steps corresponding to each edge from g on the path exist; if
// found create a step and return its index, else return no_such_index.
// pc/rule indicates the type of relation (SB or SA).
index_type InferenceTracker::find_transitive
(index_type alpha, index_type beta,
 const adjacency_list_graph& g,
 plan_constraint_type pc,
 const char* rule)
{
  assert(alpha < propagator.n_state_formulas);
  assert(beta < propagator.n_state_formulas);
  // first, if a step for (pc, alpha, beta) exists, return it:
  std::map<assertion, index_type>::iterator l =
    index.find(assertion(pc, alpha, beta));
  if (l != index.end()) return l->second;
  if (debug_mode) {
    std::cerr << "searching for " << alpha << " -> " << beta << std::endl;
  }
  // else, compute forward reachability from alpha:
  bool_vec reached(false, propagator.n_state_formulas);
  index_vec parent(no_such_index, propagator.n_state_formulas);
  index_vec step_num(no_such_index, propagator.n_state_formulas);
  index_vec q; // exploration queue
  q.append(alpha);
  index_type p = 0;
  if (alpha != beta) reached[alpha] = true;
  while ((p < q.size()) & !reached[beta]) {
    assert(q[p] < propagator.n_state_formulas);
    index_type n = q[p++];
    if (debug_mode) {
      std::cerr << "n = " << n << std::endl;
    }
    for (index_type i = 0; i < propagator.sf_implies[n].size(); i++)
      if (!reached[propagator.sf_implies[n][i]]) {
	if (debug_mode) {
	  std::cerr << " - reached " << propagator.sf_implies[n][i]
		    << " by implication" << std::endl;
	}
	reached[propagator.sf_implies[n][i]] = true;
	parent[propagator.sf_implies[n][i]] = n;
	step_num[propagator.sf_implies[n][i]] = no_such_index;
	q.append(propagator.sf_implies[n][i]);
      }
    const index_set& sn = g.successors(n);
    for (index_type i = 0; i < sn.size(); i++)
      if (!reached[sn[i]]) {
	l = index.find(assertion(pc, n, sn[i]));
	if (l != index.end()) {
	  if (debug_mode) {
	    std::cerr << " - reached " << sn[i] << " by step "
		      << l->second << std::endl;
	  }
	  reached[sn[i]] = true;
	  parent[sn[i]] = n;
	  step_num[sn[i]] = l->second;
	  q.append(sn[i]);
	}
      }
  }
  // if we didn't reach beta, fail:
  if (!reached[beta])
    return no_such_index;
  // else, make a new step and find supporting steps by tracing the path:
  step s(assertion(pc, alpha, beta), rule);
  if (alpha == beta) {
    assert(parent[beta] < propagator.n_state_formulas);
    if (step_num[beta] != no_such_index)
      s.pre.append(step_num[beta]);
    p = parent[beta];
  }
  else {
    p = beta;
  }
  while (p != alpha) {
    if (debug_mode) {
      std::cerr << "p = " << p << ", s.pre = " << s.pre << std::endl;
    }
    assert(reached[p]);
    if (step_num[p] != no_such_index)
      s.pre.append(step_num[p]);
    p = parent[p];
  }
  if (debug_mode) {
    std::cerr << "end: p = " << p << ", s.pre = " << s.pre << std::endl;
  }
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  return steps.size() - 1;
}

// N beta & beta SB alpha => N alpha
void InferenceTracker::infer_N_from_N_SB
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_never, alpha), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_never, beta));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  index_type j = find_transitive(alpha, beta, propagator.sometime_before, 
				 pc_sometime_before, "transitive SB");
  if (j != no_such_index)
    s.pre.append(j);
  else
    s.valid = false;
  // std::map<assertion, index_type>::iterator j =
  //   index.find(assertion(pc_sometime_before, alpha, beta));
  // if (j != index.end())
  //   s.pre.append(j->second);
  // else
  //   s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// E beta & alpha SB beta => E alpha
void InferenceTracker::infer_E_from_E_SB
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_sometime, alpha), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_sometime, beta));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  index_type j = find_transitive(beta, alpha, propagator.sometime_before, 
				 pc_sometime_before, "transitive SB");
  if (j != no_such_index)
    s.pre.append(j);
  else
    s.valid = false;
  // std::map<assertion, index_type>::iterator j =
  //   index.find(assertion(pc_sometime_before, beta, alpha));
  // if (j != index.end())
  //   s.pre.append(j->second);
  // else
  //   s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// N beta & beta SA alpha => N alpha
void InferenceTracker::infer_N_from_N_SA
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_never, alpha), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_never, beta));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  index_type j = find_transitive(alpha, beta, propagator.sometime_after, 
				 pc_sometime_after, "transitive SA");
  if (j != no_such_index)
    s.pre.append(j);
  else
    s.valid = false;
  // std::map<assertion, index_type>::iterator j =
  //   index.find(assertion(pc_sometime_after, beta, alpha));
  // if (j != index.end())
  //   s.pre.append(j->second);
  // else
  //   s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// E beta & alpha SA beta => E alpha
void InferenceTracker::infer_E_from_E_SA
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_sometime, alpha), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_sometime, beta));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  index_type j = find_transitive(alpha, beta, propagator.sometime_after, 
				 pc_sometime_after, "transitive SA");
  if (j != no_such_index)
    s.pre.append(j);
  else
    s.valid = false;
  // std::map<assertion, index_type>::iterator j =
  //   index.find(assertion(pc_sometime_after, beta, alpha));
  // if (j != index.end())
  //   s.pre.append(j->second);
  // else
  //   s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// N beta & (alpha -> beta) => N alpha
void InferenceTracker::infer_N_from_N_implies
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_never, alpha), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_never, beta));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  if (!propagator.sf_implies[alpha].contains(beta))
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// E beta & (beta -> alpha) => E alpha
void InferenceTracker::infer_E_from_E_implies
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_sometime, alpha), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_sometime, beta));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  if (!propagator.sf_implies[beta].contains(alpha))
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// alpha SB beta & beta SB alpha => N alpha
void InferenceTracker::infer_N_from_SB_cycle
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_never, alpha), rule);
  index_type i = find_transitive(alpha, beta, propagator.sometime_before, 
				 pc_sometime_before, "transitive SB");
  if (i != no_such_index)
    s.pre.append(i);
  else
    s.valid = false;
  // std::map<assertion, index_type>::iterator i =
  //   index.find(assertion(pc_sometime_before, alpha, beta));
  // if (i != index.end())
  //   s.pre.append(i->second);
  // else if (!propagator.sf_implies[beta].contains(alpha))
  //   s.valid = false;
  if (alpha != beta) {
    index_type j = find_transitive(beta, alpha, propagator.sometime_before, 
				   pc_sometime_before, "transitive SB");
    if (j != no_such_index)
      s.pre.append(j);
    else
      s.valid = false;
  // std::map<assertion, index_type>::iterator j =
  //   index.find(assertion(pc_sometime_before, beta, alpha));
  // if (j != index.end())
  //   s.pre.append(j->second);
  // else if (!propagator.sf_implies[alpha].contains(beta))
  //   s.valid = false;
  }
  if (s.pre.size() == 0)
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// sometime-before(alpha, beta) & never-after(beta, alpha) => N alpha
void InferenceTracker::infer_N_from_SB_NA
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_never, alpha), rule);
  index_type i = find_transitive(alpha, beta, propagator.sometime_before, 
				 pc_sometime_before, "transitive SB");
  if (i != no_such_index)
    s.pre.append(i);
  else
    s.valid = false;
  // std::map<assertion, index_type>::iterator i =
  //   index.find(assertion(pc_sometime_before, alpha, beta));
  // if (i != index.end())
  //   s.pre.append(i->second);
  // else
  //   s.valid = false;
  if (propagator.neverafter) {
    if (!propagator.neverafter->adjacent(beta, alpha))
      s.valid = false;
  }
  else
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// (alpha SA beta) & (beta SA alpha) & mutex(alpha, beta) => N alpha
void InferenceTracker::infer_N_from_SA_cycle
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_never, alpha), rule);
  index_type i = find_transitive(alpha, beta, propagator.sometime_after, 
				 pc_sometime_after, "transitive SA");
  if (i != no_such_index)
    s.pre.append(i);
  else
    s.valid = false;
  // std::map<assertion, index_type>::iterator i =
  //   index.find(assertion(pc_sometime_after, alpha, beta));
  // if (i != index.end())
  //   s.pre.append(i->second);
  // else if (!propagator.sf_implies[alpha].contains(beta))
  //   s.valid = false;
  index_type j = find_transitive(beta, alpha, propagator.sometime_after, 
				 pc_sometime_after, "transitive SA");
  if (j != no_such_index)
    s.pre.append(j);
  else
    s.valid = false;
  // std::map<assertion, index_type>::iterator j =
  //   index.find(assertion(pc_sometime_after, beta, alpha));
  // if (j != index.end())
  //   s.pre.append(j->second);
  // else if (!propagator.sf_implies[beta].contains(alpha))
  //   s.valid = false;
  if (s.pre.size() == 0)
    s.valid = false;
  if (!propagator.mutex(alpha, beta))
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// beta SA alpha & beta NA alpha => N alpha
void InferenceTracker::infer_N_from_SA_NA
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_never, alpha), rule);
  index_type i = find_transitive(alpha, beta, propagator.sometime_after, 
				 pc_sometime_after, "transitive SA");
  if (i != no_such_index)
    s.pre.append(i);
  else
    s.valid = false;
  // std::map<assertion, index_type>::iterator i =
  //   index.find(assertion(pc_sometime_after, alpha, beta));
  // if (i != index.end())
  //   s.pre.append(i->second);
  // else
  //   s.valid = false;
  if (propagator.neverafter) {
    if (!propagator.neverafter->adjacent(alpha, beta))
      s.valid = false;
  }
  else
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// N alpha & E alpha => CONTRADICTION
void InferenceTracker::contradiction_from_N_E
(index_type alpha, const char* rule)
{
  step s(assertion(pc_contradiction), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_never, alpha));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  std::map<assertion, index_type>::iterator j =
    index.find(assertion(pc_sometime, alpha));
  if (j != index.end())
    s.pre.append(j->second);
  else
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// E alpha & E beta & (alpha NA beta) & (beta NA alpha) => CONTRADICTION
void InferenceTracker::contradiction_from_NA_cycle
(index_type alpha, index_type beta, const char* rule)
{
  step s(assertion(pc_contradiction), rule);
  std::map<assertion, index_type>::iterator i =
    index.find(assertion(pc_sometime, alpha));
  if (i != index.end())
    s.pre.append(i->second);
  else
    s.valid = false;
  std::map<assertion, index_type>::iterator j =
    index.find(assertion(pc_sometime, beta));
  if (j != index.end())
    s.pre.append(j->second);
  else
    s.valid = false;
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// |R[I]| > |S| => CONTRADICTION
// R is the set of atoms (implied E's)
// I is an independent subset of R (indexes R)
// S is a set AMO-action sets (indices into act_change) that covers the
// still allowed adders of each atom in R;
//  AMO f and f not initial => (2*f) may be in S
//  AMO f => (2*f)+1 may be in S
void InferenceTracker::contradiction_from_AMO
(const index_set& R, const index_set& I,
 const index_set& S, const char* rule)
{
  step s(assertion(pc_contradiction), rule);
  // first, find support for R[I]:
  for (index_type i = 0; i < I.size(); i++) {
    // the atom is R[I[i]]
    bool found = false;
    for (index_type f = 0; (f < propagator.n_state_formulas) && !found; f++)
      if (propagator.sf_implies_atom[f].contains(R[I[i]])) {
	std::map<assertion, index_type>::iterator b =
	  index.find(assertion(pc_sometime, f));
	if (b != index.end()) {
	  s.pre.append(b->second);
	  found = true;
	}
      }
    if (!found)
      s.valid = false;
  }
  // then, for S:
  index_set amo_acts;
  for (index_type i = 0; i < S.size(); i++) {
    assert(S[i] < (2 * propagator.n_state_formulas));
    amo_acts.insert(propagator.act_change[S[i]]);
    index_type f = S[i] / 2;
    std::map<assertion, index_type>::iterator b =
      index.find(assertion(pc_at_most_once, f));
    if (b != index.end())
      s.pre.append(b->second);
    else
      s.valid = false;
    // if S[i] is a "change to true" set, must also check f is not initial
    if ((S[i] % 2) == 0)
      if (propagator.sf_initial[f])
	s.valid = false;
  }
  // finally, for each atom in R[I], every adder that is not in the
  // collected set of AMO-actions should be disallowed:
  for (index_type i = 0; i < I.size(); i++) {
    for (index_type k = 0; k < propagator.ins.atoms[R[I[i]]].add_by.size(); k++)
      if (!amo_acts.contains(propagator.ins.atoms[R[I[i]]].add_by[k])) {
	std::map<assertion, index_type>::iterator b =
	  index.find(assertion(pc_disallowed, propagator.ins.atoms[R[I[i]]].add_by[k]));
	if (b != index.end())
	  s.pre.append(b->second);
	else
	  s.valid = false;
      }
  }
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// { Da | a in acts } => N alpha
void InferenceTracker::infer_N_from_D
(index_type alpha, const index_set& acts, const char* rule)
{
  step s(assertion(pc_never, alpha), rule);
  for (index_type k = 0; k < acts.size(); k++) {
    std::map<assertion, index_type>::iterator i =
      index.find(assertion(pc_disallowed, acts[k]));
    if (i != index.end())
      s.pre.append(i->second);
    else
      s.valid = false;
  }
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}

// { Da | a in acts } => beta SB alpha
void InferenceTracker::infer_SB_from_D
(index_type alpha, index_type beta, const index_set& acts, const char* rule)
{
  step s(assertion(pc_sometime_before, alpha, beta), rule);
  for (index_type k = 0; k < acts.size(); k++) {
    std::map<assertion, index_type>::iterator i =
      index.find(assertion(pc_disallowed, acts[k]));
    if (i != index.end())
      s.pre.append(i->second);
    else
      s.valid = false;
  }
  steps.push_back(s);
  index[s.ass] = steps.size() - 1;
  if (verbose) print_step(steps.size() - 1, std::cerr);
  if (!s.valid && debug_mode) validate();
}


Propagator::Propagator(const Instance& i, Stopwatch& s,
		       const graph* lm, const graph* na)
  : ins(i), stats(s), landmarks(lm), neverafter(na), it(0)
{
  // done
}

Propagator::~Propagator()
{
  // done
}

index_type Propagator::add_constraint
(plan_constraint_type t, index_type f)
{
  assert((t == pc_always) || (t == pc_sometime) ||
	 (t == pc_at_end) || (t == pc_at_most_once));
  assert(f < n_state_formulas);
  pct.push_back(t);
  pcf.append(index_pair(f, no_such_index));
  assert(pct.size() == pcf.size());
  return pct.size();
}

index_type Propagator::add_constraint
(plan_constraint_type t, index_type ft, index_type fc)
{
  assert((t == pc_sometime_before) || (t == pc_sometime_after));
  assert((ft < n_state_formulas) && (fc < n_state_formulas));
  pct.push_back(t);
  pcf.append(index_pair(ft, fc));
  assert(pct.size() == pcf.size());
  return pct.size();
}

void Propagator::reset()
{
  // clear inference state information
  never.assign_value(false, n_state_formulas);
  sometime.assign_value(false, n_state_formulas);
  sometime_before.init(n_state_formulas);
  sometime_after.init(n_state_formulas);
  allowed.assign_value(true, ins.n_actions());
  if (prec.size() != n_state_formulas) {
    // first call, initialise prec with implications and landmarks:
    prec.init(n_state_formulas);
    for (index_type i = 0; i < n_state_formulas; i++)
      if (sf_implies[i].size() > 0)
	prec.add_edge(i, sf_implies[i]);
    if (landmarks) {
      for (index_type i = 0; i < landmarks->size(); i++) {
	const index_set& succ = landmarks->successors(i);
	for (index_type j = 0; j < succ.size(); j++) {
	  prec.add_edge_to_transitive_closure(succ[j], i);
	}
      }
    }
    for (index_type i = 0; i < n_state_formulas; i++)
      assert(!prec.adjacent(i, i));
    // save prec for later re-init
    prec0 = prec;
  }
  else {
    // on subsequent calls, copy prec0 into prec:
    assert(prec0.size() == prec.size());
    prec.clear_edges();
    for (index_type i = 0; i < n_state_formulas; i++)
      if (prec0.successors(i).size() > 0)
	prec.add_edge(i, prec0.successors(i));
  }
  if (wprec.size() != n_state_formulas)
    wprec.init(n_state_formulas);
  else
    wprec.clear_edges();
  for (index_type i = 0; i < n_state_formulas; i++)
    if (sf_implies[i].size() > 0)
     wprec.add_edge(i, sf_implies[i]);
  // call reset hook to allow subclasses to add reset steps.
  reset_hook();
}

void Propagator::reset_hook()
{
  // does nothing
}

bool Propagator::run1(const index_set& s)
{
  reset();

  allowed_actions_changed = false;
  sometime_before_changed = false;
  sometime_after_changed = false;

  // input landmarks are SB's
  if (landmarks) {
    for (index_type i = 0; i < landmarks->size(); i++) {
      const index_set& succ = landmarks->successors(i);
      for (index_type j = 0; j < succ.size(); j++) {
	if (it) it->premise_SB(succ[j], i, "landmark");
	sometime_before.add_edge(succ[j], i);
	sometime_before_changed = true;
	// if (assert_sometime_before(succ[j], i))
	//   return false;
      }
    }
  }

  // find input SB's and initialise prec/wprec
  for (index_type i = 0; i < s.size(); i++)
    if (pct[s[i]] == pc_sometime_before)
      if (!sometime_before.adjacent(pcf[s[i]].first, pcf[s[i]].second)) {
	if (it) it->premise_SB(pcf[s[i]].first, pcf[s[i]].second, "given");
	if (assert_sometime_before(pcf[s[i]].first, pcf[s[i]].second))
	  return false;
      }

  // find input SA's and add to wprec
  for (index_type i = 0; i < s.size(); i++)
    if (pct[s[i]] == pc_sometime_after) {
      if (it) it->premise_SA(pcf[s[i]].first, pcf[s[i]].second, "given");
      if (assert_sometime_after(pcf[s[i]].first, pcf[s[i]].second))
	return false;
    }

  // find input always, assert implied nevers and disalloweds
  for (index_type i = 0; i < s.size(); i++)
    if (pct[s[i]] == pc_always) {
      if (it) it->premise_A(pcf[s[i]].first, "given");
      if (it) it->infer_E_from_A(pcf[s[i]].first, "A -> E");
      if (assert_sometime(pcf[s[i]].first)) return true;
      for (index_type j = 0; j < n_state_formulas; j++)
	if (mutex(pcf[s[i]].first, j)) {
	  if (it) it->infer_N_from_A_mutex(j, pcf[s[i]].first, "A & mutex");
	  if (assert_never(j)) return true;
	}
      for (index_type k = 0; k < act_implies_not[pcf[s[i]].first].size(); k++) {
	if (it) it->infer_D_from_A(act_implies_not[pcf[s[i]].first][k],
				   pcf[s[i]].first, "A & delete");
	allowed[act_implies_not[pcf[s[i]].first][k]] = false;
	allowed_actions_changed = true;
      }
    }

  // find input AMO's, assert implied disalloweds
  for (index_type i = 0; i < s.size(); i++)
    if (pct[s[i]] == pc_at_most_once) {
      index_type f = pcf[s[i]].first;
      if (it) it->premise_AMO(f, "given");
      if (sf_initial[f])
	for (index_type k = 0; k < act_change[2*f].size(); k++) {
	  if (it) it->infer_D_from_AMO(act_change[2*f][k], f, "AMO & init");
	  allowed[act_change[2*f][k]] = false;
	  allowed_actions_changed = true;
	}
    }

  // check conditional constraints after asserting input always
  if (allowed_actions_changed)
    if (check_conditional_constraints(allowed, never, prec, it))
      return true;

  // assert input sometime's and goals; algorithm in the paper does
  // this after the main loop, but we'll do it before so that
  // sometime/never contradictions can be detected immediately during
  // the loop.
  for (index_type i = 0; i < s.size(); i++)
    if ((pct[s[i]] == pc_sometime) || (pct[s[i]] == pc_at_end)) {
      if (it) it->premise_E(pcf[s[i]].first, "given");
      if (assert_sometime(pcf[s[i]].first)) return true;
    }

  // main loop
  while (allowed_actions_changed ||
	 sometime_before_changed ||
	 sometime_after_changed) {
    allowed_actions_changed = false;

    if (sometime_before_changed) {
      sometime_before_changed = false;

      // // add missing transitive SB edges to prec
      // for (index_type i = 0; i < prec.size(); i++) {
      // 	const index_set& s1 = prec.successors(i);
      // 	for (index_type j = 0; j < s1.size(); j++)
      // 	  if (s1[j] != i) {
      // 	    const index_set& s2 = prec.successors(s1[j]);
      // 	    for (index_type k = 0; k < s2.size(); k++)
      // 	      if ((s2[k] != s1[j]) && (s2[k] != i) && !prec.adjacent(i, s2[k])) {
      // 		// std::cerr << "transitive SB: " << i << "->" << s1[j]
      // 		// 		<< "->" << s2[k] << std::endl;
      // 		if (it) it->infer_transitive_SB(i, s2[k], "transitive SB");
      // 		if (assert_sometime_before(i, s2[k])) return true;
      // 	      }
      // 	  }
      // }

      // check for cycles in prec (strict never)
      // for (index_type i = 0; i < n_state_formulas; i++)
      // 	if (!never[i]) {
      // 	  const index_set& s1 = prec.successors(i);
      // 	  for (index_type j = 0; j < s1.size(); j++)
      // 	    if (prec.adjacent(s1[j], i)) {
      // 	      if (it) it->infer_N_from_SB_cycle(i, s1[j], "SB cycle");
      // 	      if (assert_never(i)) return true;
      // 	    }
      // 	}
      for (index_type i = 0; i < n_state_formulas; i++)
	if (!never[i])
	  if (prec.adjacent(i, i)) {
	    if (it) it->infer_N_from_SB_cycle(i, i, "SB cycle");
	    if (assert_never(i)) return true;
	  }

      // check for SB conflicts with NA's
      // note: NA graph may only contain a subset of state formulas, so
      // we need to check that both indices are within neverafter.size()
      if (neverafter)
	for (index_type i = 0; i < neverafter->size(); i++) {
	  const index_set& s1 = sometime_before.successors(i);
	  for (index_type j = 0; j < s1.size(); j++)
	    if ((s1[j] < neverafter->size()) && !never[s1[j]])
	      if (neverafter->adjacent(s1[j], i)) {
		// std::cerr << "prec: " << i << " -> " << s1[j]
		// 	      << " & NA: " << i << " -> " << s1[j]
		// 	      << std::endl;
		if (it) it->infer_N_from_SB_NA(i, s1[j], "SB & NA");
		if (assert_never(i)) return true;
	      }
	}
    }

    if (sometime_after_changed) {
      sometime_after_changed = false; 

      // // add missing transitive SA edges to wprec
      // for (index_type i = 0; i < wprec.size(); i++) {
      // 	const index_set& s1 = wprec.successors(i);
      // 	for (index_type j = 0; j < s1.size(); j++)
      // 	  if (s1[j] != i) {
      // 	    const index_set& s2 = wprec.successors(s1[j]);
      // 	    for (index_type k = 0; k < s2.size(); k++)
      // 	      if ((s2[k] != s1[j]) && (s2[k] != i) && !wprec.adjacent(i, s2[k])) {
      // 		if (it) it->infer_transitive_SA(i, s2[k], "transitive SA");
      // 		if (assert_sometime_after(i, s2[k])) return true;
      // 	      }
      // 	  }
      // }

      // check for cycles in wprec
      for (index_type i = 0; i < n_state_formulas; i++)
	if (!never[i]) {
	  const index_set& s1 = wprec.successors(i);
	  for (index_type j = 0; j < s1.size(); j++)
	    if (s1[j] != i)
	      if (mutex(i, s1[j]))
		if (wprec.adjacent(s1[j], i)) {
		  if (it) it->infer_N_from_SA_cycle(i, s1[j], "SA cycle");
		  if (assert_never(i)) return true;
		}
	}

      // check for SA conflicts with NA's
      if (neverafter) {
	for (index_type i = 0; i < neverafter->size(); i++)
	  if (!never[i]) {
	    const index_set& s1 = wprec.successors(i);
	    for (index_type j = 0; j < s1.size(); j++)
	      if (s1[j] < neverafter->size())
		if (neverafter->adjacent(i, s1[j])) {
		  // std::cerr << "wprec: " << i << " -> " << s1[j]
		  // 	      << " & NA: " << i << " -> " << s1[j]
		  // 	      << std::endl;
		  if (it) it->infer_N_from_SA_NA(i, s1[j], "SA & NA");
		  if (assert_never(i)) return true;
		}
	  }
      }
    }

    // check triggers of conditional nevers and landmarks
    if (allowed_actions_changed) {
      if (check_conditional_constraints(allowed, never, prec, it))
	return true;
    }
  }

  // check NA cycles
  if (neverafter) {
    for (index_type i = 0; i < neverafter->size(); i++)
      if (sometime[i])
	for (index_type j = 0; j < neverafter->size(); j++)
	  if (sometime[j] &&
	      neverafter->adjacent(i, j) && neverafter->adjacent(j, i)) {
	    if (it) it->contradiction_from_NA_cycle(i, j, "NA cycle");
	    return true;
	  }
  }

  // finally, check AMO count.
  if (checkAMO(s)) return true;
  return false;
}

void Propagator::set_tracker(InferenceTracker* t)
{
  if (it) delete it;
  it = t;
}

bool Propagator::run(const index_set& s)
{
  stats.start();
  bool cp = run1(s);
  stats.stop();
  return cp;
}

bool Propagator::assert_never(index_type f)
{
  if (never[f]) return false; // assertion already recorded
  never[f] = true;
  if (sometime[f]) {
    if (it) it->contradiction_from_N_E(f, "E & N (assertN)");
    return true;
  }
  for (index_type k = 0; k < act_implies[f].size(); k++) {
    if (it) it->infer_D_from_N(act_implies[f][k], f, "assertN");
    allowed[act_implies[f][k]] = false;
    allowed_actions_changed = true;
  }
  for (index_type i = 0; i < n_state_formulas; i++)
    if (!never[i])
      if (sometime_before.adjacent(i, f)) {
	if (it) it->infer_N_from_N_SB(i, f, "N & SB (assertN)");
	if (assert_never(i)) return true;
      }
  for (index_type i = 0; i < sf_implied_by[f].size(); i++)
    if (!never[sf_implied_by[f][i]]) {
      if (it) it->infer_N_from_N_implies(sf_implied_by[f][i], f,
					 "implication (assertN)");
      if (assert_never(sf_implied_by[f][i])) return true;
    }
  return false;
}

bool Propagator::assert_sometime(index_type f)
{
  if (sometime[f]) return false; // assertion already recorded
  sometime[f] = true;
  if (never[f]) {
    if (it) it->contradiction_from_N_E(f, "E & N (assertE)");
    return true;
  }
  const index_set& succs = sometime_before.successors(f);
  for (index_type i = 0; i < succs.size(); i++)
    if (!sometime[succs[i]]) {
      if (it) it->infer_E_from_E_SB(succs[i], f, "E & SB (assertE)");
      if (assert_sometime(succs[i])) return true;
    }
  for (index_type i = 0; i < sf_implies[f].size(); i++)
    if (!sometime[sf_implies[f][i]]) {
      if (it) it->infer_E_from_E_implies(sf_implies[f][i], f,
					 "implication (assertE)");
      if (assert_sometime(sf_implies[f][i])) return true;
    }
  return false;
}

bool Propagator::assert_sometime_before(index_type ft, index_type fc)
{
  // assertion already recorded:
  //if (sometime_before.adjacent(ft, fc)) return false;
  // add sb assertion and edge to prec:
  sometime_before.add_edge(ft, fc);
  sometime_before_changed = true;
  //prec.add_edge(ft, fc);
  prec.add_edge_to_transitive_closure(ft, fc);
  // propagate never and sometime:
  if (never[fc] && !never[ft]) {
    if (it) it->infer_N_from_N_SB(ft, fc, "N & SB (assertSB)");
    if (assert_never(ft)) return true;
  }
  if (sometime[ft] && !sometime[fc]) {
    if (it) it->infer_E_from_E_SB(fc, ft, "E & SB (assertSB)");
    if (assert_sometime(fc)) return true;
  }
  return false;
}

bool Propagator::assert_sometime_after(index_type ft, index_type fc)
{
  // assertion already recorded:
  //if (sometime_after.adjacent(ft, fc)) return false;
  // add sa assertion and edge to wprec:
  sometime_after.add_edge(ft, fc);
  sometime_after_changed = true;
  //wprec.add_edge(ft, fc);
  wprec.add_edge_to_transitive_closure(ft, fc);
  // propagate never and sometime:
  if (never[fc] && !never[ft]) {
    if (it) it->infer_N_from_N_SA(ft, fc, "N & SA (assertSA)");
    if (assert_never(ft)) return true;
  }
  if (sometime[ft] && !sometime[fc]) {
    if (it) it->infer_E_from_E_SA(fc, ft, "E & SA (assertSA)");
    if (assert_sometime(fc)) return true;
  }
  return false;
}

bool Propagator::check_conditional_constraints
(const bool_vec& allowed_actions,
 const bool_vec& current_never,
 const adjacency_list_graph& current_prec,
 InferenceTracker* it)
{
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (!ins.atoms[i].init &&
	!ins.atoms[i].add_by.have_common_element(allowed_actions)) {
      for (index_type f = 0; f < n_state_formulas; f++)
	if (!current_never[f] && sf_implies_atom[f].contains(i)) {
	  if (it) it->infer_N_from_D(f, ins.atoms[i].add_by,
				     "conditional never (default)");
	  if (assert_never(f)) return true;
	}
    }
  return false;
}

bool Propagator::checkAMO(const index_set& s)
{
  // find input set of at-most-once constraints:
  index_set c_amo;
  for (index_type i = 0; i < s.size(); i++)
    if (pct[s[i]] == pc_at_most_once) {
      /// don't need to do this here, it's done at the start of run1
      //if (it) it->premise_AMO(pcf[s[i]].first, "given");
      c_amo.insert(pcf[s[i]].first);
    }
  // if there aren't any, we won't find any contradiction.
  if (c_amo.empty()) return false;
  // find set of actions whose number of executions are limited
  // (to one) by each at-most-once constraint
  index_set amo_act_sets;
  for (index_type i = 0; i < c_amo.size(); i++) {
    amo_act_sets.insert((2*c_amo[i])+1); // acts change c_amo[i] to false
    if (!sf_initial[c_amo[i]])
      amo_act_sets.insert(2*c_amo[i]); // acts change c_amo[i] to true
    // for (index_type k = 0; k < ins.n_actions(); k++)
    //   if (action_change_to_false(k, c_amo[i]))
    // 	amo_acts[2*i].insert(k);
    // if (initial(c_amo[i])) {
    //   for (index_type k = 0; k < ins.n_actions(); k++)
    // 	if (action_change_to_true(k, c_amo[i])) {
    // 	  if (it) it->infer_D_from_AMO(k, c_amo[i], "AMO & init");
    // 	  allowed[k] = false;
    // 	}
    // }
    // else {
    //   for (index_type k = 0; k < ins.n_actions(); k++)
    // 	if (action_change_to_true(k, c_amo[i]))
    // 	  amo_acts[(2*i)+1].insert(k);
    // }
  }
  // find set of not initially true active sometime atoms whose (still
  // allowed) adders all belong to sets of actions thus limited:
  // amo_sets and req_sets are parallel vectors (pairs):
  // req_sets[i] contains a set of atoms such that all allowed adders
  // of each one of them is contained in the union of the amo_act_sets
  // indexed by amo_sets[i].
  index_set_vec amo_sets;
  index_set_vec req_sets;
  bool_vec active_sometime(false, ins.n_atoms());
  for (index_type j = 0; j < n_state_formulas; j++)
    if (sometime[j])
      active_sometime.insert(sf_implies_atom[j]);
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (!ins.atoms[i].init && active_sometime[i]) {
      index_set s;
      bool ok = true;
      for (index_type k = 0; (k < ins.atoms[i].add_by.size()) && ok; k++)
	if (allowed[ins.atoms[i].add_by[k]]) {
	  index_type w = no_such_index;
	  for (index_type j = 0;
	       (j < amo_act_sets.size()) && (w == no_such_index); j++)
	    if (act_change[amo_act_sets[j]].contains(ins.atoms[i].add_by[k]))
	      w = amo_act_sets[j];
	  if (w == no_such_index)
	    ok = false;
	  else
	    s.insert(w);
	}
      if (ok) {
	assert(req_sets.size() == amo_sets.size());
	bool append = true;
	for (index_type j = 0; j < amo_sets.size(); j++)
	  if (amo_sets[j].contains(s)) {
	    req_sets[j].insert(i);
	    if (amo_sets[j] == s)
	      append = false;
	  }
	if (append) {
	  amo_sets.append(s);
	  req_sets.append().assign_singleton(i);
	}
      }
    }
  // if the number of active sometime atoms that require one of the
  // actions in a union of amo_act sets to make it true exceeds the
  // maximum number of actions from this set that can be executed
  // (which is the number of amo_act sets in the union), we may have
  // a contradiction, but to be sure, we have to check for actions
  // that add more than one of the atoms
  for (index_type j = 0; j < amo_sets.size(); j++)
    if (req_sets[j].size() > amo_sets[j].size()) {
      // construct a graph over atoms in req_set, with an edge between
      // any pair of atoms that are both added by the same (allowed)
      // action
      graph req_ca(req_sets[j].size());
      for (index_type i1 = 0; i1 < req_sets[j].size(); i1++)
	for (index_type i2 = i1 + 1; i2 < req_sets[j].size(); i2++) {
	  bool has_edge = false;
	  for (index_type k1 = 0;
	       (k1 < ins.atoms[req_sets[j][i1]].add_by.size()) && !has_edge;
	       k1++)
	    if (allowed[ins.atoms[req_sets[j][i1]].add_by[k1]])
	      if (ins.atoms[req_sets[j][i2]].add_by.
		  contains(ins.atoms[req_sets[j][i1]].add_by[k1]))
		has_edge = true;
	  if (has_edge)
	    req_ca.add_undirected_edge(i1, i2);
	}
      // find an independent set in this graph: if the size of this
      // independent set is greater than amo_sets[j], then we have
      // a contradiction
      index_set ind_req_set;
      req_ca.apx_independent_set(ind_req_set);
      if (ind_req_set.size() > amo_sets[j].size()) {
	if (it) it->contradiction_from_AMO(req_sets[j], ind_req_set,
					   amo_sets[j], "AMO count");
	return true;
      }
    }
  return false;
}

/// PCC implementation

PCC::PCC(Instance& i, Stopwatch& s, const ppc_vec& ppcs)
  : Propagator(i, s, 0, 0), mx(0), trlm(0), xnag(NULL)
{
  init(ppcs);
}

PCC::PCC(Instance& i, Stopwatch& s, const ppc_vec& ppcs, StaticMutex* m,
	 const graph* lmg, set_edge_vec* t, const graph* nag)
  : Propagator(i, s, lmg, nag), mx(m), trlm(t), xnag(NULL)
{
  init(ppcs);
}

PCC::~PCC()
{
  if (xnag != NULL)
    delete xnag;
}

bool PCC::conjunction_implies(const index_set& f, index_type g) const
{
  if (g >= ins.n_atoms()) { // g is non-atomic
    g = (g - ins.n_atoms()); // g now indexes atomsets!
    assert(g < atomsets.size());
    if (disjunctive[g])
      // if g is disjunctive, f implies g iff they have any common atom
      return (f.have_common_element(atomsets[g]));
    else
      // if g is conjunctive, f implies g iff f contains g
      return (f.contains(atomsets[g]));
  }
  else { // g is atomic; f implies g iff f contains g
    return (f.contains(g));
  }
}

bool PCC::disjunction_implies(const index_set& f, index_type g) const
{
  // if f is a disjunction, it can only imply a greater disjunction
  if (g >= ins.n_atoms()) {
    g = (g - ins.n_atoms()); // g now indexes atomsets!
    assert(g < atomsets.size());
    if (disjunctive[g])
      if (atomsets[g].contains(f))
	return true;
  }
  return false;
}

bool PCC::check_implies(index_type f, index_type g) const
{
  // if f is non-atomic...
  if (f >= ins.n_atoms()) {
    f = (f - ins.n_atoms()); // f now indexes atomsets!
    assert(f < atomsets.size());
    if (disjunctive[f]) {
      return disjunction_implies(atomsets[f], g);
    }
    else { // f is conjunctive
      return conjunction_implies(atomsets[f], g);
    }
  }
  // f is atomic: f implies g only iff g is disjunctive and contains f
  else {
    if (g >= ins.n_atoms()) {
      g = (g - ins.n_atoms()); // g now indexes atomsets!
      assert(g < atomsets.size());
      if (disjunctive[g])
	if (atomsets[g].contains(f))
	  return true;
    }
    return false;
  }
}

bool PCC::conjunction_is_mutex(const index_set& f, index_type g) const
{
  for (index_type i = 0; i < f.size(); i++)
    if (mutex(g, f[i])) return true;
  return false;
}

bool PCC::disjunction_is_mutex(const index_set& f, index_type g) const
{
  for (index_type i = 0; i < f.size(); i++)
    if (!mutex(g, f[i])) return false;
  return true;
}

bool PCC::mutex(index_type f, index_type g) const
{
  if (f >= ins.n_atoms()) {
    f = (f - ins.n_atoms()); // f now indexes atomsets!
    assert(f < atomsets.size());
    if (disjunctive[f]) {
      return disjunction_is_mutex(atomsets[f], g);
    }
    else { // f is conjunctive
      return conjunction_is_mutex(atomsets[f], g);
    }
  }
  else { // f is atomic
    if (g >= ins.n_atoms()) {
      return mutex(g, f); // swap arguments
    }
    else { // both are atomic
      if (ins.atoms[f].neg == g) return true;
      if (mx) return mx->mutex(f, g);
      return false;
    }
  }
}

bool PCC::check_initial(index_type f) const
{
  if (f >= ins.n_atoms()) {
    f = (f - ins.n_atoms()); // f now indexes atomsets!
    assert(f < atomsets.size());
    if (disjunctive[f]) {
      for (index_type i = 0; i < atomsets[f].size(); i++)
	if (ins.atoms[atomsets[f][i]].init) return true;
      return false;
    }
    else { // f is conjunctive
      for (index_type i = 0; i < atomsets[f].size(); i++)
	if (!ins.atoms[atomsets[f][i]].init) return false;
      return true;
    }
  }
  else {
    return ins.atoms[f].init;
  }
}

// bool PCC::action_implies_formula(index_type a, index_type f) const
// {
//   assert(a < ins.n_actions());
//   // action implies formula f iff its (conjunctive) precondition or
//   // add effect implies f:
//   if (conjunction_implies(ins.actions[a].pre, f))
//     return true;
//   if (conjunction_implies(ins.actions[a].add, f))
//     return true;
//   return false;
// }

bool PCC::action_deletes_formula(index_type a, index_type f) const
{
  if (f >= ins.n_atoms()) {
    f = (f - ins.n_atoms()); // f now indexes atomsets!
    if (disjunctive[f]) {
      return ins.actions[a].del.contains(atomsets[f]);
    }
    else {
      return ins.actions[a].del.have_common_element(atomsets[f]);
    }
  }
  else {
    return ins.actions[a].del.contains(f);
  }
}

// bool PCC::action_implies_not_formula(index_type a, index_type f) const
// {
//   assert(a < ins.n_actions());
//   if (action_deletes_formula(a, f)) return true;
//   if (conjunction_is_mutex(ins.actions[a].pre, f)) return true;
//   return false;
// }

// bool PCC::formula_implies_atom(index_type f, index_type p) const
// {
//   if (f >= ins.n_atoms()) {
//     f = (f - ins.n_atoms()); // f now indexes atomsets!
//     assert(f < atomsets.size());
//     return (!disjunctive[f] && atomsets[f].contains(p));
//   }
//   else {
//     return (f == p);
//   }
// }

bool PCC::check_action_change_to_true(index_type a, index_type f) const
{
  return (conjunction_is_mutex(ins.actions[a].pre, f) &&
	  conjunction_implies(ins.actions[a].add, f));
}

bool PCC::check_action_change_to_false(index_type a, index_type f) const
{
  return (conjunction_implies(ins.actions[a].pre, f) &&
	  action_deletes_formula(a, f));
}

void PCC::print_state_formula(index_type f, std::ostream& s) const
{
  if (f >= ins.n_atoms()) {
    f = (f - ins.n_atoms());
    assert(f < atomsets.size());
    if (disjunctive[f])
      s << "\\/{";
    else
      s << "/\\{";
    for (index_type i = 0; i < atomsets[f].size(); i++) {
      if (i > 0) s << ", ";
      s << atomsets[f][i] << "." << ins.atoms[atomsets[f][i]].name;
    }
    s << "}";
  }
  else {
    s << f << "." << ins.atoms[f].name;
  }
}

index_type PCC::map_state_formula(const index_set& s, bool dis)
{
  // the set is atomic: return the atom index.
  if (s.size() == 1) {
    assert(s[0] < ins.n_atoms());
    return s[0];
  }
  assert(atomsets.size() == disjunctive.size());
  index_type i = atomsets.first(s);
  // the set does not exist: add it and return.
  if (i == no_such_index) {
    atomsets.append(s);
    disjunctive.append(dis);
    n_state_formulas = ins.n_atoms() + atomsets.size();
    return (ins.n_atoms() + (atomsets.size() - 1));
  }
  assert(i < atomsets.size());
  // the set exists, but with opposite con/dis flag: add second variety.
  if (dis != disjunctive[i]) {
    atomsets.append(s);
    disjunctive.append(dis);
    n_state_formulas = ins.n_atoms() + atomsets.size();
    return (ins.n_atoms() + (atomsets.size() - 1));
  }
  // the set exists, with correct con/dis flag: return it.
  return (ins.n_atoms() + i);
}

void PCC::init(const ppc_vec& ppcs)
{
  n_state_formulas = ins.n_atoms();
  // find compound state formulas and store pcs in internal form
  for (index_type k = 0; k < ppcs.size(); k++) {
    switch (ppcs[k].pct) {
    case pc_at_end:
    case pc_always:
    case pc_sometime:
    case pc_at_most_once:
      {
	index_type f = map_state_formula(ppcs[k].s_c, ppcs[k].c_dis);
	add_constraint(ppcs[k].pct, f);
      }
      break;
    case pc_sometime_before:
    case pc_sometime_after:
      {
	index_type ft = map_state_formula(ppcs[k].s_t, ppcs[k].t_dis);
	index_type fc = map_state_formula(ppcs[k].s_c, ppcs[k].c_dis);
	add_constraint(ppcs[k].pct, ft, fc);
      }
      break;
    default:
      std::cerr << "PCC: ";
      ppcs[k].write_PDDL_constraint(std::cerr, ins);
      std::cerr << " not supported" << std::endl;
      exit(255);
    }
  }

  std::cerr << n_state_formulas << " state formulas ("
	    << ins.n_atoms() << " atomic)" << std::endl;
  if (landmarks) {
    std::cerr << landmarks->n_edges() << " landmark edges" << std::endl;
  }

  // make at-end pcs from problem goals
  for (index_type i = 0; i < ins.goal_atoms.size(); i++) {
    index_type j = add_constraint(pc_at_end, ins.goal_atoms[i]);
    goals.insert(j - 1);
  }

  // initialise the implication graph
  graph gimp(n_state_formulas);
  // single atoms imply any disjunction they're part of, and
  // are implied by conjunctions they're part of:
  for (index_type i = 0; i < ins.n_atoms(); i++)
    for (index_type j = 0; j < atomsets.size(); j++)
      if (atomsets[j].contains(i)) {
	if (disjunctive[j])
	  gimp.add_edge(i, ins.n_atoms() + j);
	else
	  gimp.add_edge(ins.n_atoms() + j, i);
      }
  // a compond conjunction implies any conjunction it contains:
  for (index_type i = 0; i < atomsets.size(); i++)
    if (!disjunctive[i])
      for (index_type j = i + 1; j < atomsets.size(); j++)
	if (!disjunctive[j]) {
	  if (atomsets[i].contains(atomsets[j]))
	    gimp.add_edge(ins.n_atoms() + i, ins.n_atoms() + j);
	  else if (atomsets[j].contains(atomsets[i]))
	    gimp.add_edge(ins.n_atoms() + j, ins.n_atoms() + i);
	}
  // a compond disjunction implies any disjunction it is contained in:
  for (index_type i = 0; i < atomsets.size(); i++)
    if (disjunctive[i])
      for (index_type j = i + 1; j < atomsets.size(); j++)
	if (disjunctive[j]) {
	  if (atomsets[j].contains(atomsets[i]))
	    gimp.add_edge(ins.n_atoms() + i, ins.n_atoms() + j);
	  else if (atomsets[i].contains(atomsets[j]))
	    gimp.add_edge(ins.n_atoms() + j, ins.n_atoms() + i);
	}
  // all implications from conjunctions to disjunctions are found by
  // taking the transitive closure, since they must be via an atom.
  gimp.transitive_closure();
  sf_implies.assign_value(EMPTYSET, n_state_formulas);
  sf_implied_by.assign_value(EMPTYSET, n_state_formulas);
  for (index_type i = 0; i < n_state_formulas; i++) {
    sf_implies[i] = gimp.successors(i);
    sf_implied_by[i] = gimp.predecessors(i);
  }

  // initialise sf_implies_atom
  sf_implies_atom.assign_value(EMPTYSET, n_state_formulas);
  for (index_type i = 0; i < ins.n_atoms(); i++)
    sf_implies_atom[i].insert(i);
  for (index_type i = 0; i < atomsets.size(); i++)
    if (!disjunctive[i])
      sf_implies_atom[ins.n_atoms() + i] = atomsets[i];

  // initialise other sf arrays:
  sf_initial.assign_value(false, n_state_formulas);
  for (index_type i = 0; i < n_state_formulas; i++)
    sf_initial[i] = check_initial(i);
  act_implies.assign_value(EMPTYSET, n_state_formulas);
  act_implies_not.assign_value(EMPTYSET, n_state_formulas);
  act_change.assign_value(EMPTYSET, 2*n_state_formulas);
  for (index_type i = 0; i < n_state_formulas; i++)
    for (index_type k = 0; k < ins.n_actions(); k++) {
      if (conjunction_implies(ins.actions[k].pre, i) ||
	  conjunction_implies(ins.actions[k].add, i))
	act_implies[i].insert(k);
      if (action_deletes_formula(k, i) ||
	  conjunction_is_mutex(ins.actions[k].pre, i))
	act_implies_not[i].insert(k);
      if (check_action_change_to_true(k, i))
	act_change[2*i].insert(k);
      else if (check_action_change_to_false(k, i))
	act_change[(2*i)+1].insert(k);
    }

  // create expanded never-after graph
  if ((neverafter != NULL) && (n_state_formulas > ins.n_atoms())) {
    xnag = new graph(n_state_formulas);
    for (index_type i = 0; i < ins.n_atoms(); i++)
      for (index_type j = 0; j < ins.n_atoms(); j++)
	if (neverafter->adjacent(i, j)) {
	  xnag->add_edge(i, j);
	  xnag->add_edge(sf_implied_by[i], j);
	  xnag->add_edge(i, sf_implied_by[j]);
	}
    neverafter = xnag;
  }
}

bool PCC::check_conditional_constraints
(const bool_vec& allowed_actions,
 const bool_vec& current_never,
 const adjacency_list_graph& current_prec,
 InferenceTracker* it)
{
  for (index_type i = 0; i < ins.n_atoms(); i++)
    if (!ins.atoms[i].init && !current_never[i] &&
	!ins.atoms[i].add_by.have_common_element(allowed_actions)) {
      if (it) it->infer_N_from_D(i, ins.atoms[i].add_by,
				 "conditional never");
      if (assert_never(i)) return true;
    }
  if (trlm) {
    for (index_type k = 0; k < trlm->size(); k++)
      if (rem_trlm[k])
	if (!(*trlm)[k].first.have_common_element(allowed_actions)) {
	  if (!current_prec.adjacent((*trlm)[k].second.second,
				     (*trlm)[k].second.first)) {
	    if (it) it->infer_SB_from_D((*trlm)[k].second.second,
					(*trlm)[k].second.first,
					(*trlm)[k].first,
					"conditional landmark");
	    if (assert_sometime_before((*trlm)[k].second.second,
				       (*trlm)[k].second.first))
	      return true;
	  }
	  rem_trlm[k] = false;
	}
  }
  return false;
}

void PCC::reset_hook()
{
  if (trlm) rem_trlm.assign_value(true, trlm->size());
}

bool PCC::run(const index_set& s)
{
  index_set s1(s);
  s1.insert(goals);
  return Propagator::run(s1);
}

/// ConsistencyTest implementations

ConsistencyTest::ConsistencyTest(const ppc_vec& constraints)
  : ppcs(constraints)
{
  // done
}

ConsistencyTest::~ConsistencyTest()
{
  // done
}

void ConsistencyTest::minimise_conflict()
{
  index_type i = 0;
  while (i < conflict.size()) {
    index_type e = conflict[i];
    conflict.remove(i);
    bool con = test(conflict, false, false);
    if (!con) {
      conflict.insert(e);
      assert(conflict[i] == e);
      i += 1;
    }
  }
}

PCCTest::PCCTest
(const ppc_vec& constraints,
 Instance& ins,
 Preprocessor& prep,
 Stopwatch& prep_stats,
 Stopwatch& test_stats,
 count_type& n_tests)
  : ConsistencyTest(constraints),
    pcc_test_stats(test_stats),
    n_pcc_tests(n_tests),
    verbose1(false),
    verbose2(false)
{
  prep_stats.start();
  StaticMutex* mx = prep.inconsistency();
  graph* lmg = prep.necessary_sb_graph();
  graph* nag = prep.never_after_graph();
  landmark_graph_triggered_edges(ins, *lmg, trlm);
  tester = new PCC(ins, test_stats, ppcs, mx, lmg, &trlm, nag);
  prep_stats.stop();
}

PCCTest::~PCCTest()
{
  tester->set_tracker(NULL);
  delete tester;
}

bool PCCTest::test
(const index_set& s,
 bool find_conflict = false,
 bool min_conflict = false)
{
  if (find_conflict) {
    InferenceTracker* tracker = new InferenceTracker(*tester);
    tracker->set_verbose(verbose2);
    tester->set_tracker(tracker);
    bool con = tester->run(s);
    if (con) {
      if (verbose1) tracker->print_proof(std::cerr);
      conflict.clear();
      InferenceTracker::assertion a(pc_contradiction);
      bool ok = tracker->extract_proof_premises(InferenceTracker::assertion(pc_contradiction), conflict);
      assert(ok);
      assert(conflict.size() > 0);
      conflict.remove_greater_than(ppcs.size() - 1);
      if (min_conflict) {
	tester->set_tracker(NULL);
	minimise_conflict();
      }
    }
    tester->set_tracker(NULL);
    n_pcc_tests += 1;
    return con;
  }
  else if (verbose1 || verbose2) {
    InferenceTracker* tracker = new InferenceTracker(*tester);
    tracker->set_verbose(verbose2);
    tester->set_tracker(tracker);
    bool con = tester->run(s);
    if (con && verbose1) tracker->print_proof(std::cerr);
    tester->set_tracker(NULL);
    n_pcc_tests += 1;
    return con;
  }
  else {
    bool con = tester->run(s);
    n_pcc_tests += 1;
    return con;
  }
}

hmTest::hmTest
(const ppc_vec& constraints,
 Instance& ins,
 Stopwatch& test_stats,
 count_type& n_built,
 count_type& n_tests)
  : ConsistencyTest(constraints),
    instance(ins),
    hm_test_stats(test_stats),
    n_hm_built(n_built),
    n_hm_tests(n_tests),
    verbose(false),
    opt_H2(false)
{
  // done
}

hmTest::~hmTest()
{
  // done
}

bool hmTest::test
(const index_set& s,
 bool find_conflict,
 bool min_conflict)
{
  bool con = false;
  hm_test_stats.start();
  Instance* test_ins = instance.copy();
  mapping map(test_ins->n_actions());
  for (index_type i = 0; i < s.size(); i++) {
    assert(s[i] < ppcs.size());
    ppcs[s[i]].enforce(*test_ins, map);
  }
  test_ins->cross_reference();
  Heuristic* test_h = 0;
  if (opt_H2) {
    StaticMutex* mx = new StaticMutex(*test_ins);
    test_h = mx;
  }
  else {
    test_h = new Reachability(*test_ins);
  }
  n_hm_built += 1;
  n_hm_tests += 1;
  con = INFINITE(test_h->eval(test_ins->goal_atoms));
  delete test_h;
  delete test_ins;
  if (con && find_conflict) {
    bool ok = generate_conflict(s);
    assert(ok);
    if (min_conflict) {
      //conflict.assign_copy(s);
      minimise_conflict();
    }
    // else {
    //   bool ok = generate_conflict(s);
    //   assert(ok);
    // }
  }
  hm_test_stats.stop();
  return con;
}

bool hmTest::generate_conflict(const index_set& s)
{
  conflict.clear();
  Instance* test_ins = instance.copy();
  mapping map(test_ins->n_actions());
  index_set ca;
  index_set ce;
  index_vec g_ce;
  for (index_type i = 0; i < s.size(); i++) {
    assert(s[i] < ppcs.size());
    if (ppcs[s[i]].pct == pc_sometime) {
      index_type g =
	(ppcs[s[i]].c_dis ?
	 test_ins->compile_pc_sometime_disjunction(ppcs[s[i]].s_c,
						   ppcs[s[i]].name) :
	 test_ins->compile_pc_sometime_conjunction(ppcs[s[i]].s_c,
						   ppcs[s[i]].name, &map));
      ce.append(s[i]);
      g_ce.append(g);
      assert(ce.size() == g_ce.size());
    }
    else {
      ca.append(s[i]);
    }
  }
  for (index_type i = 0; i < ca.size(); i++) {
    conflict.insert(ca[i]);
    ppcs[ca[i]].enforce(*test_ins, map);
    test_ins->cross_reference();
    Heuristic* test_h = (opt_H2 ?
			 (Heuristic*)new StaticMutex(*test_ins) :
			 (Heuristic*)new Reachability(*test_ins));
    n_hm_built += 1;
    n_hm_tests += 1;
    if (INFINITE(test_h->eval(test_ins->goal_atoms))) {
      // current conflict set is contradictory
      // std::cerr << "contradiction 1: " << conflict
      // 		<< ", goals = " << test_ins->goal_atoms
      // 		<< ", original goals = " << instance.goal_atoms
      // 		<< ", g_ce = " << g_ce
      // 		<< std::endl;
      delete test_h;
      delete test_ins;
      return true;
    }
    index_set tset(test_ins->goal_atoms);
    for (index_type j = 0; j < ce.size(); j++) {
      tset.insert(g_ce[j]);
      n_hm_tests += 1;
      if (INFINITE(test_h->eval(tset))) {
	// std::cerr << "contradiction 2: " << conflict
	// 	  << ", ce = " << ce << ", j = " << j
	// 	  << std::endl;
	// current conflict + ce[0..j] is contradictory
	for (index_type l = 0; l <= j; l++)
	  conflict.insert(ce[l]);
	delete test_h;
	delete test_ins;
	return true;
      }
    }
    delete test_h;
  }
  // no contradictory set found
  return false;
}

SerialTest::SerialTest
(const ppc_vec& constraints,
 ConsistencyTest& test1,
 ConsistencyTest& test2)
  : ConsistencyTest(constraints),
    _test1(test1),
    _test2(test2)
{
  // done
}

SerialTest::~SerialTest()
{
  // done
}

bool SerialTest::test
(const index_set& s, bool find_conflict, bool min_conflict)
{
  bool con = _test1.test(s, find_conflict, min_conflict);
  if (con) {
    if (find_conflict)
      conflict = _test1.last_conflict();
    return true;
  }
  con = _test2.test(s, find_conflict, min_conflict);
  if (con) {
    if (find_conflict)
      conflict = _test2.last_conflict();
  }
  return con;
}

END_HSPS_NAMESPACE
