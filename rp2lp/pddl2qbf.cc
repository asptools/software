
#include "parser.h"
#include "preprocess.h"
#include "exec.h"
#include "hypergraph.h"
#include "cost_table.h"
#include "enumerators.h"

#include <sstream>
#include <fstream>

typedef long int literal;
typedef std::vector<literal> literal_vec;
typedef literal_vec clause;
typedef std::vector<clause> clause_vec;

class QBF {
  long int v_next;
  HSPS::index_vec block;
  bool first_block_is_universal;
  bool last_block_is_universal;
  clause_vec clauses;

public:
  QBF();

  literal new_var(bool u, unsigned int k = 1);
  void new_var_vec(bool u, unsigned int k, literal_vec& v);

  unsigned int n_variables() const;
  unsigned int n_clauses() const;

  literal negate(literal x) const;

  void post(literal l1);
  void post(literal l1, literal l2);
  void post(literal l1, literal l2, literal l3);
  void post(const clause& c);

  void post_reified_or(literal l1, literal l2, literal r);
  void post_reified_or(const literal_vec& l, literal r);
  void post_reified_and(literal l1, literal l2, literal r);
  void post_reified_and(const literal_vec& l, literal r);

  void write_qdimacs(std::ostream& s);
};

QBF::QBF()
  : v_next(1), first_block_is_universal(false), last_block_is_universal(false)
{
  // done
}

unsigned int QBF::n_variables() const
{
  return v_next;
}

unsigned int QBF::n_clauses() const
{
  return clauses.size();
}

literal QBF::new_var(bool u, unsigned int k)
{
  assert(k > 0);
  literal first = v_next;
  v_next += k;
  if (block.size() == 0) {
    first_block_is_universal = u;
    block.append(first);
    last_block_is_universal = u;
  }
  else if (last_block_is_universal != u) {
    block.append(first);
    last_block_is_universal = u;
  }
  return first;
}

void QBF::new_var_vec(bool u, unsigned int k, literal_vec& v)
{
  literal first = new_var(u, k);
  v.resize(k);
  for (unsigned int i = 0; i < k; i++)
    v[i] = first + i;
}

literal QBF::negate(literal x) const
{
  assert(x != 0);
  return -x;
}

void QBF::post(literal l1)
{
  clause c;
  c.push_back(l1);
  post(c);
}

void QBF::post(literal l1, literal l2)
{
  clause c;
  c.push_back(l1);
  c.push_back(l2);
  post(c);
}

void QBF::post(literal l1, literal l2, literal l3)
{
  clause c;
  c.push_back(l1);
  c.push_back(l2);
  c.push_back(l3);
  post(c);
}

void QBF::post(const clause& c)
{
  clauses.push_back(c);
}

void QBF::post_reified_or(literal l1, literal l2, literal r)
{
  // l1 -> r
  post(negate(l1), r);
  // l2 -> r
  post(negate(l2), r);
  // r -> (l1 \/ l2)
  post(negate(r), l1, l2);
}

void QBF::post_reified_or(const literal_vec& l, literal r)
{
  // l[k] -> r, for all k
  for (unsigned int k = 0; k < l.size(); k++)
    post(negate(l[k]), r);
  // r -> \/_k l[k]
  literal_vec r_imp_l(l.size() + 1, 0);
  r_imp_l[0] = negate(r);
  for (unsigned int k = 0; k < l.size(); k++)
    r_imp_l[k + 1] = l[k];
  post(r_imp_l);
}

void QBF::post_reified_and(literal l1, literal l2, literal r)
{
  // r -> l1
  post(negate(r), l1);
  // r -> l2
  post(negate(r), l2);
  // (l1 /\ l2) -> r
  post(negate(l1), negate(l2), r);
}

void QBF::post_reified_and(const literal_vec& l, literal r)
{
  // r -> l[k], for all k
  for (unsigned int k = 0; k < l.size(); k++)
    post(negate(r), l[k]);
  // /\_k l[k] -> r
  literal_vec l_imp_r(l.size() + 1, 0);
  for (unsigned int k = 0; k < l.size(); k++)
    l_imp_r[k] = negate(l[k]);
  l_imp_r[l.size()] = r;
  post(l_imp_r);
}

void QBF::write_qdimacs(std::ostream& s)
{
  s << "p cnf " << v_next << " " << clauses.size() << std::endl;
  bool block_is_universal = first_block_is_universal;
  block.append(v_next);
  for (HSPS::index_type k = 0; k < block.size() - 1; k++) {
    s << (block_is_universal ? "a" : "e");
    for (HSPS::index_type i = block[k]; i < block[k + 1]; i++)
      s << " " << i;
    s << " 0" << std::endl;
    block_is_universal = !block_is_universal;
  }
  for (HSPS::index_type k = 0; k < clauses.size(); k++) {
    for (HSPS::index_type i = 0; i < clauses[k].size(); i++)
      s << clauses[k][i] << " ";
    s << "0" << std::endl;
  }
}


int main(int argc, char *argv[]) {
  bool     opt_preprocess = true;
  bool     opt_verbose = false;
  HSPS::index_type alpha_m = 5;
  unsigned int random_seed = 0;

  HSPS::PDDL_Base::warning_level = 0;
  HSPS::PDDL_Base::compile_away_plan_constraints = false;

  HSPS::Instance::write_PDDL2 = false;
  HSPS::Instance::write_PDDL3 = false;
  HSPS::Instance::write_metric = false;
  HSPS::Instance::write_time = false;

  HSPS::Statistics stats;
  stats.enable_interrupt();
  stats.start();

  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  HSPS::Parser* reader = new HSPS::Parser(symbols);

  for (int k = 1; k < argc; k++) {
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      HSPS::Instance::default_trace_level = atoi(argv[++k]);
      opt_verbose = true;
    }
    else if (strcmp(argv[k],"-no-prep") == 0) {
      opt_preprocess = false;
    }
    else if ((strcmp(argv[k],"-catc") == 0) && (k < argc - 1)) {
      HSPS::PDDL_Name::catc = *argv[++k];
    }
    else if (strcmp(argv[k],"-silly-cap") == 0) {
      HSPS::PDDL_Name::silly_cap = true;
    }
    else if ((strcmp(argv[k],"-m") == 0) && (k < argc - 1)) {
      alpha_m = atoi(argv[++k]);
      assert(alpha_m > 0);
    }
    else if (((strcmp(argv[k],"-rnd") == 0) ||
	      (strcmp(argv[k],"-r") == 0)) &&
	     (k < argc - 1)) {
      random_seed = atoi(argv[++k]);
    }
    else if (*argv[k] != '-') {
      reader->read(argv[k], false);
    }
  }

  HSPS::Instance instance;
  reader->post_process();
  stats.start();
  reader->instantiate(instance);
  HSPS::Preprocessor prep(instance, stats);
  if (opt_preprocess) {
    prep.preprocess();
  }
  else {
    instance.cross_reference();
  }
  stats.stop();
  std::cerr << "instantiation and preprocessing finished in "
	    << stats.time() << " seconds" << std::endl;
  std::cerr << instance.name << ": " << instance.n_atoms() << " atoms, "
	    << instance.n_actions() << " actions" << std::endl;

  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++) {
    std::cout << "c " << i << " = " << instance.atoms[i].name << std::endl;
  }

  QBF theory;

  // outermost existential
  std::cout << "c var " << theory.n_variables() << ": alpha ("
	    << instance.n_atoms() * 2 * alpha_m << ")" << std::endl;
  literal_vec alpha;
  theory.new_var_vec(false, instance.n_atoms() * 2 * alpha_m, alpha);

#define Pj_IN_ALPHAi(i,j) (alpha[(instance.n_atoms() * 2 * (i)) + (2 * (j)) + 1])
#define NOT_Pj_IN_ALPHAi(i,j) (alpha[(instance.n_atoms() * 2 * (i)) + (2 * (j))])

  for (HSPS::index_type i = 0; i < alpha_m; i++) {
    for (HSPS::index_type j = 0; j < instance.n_atoms(); j++) {
      std::cout << "c map " << NOT_Pj_IN_ALPHAi(i,j) << " = (not "
		<< instance.atoms[j].name << ") in beta" << i
		<< std::endl;
      std::cout << "c map " << Pj_IN_ALPHAi(i,j) << " = "
		<< instance.atoms[j].name << " in beta" << i
		<< std::endl;
    }
  }

  // next, universal
  std::cout << "c var " << theory.n_variables() << ": x0 ("
	    << instance.n_atoms() << ")" << std::endl;
  literal_vec x0;
  theory.new_var_vec(true, instance.n_atoms(), x0);
  std::cout << "c var " << theory.n_variables() << ": x1 ("
	    << instance.n_atoms() << ")" << std::endl;
  literal_vec x1;
  theory.new_var_vec(true, instance.n_atoms(), x1);
  std::cout << "c var " << theory.n_variables() << ": a0 ("
	    << instance.n_actions() << ")" << std::endl;
  literal_vec a0;
  theory.new_var_vec(true, instance.n_actions(), a0);

  // all following allocations belong to the innermost existential

  // construct alpha(x0)
  std::cout << "c cl " << theory.n_clauses() << ": alpha(x0)" << std::endl;
  std::cout << "c var " << theory.n_variables() << ": beta0 ("
	    << alpha_m << ")" << std::endl;
  literal_vec beta0;
  theory.new_var_vec(false, alpha_m, beta0);
  for (HSPS::index_type i = 0; i < alpha_m; i++) {
    std::cout << "c var " << theory.n_variables() << ": gamma ("
	      << instance.n_atoms() * 2 << ")" << std::endl;
    literal_vec gamma;
    theory.new_var_vec(false, instance.n_atoms() * 2, gamma);
    for (HSPS::index_type j = 0; j < instance.n_atoms(); j++) {
      theory.post_reified_or(theory.negate(Pj_IN_ALPHAi(i,j)), x0[j], gamma[2 * j]);
      theory.post_reified_or(theory.negate(NOT_Pj_IN_ALPHAi(i,j)), theory.negate(x0[j]), gamma[(2 * j) + 1]);
    }
    theory.post_reified_and(gamma, beta0[i]);
  }
  std::cout << "c var " << theory.n_variables() << ": alpha0" << std::endl;
  literal alpha0 = theory.new_var(false);
  theory.post_reified_or(beta0, alpha0);

  // construct alpha(x1)
  std::cout << "c cl " << theory.n_clauses() << ": alpha(x1)" << std::endl;
  std::cout << "c var " << theory.n_variables() << ": beta1 ("
	    << alpha_m << ")" << std::endl;
  literal_vec beta1;
  theory.new_var_vec(false, alpha_m, beta1);
  for (HSPS::index_type i = 0; i < alpha_m; i++) {
    std::cout << "c var " << theory.n_variables() << ": gamma ("
	      << instance.n_atoms() * 2 << ")" << std::endl;
    literal_vec gamma;
    theory.new_var_vec(false, instance.n_atoms() * 2, gamma);
    for (HSPS::index_type j = 0; j < instance.n_atoms(); j++) {
      theory.post_reified_or(theory.negate(Pj_IN_ALPHAi(i,j)), x1[j], gamma[2 * j]);
      theory.post_reified_or(theory.negate(NOT_Pj_IN_ALPHAi(i,j)), theory.negate(x1[j]), gamma[(2 * j) + 1]);
    }
    theory.post_reified_and(gamma, beta1[i]);
  }
  std::cout << "c var " << theory.n_variables() << ": alpha1" << std::endl;
  literal alpha1 = theory.new_var(false);
  theory.post_reified_or(beta1, alpha1);

  // construct T(x0,x1)
  std::cout << "c cl " << theory.n_clauses() << ": T(x0,x1)" << std::endl;
  literal_vec T_vec;

  // action implies preconds
  std::cout << "c cl " << theory.n_clauses() << ": act -> pre" << std::endl;
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
    for (HSPS::index_type i = 0; i < instance.actions[k].pre.size(); i++) {
      literal c = theory.new_var(false);
      theory.post_reified_or(theory.negate(a0[k]), x0[instance.actions[k].pre[i]], c);
      T_vec.push_back(c);
    }

  // action implies direct effects (adds)
  std::cout << "c cl " << theory.n_clauses() << ": act -> X add" << std::endl;
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
    for (HSPS::index_type i = 0; i < instance.actions[k].add.size(); i++) {
      literal c = theory.new_var(false);
      theory.post_reified_or(theory.negate(a0[k]), x1[instance.actions[k].add[i]], c);
      T_vec.push_back(c);
    }

  // action implies direct effects (dels)
  std::cout << "c cl " << theory.n_clauses() << ": act -> X ~del" << std::endl;
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
    for (HSPS::index_type i = 0; i < instance.actions[k].del.size(); i++) {
      literal c = theory.new_var(false);
      theory.post_reified_or(theory.negate(a0[k]), theory.negate(x1[instance.actions[k].del[i]]), c);
      T_vec.push_back(c);
    }

  // action mutex constraints
  std::cout << "c cl " << theory.n_clauses() << ": action exclusion"
	    << std::endl;
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
    for (HSPS::index_type j = k + 1; j < instance.n_actions(); j++)
      if (!instance.non_interfering(k, j)) {
	literal c = theory.new_var(false);
	theory.post_reified_or(theory.negate(a0[k]), theory.negate(a0[j]), c);
	T_vec.push_back(c);
      }

  // positive persistence
  std::cout << "c cl " << theory.n_clauses() << ": positive persistence"
	    << std::endl;
  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++) {
    literal_vec c1(instance.atoms[i].del_by.size() + 2, 0);
    c1[0] = theory.negate(x0[i]);
    c1[1] = x1[i];
    for (HSPS::index_type j = 0; j < instance.atoms[i].del_by.size(); j++)
      c1[j + 2] = a0[instance.atoms[i].del_by[j]];
    literal c = theory.new_var(false);
    theory.post_reified_or(c1, c);
    T_vec.push_back(c);
  }

  // negative persistence
  std::cout << "c cl " << theory.n_clauses() << ": negative persistence"
	    << std::endl;
  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++) {
    literal_vec c1(instance.atoms[i].add_by.size() + 2, 0);
    c1[0] = x0[i];
    c1[1] = theory.negate(x1[i]);
    for (HSPS::index_type j = 0; j < instance.atoms[i].add_by.size(); j++)
      c1[j + 2] = a0[instance.atoms[i].add_by[j]];
    literal c = theory.new_var(false);
    theory.post_reified_or(c1, c);
    T_vec.push_back(c);
  }

  std::cout << "c var " << theory.n_variables() << ": T01" << std::endl;
  literal T01 = theory.new_var(false);
  theory.post_reified_and(T_vec, T01);

  // main loop clause
  std::cout << "c cl " << theory.n_clauses() << ": loop clause" << std::endl;
  theory.post(theory.negate(alpha0), theory.negate(T01), alpha1);

  // alpha(x1) -> not G(x1)
  std::cout << "c cl " << theory.n_clauses() << ": alpha(x1) -> ~goal"
	    << std::endl;
  literal_vec notg1;
  notg1.push_back(theory.negate(alpha1));
  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++)
    if (instance.atoms[i].goal)
      notg1.push_back(theory.negate(x1[i]));
  theory.post(notg1);

  // filter alpha term 0 with init state
  std::cout << "c cl " << theory.n_clauses() << ": alpha(x0) true in init"
	    << std::endl;
  for (HSPS::index_type i = 0; i < instance.n_atoms(); i++) {
    if (instance.atoms[i].init) {
      theory.post(theory.negate(NOT_Pj_IN_ALPHAi(0, i)));
    }
    else {
      theory.post(theory.negate(Pj_IN_ALPHAi(0, i)));
    }
  }

  // all done!
  std::cout << "c cl " << theory.n_clauses() << ": end" << std::endl;
  theory.write_qdimacs(std::cout);
}
