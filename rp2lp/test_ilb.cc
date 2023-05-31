
#include "problem.h"
#include "preprocess.h"
#include "parser.h"
#include "cost_table.h"
#include "ilb.h"
#include "plans.h"
#include "enumerators.h"


//#define MEASURE_DENSITY
//#define UIPC_OUTPUT

// void compile_cadd_action
// (HSPS::Instance& instance,
//  HSPS::index_type act,
//  const HSPS::index_set_vec& ces,
//  const HSPS::bool_vec& sel)
// {
//   HSPS::Instance::Action& a = instance.copy_action(act);
//   for (HSPS::index_type i = 0; i < ces.size(); i++)
//     if (sel[i]) {
//       for (HSPS::index_type j = 0; j < ces[i].size(); j++) {
// 	a.pre.insert(instance.actions[act].cadd[ces[i][j]].antecedent);
// 	a.add.insert(instance.actions[act].cadd[ces[i][j]].consequent);
//       }
//     }
// }
//
// void compile_cadd_rec
// (HSPS::Instance& instance,
//  HSPS::index_type act,
//  const HSPS::graph& imp_tree,
//  const HSPS::index_set_vec& ces,
//  const HSPS::index_vec& node_level,
//  const HSPS::index_type max_level,
//  const HSPS::index_type cur_level,
//  HSPS::bool_vec& sel,
//  HSPS::Stopwatch& stats)
// {
//   HSPS::index_set choices_at_cur_level;
//   for (HSPS::index_type k = 0; k < imp_tree.size(); k++)
//     if ((node_level[k] == cur_level) && !sel[k])
//       choices_at_cur_level.insert(k);
//   if (choices_at_cur_level.empty()) {
//     assert(sel.count(true) > 0);
//     compile_cadd_action(instance, act, ces, sel);
//   }
//   else {
//     HSPS::SubsetEnumerator se(choices_at_cur_level.size());
//     HSPS::bool_vec new_sel(sel);
//     bool more = se.first();
//     while (more && !stats.break_signal_raised()) {
//       new_sel.assign_copy(sel);
//       for (HSPS::index_type i = 0; i < choices_at_cur_level.size(); i++)
// 	if (se.current_set()[i]) {
// 	  new_sel[choices_at_cur_level[i]] = true;
// 	  new_sel.insert(imp_tree.successors(choices_at_cur_level[i]));
// 	}
//       if (cur_level == max_level) {
// 	if (new_sel.count(true) > 0) {
// 	  compile_cadd_action(instance, act, ces, new_sel);
// 	}
//       }
//       else {
// 	compile_cadd_rec(instance, act, imp_tree, ces, node_level, max_level,
// 			 cur_level + 1, new_sel, stats);
//       }
//       more = se.next();
//     }
//   }
// }
//
// void compile_cadd
// (HSPS::Instance& instance, HSPS::Stopwatch& stats)
// {
//   HSPS::index_type n = instance.n_actions();
//   for (HSPS::index_type k = 0; k < n; k++)
//     if (!instance.actions[k].cadd.empty()) {
//       std::cerr << "compiling " << instance.actions[k].cadd.size()
// 		<< " ce's of action " << instance.actions[k].name
// 		<< std::endl;
//       HSPS::index_type m = instance.actions[k].cadd.size();
//       HSPS::graph imp_graph(m);
//       for (HSPS::index_type i = 0; i < m; i++)
// 	for (HSPS::index_type j = i + 1; j < m; j++) {
// 	  if (instance.actions[k].cadd[i].antecedent.
// 	      contains(instance.actions[k].cadd[j].antecedent))
// 	    imp_graph.add_edge(i, j);
// 	  if (instance.actions[k].cadd[j].antecedent.
// 	      contains(instance.actions[k].cadd[i].antecedent))
// 	    imp_graph.add_edge(j, i);
// 	}
//       // std::cerr << "imp. graph = " << imp_graph << std::endl;
//       imp_graph.strongly_connected_components();
//       std::cerr << imp_graph.n_components() << " components" << std::endl;
//       HSPS::graph imp_tree;
//       imp_graph.component_tree(imp_tree);
//       // std::cerr << "tree = " << imp_tree << std::endl;
//       HSPS::index_set_vec ces;
//       imp_graph.component_node_sets(ces);
//       // std::cerr << "ce's = " << ces << std::endl;
//       assert(ces.size() == imp_tree.size());
//       HSPS::index_vec node_level;
//       bool ok = imp_tree.assign_node_level_top_down(node_level);
//       assert(ok);
//       // std::cerr << "levels = " << node_level << std::endl;
//       HSPS::index_type max_level = HSPS::index_vec_util::max(node_level);
//       assert(max_level != HSPS::no_such_index);
//       imp_tree.transitive_closure();
//       HSPS::bool_vec sel(false, imp_tree.size());
//       compile_cadd_rec(instance, k, imp_tree, ces, node_level, max_level,
// 		       0, sel, stats);
//       if (stats.break_signal_raised()) {
// 	return;
//       }
//       std::cerr << "now " << instance.n_actions() << " actions" << std::endl;
//     }
//   std::cerr << "re-cross referencing..." << std::endl;
//   instance.clear_cross_reference();
//   instance.cross_reference();
//   std::cout << ";; full conditional effects compilation: "
// 	    << n << " => " << instance.n_actions()
// 	    << " actions" << std::endl;
// }

void addedge(int **g, int *gs, int *gn, int s, int t, int check){
  if (check==1)
    for (int i=0; i<gn[s]; i++)
      if (g[s][i]==t)
        return;
      if (gs[s]==gn[s]){
        g[s]= (int*) realloc(g[s],gs[s]*2*sizeof(int) );
        gs[s]=gs[s]*2;
        
      }
      g[s][gn[s]]=t;
      gn[s]++;
}



class ExpandWithSymmetry : public HSPS::ILB::ConflictModifier
{
  const HSPS::lvector<HSPS::mapping>& symgen;
public:
  ExpandWithSymmetry(const HSPS::lvector<HSPS::mapping>& sg);
  virtual ~ExpandWithSymmetry();
  virtual void apply(HSPS::ILB& ilb, HSPS::index_set_vec& cs);
};

ExpandWithSymmetry::ExpandWithSymmetry
  (const HSPS::lvector<HSPS::mapping>& sg)
  : symgen(sg)
{
  // done
}

ExpandWithSymmetry::~ExpandWithSymmetry()
{
  // done
}

void ExpandWithSymmetry::apply
  (HSPS::ILB& ilb, HSPS::index_set_vec& cs)
{
  std::cerr << "EWS: cs in = " << cs << std::endl;
  HSPS::index_type p = 0;
  while (p < cs.size()) {
    for (HSPS::index_type k = 0; k < symgen.size(); k++) {
      HSPS::index_set s(cs[p]);
      s.remap(symgen[k]);
      HSPS::index_type i = cs.first(s);
      if (i == HSPS::no_such_index)
        cs.append(s);
    }
    p += 1;
  }
  std::cerr << "EWS: cs out = " << cs << std::endl;
}
int probnum=0;

int main(int argc, char *argv[]) {
  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  int  verbose_level = 1;
  int  crp_strategy = 2;
  bool opt_rdc = false;
  bool opt_prep = true;
  bool opt_symmetry = false;
  bool opt_unit_cost = false;
  bool opt_max_with_h1 = true;
  bool opt_hplus = false;
  bool opt_hplus_with_selected_actions = false;
  // bool opt_one = false;
  bool opt_ce = false;
  bool opt_ce_ti = false;
  bool opt_cc = false;
  bool opt_print_plan = true;
  bool opt_cbd = false;
  bool opt_no_R1 = false;
  bool opt_no_R2 = true;
  bool opt_ila_apx = false;
  bool opt_ila_sat = false;
  bool opt_ila_lmc = false;
  bool opt_zero_cost_fill = false;
  bool opt_hs_split = false;
  bool opt_hs_dom = false;
  bool opt_use_inc = false;
  bool opt_rp_bb = false;
  bool opt_compare_with_lmcut = false;
  bool opt_val = false;
  bool opt_test = false;
  bool opt_wsat = false;
  HSPS::count_type wsat_limit = 0;
  bool opt_t0_hack = false;
  HSPS::index_type opt_add_all = 0;
  
  NTYPE ub = POS_INF;
  
  long          time_limit = 0;
  long          ILB_time_limit = 0;
  unsigned long memory_limit = 0;
  
  HSPS::Statistics  stats;
  stats.enable_interrupt();
  stats.start();
  HSPS::Statistics  read_stats(&stats);
  HSPS::Statistics  prep_stats(&stats);
  HSPS::Statistics  h_stats(&stats);
  HSPS::Statistics  ilb_stats(&stats);
  HSPS::Statistics  lmc_stats(&stats);
  
  HSPS::Parser* reader = new HSPS::Parser(symbols);
  
  // set dba semantics on by default:
  HSPS::PDDL_Base::del_before_add_semantics = true;
  // ignore undefined fluent error, assume them to be zero
  HSPS::PDDL_Base::use_default_function_value = true;
  HSPS::PDDL_Base::default_function_value = ZERO;
  
  // set make_types off by default (it's buggy):
  HSPS::PDDL_Base::make_types_from_static_predicates = false;
  
  for (int k = 1; k < argc; k++) {
	if ((strcmp(argv[k],"-n") == 0) && (k < argc - 1)) {
      probnum = atoi(argv[++k]);
    } 
    else if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      verbose_level = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-no-warnings") == 0) {
      HSPS::PDDL_Base::warning_level = 0;
    }
    else if ((strcmp(argv[k],"-c") == 0) && (k < argc - 1)) {
      crp_strategy = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-u") == 0) {
      opt_unit_cost = true;
    }
    else if (strcmp(argv[k],"-t0") == 0) {
      opt_t0_hack = true;
    }
    else if ((strcmp(argv[k],"-all") == 0) && (k < argc - 1)) {
      opt_add_all = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-h") == 0) {
      opt_hplus = true;
      opt_print_plan = false;
    }
    else if (strcmp(argv[k],"-l") == 0) {
      opt_hplus = true;
      opt_hplus_with_selected_actions = true;
      opt_print_plan = false;
    }
    // else if (strcmp(argv[k],"-1") == 0) {
    //   opt_one = true;
    // }
    else if (strcmp(argv[k],"-ce") == 0) {
      opt_ce = true;
      HSPS::PDDL_Base::compile_away_conditional_effects = false;
    }
    else if (strcmp(argv[k],"-ceti") == 0) {
      opt_ce = true;
      HSPS::PDDL_Base::compile_away_conditional_effects = false;
      opt_ce_ti = true;
    }
    else if (strcmp(argv[k],"-cc") == 0) {
      opt_cc = true;
      HSPS::PDDL_Base::compile_away_conditional_effects = false;
    }
    else if (strcmp(argv[k],"-no-h1") == 0) {
      opt_max_with_h1 = false;
    }
    else if (strcmp(argv[k],"-no-R1") == 0) {
      opt_no_R1 = true;
    }
    else if (strcmp(argv[k],"-R2") == 0) {
      opt_no_R2 = false;
    }
    else if (strcmp(argv[k],"-no-prep") == 0) {
      opt_prep = false;
    }
    else if (strcmp(argv[k],"-sym") == 0) {
      opt_symmetry = true;
    }
    else if (strcmp(argv[k],"-rdc") == 0) {
      opt_rdc = true;
    }
    else if (strcmp(argv[k],"-p") == 0) {
      opt_print_plan = true;
    }
    else if (strcmp(argv[k],"-rpb") == 0) {
      opt_rp_bb = true;
    }
    else if (strcmp(argv[k],"-compare-with-lmcut") == 0) {
      opt_compare_with_lmcut = true;
    }
    else if (strcmp(argv[k],"-cbd") == 0) {
      opt_cbd = true;
      ub = 1;
    }
    else if (strcmp(argv[k],"-val") == 0) {
      opt_val = true;
    }
    else if (strcmp(argv[k],"-test") == 0) {
      opt_test = true;
    }
    else if (strcmp(argv[k],"-wsat") == 0) {
      opt_wsat = true;
    }
    else if ((strcmp(argv[k],"-wsat-limit") == 0) && (k < argc - 1)) {
      opt_wsat = true;
      wsat_limit = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-a") == 0) {
      opt_ila_apx = true;
    }
    else if (strcmp(argv[k],"-s") == 0) {
      opt_ila_sat = true;
    }
    else if (strcmp(argv[k],"-lmcut") == 0) {
      opt_ila_lmc = true;
    }
    else if (strcmp(argv[k],"-z") == 0) {
      opt_zero_cost_fill = true;
    }
    else if (strcmp(argv[k],"-x") == 0) {
      opt_use_inc = true;
    }
    else if ((strcmp(argv[k],"-b") == 0) && (k < argc - 1)) {
      ub = A_TO_N(argv[++k]);
    }
    else if ((strcmp(argv[k],"-t") == 0) && (k < argc - 1)) {
      time_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-tilb") == 0) && (k < argc - 1)) {
      ILB_time_limit = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-y") == 0) && (k < argc - 1)) {
      memory_limit = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-split") == 0) {
      opt_hs_split = true;
    }
    else if (strcmp(argv[k],"-nosplit") == 0) {
      opt_hs_split = false;
    }
    else if (strcmp(argv[k],"-dom") == 0) {
      opt_hs_dom = true;
    }
    else if (strcmp(argv[k],"-dba-semantics") == 0) {
      HSPS::PDDL_Base::del_before_add_semantics = true;
    }
    else if (strcmp(argv[k],"-no-make-types") == 0) {
      HSPS::PDDL_Base::make_types_from_static_predicates = false;
    }
    else if (*argv[k] != '-') {
      read_stats.start();
      reader->read(argv[k], false);
      read_stats.stop();
    }
  }
  
  if (time_limit > 0) stats.enable_time_out(time_limit);
  if (ILB_time_limit > 0) ilb_stats.enable_time_out(ILB_time_limit);
  if (memory_limit > 0) stats.enable_memory_limit(memory_limit);
  
  //HSPS::Heuristic::default_trace_level = verbose_level;
  HSPS::Instance::default_trace_level = verbose_level - 2;
  //HSPS::Preprocessor::default_trace_level = verbose_level;
  if (verbose_level <= 0) HSPS::PDDL_Base::warning_level = 0;
  //if (verbose_level > 1) HSPS::PDDL_Base::write_info = true;
  
  prep_stats.start();
  std::cerr << "instantiating..." << std::endl;
  HSPS::PDDL_Base::name_instance_by_problem_file = true;
  HSPS::Instance instance;
  reader->instantiate(instance);
  // special fix if we're trying to be PDDL compliant: we don't want
  // actions to add anything that is in their preconditions (because
  // we're only considering sequential plans, so "temporary deletes"
  // don't matter)
  if (HSPS::PDDL_Base::del_before_add_semantics) {
    for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
      instance.actions[k].add.subtract(instance.actions[k].pre);
  }
  else {
    for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
      if (instance.actions[k].add.
            have_common_element(instance.actions[k].del)) {
        std::cerr << "serious warning: non-empty add-delete intersection"
                  << std::endl;
        instance.print_action(std::cerr, instance.actions[k]);
        exit(1);
      }
  }
  
  HSPS::lvector<HSPS::mapping> symgen;
  if (opt_symmetry) {
    reader->export_generators(instance, symgen);
    if (verbose_level > 0) {
      std::cerr << symgen.size() << " symmetry generators found in input"
                << std::endl;
    }
  }
  
  HSPS::Preprocessor prep(instance, prep_stats);
  if (opt_ce || opt_cc) {
    std::cerr << "preprocessing (with ce)..." << std::endl;
    if (verbose_level > 1) {
      std::cerr << instance.n_atoms() << " atoms and "
                << instance.n_actions() << " actions"
                << std::endl;
      if (verbose_level > 3) {
        instance.print(std::cerr);
        prep.trace_level = 2;
      }
    }
    if (opt_prep) {
      prep.preprocess_with_ce(false, true);
      if (verbose_level > 3) {
        instance.print(std::cerr);
      }
    }
    if (opt_cc) {
      // full compilation of conditional effects...
      std::cerr << "compiling all conditional effects..." << std::endl;
      instance.compile_conditional_effects(stats, !opt_hplus);
      if (stats.break_signal_raised()) {
        std::cout << ";; compilation interrupted! "
                  << instance.n_actions() << " actions, "
                  << stats.peak_memory() << "k, "
                  << stats.total_time() << " seconds"
                  << std::endl;
        exit(0);
      }
      if (verbose_level > 3) {
        instance.print(std::cerr);
      }
    }
  }
  else if (opt_prep) {
    std::cerr << "preprocessing..." << std::endl;
    prep.preprocess(false);
    prep.compute_irrelevant_atoms();
    prep.remove_irrelevant_atoms();
  }
  if (!instance.cross_referenced()) {
    std::cerr << "re-cross referencing..." << std::endl;
    instance.cross_reference();
  }
  
  if (opt_prep && opt_symmetry) {
    for (HSPS::index_type k = 0; k < symgen.size(); k++) {
      HSPS::mapping new_sg;
      bool ok = HSPS::mapping::remap(symgen[k], prep.atom_map, new_sg);
      assert(ok);
      symgen[k] = new_sg;
    }
    if (verbose_level > 1) {
      for (HSPS::index_type k = 0; k < symgen.size(); k++) {
        std::cerr << "generator #" << k << " (post-prep):" << std::endl;
        for (HSPS::index_type i = 0; i < instance.n_atoms(); i++)
          if (symgen[k][i] != i) {
            std::cerr << "   " << i << ".";
            instance.atoms[i].name->write(std::cerr);
            std::cerr << " -> " << symgen[k][i] << ".";
            instance.atoms[symgen[k][i]].name->write(std::cerr);
            std::cerr << std::endl;
          }
      }
    }
  }
  
  if (opt_add_all > 1) {
    if (verbose_level > 3) {
      instance.print(std::cerr);
    }
    HSPS::rule_set ma_map;
    HSPS::StaticMutex* mx = (opt_use_inc ? prep.inconsistency() : 0);
    for (HSPS::index_type s = 2; s <= opt_add_all; s++) {
      std::cerr << "adding all " << s << "-conjunctions..." << std::endl;
      HSPS::mSubsetEnumerator e(instance.n_atoms(), s);
      bool more = e.first();
      while (more) {
        HSPS::index_set c;
        e.current_set(c);
        std::cerr << c << std::endl;
        if (opt_ce) {
          instance.create_meta_atom_with_ce(c, ma_map, mx);
        }
        else {
          instance.create_meta_atom(c, ma_map, mx);
        }
        if (stats.break_signal_raised()) exit(1);
        more = e.next();
      }
    }
    if (!instance.cross_referenced()) {
      std::cerr << "re-cross referencing..." << std::endl;
      instance.cross_reference();
    }
    if (verbose_level > 3) {
      instance.print(std::cerr);
    }
  }
  
  if (opt_t0_hack) {
    for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
      HSPS::PDDL_Name* n =
        (HSPS::PDDL_Name*)instance.actions[k].name->cast_to("PDDL_Name");
      if (n) {
        if (strncmp(n->symbol()->print_name, "merge_", 6) == 0) {
          instance.actions[k].cost = 0;
          std::cerr << "set cost of " << instance.actions[k].name
                    << " to zero" << std::endl;
        }
        else if (strcmp(n->symbol()->print_name, "make_end_disj_goal") == 0) {
          instance.actions[k].cost = 0;
          std::cerr << "set cost of " << instance.actions[k].name
                    << " to zero" << std::endl;
        }
      }
    }
  }
  
  prep_stats.stop();
  std::cerr << "instance " << instance.name << " built in "
            << prep_stats.total_time() << " seconds" << std::endl;
  HSPS::index_type max_ces = 0;
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
    if (instance.actions[k].cadd.size() > max_ces)
      max_ces = instance.actions[k].cadd.size();
    std::cout << ";; "<< instance.name << ": "
              << instance.n_atoms() << " atoms, "
              << instance.n_actions() << " actions, "
              << instance.n_conditional_effects() << " ce's (max "
              << max_ces << " / action)"
              << std::endl;
    
    if (verbose_level > 3) {
      instance.print(std::cerr);
    }
    
    if (opt_val) {
      for (HSPS::index_type p = 0; p < reader->n_plans(); p++) {
        std::cerr << "validating plan #" << p << std::endl;
        HSPS::ActionSequence plan;
        reader->export_plan(p, instance, plan);
        HSPS::bool_vec s(instance.init_atoms, instance.n_atoms());
        bool ok = true;
        for (HSPS::index_type k = 0; (k < plan.length()) && ok; k++) {
          assert(plan[k] < instance.n_actions());
          const HSPS::Instance::Action& act = instance.actions[plan[k]];
          std::cerr << "next action: " << act.name << std::endl;
          if (!s.contains(act.pre)) {
            std::cerr << "plan failed! unsatisfied preconditions:" << std::endl;
            for (HSPS::index_type i = 0; i < act.pre.size(); i++)
              if (!s[act.pre[i]])
                std::cerr << instance.atoms[act.pre[i]].name << std::endl;
              ok = false;
          }
          else {
            if (!opt_hplus) {
              for (HSPS::index_type i = 0; i < act.del.size(); i++) {
                s[act.del[i]] = false;
                std::cerr << " del: " << instance.atoms[act.del[i]].name
                          << std::endl;
              }
            }
            for (HSPS::index_type i = 0; i < act.add.size(); i++) {
              s[act.add[i]] = true;
              std::cerr << " add: " << instance.atoms[act.add[i]].name
                        << std::endl;
            }
            if (!opt_hplus) {
              for (HSPS::index_type i = 0; i < act.cdel.size(); i++)
                if (s.contains(act.cdel[i].antecedent)) {
                  HSPS::index_type j = act.cdel[i].consequent;
                  s[j] = false;
                  std::cerr << " del: " << instance.atoms[j].name
                            << " because ";
                  instance.write_atom_set(std::cerr, act.cdel[i].antecedent);
                  std::cerr << std::endl;
                }
            }
            for (HSPS::index_type i = 0; i < act.cadd.size(); i++)
              if (s.contains(act.cadd[i].antecedent)) {
                HSPS::index_type j = act.cadd[i].consequent;
                s[j] = true;
                std::cerr << " add: " << instance.atoms[j].name
                          << " because ";
                instance.write_atom_set(std::cerr, act.cadd[i].antecedent);
                std::cerr << std::endl;
              }
          }
        }
        if (ok) {
          if (!s.contains(instance.goal_atoms)) {
            std::cerr << "goals not achieved!" << std::endl;
            for (HSPS::index_type i = 0; i < instance.goal_atoms.size(); i++)
              if (!s[instance.goal_atoms[i]])
                std::cerr << instance.atoms[instance.goal_atoms[i]].name
                          << std::endl;
                ok = false;
          }
        }
        if (ok) {
          std::cerr << "plan successful!" << std::endl;
        }
      }
      exit(0);
    }
    
    HSPS::ACF* costf = (opt_unit_cost ?
                          (HSPS::ACF*)new HSPS::UnitACF() :
                          (HSPS::ACF*)new HSPS::CostACF(instance));
    NTYPE v_max = ZERO;
    if (opt_max_with_h1) {
      if (opt_ce) {
        HSPS::Instance irce(instance);
        irce.relax_conditional_effects(false);
        HSPS::CostTable* h1 = new HSPS::CostTable(irce, h_stats);
        if (opt_unit_cost)
          h1->compute_H1(HSPS::UnitACF());
        else
          h1->compute_H1(HSPS::CostACF(irce));
        //h1->write_pddl(std::cerr, irce);
        v_max = h1->eval(irce.goal_atoms);
        delete h1;
      }
      else {
        HSPS::CostTable* h1 = new HSPS::CostTable(instance, h_stats);
        h1->compute_H1(*costf);
        v_max = h1->eval(instance.goal_atoms);
        delete h1;
      }
      std::cerr << "h^1 = " << PRINT_NTYPE(v_max) << std::endl;
    }
    
    prep_stats.start();
    HSPS::StaticMutex* mx = (opt_use_inc ? prep.inconsistency() : 0);
    prep_stats.stop();
    HSPS::ILB ilb(instance, *costf, mx, ilb_stats);
    
    
    
    ilb.SstepOrder =(int **)malloc(instance.atoms.size() * sizeof(int *));
	ilb.probnum = probnum;
    for (int prop1=0;prop1<instance.atoms.size();prop1++){
      ilb.SstepOrder[prop1]=(int *)malloc(instance.atoms.size() * sizeof(int ));
      for (int prop2=0;prop2<instance.atoms.size();prop2++)
        ilb.SstepOrder[prop1][prop2]=-1;
    }
    
    for (int a=0; a<instance.actions.size(); a++){
      for (int i = 0; i < instance.actions[a].pre.size(); i++){
        
        for (int j = 0; j < instance.actions[a].add.size(); j++){
          
          if (instance.actions[a].pre[i]!=instance.actions[a].add[j]){
            ilb.SstepOrder[instance.actions[a].pre[i]][instance.actions[a].add[j]]=1;
            
          }
        }
        
      }
      
    }
    
    if(0)
    {
      ilb.adjorddes =(int **)malloc(instance.atoms.size() * sizeof(int *));
      ilb.nadjorddes = (int *)malloc(instance.atoms.size() * sizeof(int ));
      ilb.sadjorddes = (int *)malloc(instance.atoms.size() * sizeof(int ));
      ilb.adjorddesrev =(int **)malloc(instance.atoms.size() * sizeof(int *));;
      ilb.nadjorddesrev = (int *)malloc(instance.atoms.size() * sizeof(int ));
      ilb.sadjorddesrev = (int *)malloc(instance.atoms.size() * sizeof(int ));
      for (int i=0; i<instance.atoms.size(); i++){
        ilb.adjorddes[i] = (int *)malloc(100*sizeof(int ));
        ilb.adjorddesrev[i] = (int *)malloc(100*sizeof(int ));
        ilb.nadjorddes[i]=0;
        ilb.sadjorddes[i]=100;
        ilb.nadjorddesrev[i]=0;
        ilb.sadjorddesrev[i]=100;            
      }
      
      
      for (int i=0; i<instance.atoms.size(); i++)
        for (int j=0; j<instance.atoms.size(); j++){
          if (i==j)
            continue;
          if (ilb.SstepOrder[i][j]==1) {
            //exit(0);
            addedge(ilb.adjorddes,ilb.sadjorddes,ilb.nadjorddes,i,j,1);
            addedge(ilb.adjorddesrev,ilb.sadjorddesrev,ilb.nadjorddesrev,j,i,1);
          }
          
          
          
        }
        
        ilb.triangles = (int **)malloc(ilb.striangles * sizeof(int *));
      for (int i=0; i<ilb.striangles; i++){
        ilb.triangles[i] = (int *)malloc(3*sizeof(int ));
      }
      
      int degreesOut[instance.atoms.size()];
      int degreesIn[instance.atoms.size()];
      for (int i=0; i<instance.atoms.size(); i++){
        degreesOut[i]=0;
        degreesIn[i]=0;
      }
      for (int i=0; i<instance.atoms.size(); i++){
        
        degreesOut[i]=ilb.nadjorddes[i];
        degreesIn[i]=ilb.nadjorddesrev[i];
        
      }  
      int delta=-1;
      int ts=0;
      
      
      //RUN(example1);

      for (int i=0; i<instance.atoms.size(); i++){
        
        long int mindegree=1000000000;
        int minv =-1;
        for (int j=0; j<instance.atoms.size(); j++){
          if (degreesOut[j]==-1)
            continue;
          if (degreesIn[j]*degreesOut[j]<mindegree){
            mindegree=degreesIn[j]*degreesOut[j];
            minv=j;
          }
          //minv=i;
        }
        ts+=mindegree;
        //printf("\n %d",  degreesOut[minv]);
        if (degreesOut[minv]>delta)
          delta=degreesOut[minv];
        
        
        
        for (int j=0; j<ilb.nadjorddesrev[minv]; j++){
          if ((ilb.SstepOrder[ilb.adjorddesrev[minv][j]][minv]!=1)) 
            continue;
          if (ilb.adjorddesrev[minv][j]==minv)
            continue;
          
          for (int k=0; k<ilb.nadjorddes[minv]; k++){
            if (ilb.adjorddes[minv][k]==ilb.adjorddesrev[minv][j])
              continue;
            if ((ilb.SstepOrder[minv][ilb.adjorddes[minv][k]]!=1) )
              continue;
            if (ilb.adjorddes[minv][k]==minv)
              continue;
            
            
            
            //exit(0);
            //if (k==j)
            //continue;
            if ((ilb.SstepOrder[ilb.adjorddesrev[minv][j]][minv]==1) && (ilb.SstepOrder[minv][ilb.adjorddes[minv][k]]==1)){
              //exit(0);
              
              
              ilb.triangles[ilb.ntriangles][0] = minv;
              ilb.triangles[ilb.ntriangles][1] = ilb.adjorddesrev[minv][j];
              ilb.triangles[ilb.ntriangles][2] = ilb.adjorddes[minv][k];
              ilb.ntriangles++;
              if (ilb.ntriangles==ilb.striangles){
                
                int **triangles2 = (int **)malloc(ilb.striangles * 2 * sizeof(int *));
                for (int i=0; i<ilb.striangles * 2; i++){
                  triangles2[i] = (int *)malloc(3*sizeof(int ));
                  
                }
                for (int i=0; i<ilb.striangles ; i++){
                  triangles2[i][0] =  ilb.triangles[i][0];
                  triangles2[i][1] =  ilb.triangles[i][1];
                  triangles2[i][2] =  ilb.triangles[i][2];
                  free(ilb.triangles[i]);
                }
                free(ilb.triangles);
                ilb.striangles *= 2;
                ilb.triangles = triangles2;
                
              }
              
              
              if ((ilb.SstepOrder[ilb.adjorddesrev[minv][j]][ilb.adjorddes[minv][k]]!=1)){
                
                ilb.SstepOrder[ilb.adjorddesrev[minv][j]][ilb.adjorddes[minv][k]]=1;
                addedge(ilb.adjorddesrev,ilb.sadjorddesrev,ilb.nadjorddesrev,ilb.adjorddes[minv][k],ilb.adjorddesrev[minv][j],1);
                addedge(ilb.adjorddes,ilb.sadjorddes,ilb.nadjorddes,ilb.adjorddesrev[minv][j],ilb.adjorddes[minv][k],1);
                degreesOut[ilb.adjorddesrev[minv][j]]++;
                degreesIn[ilb.adjorddes[minv][k]]++;
              }
              continue;
            } 
            
          }
        }                      
        
        
        
        
        
        for (int k=0; k<instance.atoms.size(); k++){
          // if ((SstepOrder[minv][k]!=1) && (SdestroyedTemp[minv][k]!=1) )
          // continue;
          if (k==minv)
            continue;
          
          if ((ilb.SstepOrder[minv][k]==1) ){
            ilb.SstepOrder[minv][k]=-1;
            degreesOut[minv]--;
            degreesIn[k]--;
          }
          
          if ((ilb.SstepOrder[k][minv]==1)){
            ilb.SstepOrder[k][minv]=-1;
            degreesOut[k]--;
            degreesIn[minv]--;
          }
          
        }
        if ((degreesOut[minv]!=0) || (degreesIn[minv]!=0))
        {
          printf("\n %d %d %d DEGREES ERROR2\n",degreesOut[minv], degreesIn[minv],i);
          exit(0);
        }          
        degreesOut[minv]=-1;
        degreesIn[minv]=-1;  
        
      }
      
      
      
      
      
      
    }
    
    
    
    
    
    
    ilb.verbose_level = verbose_level;
    ilb.check_rp_strategy = crp_strategy;
    ilb.remove_dominated_conditions = opt_rdc;
    ilb.prune_relaxed_irrelevant = !opt_no_R1;
    ilb.prune_relaxed_dominated = !opt_no_R2;
    ilb.ILA_use_approximate = opt_ila_apx;
    ilb.ILA_use_saturation = opt_ila_sat;
    ilb.ILA_use_lmcut = opt_ila_lmc;
    ilb.zero_cost_fill = opt_zero_cost_fill;
    ilb.hitting_set_use_split = opt_hs_split;
    ilb.hitting_set_use_dominance = opt_hs_dom;
    if (opt_symmetry) {
      ExpandWithSymmetry* ews = new ExpandWithSymmetry(symgen);
      ilb.conflict_mod = ews;
    }
    HSPS::Plan* sol = NULL;
    
    if (opt_test) {
      assert(reader->n_plans() > 0);
      HSPS::ActionSequence plan;
      reader->export_plan(0, instance, plan);
      ilb.test(plan);
      exit(0);
    }
    
    if (opt_wsat) {
      ilb.save_hplus_wcnf(instance.init_atoms, instance.goal_atoms,
                          reader->problem_file_basename(), wsat_limit);
      exit(0);
    }
    
    if (opt_rp_bb) {
      std::cerr << "computing h^+..." << std::endl;
      NTYPE v = ilb.compute_relaxed_plan_BDGBT(instance.init_atoms,
                                               instance.goal_atoms);
      // NTYPE v = ilb.compute_relaxed_plan_BB(instance.init_atoms,
      // 					  instance.goal_atoms);
      std::cerr << "h^+ = " << PRINT_NTYPE(v) << std::endl;
      v_max = MAX(v_max, v);
    }
    else if (opt_compare_with_lmcut) {
      HSPS::ForwardLMCut lmc(instance, instance.goal_atoms, *costf, lmc_stats);
      lmc_stats.start();
      NTYPE v_lmc = lmc.eval(instance.init_atoms);
      lmc_stats.stop();
      NTYPE v = ilb.hplus(instance.init_atoms, instance.goal_atoms,
                          v_lmc + 1, NULL);
      std::cout << "lmcut = " << PRINT_NTYPE(v_lmc)
                << ", " << lmc_stats.total_time() << " seconds" << std::endl;
      std::cout << "h^+ >= " << PRINT_NTYPE(v)
                << ", " << ilb_stats.total_time() << " seconds" << std::endl;
      exit(0);
    }
    else if (opt_hplus) {
      std::cerr << "computing h^+..." << std::endl;
      NTYPE v = 0;
      if (opt_hplus_with_selected_actions) {
        HSPS::name_vec action_names(0, 0);
        instance.action_names(action_names);
        HSPS::index_set_vec sets;
        HSPS::PDDL_Base::strict_set_export = true;
        reader->export_sets(action_names, sets);
        assert(sets.length() == 1);
        v = ilb.hplus_with_selected_actions(instance.init_atoms,
                                            instance.goal_atoms,
                                            sets[0],
                                                ub, (opt_print_plan ? &sol : NULL));
        if (verbose_level > 1) {
          ilb.validate_rp(ilb.relaxed_plan());
        }
      }
      else if (opt_ce_ti) {
        HSPS::index_vec rp;
        HSPS::index_set_vec rp_ce;
        v = ilb.hplus_with_ce_time_indexed(instance.init_atoms,
                                           instance.goal_atoms,
                                           ub, rp, rp_ce);
      }
      else if (opt_ce) {
        v = ilb.hplus_with_ce(instance.init_atoms, instance.goal_atoms,
                              ub, (opt_print_plan ? &sol : NULL));
      }
      else {
        v = ilb.hplus(instance.init_atoms, instance.goal_atoms,
                      ub, (opt_print_plan ? &sol : NULL));
        if (verbose_level > 1) {
          ilb.validate_rp(ilb.relaxed_plan());
        }
      }
      std::cerr << "h^+ = " << PRINT_NTYPE(v) << std::endl;
      v_max = MAX(v_max, v);
    }
    // else if (opt_one) {
    //   HSPS::name_vec atom_names(0, 0);
    //   instance.atom_names(atom_names);
    //   HSPS::index_set_vec cs_in;
    //   reader->export_sets(atom_names, cs_in);
    //   HSPS::index_set_vec cs_out;
    //   bool completed;
    //   v_max = ilb.hplusplus1(cs_in, cs_out, &sol, completed, true);
    //   if (completed && !cs_out.empty()) {
    //     std::cout << ";; new conflicts:" << std::endl;
    //     const HSPS::rule_set& ma_map = ilb.get_meta_atom_map();
    //     for (HSPS::index_type k = 0; k < cs_out.size(); k++) {
    // 	HSPS::index_set s(cs_out[k]);
    // 	ma_map.backchain_to_fixpoint(s);
    // 	std::cout << "(:set";
    // 	for (HSPS::index_type i = 0; i < s.size(); i++) {
    // 	  assert(s[i] < instance.n_atoms());
    // 	  std::cout << " " << instance.atoms[s[i]].name;
    // 	}
    // 	std::cout << ")" << std::endl;
    //     }
    //   }
    // }
    else {
      std::cerr << "computing h^++..." << std::endl;
      NTYPE v = 0;
      if (opt_ce) {
        v = ilb.hplusplus_with_ce(ub, &sol);
      }
      else {
        v = ilb.hplusplus(ub, &sol);
      }
      v_max = MAX(v_max, v);
    }
    
    std::cout << ";; highest lower bound = " << PRINT_NTYPE(v_max) << std::endl;
    if (opt_ce) {
      std::cout << ";; h+/ce calls = " << ilb.calls_to_hplus_ce
                << " (avg " << ilb.calls_to_hplus/(double)ilb.calls_to_hplus_ce
                << " iterations/call)"
                << std::endl;
      std::cout << ";; average ce's compiled = "
                << ilb.n_ce_compiled/(double)ilb.calls_to_hplus_ce
                << std::endl;
      std::cout << ";; average domain size increase (ce) = "
                << ilb.cce_actions_ratio/(double)ilb.calls_to_hplus_ce
                << std::endl;
      std::cout << ";; max ce's compiled / action = "
                << ilb.max_ce_compiled << std::endl;
    }
    std::cout << ";; h+ calls = " << ilb.calls_to_hplus << std::endl;
    std::cout << ";; average relevant actions = "
              << ilb.relevant_actions_ratio/(double)ilb.calls_to_hplus
              << std::endl;
    std::cout << ";; lms generated = " << ilb.calls_to_newlm <<" (avg "
              << PRINT_NTYPE(R_TO_N(ilb.calls_to_newlm,
    ilb.calls_to_hplus))
      << " lms/h+ call)" << std::endl;
    std::cout << ";; hs(opt) calls = " << ilb.calls_to_hs_opt
              << ", nodes = " << ilb.hs_nodes << " (avg "
              << PRINT_NTYPE(R_TO_N(ilb.hs_nodes, ilb.calls_to_hs_opt))
              << " nodes/call)" << std::endl;
    std::cout << ";; average width = "
              << PRINT_NTYPE(R_TO_N(ilb.sum_width, ilb.n_width))
              << std::endl;
    std::cout << ";; hs(apx) calls = " << ilb.calls_to_hs_apx << " ("
              << ilb.hs_apx_improve << " improved)" << std::endl;
#ifdef MEASURE_DENSITY
    std::cout << ";; average density = "
              << (ilb.sum_density/ilb.calls_to_hs_opt)*100
              << "%" << std::endl;
#endif
    std::cout << ";; hs lb1 = "
              << PRINT_NTYPE(R_TO_N(ilb.hs_lb1_max, ilb.hs_lb_calls)*100)
              << ", lb2 = "
              << PRINT_NTYPE(R_TO_N(ilb.hs_lb2_max, ilb.hs_lb_calls)*100)
              << ", lb3 = "
              << PRINT_NTYPE(R_TO_N(ilb.hs_lb3_max, ilb.hs_lb_calls)*100)
              << ", lb4 = "
              << PRINT_NTYPE(R_TO_N(ilb.hs_lb4_max, ilb.hs_lb_calls)*100)
              << std::endl;
    if (opt_hs_split) {
      std::cout << ";; " << PRINT_NTYPE(R_TO_N(ilb.hs_split1, ilb.hs_nodes)*100)
                << "% hs nodes split" << std::endl;
      std::cout << ";; " << PRINT_NTYPE(R_TO_N(ilb.hs_splits, ilb.hs_nodes)*100)
                << "% hs nodes properly split" << std::endl;
    }
    std::cout << ";; " << PRINT_NTYPE(R_TO_N(ilb.hs_branch, ilb.hs_nodes)*100)
              << "% hs nodes branched" << std::endl;
    std::cout << ";; " << PRINT_NTYPE(R_TO_N(ilb.hs_cache_hits, ilb.hs_cache_hits + ilb.hs_cache_miss)*100)
              << "% hs cache hits" << std::endl;
    std::cout << ";; average conflict set size = "
              << PRINT_NTYPE(R_TO_N(ilb.conflict_set_size, ilb.hpp_iterations))
              << std::endl;
    std::cout << ";; average domain size increase (pc) = x"
              << ilb.pc_actions_ratio/ilb.hpp_iterations
              << std::endl;
    std::cout << ";; time (reading input) = " << read_stats.total_time()
              << std::endl;
    std::cout << ";; time (preprocessing) = " << prep_stats.total_time()
              << std::endl;
    std::cout << ";; time (ILB) = " << ilb_stats.total_time() << std::endl;
    std::cout << ";; time (relevance analysis) = "
              << ilb.ra_stats.total_time() << std::endl;
    std::cout << ";; time (approx hitting set) = "
              << ilb.hsa_stats.total_time()
              << " (" << ilb.hsa_stats.total_time() / ilb.calls_to_hs_apx
              << " sec/call)" << std::endl;
    std::cout << ";; time (optimal hitting set) = "
              << ilb.hss_stats.total_time()
              << " (" << ilb.hss_stats.total_time() / ilb.calls_to_hs_opt
              << " sec/call)" << std::endl;
    std::cout << ";; time (lm generation) = " << ilb.lm_stats.total_time()
              << " (" << ilb.lm_stats.total_time() / ilb.calls_to_newlm
              << " sec/call)" << std::endl;
    if (opt_ce) {
      std::cout << ";; time (scheduling ce) = " << ilb.sce_stats.total_time()
                << " (" << ilb.sce_stats.total_time() / ilb.calls_to_hplus
                << " sec/iteration)" << std::endl;
    }
    std::cout << ";; time (rp normalisation) = "
              << ilb.rpn_stats.total_time() << std::endl;
    std::cout << ";; time (conflict extraction) = "
              << ilb.ce_stats.total_time() << std::endl;
    std::cout << ";; time (total) = " << stats.total_time()
              << ", peak memory = " << stats.peak_memory()
              << "k, flags = " << stats.flags() << " / " << ilb_stats.flags()
              << std::endl;
    
    if (sol != NULL) {
      HSPS::PrintActions ap(instance, std::cout, '\n', '\n');
      sol->output(ap);
      std::cout << std::endl;
    }
    
    else if (opt_cbd) {
      std::cout << "(define (problem conflict)" << std::endl;
      const HSPS::index_set_vec& lms = ilb.action_landmarks();
      for (HSPS::index_type k = 0; k < lms.size(); k++)
        if (!IS_ZERO(costf->min_cost(lms[k]))) {
          std::cout << "  (:set" << std::endl;
          for (HSPS::index_type i = 0; i < lms[k].size(); i++) {
            HSPS::index_type a = lms[k][i];
            if (instance.actions[a].src)
              a = *((HSPS::index_type*)instance.actions[a].src);
            std::cout << "    " << instance.actions[a].name << std::endl;
          }
          std::cout << "   )" << std::endl;
        }
        std::cout << " )" << std::endl;
    }
    
#ifdef UIPC_OUTPUT
    if (INFINITE(v_max) && (stats.flags() == 0)) {
      std::cout << "unsolvable" << std::endl;
    }
    else if (sol != NULL) {
      std::cout << "solvable" << std::endl;
    }
    else {
      std::cout << "unknown" << std::endl;
    }
#endif
    
    return 0;
}
