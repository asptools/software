
#include "problem.h"
#include "parser.h"
#include "exec.h"
#include "pop.h"

#include <fstream>
#include <sstream>

#define PDDL21_EPSILON R_TO_N(1, 100)
#define VERBOSE
// #define VERY_VERBOSE

int main(int argc, char *argv[]) {
  HSPS::Statistics stats;
  stats.enable_interrupt();
  stats.start();

  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  HSPS::Parser* reader = new HSPS::Parser(symbols);

  if (argc < 4) {
    std::cerr << argv[0] << " <domain> <problem> <plan>+" << std::endl;
    exit(1);
  }

  HSPS::PDDL_Base::action_vec ea;
  HSPS::cost_vec tp;

#ifndef VERY_VERBOSE
  HSPS::PDDL_Base::warning_level = 0;
#endif
  reader->read(argv[1], false);
  reader->read(argv[2], false);
  reader->detil(false, 0, ea, tp);
  for (HSPS::index_type k = 0; k < ea.size(); k++) {
    HSPS::StringTable::Cell* c = symbols.inserta(ea[k]->print_name);
    c->val = ea[k];
  }
  reader->post_process();

  for (HSPS::index_type k = 3; k < argc; k++) {
    reader->read(argv[k], false);
  }

  HSPS::Instance instance;
  reader->instantiate(instance);
  instance.cross_reference();

  HSPS::ScheduleSet plans(instance);
  for (HSPS::index_type k = 0; k < reader->n_plans(); k++) {
    HSPS::Schedule* s = new HSPS::Schedule(instance);
    bool ok = reader->export_plan(k, instance, *s);
    if (!ok) {
      std::cerr << "error: export plan " << k << " failed" << std::endl;
    }
    else {
      plans.add_schedule(s);
    }
  }

  for (HSPS::index_type k = 0; k < plans.size(); k++) {
#ifdef VERY_VERBOSE
    std::cerr << "plan " << k << ":" << std::endl;
    plans[k]->write(std::cerr);
#endif

    HSPS::graph prec;
    plans[k]->deorder_non_temporal(prec);
    assert(prec.acyclic());

#ifdef VERY_VERBOSE
    HSPS::graph p1(prec);
    p1.transitive_reduction();
    HSPS::name_vec step_act_names;
    plans[k]->step_action_names(step_act_names);
    HSPS::write_labeled_digraph<HSPS::name_vec>
      (std::cerr, p1, step_act_names, false, "Deordered Plan");
#endif

    const HSPS::Schedule::step_vec& steps = plans[k]->plan_steps();
    HSPS::bool_vec is_ea(false, steps.size());
    HSPS::STN tcn((2 * steps.size()) + 1);
    for (HSPS::index_type i = 0; i < steps.size(); i++) {
      HSPS::index_type ea_index = HSPS::no_such_index;
      HSPS::index_type step_act = steps[i].act;
      if (instance.actions[step_act].src) {
	HSPS::ptr_pair* pp =
	  (HSPS::ptr_pair*)instance.actions[step_act].src;
	HSPS::PDDL_Base::ActionSymbol* src_act = 
	  (HSPS::PDDL_Base::ActionSymbol*)pp->first;
	ea_index = ea.first(src_act);
      }
      if (ea_index != HSPS::no_such_index) {
	assert(ea_index < tp.size());
#ifdef VERY_VERBOSE
	std::cerr << "step " << i << " = "
		  << instance.actions[step_act].name
		  << " fixed at " << tp[ea_index]
		  << std::endl;
#endif
	tcn.set_min(0, (2*i)+1, tp[ea_index]);
	tcn.set_max(0, (2*i)+2, tp[ea_index]);
	tcn.set_min((2*i)+1, (2*i)+2, 0);
	is_ea[i] = true;
      }
      else {
	// step i starts after init
	tcn.set_min(0, (2*i)+1, 0);
	// duration of action
	tcn.set_min((2*i)+1, (2*i)+2, instance.actions[step_act].dur);
      }
      for (HSPS::index_type j = 0; j < steps.size(); j++)
	if (prec.adjacent(j, i)) {
	  tcn.set_min((2*j)+2, (2*i)+1, PDDL21_EPSILON);
	}
    }
    // std::cerr << "TCN before 1st compute_min: " << tcn << std::endl;
    tcn.compute_minimal();
    std::cerr << "plan " << k << ": initially consistent = "
	      << tcn.consistent() << std::endl;

    if (tcn.consistent()) {
      bool ok = true;
      for (HSPS::index_type i = 0; (i < steps.size()) && ok; i++)
	for (HSPS::index_type j = i + 1; (j < steps.size()) && ok; j++)
	  if ((steps[i].at == steps[j].at) &&
	      !instance.lock_compatible(steps[i].act, steps[j].act)) {
#ifdef VERBOSE
	    std::cerr << "separating steps " << i << " = "
		      << instance.actions[steps[i].act].name
		      << " and " << j << " = "
		      << instance.actions[steps[j].act].name
		      << std::endl;
#endif
	    if (tcn.admits_min((2*i)+2, (2*j)+1, PDDL21_EPSILON)) {
	      tcn.set_min((2*i)+2, (2*j)+1, PDDL21_EPSILON);
	    }
	    else {
	      tcn.set_min((2*j)+2, (2*i)+1, PDDL21_EPSILON);
	    }
	    tcn.compute_minimal();
	    ok = tcn.consistent();
	  }
    }
    std::cerr << "plan " << k << ": finally consistent = "
	      << tcn.consistent() << std::endl;

#ifdef VERY_VERBOSE
    for (HSPS::index_type i = 0; i < steps.size(); i++) {
      std::cerr << "step " << i << " = "
		<< instance.actions[steps[i].act].name
		<< ": start in [" << tcn.min_distance(0, (2*i)+1)
		<< "," << tcn.max_distance(0, (2*i)+1)
		<< "], end in [" << tcn.min_distance(0, (2*i)+2)
		<< "," << tcn.max_distance(0, (2*i)+2)
		<< "]" << std::endl;
    }
#endif

    if (tcn.consistent()) {
#ifdef IPC_PLAN_FORMAT
      std::cout << ";; plan " << k << std::endl;
      for (HSPS::index_type i = 0; i < steps.size(); i++)
	if (!is_ea[i]) {
	  std::cout << PRINT_NTYPE(tcn.min_distance(0, (2*i)+1))
		    << " : " << instance.actions[steps[i].act].name
		    << " [" << PRINT_NTYPE(instance.actions[steps[i].act].dur)
		    << "]" << std::endl;
	}
#else
      std::cout << "(:plan  ;; #" << k << std::endl;
      for (HSPS::index_type i = 0; i < steps.size(); i++)
	if (!is_ea[i]) {
	  std::cout << PRINT_NTYPE(tcn.min_distance(0, (2*i)+1))
		    << " : " << instance.actions[steps[i].act].name
		    << std::endl;
	}
      std::cout << " )" << std::endl;
#endif
    }
  }

//   for (HSPS::index_type k = 0; k < plans.size(); k++) {
//     HSPS::SafePOP* pop = new HSPS::SafePOP(instance);
//     pop->construct(plans[k]->plan_steps(), true, true);
//     HSPS::index_type t_init = pop->steps[HSPS::SafePOP::INIT_STEP].t_end;
//     HSPS::index_type t_goal = pop->steps[HSPS::SafePOP::GOAL_STEP].t_start;
//     pop->enforce_min_durations();
//     pop->enforce_max_durations();
//     for (HSPS::index_type i = 0; i < pop->steps.size(); i++) {
//       HSPS::index_type step_act = pop->steps[i].act;
//       if (instance.actions[step_act].src) {
// 	HSPS::ptr_pair* pp =
// 	  (HSPS::ptr_pair*)instance.actions[step_act].src;
// 	HSPS::PDDL_Base::ActionSymbol* src_act = 
// 	  (HSPS::PDDL_Base::ActionSymbol*)pp->first;
// 	HSPS::index_type ea_index = ea.first(src_act);
// 	if (ea_index != HSPS::no_such_index) {
// 	  assert(ea_index < tp.size());
// 	  std::cerr << "step " << i << " = "
// 		    << instance.actions[step_act].name
// 		    << " fixed at " << tp[ea_index]
// 		    << std::endl;
// 	  pop->tcn.set_min(t_init, pop->steps[i].t_start, tp[ea_index]);
// 	  pop->tcn.set_max(t_init, pop->steps[i].t_start, tp[ea_index]);
// 	}
//       }
//     }
//     // pop->find_safe_causal_links();
//     std::cerr << "plan " << k
// 	      << ": consistent = " << pop->tcn.consistent()
// 	      << ", makespan = " << pop->tcn.min_distance(t_init, t_goal)
// 	      << std::endl;
//     /// pop->write(std::cout);
//   }

}
