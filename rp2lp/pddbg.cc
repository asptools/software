
#include "exec.h"
#include "preprocess.h"
#include "plans.h"
#include "parser.h"

#include <readline/readline.h>
#include <readline/history.h>

#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

#include <iostream>
#include <fstream>

void print_state
(HSPS::Instance& instance,
 std::vector<std::string> atom_names,
 std::vector<std::string> resource_names,
 std::vector<std::string> action_names,
 HSPS::index_type depth,
 HSPS::ExecState& state,
 const std::string& args,
 bool print_true_atoms_only)
{
  HSPS::bool_vec atoms_to_print(true, instance.n_atoms());
  HSPS::bool_vec resources_to_print(true, instance.n_resources());
  const HSPS::ExecState::exec_act_vec& on_going_acts = state.current_actions();
  HSPS::bool_vec actions_to_print(true, on_going_acts.size());
  if (args.length() > 0) {
    boost::regex sel(args, boost::regex::emacs);
    boost::smatch m;
    for (HSPS::index_type k = 0; k < instance.n_atoms(); k++) {
      if (!boost::regex_search(atom_names[k], m, sel, boost::match_any))
	atoms_to_print[k] = false;
    }
    for (HSPS::index_type k = 0; k < instance.n_resources(); k++) {
      if (!boost::regex_search(resource_names[k], m, sel, boost::match_any))
	resources_to_print[k] = false;
    }
    for (HSPS::index_type k = 0; k < on_going_acts.size(); k++) {
      if (!boost::regex_search(action_names[on_going_acts[k].act], m, sel,
			      boost::match_any))
	actions_to_print[k] = false;
    }
  }
  std::cout << "[" << depth << "] current state:" << std::endl;
  if (atoms_to_print.count(false) == 0) {
    std::cout << "atoms: ";
    instance.write_atom_set(std::cout, state.current_atoms());
    std::cout << std::endl;
  }
  else if (atoms_to_print.count(true) > 0) {
    for (HSPS::index_type k = 0; k < instance.n_atoms(); k++)
      if (atoms_to_print[k]) {
	if (print_true_atoms_only) {
	  if (state.current_atoms()[k]) {
	    std::cout << instance.atoms[k].name << std::endl;
	  }
	}
	else {
	  std::cout << instance.atoms[k].name << " is "
		    << (state.current_atoms()[k] ? "true" : "false")
		    << std::endl;
	}
      }
  }
  if (resources_to_print.count(true) > 0) {
    std::cout << "resource levels:" << std::endl;
    for (HSPS::index_type k = 0; k < instance.n_resources(); k++)
      if (resources_to_print[k]) {
	NTYPE avail;
	NTYPE inuse;
	state.current_resource_levels(k, avail, inuse);
	std::cout << instance.resources[k].name << ": "
		  << PRINT_NTYPE(avail) << " total, "
		  << PRINT_NTYPE(avail - inuse) << " free"
		  << std::endl;
      }
  }
  if (actions_to_print.count(true) > 0) {
    std::cout << "on-going actions:" << std::endl;
    for (HSPS::index_type k = 0; k < on_going_acts.size(); k++)
      if (actions_to_print[k]) {
	HSPS::index_type i = on_going_acts[k].act;
	std::cout << instance.actions[i].name << ", with "
		  << PRINT_NTYPE(on_going_acts[k].rem) << " remaining"
		  << std::endl;
      }
  }
}

void find_applicable
(HSPS::Instance& instance,
 HSPS::ExecState& state,
 HSPS::bool_vec& app)
{
  app.assign_value(false, instance.n_actions());
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
    app[k] = state.applicable(k, 0, HSPS::no_such_index);
}

void find_auto_applicable
(HSPS::Instance& instance,
 std::vector<std::string> action_names,
 std::string pattern,
 HSPS::ExecState& state,
 HSPS::bool_vec& app)
{
  //std::cerr << "find_auto_applicable:" << std::endl;
  app.assign_value(false, instance.n_actions());
  if (pattern.length() > 0) {
    boost::regex sel(pattern, boost::regex::emacs);
    boost::smatch m;
    for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
      if (boost::regex_search(action_names[k], m, sel, boost::match_any)) {
	app[k] = state.applicable(k, 0, HSPS::no_such_index);
	//std::cerr << "match: " << action_names[k] << " / " << app[k]
	//	  << std::endl;
      }
    }
  }
  else {
    for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
      app[k] = state.applicable(k, 0, HSPS::no_such_index);
    }
    if (app.count(true) > 1)
      app.assign_value(false, instance.n_actions());
  }
  //std::cerr << "#auto applicable = " << app.count(true) << std::endl;
}

void print_applicable
(HSPS::Instance& instance,
 std::vector<std::string> action_names,
 HSPS::index_type depth,
 HSPS::ExecState& state,
 const std::string& args)
{
  std::cout << "[" << depth << "] applicable actions:" << std::endl;
  HSPS::bool_vec app;
  find_applicable(instance, state, app);
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
    if (app[k]) {
      bool print_app = true;
      if (args.length() > 0) {
	boost::regex sel(args, boost::regex::emacs);
	boost::smatch m;
	if (!boost::regex_search(action_names[k], m, sel, boost::match_any))
	  print_app = false;
      }
      if (print_app) {
	std::cout << " " << k << ". " << instance.actions[k].name
		  << "[c=" << PRINT_NTYPE(instance.actions[k].cost)
		  << ",d=" << PRINT_NTYPE(instance.actions[k].dur) << "]"
		  << std::endl;
      }
    }
  }
  std::cout << app.count(true) << " total" << std::endl;
}

void print_action_details
(HSPS::Instance& instance,
 std::vector<std::string> action_names,
 HSPS::index_type depth,
 HSPS::ExecState& state,
 const std::string& args)
{
  boost::regex sel(args, boost::regex::emacs);
  boost::smatch m;
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++) {
    if (boost::regex_search(action_names[k], m, sel, boost::match_any)) {
      std::cout << k << ". " << instance.actions[k].name
		<< "[c=" << PRINT_NTYPE(instance.actions[k].cost)
		<< ",d=" << PRINT_NTYPE(instance.actions[k].dur) << "]"
		<< std::endl;
      HSPS::ExecErrorSet* errors = new HSPS::ExecErrorSet();
      bool is_app =
	state.applicable(k, errors, HSPS::no_such_index);
      if (is_app) {
	std::cout << " applicable" << std::endl;
      }
      else {
	for (HSPS::index_type i = 0; i < errors->size(); i++) {
	  std::cout << " ";
	  (*errors)[i]->write(std::cout);
	  std::cout << std::endl;
	}
      }
      delete errors;
    }
  }
}

int main(int argc, char *argv[]) {
  HSPS::StringTable symbols(50, HSPS::lowercase_map);
  HSPS::Statistics stats;
  HSPS::Parser* reader = new HSPS::Parser(symbols);
  HSPS::PDDL_Base::make_types_from_static_predicates = false;

  bool opt_preprocess = true;
  bool opt_remove_irrelevant = false;
  bool opt_apply_plan = false;
  int verbose_level = 0;

  for (int k = 1; k < argc; k++) {
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      verbose_level = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-dba-semantics") == 0) {
      HSPS::PDDL_Base::del_before_add_semantics = true;
    }
    else if (strcmp(argv[k],"-create-all") == 0) {
      HSPS::PDDL_Base::create_all_atoms = true;
      HSPS::PDDL_Base::create_all_actions = true;
      HSPS::PDDL_Base::check_precondition_consistency = false;
      opt_preprocess = false;
    }
    else if (strcmp(argv[k],"-no-prep") == 0) {
      opt_preprocess = false;
    }
    else if (strcmp(argv[k],"-remove") == 0) {
      opt_remove_irrelevant = true;
    }
    else if (strcmp(argv[k],"-apply-plan") == 0) {
      opt_apply_plan = true;
    }
    else if (*argv[k] != '-') {
      reader->read(argv[k], false);
    }
  }

  HSPS::Heuristic::default_trace_level = verbose_level;
  HSPS::Instance::default_trace_level = verbose_level;
  HSPS::Preprocessor::default_trace_level = verbose_level;
  if (verbose_level <= 0) HSPS::PDDL_Base::warning_level = 0;
  if (verbose_level > 1) HSPS::PDDL_Base::write_info = true;

  stats.start();
  std::cout << "instantiating..." << std::endl;
  HSPS::Instance instance;
  reader->instantiate(instance);

  HSPS::Preprocessor prep(instance, stats);
  if (opt_preprocess) {
    std::cout << "preprocessing..." << std::endl;
    prep.preprocess(false);
    if (opt_remove_irrelevant) {
      prep.compute_irrelevant_atoms();
      prep.remove_irrelevant_atoms();
      if (!instance.cross_referenced()) {
	std::cout << "re-cross referencing..." << std::endl;
	instance.cross_reference();
      }
    }
  }
  else {
    instance.cross_reference();
  }

  stats.stop();
  std::cout << "instance " << instance.name << " built in "
	    << stats.total_time() << " seconds" << std::endl;
  std::cout << instance.n_atoms() << " atoms, "
	    << instance.n_resources() << " resources ("
	    << instance.n_reusable_resources() << " reusable, "
	    << instance.n_consumable_resources() << " consumable), "
	    << instance.n_actions() << " actions, "
	    << instance.n_invariants() << " invariants"
	    << std::endl;

  std::vector<std::string> atom_names(instance.n_atoms(), "");
  for (HSPS::index_type k = 0; k < instance.n_atoms(); k++)
    atom_names[k] = instance.atoms[k].name->to_string();
  std::vector<std::string> resource_names(instance.n_resources(), "");
  for (HSPS::index_type k = 0; k < instance.n_resources(); k++)
    resource_names[k] = instance.resources[k].name->to_string();
  std::vector<std::string> action_names(instance.n_actions(), "");
  for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
    action_names[k] = instance.actions[k].name->to_string();

  // void using_history();  // is this necessary?

  bool sequential = true;
  bool autoapply = false;
  std::string autoapply_pattern = "";
  bool watching_state = false;
  bool watching_applic = true;
  std::string watched_state = "";
  std::string watched_applic = "";
  bool triggered_watch = true;
  std::vector<HSPS::ExecState*> state_stack;
  std::vector<HSPS::Schedule*> plan_stack;
  HSPS::ExecState* current_state = NULL;
  HSPS::Schedule* current_plan = NULL;
  if (opt_apply_plan) {
    if (reader->n_plans() < 1) {
      std::cout << "error: no initial plan provided" << std::endl;
      exit(0);
    }
    if (reader->n_plans() > 1) {
      std::cout << "error: more than one initial plan provided" << std::endl;
      exit(0);
    }
    HSPS::Schedule* initial_plan = new HSPS::Schedule(instance);
    bool ok = reader->export_plan(0, instance, prep.action_map, *initial_plan);
    if (!ok) {
      std::cout << "error: failed to export initial plan" << std::endl;
      exit(0);
    }
    HSPS::ExecTrace* trace = new HSPS::ExecTrace(instance);
    HSPS::ExecErrorSet* errors = new HSPS::ExecErrorSet();
    ok = initial_plan->simulate(trace, errors, true);
    if (!errors->executable()) {
      std::cout << "error: initial plan failed to execute" << std::endl;
      errors->write(std::cout);
      std::cout << std::endl;
      exit(0);
    }
    delete errors;
    current_state = (HSPS::ExecState*)trace->final_state()->copy();
    delete trace;
    current_plan = new HSPS::Schedule(instance);
    initial_plan->output(*current_plan, false);
    delete initial_plan;
  }

  // not applying initial plan, or failed to apply it
  if (!opt_apply_plan) {
    current_state = new HSPS::ExecState(instance, instance.init_atoms);
    current_plan = new HSPS::Schedule(instance);
  }

  boost::regex split_cmd_and_args
    ("[[:space:]]*\\([[:alnum:]]+\\)\\([[:space:]]+.*\\)?",
     boost::regex::emacs);

  boost::regex index_type_exp("[0-9]+", boost::regex::emacs);

  bool done = false;
  while (!done) {
    if ((watching_state || watching_applic) && triggered_watch) {
      if (watching_state)
	print_state(instance, atom_names, resource_names, action_names,
		    state_stack.size(), *current_state, watched_state, true);
      if (watching_applic)
	print_applicable(instance, action_names,
			 state_stack.size(), *current_state, watched_applic);
    }
    triggered_watch = false;
    char* _cmd = readline("> ");
    if (_cmd && *_cmd) {
      add_history(_cmd);
      boost::cmatch m;
      if (boost::regex_match(_cmd, m, split_cmd_and_args)) {
	assert(m.size() >= 2);
	std::string cmd = m.str(1);
	boost::to_lower(cmd);
	std::string args = (m.size() >= 3 ? m.str(2) : "");
	boost::trim(args);
	boost::smatch m2;
	bool cmd_is_apply = false;
	HSPS::index_type act = HSPS::no_such_index;

	if (cmd == std::string("quit")) {
	  done = true;
	}

	else if (cmd == std::string("state")) {
	  print_state(instance, atom_names, resource_names, action_names,
		      state_stack.size(), *current_state, args, false);
	}

	else if (cmd == std::string("true")) {
	  print_state(instance, atom_names, resource_names, action_names,
		      state_stack.size(), *current_state, args, true);
	}

	else if (cmd == std::string("applicable")) {
	  print_applicable(instance, action_names,
			   state_stack.size(), *current_state, args);
	}

	else if (cmd == std::string("action")) {
	  if (args.length() == 0) {
	    std::cout << "error:  \"action\" without argument" << std::endl;
	  }
	  else {
	    print_action_details(instance, action_names,
				 state_stack.size(), *current_state, args);
	  }
	}

	else if (cmd == std::string("back")) {
	  bool valid_arg = true;
	  HSPS::index_type n = 1;
	  boost::smatch m2;
	  if (boost::regex_match(args, m2, index_type_exp)) {
	    n = atoi(args.c_str());
	    if (n > state_stack.size()) {
	      std::cout << "error: " << n << " is not a valid index"
			<< std::endl;
	      valid_arg = false;
	    }
	  }
	  if (valid_arg) {
	    HSPS::index_type p = (state_stack.size() - n);
	    for (HSPS::index_type k = p + 1; k < state_stack.size(); k++) {
	      delete state_stack[k];
	      delete plan_stack[k];
	    }
	    delete current_state;
	    delete current_plan;
	    current_state = state_stack[p];
	    current_plan = plan_stack[p];
	    state_stack.resize(p);
	    plan_stack.resize(p);
	    triggered_watch = true;
	  }
	}

	else if (cmd == std::string("backto")) {
	  bool valid_arg = true;
	  HSPS::index_type p = HSPS::no_such_index;
	  if (args.length() == 0) {
	    std::cout << "error: \"backto\" without argument!" << std::endl;
	    valid_arg = false;
	  }
	  else if (boost::regex_match(args, m2, index_type_exp)) {
	    p = atoi(args.c_str());
	    if (p >= state_stack.size()) {
	      std::cout << "error: " << p << " is not a valid index"
			<< std::endl;
	      valid_arg = false;
	    }
	  }
	  if (valid_arg) {
	    for (HSPS::index_type k = p + 1; k < state_stack.size(); k++) {
	      delete state_stack[k];
	      delete plan_stack[k];
	    }
	    delete current_state;
	    delete current_plan;
	    current_state = state_stack[p];
	    current_plan = plan_stack[p];
	    state_stack.resize(p);
	    plan_stack.resize(p);
	    triggered_watch = true;
	  }
	}

	else if (cmd == std::string("dumpplan")) {
	  current_plan->write_steps(std::cout);
	}

	else if (cmd == std::string("printplan")) {
	  HSPS::PrintIPC planprinter(instance, std::cout);
	  current_plan->output(planprinter);
	}

	else if (cmd == std::string("saveplan")) {
	  if (args.length() > 0) {
	    std::cout << "saving current plan to " << args << std::endl;
	    std::ofstream save_to(args.c_str());
	    HSPS::PrintIPC planprinter(instance, save_to);
	    current_plan->output(planprinter);
	    save_to.close();
	  }
	  else {
	    std::cerr << "error: \"saveplan\" without file name" << std::endl;
	  }
	}

	else if (cmd == std::string("watch")) {
	  if (args.length() > 0) {
	    boost::smatch m2;
	    if (boost::regex_match(args, m2, split_cmd_and_args)) {
	      assert(m2.size() >= 2);
	      if (m2.str(1) == "state") {
		watching_state = true;
		if (m2.size() >= 3) {
		  watched_state = m2.str(2);
		  boost::trim(watched_state);
		}
		else {
		  watched_state = "";
		}
		watching_applic = false;
	      }
	      else if (m2.str(1) == "applicable") {
		watching_state = false;
		watching_applic = true;
		if (m2.size() >= 3) {
		  watched_applic = m2.str(2);
		  boost::trim(watched_applic);
		}
		else {
		  watched_applic = "";
		}
	      }
	    }
	    else {
	      watching_state = true;
	      watched_state = args;
	      watching_applic = true;
	      watched_applic = args;
	    }
	  }
	  else {
	    watching_state = true;
	    watched_state = "";
	    watching_applic = true;
	    watched_applic = "";
	  }
	  triggered_watch = true;
	}

	else if (cmd == std::string("nowatch")) {
	  watching_state = false;
	  watching_applic = false;
	}

	else if (cmd == std::string("autoapply")) {
	  autoapply = true;
	  autoapply_pattern = args;
	}

	else if (cmd == std::string("noautoapply")) {
	  autoapply = false;
	}

	else if (boost::regex_match(cmd, m2, index_type_exp)) {
	  act = atoi(cmd.c_str());
	  if (act >= instance.n_actions()) {
	    std::cout << "error: " << act << " is not a valid action index"
		      << std::endl;
	  }
	  else {
	    cmd_is_apply = true;
	  }
	}

	else {
	  for (HSPS::index_type k = 0; k < instance.n_actions(); k++)
	    if (action_names[k] == cmd)
	      act = k;
	  if (act == HSPS::no_such_index) {
	    std::cout << "error: " << cmd << " is not a valid action name"
		      << std::endl;
	  }
	  else {
	    cmd_is_apply = true;
	  }
	}

	if (cmd_is_apply) {
	  assert(act < instance.n_actions());
	  HSPS::ExecErrorSet* errors = new HSPS::ExecErrorSet();
	  bool is_app =
	    current_state->applicable(act, errors, HSPS::no_such_index);
	  if (!is_app) {
	    std::cout << "error: " << instance.actions[act].name
		      << " is not applicable in current state"
		      << std::endl;
	    for (HSPS::index_type i = 0; i < errors->size(); i++) {
	      (*errors)[i]->write(std::cout);
	      std::cout << std::endl;
	    }
	  }
	  delete errors;
	  if (is_app) {
	    std::cout << "[" << state_stack.size() << "] applying "
		      << instance.actions[act].name << "..." << std::endl;
	    plan_stack.push_back(current_plan);
	    HSPS::Schedule* new_plan =
	      new HSPS::Schedule(*current_plan);
	    state_stack.push_back(current_state);
	    HSPS::ExecState* new_state =
	      (HSPS::ExecState*)current_state->copy();
	    new_plan->insert(act);
	    new_state->apply(act, NULL, HSPS::no_such_index);
	    if (sequential) {
	      NTYPE d = new_state->max_delta();
	      new_plan->advance(d);
	      new_state->advance(d, NULL);
	    }
	    if (autoapply) {
	      HSPS::bool_vec app;
	      find_auto_applicable(instance, action_names,
				   autoapply_pattern, *new_state, app);
	      while (app.count(true) > 0) {
		act = app.first(true);
		assert(act < instance.n_actions());
		std::cout << "[" << state_stack.size() << "] auto-applying "
			  << act << ". " << instance.actions[act].name
			  << "..." << std::endl;
		new_plan->insert(act);
		new_state->apply(act, NULL, HSPS::no_such_index);
		if (sequential) {
		  NTYPE d = new_state->max_delta();
		  new_plan->advance(d);
		  new_state->advance(d, NULL);
		}
		find_auto_applicable(instance, action_names,
				     autoapply_pattern, *new_state, app);
	      }
	    }
	    current_plan = new_plan;
	    current_state = new_state;
	    triggered_watch = true;
	  }
	}
      }
      else {
	std::cout << "can't parse that" << std::endl;
      }
      free(_cmd);
    }
  }

  return 0;
}
