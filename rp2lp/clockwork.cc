
#include "sas.h"
#include "preprocess.h"
#include "cost_table.h"
#include "stats.h"
#include "base.h"
#include <sstream>
#include <fstream>
#include <math.h>

HSPS::index_type choose_prevail
(HSPS::index_type layer,
 HSPS::index_type pos,
 HSPS::index_type  width,
 HSPS::index_type  conc,
 bool asym,
 HSPS::RNG& rng)
{
  HSPS::index_type r = width + conc - 1;
  HSPS::index_type s1 = rng.random_in_range(r);
  HSPS::index_type s = (s1 < width) ? s1 : 0;
  if (!asym) {
    s = (s + pos) % width;
  }
  return (((layer - 1) * width) + s);
}

int main(int argc, char *argv[]) {
  HSPS::index_type depth   = 2;
  HSPS::index_type width   = 2;
  HSPS::index_type n_val   = 10;
  HSPS::index_type conc = 1;
  bool  asym = false;

  unsigned long rnd_seed = 0;
  HSPS::index_type n_instances = 1;

  bool opt_save = false;
  bool opt_single_file = false;
  const char* save_dir = "";

  unsigned int time_limit = 0;
  int verbose_level = 1;

  for (int k = 1; k < argc; k++) {
    if ((strcmp(argv[k],"-v") == 0) && (k < argc - 1)) {
      verbose_level = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-d") == 0) && (k < argc - 1)) {
      depth = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-w") == 0) && (k < argc - 1)) {
      width = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-m") == 0) && (k < argc - 1)) {
      n_val = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-c") == 0) && (k < argc - 1)) {
      conc = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-seq") == 0) {
      asym = true;
    }
    else if (((strcmp(argv[k],"-rnd") == 0) ||
	      (strcmp(argv[k],"-r") == 0)) &&
	     (k < argc - 1)) {
      rnd_seed = atoi(argv[++k]);
    }
    else if ((strcmp(argv[k],"-i") == 0) && (k < argc - 1)) {
      n_instances = atoi(argv[++k]);
    }
    else if (strcmp(argv[k],"-save") == 0) {
      opt_save = true;
    }
    else if (strcmp(argv[k],"-save-1") == 0) {
      opt_save = true;
      opt_single_file = true;
    }
    else if (strcmp(argv[k],"-save-2") == 0) {
      opt_save = true;
      opt_single_file = false;
    }
    else if ((strcmp(argv[k],"-dir") == 0) && (k < argc - 1)) {
      save_dir = argv[++k];
    }
    else if (strcmp(argv[k],"-pq-param") == 0)
      HSPS::Instance::always_write_parameters = true;
    else if (strcmp(argv[k],"-pq-req") == 0)
      HSPS::Instance::always_write_requirements = true;
  }

  HSPS::Instance::write_DKEL = false;
  HSPS::Instance::write_PDDL2 = false;
  HSPS::Instance::write_time = false;
  HSPS::Instance::write_metric = false;
  HSPS::SASInstance::STRIPS_eq_sign = '-';

  HSPS::Instance::always_write_parameters = true;
  HSPS::Instance::always_write_requirements = true;

  HSPS::Statistics stats;

  HSPS::LC_RNG rng;
  if (rnd_seed > 0) rng.seed(rnd_seed);

  std::cout << "rnd " << rng.seed_value() << std::endl;

  std::ostringstream d_name;
  d_name << "CW" << width << "X" << depth << "X" << n_val;
  if (conc == 1) {
    d_name << "U";
  }
  else if (asym) {
    d_name << "S" << conc;
  }
  else {
    d_name << "P" << conc;
  }

  //HSPS::index_type n_var = (width * depth) + 1;
  HSPS::index_type n_var = (width * depth);
  HSPS::index_type n_top_val = 2*width + 1;

  HSPS::index_type ins_no = 0;
  while (ins_no < n_instances) {
    unsigned long instance_id = rng.seed_value();
    std::cerr << "generating instance " << instance_id << "..." << std::endl;

    std::ostringstream i_name;
    i_name << "I" << instance_id;
    HSPS::Name* n = new HSPS::InstanceName(strdup(d_name.str().c_str()),
					   strdup(i_name.str().c_str()));
    HSPS::SASInstance* sas_ins = new HSPS::SASInstance(n);

    // create variables
    for (HSPS::index_type l = 0; l < depth; l++)
      for (HSPS::index_type p = 0; p < width; p++) {
	n = new HSPS::EnumName("var", l, p);
	HSPS::SASInstance::Variable& v = sas_ins->new_variable(n);
	for (HSPS::index_type i = 0; i < n_val; i++)
	  v.domain.append(new HSPS::EnumName("val", i));
      }
    //n = new HSPS::StringName("top");
    //HSPS::SASInstance::Variable& v = sas_ins->new_variable(n);
    //for (HSPS::index_type i = 0; i < n_top_val + 1; i++)
    //  v.domain.append(new HSPS::EnumName("val", i));

    // create clockwork actions
    for (HSPS::index_type l = 0; l < depth; l++)
      for (HSPS::index_type p = 0; p < width; p++) {
	HSPS::SASInstance::Variable& vlp =
	  sas_ins->variables[(l*width)+p];
	for (HSPS::index_type i = 0; i < n_val; i++) {
	  // clockwise (increasing)
	  n = new HSPS::EnumName("cw", l, p, i);
	  HSPS::SASInstance::Action& a1 = sas_ins->new_action(n);
	  a1.pre.assign(vlp.index, i);
	  a1.post.assign(vlp.index, ((i + 1) < n_val ? i + 1 : 0));
	  if (l > 0) {
	    HSPS::index_type pi = choose_prevail(l, p, width, conc, asym, rng);
	    HSPS::index_type pv = rng.random_in_range(n_val);
	    a1.prv.assign(pi, pv);
	  }
	  // anti-clockwise (decreasing)
	  n = new HSPS::EnumName("acw", l, p, i);
	  HSPS::SASInstance::Action& a2 = sas_ins->new_action(n);
	  a2.pre.assign(vlp.index, i);
	  a2.post.assign(vlp.index, (i > 0 ? i - 1 : n_val - 1));
	  if (l > 0) {
	    HSPS::index_type pi = choose_prevail(l, p, width, conc, asym, rng);
	    HSPS::index_type pv = rng.random_in_range(n_val);
	    a2.prv.assign(pi, pv);
	  }
	}
      }

    // top level actions
    //for (HSPS::index_type i = 0; i < n_top_val; i++) {
    //  n = new HSPS::EnumName("step", i);
    //  HSPS::SASInstance::Action& a = sas_ins->new_action(n);
    //  a.pre.assign(n_var - 1, i);
    //  a.post.assign(n_var - 1, i + 1);
    //  if (depth > 0) {
    //	HSPS::index_type pi = choose_prevail(depth, 0, width, conc, asym, rng);
    //	HSPS::index_type pv = rng.random_in_range(n_val);
    //	a.prv.assign(pi, pv);
    //  }
    //}

    // assign initial values & goal values
    for (HSPS::index_type l = 0; l < depth; l++)
      for (HSPS::index_type p = 0; p < width; p++) {
	HSPS::index_type iv = rng.random_in_range(n_val);
	sas_ins->init_state.assign((l * width) + p, iv);
	if ((l + 1) == depth) {
	  HSPS::index_type gv = rng.random_in_range(n_val, iv);
	  sas_ins->goal_state.assign((l * width) + p, gv);
	}
      }
    //sas_ins->init_state.assign(n_var - 1, 0);
    //sas_ins->goal_state.assign(n_var - 1, n_top_val);

    if (verbose_level > 1) {
      sas_ins->write_domain(std::cout);
    }

    if (opt_save) {
      std::cerr << "converting instance to STRIPS..." << std::endl;
      HSPS::Instance* strips_ins = sas_ins->convert_to_STRIPS();
      strips_ins->cross_reference();

      std::string dir_name(save_dir);
      if (strlen(save_dir) > 0) {
	if (save_dir[strlen(save_dir) - 1] != '/') dir_name += "/";
      }

      std::ostringstream base_file_name;
      base_file_name << d_name.str() << "_" << i_name.str();

      if (opt_single_file) {
	std::string file_name =
	  dir_name + base_file_name.str() + std::string(".pddl");
	if (verbose_level > 0)
	  std::cerr << "writing domain/problem "
		    << base_file_name.str() << "..."
		    << std::endl;
	std::ofstream save_file(base_file_name.str().c_str());
	strips_ins->write_domain(save_file);
	strips_ins->write_problem(save_file);
	save_file.close();
      }

      else {
	std::string domain_file_name =
	  dir_name + base_file_name.str() + std::string(".domain.pddl");
	if (verbose_level > 0)
	  std::cerr << "writing domain " << domain_file_name << "..."
		  << std::endl;
	std::ofstream domain_file(domain_file_name.c_str());
	strips_ins->write_domain(domain_file);
	domain_file.close();

	std::string problem_file_name =
	  dir_name + base_file_name.str() + std::string(".pddl");
	if (verbose_level > 0)
	  std::cerr << "writing problem " << problem_file_name << "..."
		    << std::endl;
	std::ofstream problem_file(problem_file_name.c_str());
	strips_ins->write_problem(problem_file);
	problem_file.close();
      }
      delete strips_ins;
    }

    delete sas_ins;
    ins_no += 1;
  }

  std::cout << "next " << rng.seed_value() << std::endl;

  return 0;
}
