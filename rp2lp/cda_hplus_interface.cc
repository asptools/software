
#include "cda_hplus_interface.h"
#include "preprocess.h"
#include "parser.h"
#include "ilb.h"

extern "C" {

  void* init(char* files[], int n_files)
  {
    HSPS::Stopwatch* stats = new HSPS::Stopwatch;
    stats->start();
    HSPS::StringTable symbols(50, HSPS::lowercase_map);
    HSPS::Parser* reader = new HSPS::Parser(symbols);
    HSPS::PDDL_Base::del_before_add_semantics = true;
    HSPS::PDDL_Base::make_types_from_static_predicates = false;
    for (unsigned int k = 0; k < n_files; k++) {
      std::cerr << "reading " << files[k] << "..." << std::endl;
      reader->read(files[k], false);
    }
    std::cerr << "instantiating..." << std::endl;
    HSPS::Instance* instance = new HSPS::Instance;
    reader->instantiate(*instance);
    if (HSPS::PDDL_Base::del_before_add_semantics) {
      for (HSPS::index_type k = 0; k < instance->n_actions(); k++)
	instance->actions[k].add.subtract(instance->actions[k].pre);
    }
    HSPS::Preprocessor prep(*instance, *stats);
    std::cerr << "preprocessing..." << std::endl;
    prep.preprocess(false);
    prep.compute_irrelevant_atoms();
    prep.remove_irrelevant_atoms();
    if (!instance->cross_referenced()) {
      std::cerr << "re-cross referencing..." << std::endl;
      instance->cross_reference();
    }
    HSPS::StaticMutex* mx = prep.inconsistency();
    HSPS::CostACF* cost = new HSPS::CostACF(*instance);
    stats->stop();
    return new HSPS::ILB(*instance, *cost, mx, *stats);
  }

  unsigned int number_of_actions(void* obj)
  {
    assert(obj != NULL);
    HSPS::ILB* ilb = (HSPS::ILB*)obj;
    return ilb->instance().n_actions();
  }

  void costs(void* obj, int num[], int div[])
  {
    assert(obj != NULL);
    HSPS::ILB* ilb = (HSPS::ILB*)obj;
    for (unsigned int i = 0; i < ilb->instance().n_actions(); i++) {
      num[i] = N_TO_R(ilb->instance().actions[i].cost).numerator();
      div[i] = N_TO_R(ilb->instance().actions[i].cost).divisor();
    }
  }

  bool test(void* obj,
	    unsigned int acts[], unsigned int n_acts,
	    unsigned int conflict[], unsigned int n_conflict)
  {
    assert(obj != NULL);
    HSPS::ILB* ilb = (HSPS::ILB*)obj;
    HSPS::bool_vec act_set_in(false, ilb->instance().n_actions());
    for (unsigned int i = 0; i < n_acts; i++) {
      assert(acts[i] < ilb->instance().n_actions());
      act_set_in[acts[i]] = true;
    }
    bool pass = ilb->is_relaxed_plan(act_set_in);
    if (pass) return true;
    HSPS::index_set as2(act_set_in);
    const HSPS::index_set& new_lm = ilb->make_new_landmark(as2);
    for (unsigned int i = 0; i < new_lm.size(); i++)
      conflict[i] = new_lm[i];
    n_conflict = new_lm.size();
    return false;
  }

}
