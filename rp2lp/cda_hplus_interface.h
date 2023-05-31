#ifndef CDASTAR_HPLUS_INTERFACE
#define CDASTAR_HPLUS_INTERFACE

extern "C" {

  // Read PDDL input (domain and problem), create ILB object.
  // Note that input files will be read in the order they are
  // given; the domain must precede the problem.
  void* init(char* files[], int n_files);

  // Get number of (ground) actions in the domain from the ILB object.
  unsigned int number_of_actions(void* obj);

  // Get action costs.
  // On call: num and div must be allocated arrays of #actions ints.
  // On return: num[i]/div[i] is the rational cost of action i.
  void costs(void* obj, int num[], int div[]);

  // Perform a relaxed reachability test.
  // On call:
  // - set is the input (set of actions to be tested); it does not
  //   need to be sorted;
  // - conflict must be allocated an array of #actions ints.
  // On return, if return value == true:
  // - the set is a relaxed plan; the conflict array is untouched.
  // On return, if return value == false:
  // - the set is not a relaxed plan;
  // - the first n_conflict entries of the conflict array contains
  //   a new disjunctive action landmark (which is sorted).
  bool test(void* obj,
	    unsigned int set[], unsigned int n_set,
	    unsigned int conflict[], unsigned int n_conflict);
}

#endif
