
#include "factor.h"
#include <stdlib.h>

int main(int argc, char *argv[])
{
  if (argc < 2) {
    std::cerr << argv[0] << " <positive number>" << std::endl;
    exit(0);
  }

  if (argc > 2) {
    HSPS::index_type n = atoi(argv[1]);
    HSPS::index_type k = atoi(argv[2]);
    if (k > n) {
      std::cerr << "BAD USER!" << std::endl;
      exit(1);
    }
    HSPS::index_vec cf;
    HSPS::fcombinations(n, k, cf);
    HSPS::index_type cn = HSPS::number(cf);
    std::cout << "c (in factor form): " << cf << std::endl;
    std::cout << "c (in number form): " << cn << std::endl;
  }

  else {
    HSPS::index_type t = atoi(argv[1]);
    if (t == 0) {
      std::cerr << "BAD USER! target number must be POSITIVE!" << std::endl;
      exit(1);
    }

    HSPS::index_vec f(0, 0);
    HSPS::factors(t, f);
    if (f.length() == 1) {
      std::cout << t << " is prime" << std::endl;
      exit(0);
    }

    HSPS::index_type d = 1;
    bool done = false;

    while (!done) {
      if (t > d) {
	HSPS::factors(t - d, f);
	if (f.length() == 1) {
	  std::cerr << t - d << " is prime" << std::endl;
	  done = true;
	}
      }
      if (!done) {
	HSPS::factors(t + d, f);
	if (f.length() == 1) {
	  std::cerr << t + d << " is prime" << std::endl;
	  done = true;
	}
      }
      d += 1;
    }
  }

  return 0;
}
