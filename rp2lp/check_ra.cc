
// #define CHECK_RATIONAL_ARITHMETIC
// #define TRACE_PRINT_RATIONAL_ARITHMETIC

// #define VECTOR_ASSIGNMENT_BUG_TEST
// #define ROUNDING_TEST
// #define PRINT_TEST
#define SUM_TEST
// #define LTGT_TEST

#include "rational.h"
#ifdef SUM_TEST
#include "stats.h"
#endif
#ifdef LTGT_TEST
#include "rng.h"
#endif
#ifdef VECTOR_ASSIGNMENT_BUG_TEST
#include "vector.h"
#endif

int main(int argc, char *argv[])
{
  // HSPS::rational r(377881733, 707099250);
  // std::cerr << r << std::endl;
  // std::cerr << HSPS::rational::printdecimal(r) << std::endl;
  // std::cerr << "done" << std::endl;

#ifdef ROUNDING_TEST
  for (long i = 1; i < 20; i++) {
    HSPS::rational r1(i, i + 1);
    std::cerr << "r = << " << r1 << ", " << -r1 << std::endl;
    std::cerr << " floor to 1 = " << r1.floor_to(1)
	      << ", " << (-r1).floor_to(1) << std::endl;
    std::cerr << " floor to 2 = " << r1.floor_to(2)
	      << ", " << (-r1).floor_to(2) << std::endl;
    std::cerr << " floor to 10 = " << r1.floor_to(10)
	      << ", " << (-r1).floor_to(10) << std::endl;
    std::cerr << " ceil to 1 = " << r1.ceil_to(1)
	      << ", " << (-r1).ceil_to(1) << std::endl;
    std::cerr << " ceil to 2 = " << r1.ceil_to(2)
	      << ", " << (-r1).ceil_to(2) << std::endl;
    std::cerr << " ceil to 10 = " << r1.ceil_to(10)
	      << ", " << (-r1).ceil_to(10) << std::endl;
    HSPS::rational r2(i + 1, i);
    std::cerr << "r = << " << r2 << ", " << -r2 << std::endl;
    std::cerr << " floor to 1 = " << r2.floor_to(1)
	      << ", " << (-r2).floor_to(1) << std::endl;
    std::cerr << " floor to 2 = " << r2.floor_to(2)
	      << ", " << (-r2).floor_to(2) << std::endl;
    std::cerr << " floor to 10 = " << r2.floor_to(10)
	      << ", " << (-r2).floor_to(10) << std::endl;
    std::cerr << " ceil to 1 = " << r2.ceil_to(1)
	      << ", " << (-r2).ceil_to(1) << std::endl;
    std::cerr << " ceil to 2 = " << r2.ceil_to(2)
	      << ", " << (-r2).ceil_to(2) << std::endl;
    std::cerr << " ceil to 10 = " << r2.ceil_to(10)
	      << ", " << (-r2).ceil_to(10) << std::endl;
  }
#endif

#ifdef VECTOR_ASSIGNMENT_BUG_TEST
  cost_vec v0(0, 0);
  cost_vec v1(0, 0);

  v0.assign_value(rational(7,11), 5);
  // std::cerr << "v0 = " << v0 << std::endl;
  v1.set_length(8);
  // v1.assign_value(0, 8);
  for (index_type k = 0; k < 7; k++) {
    // std::cerr << "v0[" << k << "] = " << v0.read_only(k) << std::endl;
    //     rational d = v0[k];
    //     v1[k] = d;
    v1[k] = v0[k];
    // std::cerr << "v1[" << k << "] = " << v1[k] << std::endl;
  }

  std::cerr << "v1 = " << v1 << std::endl;
#endif

#ifdef PRINT_TEST
  for (int i = 1; i < 100; i++) {
    std::cerr << HSPS::rational(i, i+1) << " = "
	      << HSPS::rational::printdecimal(HSPS::rational(i, i + 1))
	      << std::endl;
  }
  int d = 1;
  for (int i = 1; i < 10; i++) {
    std::cerr << HSPS::rational(1, d) << " = "
	      << HSPS::rational::printdecimal(HSPS::rational(1, d))
	      << std::endl;
    d = d * 10;
  }
#endif

#ifdef SUM_TEST
  if (argc < 2) {
    std::cerr << argv[0] << " <n> [<N>]" << std::endl;
    exit(0);
  }

  long n = atoi(argv[1]);
  long N = 1;
  if (argc > 2) N = atoi(argv[2]);

  HSPS::Stopwatch t;
  t.enable_interrupt();

  t.start();
  HSPS::rational s1;
  HSPS::rational s2;
  for (long k = 0; k < N; k++) {
    s1 = 0;
    s2 = 0;
    for (long i = 1; i <= n; i++) {
      HSPS::rational r(1, i);
      bool t1 = (s1 < r);
      bool t2 = (r < s1);
      assert(!t1 || !t2);
      assert((s1 == 0) || t2);
      s1 += r;
      s2 += i;
    }
  }
  std::cerr << "sum i=1.." << n << " 1/i = " << s1 << std::endl;
  std::cerr << "sum i=1.." << n << " i = " << s2 << std::endl;

  //std::cerr << "before stop: time = " << t.time()
  //	    << ", total time = " << t.total_time() << std::endl;
  t.stop();
  //std::cerr << "after stop: time = " << t.time()
  //	    << ", total time = " << t.total_time() << std::endl;
  std::cerr << 2 * n * N << " additions in " << t.time()
	    << " seconds (" << (2 * n * N)/t.time() << " add/sec)"
	    << std::endl;
#endif

#ifdef LTGT_TEST
  HSPS::LC_RNG rng;
  rng.seed_with_time();

  while (true) {
    HSPS::rational r1 =
      HSPS::rational(rng.random() % LONG_MAX,
		     rng.random() % LONG_MAX).reduce();
    HSPS::rational r2 =
      HSPS::rational(rng.random() % LONG_MAX,
		     rng.random() % LONG_MAX).reduce();
    std::cerr << "r1 = " << r1 << ", r2 = " << r2;
    bool lt12 = (r1 < r2);
    std::cerr << ", r1 < r2? " << lt12;
    bool lt21 = (r2 < r1);
    std::cerr << ", r2 < r1? " << lt21;
    assert(!lt12 || !lt21);
    //HSPS::rational d = r1 - r2;
    //std::cerr << "d = " << d << std::endl;
    //if (d.sign() < 0)
    //  assert(lt12);
    //else if (d.sign() > 0)
    //  assert(lt21);
    std::cerr << std::endl;
  }

#endif

  return 0;
}
