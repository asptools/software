#ifndef NUMERIC_TYPE_H
#define NUMERIC_TYPE_H

#include "config.h"
#include "rational.h"
#include "index_type.h"
#include "rng.h"

#ifdef NTYPE_RATIONAL
#include <iostream>
#include <iomanip>

#define NTYPE  HSPS::rational
#define NNTYPE HSPS::rational::XR

#define N_TO_NN(x) x
#define NN_TO_N(x) HSPS::rational(x)

#define A_TO_N(x)   HSPS::rational::ator(x)
#define I_TO_N(x)   HSPS::rational(x)
#define R_TO_N(x,y) (HSPS::rational(x,y).reduce())
#define R1_TO_N(x)  (x)
#define D_TO_N(x)   HSPS::rational::dtor(x)
#define N_TO_D(x)   ((x).decimal())
#define N_TO_R(x)   (x)

#define INFINITE(x) (x).infinite()
#define FINITE(x)   (x).finite()
const HSPS::rational POS_INF(1,0);
const HSPS::rational NEG_INF(-1,0);
const HSPS::rational ZERO(0,1);
#define INTEGRAL(x) (x).integral()
#define FLOOR(x) (x).floor()
#define FLOOR_TO_INT(x) ((x).floor().numerator())
#define CEIL(x) (x).ceil()
#define CEIL_TO_INT(x) ((x).ceil().numerator())
#define FRAC(x) (x).frac()
#define FLOOR_TO(x,p) rational::floor_to(x, p)
#define CEIL_TO(x,p) rational::floor_to(x, p)
#define IS_ZERO(x) (x).zero()
#define ABS(x)     ((x) * (x).sign())
#define SAFE_GT(x,y) ((x) > (y))

#define MIN(x,y)    HSPS::rational::min(x,y)
#define MAX(x,y)    HSPS::rational::max(x,y)

#define HASH1(x)     ((index_type)((x).numerator() - (x).divisor()))
// #define HASH2(x,y)   HASH1(((x) + 1) * ((y) + 1))
#define HASH2(x,y)   (HASH1(x) + (y * 13))

#ifdef USE_HSPS_NAMESPACE
#define PRINT_NTYPE(x) ::HSPS::rational::printdecimal(x)
#else
#define PRINT_NTYPE(x) rational::printdecimal(x)
#endif

#endif /* NTYPE_RATIONAL */

#ifdef NTYPE_FLOAT

#include <limits>
#include <cmath>

#define NTYPE  double
#define NNTYPE double

#define N_TO_NN(x) x
#define NN_TO_N(x) x
#define A_TO_NN(x) atof(x)

#define A_TO_N(x)   atof(x)
#define I_TO_N(x)   ((double)(x))
#define R_TO_N(x,y) (((double)(x))/((double)(y)))
#define R1_TO_N(x)  (((double)(x.numerator()))/((double)(x.divisor())))
#define D_TO_N(x)   (x)
#define N_TO_D(x)   (x)
#define N_TO_R(x)   HSPS::rational::dtor(x)

#define INFINITE(x) (!std::isfinite(x))
#define FINITE(x)   std::isfinite(x)

/* definitions of isinf/finite from cygwin ieeefp.h (should not be
   needed anymore since isfinite is in c++ std):
#define INFINITE(x) (((*(long *)&(x) & 0x7f800000L)==0x7f800000L) && \
		    ((*(long *)&(x) & 0x007fffffL)==0000000000L))
#define FINITE(x)   (((*(long *)&(x) & 0x7f800000L)!=0x7f800000L))
*/

//const double POS_INF = ((double)(1.0/0.0));
//const double NEG_INF = ((double)(-1.0/0.0));
const double POS_INF = std::numeric_limits<double>::infinity();
const double NEG_INF = -POS_INF;
const double ZERO = ((double)0);

#define INTEGRAL(x) (std::floor(x) == (x))
#define FLOOR(x) std::floor(x)
#define FLOOR_TO_INT(x) ((long)std::floor(x))
#define CEIL(x) std::ceil(x)
#define CEIL_TO_INT(x) ((long)std::ceil(x))
#define FRAC(x)  (x - std::floor(x))
#define FLOOR_TO(x,p) (std::floor((x) * (p)) / (p))
#define CEIL_TO(x,p) (std::ceil((x) * (p)) / (p))
#define ABS(x) (std::abs(x))

inline bool SAFE_EQ(double x, double y) {
  if (x == y) {
    return true;
  }
  else {
    double diff = std::abs(x - y);
    double abs_x = std::abs(x);
    double abs_y = std::abs(y);
    if ((abs_x < std::numeric_limits<double>::epsilon()) ||
	(abs_y < std::numeric_limits<double>::epsilon())) {
      return (diff < (2 * std::numeric_limits<double>::epsilon()));
    }
    else {
      double m = (abs_x > abs_y ? abs_x : abs_y);
      return ((diff / m) < (2 * std::numeric_limits<double>::epsilon()));
    }
  }
}

#define IS_ZERO(x) SAFE_EQ(x, ZERO)

inline bool SAFE_GT(double x, double y) {
  // if x is +INF then x > y unless y is also +INF:
  if (!std::isfinite(x)) {
    return std::isfinite(y);
  }
  // else if y is +INF, x > y is FALSE:
  else if (!std::isfinite(y)) {
    return false;
  }
  // else, if x - y is positive and SAFE_EQ(x, y) fails, x > y is TRUE:
  else {
    double diff = (x - y);
    return ((diff > ZERO) && !SAFE_EQ(x, y));
  }
}

inline double MIN(double x, double y) {
  if (x < y) return x; else return y;
}

inline double MAX(double x, double y) {
  if (x > y) return x; else return y;
}

#define HASH1(x)     ((index_type)(x))
#define HASH2(x,y)   HASH1((x) * (y))

#define PRINT_NTYPE(x) x

#endif /* NTYPE_FLOAT */

BEGIN_HSPS_NAMESPACE

inline NTYPE random_numeric
(NTYPE min, NTYPE max, unsigned long prec, RNG& rng)
{
  NTYPE d = (max - min);
  NTYPE s = (d / prec);
  unsigned long r = rng.random_in_range(prec + 1);
  return ((r*s) + min);
}

class amt_vec : public auto_expanding_vector<NTYPE> {
 public:
  amt_vec()
    : auto_expanding_vector<NTYPE>() { };
  amt_vec(const NTYPE& v, index_type l)
    : auto_expanding_vector<NTYPE>(v, l) { };
  amt_vec(const amt_vec& vec)
    : auto_expanding_vector<NTYPE>(vec) { };

  int compare(const amt_vec& vec, index_type n);
  //  1 if this dominates vec
  // -1 if vec dominates this
  //  0 if incomparable
  int dcompare(const amt_vec& vec, index_type n);
  index_type hash(index_type n);
  void write(std::ostream& s, index_type n);
};

inline int amt_vec::compare(const amt_vec& vec, index_type n)
{
  for (index_type k = 0; k < n; k++) {
    if ((*this)[k] < vec[k]) return -1;
    else if ((*this)[k] > vec[k]) return 1;
  }
  return 0;
}

inline int amt_vec::dcompare(const amt_vec& vec, index_type n)
{
  bool this_less_than_vec = false;
  bool vec_less_than_this = false;
  for (index_type k = 0; k < n; k++) {
    if ((*this)[k] < vec[k]) this_less_than_vec = true;
    else if ((*this)[k] > vec[k]) vec_less_than_this = true;
  }
  if (this_less_than_vec && !vec_less_than_this) return -1;
  else if (!this_less_than_vec && vec_less_than_this) return 1;
  else return 0;
}

inline index_type amt_vec::hash(index_type n)
{
  if (n == 0) return 0;
  if (n == 1) return HASH1((*this)[0]);
  index_type h = 0;
  // for (index_type k = 0; k < n - 1; k++) {
  //   h += HASH2((*this)[k], (*this)[k + 1]);
  // }
  for (index_type k = 0; k < n; k++) {
    h = HASH2((*this)[k], h);
  }
  return h;
}

inline void amt_vec::write(std::ostream& s, index_type n)
{
  s << '[';
  for (index_type k = 0; k < n; k++) {
    if (k > 0) s << ',';
    s << PRINT_NTYPE((*this)[k]);
  }
  s << ']';
}

typedef lvector<NTYPE> cost_vec;
typedef svector<NTYPE> cost_set;

class cost_vec_util : public cost_vec
{
 public:

  class decreasing_cost_order : public cost_vec::order {
  public:
    virtual bool operator()
      (const NTYPE& v0, const NTYPE& v1) const
      { return (v0 > v1); };
  };

  class increasing_cost_order : public cost_vec::order {
  public:
    virtual bool operator()
      (const NTYPE& v0, const NTYPE& v1) const
      { return (v0 < v1); };
  };

  static class decreasing_cost_order decreasing;
  static class increasing_cost_order increasing;

  static NTYPE max(const cost_vec& v);
  static NTYPE min(const cost_vec& v);
  static NTYPE sum(const cost_vec& v);

  // min/max/sum of a subset of elements:
  static NTYPE max(const cost_vec& v, const index_set& s);
  static NTYPE max(const cost_vec& v, const bool_vec& s);
  static NTYPE min(const cost_vec& v, const index_set& s);
  static NTYPE min(const cost_vec& v, const bool_vec& s);
  static NTYPE sum(const cost_vec& v, const index_set& s);
  static NTYPE sum(const cost_vec& v, const bool_vec& s);

  NTYPE max() const { return max(*this); };
  NTYPE min() const { return min(*this); };
  NTYPE sum() const { return sum(*this); };
};

class cost_matrix : public matrix<NTYPE> {
 public:
  cost_matrix()
    : matrix<NTYPE>() { };
  cost_matrix(const NTYPE _val, index_type r, index_type c)
    : matrix<NTYPE>(_val, c, r) { };
  cost_matrix(const cost_matrix& _mat)
    : matrix<NTYPE>(_mat) { };
  cost_matrix(const cost_matrix& _mat,
	      const bool_vec& rs,
	      const bool_vec& cs)
    : matrix<NTYPE>(_mat, rs, cs) { };
  cost_matrix(const cost_matrix& _mat,
	      const index_set& rs,
	      const index_set& cs)
    : matrix<NTYPE>(_mat, rs, cs) { };

  void row_max(cost_vec& m);
  void row_min(cost_vec& m);
  void column_max(cost_vec& m);
  void column_min(cost_vec& m);

  void transitive_closure();
};

struct interval : public comparable_pair<NTYPE> {
  interval(const NTYPE& v1, const NTYPE& v2) :
    comparable_pair<NTYPE>(v1, v2) { };
  interval(const NTYPE& v) :
    comparable_pair<NTYPE>(v) { };
  interval(const interval& p) :
    comparable_pair<NTYPE>(p) { };
  interval() :
    comparable_pair<NTYPE>(NEG_INF, POS_INF) { };
};

typedef std::pair<index_type, NTYPE> index_cost_pair;
typedef lvector<index_cost_pair> index_cost_vec;

class index_cost_vec_util : public index_cost_vec
{
 public:

  class decreasing_cost_order : public index_cost_vec::order {
  public:
    virtual bool operator()
      (const index_cost_pair& v0, const index_cost_pair& v1) const
      { return (v0.second > v1.second); };
  };

  class increasing_cost_order : public index_cost_vec::order {
  public:
    virtual bool operator()
      (const index_cost_pair& v0, const index_cost_pair& v1) const
      { return (v0.second < v1.second); };
  };

  static class decreasing_cost_order decreasing_cost;
  static class increasing_cost_order increasing_cost;
};

inline std::ostream& operator<<(std::ostream& s, const index_cost_pair& p)
{
  s << '(' << p.first << ',' << p.second << ')';
}

inline std::ostream& operator<<(std::ostream& s, const interval& i)
{
  s << '[' << i.first << ',' << i.second << ']';
}

END_HSPS_NAMESPACE

#endif
