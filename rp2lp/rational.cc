
#include <ctype.h>
#include <math.h>
#include "rational.h"

BEGIN_HSPS_NAMESPACE

const long CONVERSION_MAX_DIVISOR = 100000;
bool rational::overflow = false;

inline double floor_d(double d) {
  return floor(d);
}

long euclid(long n, long k, long& a, long& b) {
  long q = n / k;
  long r = n % k;
  /* n = q*k + r */
  if (r == 0) {
    /* 1*n + (1-q)*k = k */
    a = 1;
    b = 1-q;
    return k;
  }
  else {
    long a1, b1;
    long c = euclid(k, r, a1, b1);
    /* a1*k + b1*r = c, n - q*k = r => a1*k + b1*(n - q*k) = c */
    a = b1;
    b = a1 - (b1*q);
    return c;
  }
}

long gcd(long n, long k) {
  long q = n / k;
  long r = n % k;
  /* n = q*k + r */
  if (r == 0) return k;
  else return gcd(k, r);
}

long lcm(long n, long k) {
  long d = gcd(n, k);
  return ((n / d) * (k / d) * d);
}

unsigned long ilog(unsigned long n)
{
  unsigned long k = 0;
  while (n > 0) {
    k += 1;
    n = (n / 2);
  }
  return k;
}

long imag(long n)
{
  return (n < 0 ? n * -1 : n);
}

rational rational::dtor(double v) {
  long d = 1;
  while (((v - floor_d(v)) > 0.0) && (d < CONVERSION_MAX_DIVISOR)) {
    v *= 10.0;
    d *= 10;
  }
  return reduce(rational((long)floor_d(v), d));
}

rational rational::ator(const char* s) {
  long n = 0;
  long d = 1;
  if (*s == '-') {
    d = -1;
    s += 1;
  }
  while (isdigit(*s)) {
    n = (n*10) + (*s - 48);
    s += 1;
  }
  if (*s == '.') {
    s += 1;
    while (isdigit(*s)) {
      n = (n*10) + (*s - 48);
      d *= 10;
      s += 1;
    }
  }
  return reduce(rational(n, d));
}

rational rational::floor_to(const rational r, long d)
{
  assert(d > 0);
  if (r.dv < d) return r;
  if (r.sign() < 0) return -ceil_to(-r, d);
  long c = gcd(r.dv, d);
  long a = (r.nm * (d/c));
  long b = (r.dv/c);
  return reduce(rational(a/b, d));
}

rational rational::ceil_to(const rational r, long d)
{
  assert(d > 0);
  if (r.dv < d) return r;
  if (r.sign() < 0) return -floor_to(-r, d);
  long c = gcd(r.dv, d);
  long a = (r.nm * (d/c));
  long b = (r.dv/c);
  if ((a % b) == 0) {
    return reduce(rational(a/b, d));
  }
  else {
    return reduce(rational((a/b)+1, d));
  }
}

rational rational::frac(const rational r) {
  if (r.infinite()) return r;
  return reduce(rational(r.nm % r.dv, r.dv));
}

rational rational::round(const rational r, long d_max)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r << ") round " << d_max << ::std::endl;
#endif
  if (r.infinite()) return r;
  rational s(r);
  while (s.dv > d_max) {
    s.dv = (s.dv / 2);
    s.nm = (s.nm / 2);
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
    ::std::cerr << "s = " << s.nm << " / " << s.dv << ::std::endl;
#endif
    s.reduce();
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
    ::std::cerr << "s = " << s << ::std::endl;
#endif
  }
  return s;
}

rational rational::rgcd(const rational r0, const rational r1)
{
  long c = gcd(r0.divisor(), r1.divisor());
  long a0 = r0.numerator() * (r1.divisor() / c);
  long a1 = r1.numerator() * (r0.divisor() / c);
  long d = gcd(a0, a1);
  return rational(d, (r0.divisor() / c) * (r1.divisor() / c) * c).reduce();
}

// subroutine used in comparison of rationals: test if m0/d0 < m1/d1
// assumptions: m0/d0 and m1/d1 are both positive, finite and
// fractional (i.e., m0 < d0 and m1 < d1)
bool rlt(long m0, long d0, long m1, long d1)
{
  // m0/d0 < m1/d1 iff m0 < m1*(d0/d1) iff m0*(d1/d0) < m1
  if (d0 > d1) {
    // d0/d1 = k + m/d1 (> 1): check m0 < m1*k + m1*m/d1
    long k = d0/d1;
    long m = d0 - (k * d1);
    // if m0 < m1*k, return TRUE
    // if m1*k overflows, it's greater than m0 for sure
    if (k > 0)
      if ((m1 - (LONG_MAX % k)) > (LONG_MAX/k))
	return true;
    if (m0 < m1*k) return true;
    // if m1 or m is 0, so is m1*m/d1; return false
    if ((m1 == 0) || (m == 0)) return false;
    // can we compute m1*m without overflow?
    if (LONG_MAX/(m+1) <= m1) {
      return ((m0 - m1*k) < (m1*m)/d1);
    }
    else {
      m0 = (m0 - m1*k);
      long p = 0;
      while ((m0 > 0) && (m1 > 0)) {
	if ((LONG_MAX - p) >= m) {
	  p += m;
	  m1 -= 1;
	  if (p >= d1) {
	    p -= d1;
	    m0 -= 1;
	  }
	}
	else {
	  p -= d1;
	  p += m;
	  m0 -= 1;
	}
      }
      return (m1 > 0);
    }
  }
  else {
    // d1/d0 = k + m/d0 (> 1): check m0*k + m0*m/d0 < m1
    long k = d1/d0;
    long m = d1 - (k * d0);
    // if m0*k >= m1, return FALSE
    // if m0*k overflows, it's greater than m1 for sure
    if (k > 0)
      if ((m0 - (LONG_MAX % k)) > (LONG_MAX/k))
	return false;
    if (m0*k >= m1) return false;
    // if m0 or m is 0, so is m0*m/d0; return TRUE
    if ((m1 == 0) || (m == 0)) return true;
    // can we compute m0*m without overflow?
    if (LONG_MAX/(m+1) <= m0) {
      return ((m0*m/d0) < (m1 - m0*k));
    }
    else {
      m1 = (m1 - m0*k);
      long p = 0;
      while ((m0 > 0) && (m1 > 0)) {
	if ((LONG_MAX - p) >= m) {
	  p += m;
	  m0 -= 1;
	  if (p >= d0) {
	    p -= d0;
	    m1 -= 1;
	  }
	}
	else {
	  p -= d0;
	  p += m;
	  m1 -= 1;
	}
      }
      return (m1 > 0);
    }
  }
}

bool operator<(const rational r0, const rational r1)
{
  if (r0.infinite()) {
    // r0 == -INF: r0 < r1 iff r1 != -INF
    if (r0.sign() < 0) {
      if (r1.infinite() && (r1.sign() < 0))
	return false;
      else
	return true;
    }
    // r0 == +INF: r0 < r1 is FALSE
    return false;
  }
  else if (r1.infinite()) {
    // r0 is finite and r1 == -INF: r0 < r1 is FALSE
    if (r1.sign() < 0)
      return false;
    // r0 is finite and r1 == +INF: r0 < r1 is TRUE
    return true;
  }
  // both are finite
  else {
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
    if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
      return (r0.numerator() < r1.numerator());
    }
#endif
    long s0, k0, m0, s1, k1, m1;
    r0.split(s0, k0, m0);
    r1.split(s1, k1, m1);
    // if r0 is negative
    if (s0 < 0) {
      // if r0 is negative and r1 zero or positive, r0 < r1 is TRUE
      if (s1 >= 0)
	return true;
      // both are negative: r0 < r1 iff k0 + (m0/r0.dv) > k1 + (m1/r1.dv)
      else {
	// if r0 (neg) integral part > r1 (neg) integral, r0 < r1 is TRUE
	if (k0 > k1) return true;
	// if r0 (neg) integral part < r1 (neg) integral, r0 < r1 is FALSE
	if (k0 < k1) return false;
	// integral parts and signs are same, return m0/r0.dv > m1/r1.dv
	return rlt(m1, r1.divisor(), m0, r0.divisor());
	//return ((rational(m0, r0.divisor()) -
	//	 rational(m1, r1.divisor())).sign() > 0);
      }
    }
    // r0 is zero or positive
    else {
      // if r1.sign < r0.sign, r0 < r1 is FALSE
      if (s1 < s0)
	return false;
      // both are zero or pos: r0 < r1 iff k0 + (m0/r1.dv) < k1 + (m1/r1.dv)
      else {
	// if r0 integral part < r1 integral, r0 < r1 is TRUE
	if (k0 < k1) return true;
	// if r0 integral part > r1 integral, r0 < r1 is FALSE
	if (k0 > k1) return false;
	// integral parts and signs are same, return m0/r1.dv < m1/r1.dv
	return rlt(m0, r0.divisor(), m1, r1.divisor());
	//return ((rational(m0, r0.divisor()) -
	//	 rational(m1, r1.divisor())).sign() < 0);
      }
    }
  }
}

bool operator>(const rational r0, const rational r1)
{
  if (r0.infinite()) {
    // r0 == +INF: r0 > r1 iff r1 != +INF
    if (r0.sign() >= 0) {
      if (r1.infinite() && (r1.sign() >= 0))
	return false;
      else
	return true;
    }
    // r0 == -INF: r0 > r1 is FALSE
    return false;
  }
  else if (r1.infinite()) {
    // r0 is finite and r1 == +INF: r0 > r1 is FALSE
    if (r1.sign() >= 0)
      return false;
    // r0 is finite and r1 == -INF: r0 > r1 is TRUE
    return true;
  }
  // both are finite
  else {
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
    if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
      return (r0.numerator() > r1.numerator());
    }
#endif
    long s0, k0, m0, s1, k1, m1;
    r0.split(s0, k0, m0);
    r1.split(s1, k1, m1);
    // if r0 is negative
    if (s0 < 0) {
      // if r0 is negative and r1 zero or positive, r0 > r1 is FALSE
      if (s1 >= 0)
	return false;
      // both are negative: r0 > r1 iff k0 + (m0/r1.dv) < k1 + (m1/r1.dv)
      else {
	// if r0 (neg) integral part < r1 (neg) integral, r0 > r1 is TRUE
	if (k0 < k1) return true;
	// if r0 (neg) integral part > r1 (neg) integral, r0 > r1 is FALSE
	if (k0 > k1) return false;
	// integral parts and signs are same, return m0/r1.dv < m1/r1.dv
	return rlt(m0, r0.divisor(), m1, r1.divisor());
	//return ((rational(m0, r0.divisor()) -
	//	 rational(m1, r1.divisor())).sign() < 0);
      }
    }
    // r0 is zero or positive
    else {
      // if r1.sign < r0.sign, r0 > r1 is TRUE
      if (s1 < s0)
	return true;
      // both are zero or pos: r0 > r1 iff k0 + (m0/r1.dv) > k1 + (m1/r1.dv)
      else {
	// if r0 integral part > r1 integral, r0 > r1 is TRUE
	if (k0 > k1) return true;
	// if r0 integral part < r1 integral, r0 > r1 is FALSE
	if (k0 < k1) return false;
	// integral parts and signs are same, return m0/r1.dv > m1/r1.dv
	return rlt(m1, r1.divisor(), m0, r0.divisor());
	//return ((rational(m0, r0.divisor()) -
	//	 rational(m1, r1.divisor())).sign() > 0);
      }
    }
  }
}

rational safeadd(const rational r0, const rational r1)
{
  if (r1.infinite()) {
    if (r0.infinite()) {
      if (r1.sign() == r0.sign()) return rational::infinity(r0);
      else {
	::std::cerr << "error: " << r0 << " + " << r1 << " not defined"
		    << ::std::endl;
	abort();
      }
    }
    else {
      return rational::infinity(r1);
    }
  }
  else if (r0.infinite()) {
    return rational::infinity(r0);
  }
  else {
    long c = gcd(r0.divisor(), r1.divisor());
    long n0 = r0.numerator();
    long d0 = r0.divisor() / c;
    long n1 = r1.numerator();
    long d1 = r1.divisor() / c;
    assert(d0 > 0);
    assert(d1 > 0);
    while (((LONG_MAX / (2*d1)) < (imag(n0) + 1)) ||
	   ((LONG_MAX / (2*d0)) < (imag(n1) + 1)) ||
	   ((LONG_MAX / (d0 * c)) < (d1 + 1))) {
      if ((d1 < 2) && (d0 < 2)) {
	std::cerr << "error: overflow in safeadd(" << r0 << ", " << r1 << ")"
		  << std::endl;
	abort();
      }
      if (d0 < 2) {
	n1 = (n1 / 2);
	d1 = (d1 / 2);
      }
      else if (d1 < 2) {
	n0 = (n0 / 2);
	d0 = (d0 / 2);
      }
      else if ((imag(n0) - d0) > (imag(n1) - d1)) {
	n0 = (n0 / 2);
	d0 = (d0 / 2);
      }
      else {
	n1 = (n1 / 2);
	d1 = (d1 / 2);
      }
    }
    long f0 = n0 * d1;
    long f1 = n1 * d0;
    long n = f0 + f1;
    long d = d0 * d1 * c;
    return rational(n, d).reduce();
  }
}

rational safemul(const rational r0, const rational r1)
{
  if (r0.infinite() || r1.infinite()) {
    return rational::infinity(r0.sign() * r1.sign());
  }
  else {
    long c0 = gcd(r0.numerator(), r1.divisor());
    long c1 = gcd(r1.numerator(), r0.divisor());
    long n0 = r0.numerator() / c0;
    long n1 = r1.numerator() / c1;
    long d0 = r0.divisor() / c1;
    long d1 = r1.divisor() / c0;
    while (((ilog(n0) + ilog(n1)) > LONG_BITS) ||
	   ((ilog(d0) + ilog(d1)) > LONG_BITS)) {
      if ((d1 < 2) && (d0 < 2)) {
	std::cerr << "error: overflow in safemul(" << r0 << ", " << r1 << ")"
		  << std::endl;
	abort();
      }
      if (d0 < 2) {
	n1 = (n1 / 2);
	d1 = (d1 / 2);
      }
      else if (d1 < 2) {
	n0 = (n0 / 2);
	d0 = (d0 / 2);
      }
      else if (ilog(d0) > ilog(d1)) {
	n0 = (n0 / 2);
	d0 = (d0 / 2);
      }
      else {
	n1 = (n1 / 2);
	d1 = (d1 / 2);
      }
    }
    return rational(n0 * n1, d0 * d1).reduce();
  }
}

rational operator+(const rational r0, const rational r1)
{
  bool overflow = false;
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") + rational(" << r1 << ")"
	    << ::std::endl;
#endif
  if (r1.infinite()) {
    if (r0.infinite()) {
      if (r1.sign() == r0.sign()) return rational::infinity(r0);
      else {
	::std::cerr << "error: " << r0 << " + " << r1 << " not defined"
		    << ::std::endl;
	abort();
      }
    }
    else {
      return rational::infinity(r1);
    }
  }
  else if (r0.infinite()) {
    return rational::infinity(r0);
  }
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
  else if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
    return rational(r0.numerator() + r1.numerator(), 1);
  }
#endif
  else {
    long c = gcd(r0.divisor(), r1.divisor());
    long d0 = (r0.divisor() / c);
    long d1 = (r1.divisor() / c);
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
    long a0 = (r0.numerator() < 0 ? -r0.numerator() : r0.numerator());
    long a1 = (r1.numerator() < 0 ? -r1.numerator() : r1.numerator());
    assert(d1 > 0);
    if ((a0 - (LONG_MAX % d1)) > (LONG_MAX/d1)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "a0 = " << a0 << ", d1 = " << d1
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
      overflow = true;
    }
    assert(d0 > 0);
    if ((a1 - (LONG_MAX % d0)) > (LONG_MAX/d0)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "a1 = " << a1 << ", d0 = " << d0
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
      overflow = true;
    }
    if (LONG_MAX/d0 < r1.divisor()) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "d0 = " << d0 << ", r1.dv = " << r1.divisor()
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
      overflow = true;
    }
#endif // check
    if (!overflow) {
      long n0 = (r0.numerator() * d1);
      long n1 = (r1.numerator() * d0);
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
      if ((n0 > 0) && (n1 > 0)) {
	if ((LONG_MAX - n0) < n1) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
	  std::cerr << "n0 = " << n0 << ", n1 = " << n1
		    << ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
	  overflow = true;
	}
      }
      else if  ((n0 < 0) && (n1 < 0)) {
	if ((LONG_MIN - n1) > n0) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
	  std::cerr << "n0 = " << n0 << ", n1 = " << n1
		    << ", LONG_MIN = " << LONG_MIN << std::endl;
#endif
	  overflow = true;
	}
      }
#endif // check
      if (!overflow) {
	long n = (n0 + n1);
	long d = (d0 * r1.divisor());
	return rational(n, d).reduce();
      }
    }
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
    ::std::cerr << "error: numeric overflow  in " << r0 << " + " << r1
		<< ::std::endl;
#endif
#ifdef RATIONAL_ARITHMETIC_ABORT_ON_OVERFLOW
    abort();
#else
    rational::overflow = true;
#endif
  }
}

rational operator-(const rational r0, const rational r1)
{
  bool overflow = false;
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") - rational(" << r1 << ")"
	      << ::std::endl;
#endif
  if (r1.infinite()) {
    if (r0.infinite()) {
      if (r1.sign() != r0.sign()) return rational::infinity(r0);
      else {
	::std::cerr << "error: " << r0 << " - " << r1 << " not defined"
		    << ::std::endl;
	abort();
      }
    }
    else {
      return rational::infinity(r1.sign() * -1);
    }
  }
  else if (r0.infinite()) {
    return rational::infinity(r0);
  }
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
  else if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
    return rational(r0.numerator() - r1.numerator(), 1);
  }
#endif
  else {
    long c = gcd(r0.divisor(), r1.divisor());
    long d0 = (r0.divisor() / c);
    long d1 = (r1.divisor() / c);
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
    long a0 = (r0.numerator() < 0 ? -r0.numerator() : r0.numerator());
    long a1 = (r1.numerator() < 0 ? -r1.numerator() : r1.numerator());
    assert(d1 > 0);
    if ((a0 - (LONG_MAX % d1)) > (LONG_MAX/d1)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "a0 = " << a0 << ", d1 = " << d1
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
      std::cerr << "lhs = " << (a0 - (LONG_MAX % d1))
		<< ", rhs = " << (LONG_MAX/d1) << std::endl;
#endif
      overflow = true;
    }
    assert(d0 > 0);
    if ((a1 - (LONG_MAX % d0)) > (LONG_MAX/d0)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "a1 = " << a1 << ", d0 = " << d0
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
      overflow = true;
    }
    if (LONG_MAX/d0 < r1.divisor()) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "d0 = " << d0 << ", r1.dv = " << r1.divisor()
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
      overflow = true;
    }
#endif // check
    if (!overflow) {
      long n0 = (r0.numerator() * d1);
      long n1 = (r1.numerator() * d0);
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
      if ((n0 > 0) && (n1 < 0)) {
	if (n0 > (LONG_MAX + n1)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
	  std::cerr << "n0 = " << n0 << ", n1 = " << n1
		    << ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
	  overflow = true;
	}
      }
      else if  ((n0 < 0) && (n1 > 0)) {
	if (n0 < (LONG_MIN + n1)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
	  std::cerr << "n0 = " << n0 << ", n1 = " << n1
		    << ", LONG_MIN = " << LONG_MIN << std::endl;
#endif
	  overflow = true;
	}
      }
#endif // check
      if (!overflow) {
	long n = (n0 - n1);
	long d = (d0 * r1.divisor());
	return rational(n, d).reduce();
      }
    }
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
    ::std::cerr << "error: numeric overflow  in " << r0 << " - " << r1
		<< ::std::endl;
#endif
#ifdef RATIONAL_ARITHMETIC_ABORT_ON_OVERFLOW
    abort();
#else
    rational::overflow = true;
#endif
  }
}

rational operator*(const rational r0, const rational r1)
{
  bool overflow = false;
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") * rational(" << r1 << ")"
	      << ::std::endl;
#endif
  if (r0.infinite()) {
    return rational::infinity(r0.sign() * r1.sign());
  }
  else if (r1.infinite()) {
    return rational::infinity(r0.sign() * r1.sign());
  }
#ifdef RATIONAL_ARITHMETIC_INTEGER_SPECIAL_CASE
  else if ((r0.divisor() == 1) && (r1.divisor() == 1)) {
    return rational(r0.numerator() * r1.numerator(), 1);
  }
#endif
  else {
    long c0 = gcd(r0.numerator(), r1.divisor());
    long c1 = gcd(r1.numerator(), r0.divisor());
    long n0 = (r0.numerator() / c0);
    long n1 = (r1.numerator() / c1);
    long d0 = (r0.divisor() / c1);
    long d1 = (r1.divisor() / c0);
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
    long a0 = (n0 < 0 ? -n0 : n0);
    long a1 = (n1 < 0 ? -n1 : n1);
    if (a1 > 0) {
      if ((a0 - (LONG_MAX % a1)) > (LONG_MAX/a1)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
	std::cerr << "a0 = " << a0 << ", a1 = " << a1
		  << ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
	overflow = true;
      }
    }
    assert((d0 > 0) && (d1 > 0));
    if ((d0 - (LONG_MAX % d1)) > (LONG_MAX/d1)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "d0 = " << d0 << ", d1 = " << d1
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
      overflow = true;
    }
#endif // check
    if (!overflow)
      return rational(n0 * n1, d0 * d1).reduce();
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
    ::std::cerr << "error: numeric overflow  in " << r0 << " * " << r1
		<< ::std::endl;
#endif
#ifdef RATIONAL_ARITHMETIC_ABORT_ON_OVERFLOW
    abort();
#else
    rational::overflow = true;
#endif
  }
}

rational operator+(const rational r0, long n1)
{
  bool overflow = false;
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") + long(" << n1 << ")"
	      << ::std::endl;
#endif
  if (r0.infinite()) {
    return r0;
  }
  else {
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
    long a1 = (n1 < 0 ? -n1 : n1);
    assert(r0.divisor() > 0);
    if ((a1 - (LONG_MAX % r0.divisor())) > (LONG_MAX/r0.divisor())) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "a1 = " << a1 << ", r0.dv = " << r0.divisor()
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
      overflow = true;
    }
#endif // check
    if (!overflow) {
      long m1 = (n1 * r0.divisor());
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
      long m0 = r0.numerator();
      if ((m0 > 0) && (m1 > 0)) {
	if ((LONG_MAX - m0) < m1) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
	  std::cerr << "m0 = " << m0 << ", m1 = " << m1
		    << ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
	  overflow = true;
	}
      }
      else if  ((m0 < 0) && (m1 < 0)) {
	if ((LONG_MIN - m1) > m0) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
	  std::cerr << "m0 = " << m0 << ", m1 = " << m1
		    << ", LONG_MIN = " << LONG_MIN << std::endl;
#endif
	  overflow = true;
	}
      }
#endif // check
      if (!overflow)
	return rational(r0.numerator() + m1, r0.divisor()).reduce();
    }
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
    ::std::cerr << "error: numeric overflow  in " << r0 << " + " << n1
		<< ::std::endl;
#endif
#ifdef RATIONAL_ARITHMETIC_ABORT_ON_OVERFLOW
    abort();
#else
    rational::overflow = true;
#endif
  }
}

rational operator-(const rational r0, long n1)
{
  bool overflow = false;
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") - long(" << n1 << ")"
	      << ::std::endl;
#endif
  if (r0.infinite()) {
    return r0;
  }
  else {
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
    long a1 = (n1 < 0 ? -n1 : n1);
    assert(r0.divisor() > 0);
    if ((a1 - (LONG_MAX % r0.divisor())) > (LONG_MAX/r0.divisor())) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "a1 = " << a1 << ", r0.dv = " << r0.divisor()
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
      overflow = true;
    }
#endif // check
    if (!overflow) {
      long m1 = (n1 * r0.divisor());
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
      long m0 = r0.numerator();
      if ((m0 > 0) && (m1 < 0)) {
	if (m0 > (LONG_MAX + m1)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
	  std::cerr << "m0 = " << m0 << ", m1 = " << m1
		    << ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
	  overflow = true;
	}
      }
      else if  ((m0 < 0) && (m1 > 0)) {
	if (m0 < (LONG_MIN + m1)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
	  std::cerr << "m0 = " << m0 << ", m1 = " << m1
		    << ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
	  overflow = true;
	}
      }
#endif // check
      if (!overflow)
	return rational(r0.numerator() - m1, r0.divisor()).reduce();
    }
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
    ::std::cerr << "error: numeric overflow  in " << r0 << " - " << n1
		<< ::std::endl;
#endif
#ifdef RATIONAL_ARITHMETIC_ABORT_ON_OVERFLOW
    abort();
#else
    rational::overflow = true;
#endif
  }
}

rational operator*(const rational r0, long n1)
{
  bool overflow = false;
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") * long(" << n1 << ")"
	    << ::std::endl;
#endif
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
  long a0 = (r0.numerator() < 0 ? -r0.numerator() : r0.numerator());
  long a1 = (n1 < 0 ? -n1 : n1);
  if (a1 > 0) {
    if ((a0 - (LONG_MAX % a1)) > (LONG_MAX/a1)) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
      std::cerr << "a0 = " << a0 << ", a1 = " << a1
		<< ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
      overflow = true;
    }
  }
#endif // check
  if (!overflow)
    return rational(r0.numerator() * n1, r0.divisor()).reduce();
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
  ::std::cerr << "error: numeric overflow  in " << r0 << " * " << n1
	      << ::std::endl;
#endif
#ifdef RATIONAL_ARITHMETIC_ABORT_ON_OVERFLOW
    abort();
#else
    rational::overflow = true;
#endif
}

rational operator/(const rational r0, long n1)
{
  bool overflow = false;
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") / long(" << n1 << ")"
	    << ::std::endl;
#endif
#ifdef RATIONAL_ARITHMETIC_CHECK_OVERFLOW
  long a1 = (n1 < 0 ? -n1 : n1);
  assert(r0.divisor() > 0);
  if ((a1 - (LONG_MAX % r0.divisor())) > (LONG_MAX/r0.divisor())) {
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
    std::cerr << "a1 = " << a1 << ", r0.dv = " << r0.divisor()
	      << ", LONG_MAX = " << LONG_MAX << std::endl;
#endif
    overflow = true;
  }
#endif
  if (!overflow)
    return rational((n1 < 0 ? -r0.numerator() : r0.numerator()),
		    r0.divisor() * (n1 < 0 ? -n1 : n1)).reduce();
#ifdef RATIONAL_ARITHMETIC_PRINT_OVERFLOW_WARNINGS
  ::std::cerr << "error: numeric overflow  in " << r0 << " / " << n1
	      << ::std::endl;
#endif
#ifdef RATIONAL_ARITHMETIC_ABORT_ON_OVERFLOW
    abort();
#else
    rational::overflow = true;
#endif
}

END_HSPS_NAMESPACE
