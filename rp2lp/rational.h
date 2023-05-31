#ifndef RATIONAL_H
#define RATIONAL_H

#include "config.h"
#include <iostream>
#include <limits.h>
#include <stdlib.h>

#ifndef LONG_BITS
#define LONG_BITS ilog(LONG_MAX)
#endif

#ifndef RATIONAL_PRINT_MAX_DECIMALS
#define RATIONAL_PRINT_MAX_DECIMALS 9
#endif

BEGIN_HSPS_NAMESPACE

long euclid(long n, long k, long& a, long& b);
long gcd(long n, long k);
long lcm(long n, long k);
unsigned long ilog(unsigned long n);
long imag(long n);

class rational {
  long nm;
  long dv;

 public:
  rational();
  rational(long n);
  rational(long n, long d);
  rational(const rational& r);

  struct XR {
    long x_nm;
    long x_dv;

    XR& operator=(const rational r);
  };
  rational(const XR& x);

  long numerator() const;
  long divisor() const;
  long sign() const;
  bool zero() const;
  bool finite() const;
  bool infinite() const;
  bool integral() const;

  static rational reduce(rational r);
  static void     split(const rational r, long& s, long& k, long& m);
  static rational invert(const rational r);
  static rational infinity(const rational r);
  static rational infinity(const long s);
  static rational floor_to(const rational r, long div);
  static rational ceil_to(const rational r, long div);
  static rational frac(const rational r);
  static rational round(const rational r, long div);
  static rational min(const rational r0, const rational r1);
  static rational max(const rational r0, const rational r1);
  static rational rgcd(const rational r0, const rational r1);
  static rational dtor(double v);
  static rational ator(const char* s);

  rational reduce() const;
  void     split(long& s, long& k, long& m) const;
  rational invert() const;
  rational floor() const;
  rational floor_to(long d) const;
  rational ceil() const;
  rational ceil_to(long d) const;
  rational frac() const;
  rational round(long d) const;
  rational round() const;

  rational operator=(const rational r);
  rational operator=(long n);

  rational operator+=(const rational r);
  rational operator-=(const rational r);
  rational operator*=(const rational r);
  rational operator/=(const rational r);
  rational operator+=(long n);
  rational operator-=(long n);
  rational operator*=(long n);
  rational operator/=(long n);

  double decimal() const;

  struct printdecimal {
    const rational& r;
    printdecimal(const rational& _r) : r(_r) { };
  };

  static bool overflow;
};

bool operator==(const rational r0, const rational r1);
bool operator==(const rational r0, long n1);
bool operator==(long n0, const rational r1);

bool operator!=(const rational r0, const rational r1);
bool operator!=(const rational r0, long n1);
bool operator!=(long n0, const rational r1);

bool rlt(long m0, long d0, long m1, long d1);
bool operator<(const rational r0, const rational r1);
bool operator<=(const rational r0, const rational r1);
bool operator>(const rational r0, const rational r1);
bool operator>=(const rational r0, const rational r1);
rational operator-(const rational r0);
rational operator+(const rational r0, const rational r1);
rational operator-(const rational r0, const rational r1);
rational operator*(const rational r0, const rational r1);
rational operator/(const rational r0, const rational r1);
rational operator+(const rational r0, long n1);
rational operator-(const rational r0, long n1);
rational operator*(const rational r0, long n1);
rational operator/(const rational r0, long n1);
rational operator+(long n0, const rational r1);
rational operator-(long n0, const rational r1);
rational operator*(long n0, const rational r1);
rational operator/(long n0, const rational r1);

rational safeadd(const rational r0, const rational r1);
rational safemul(const rational r0, const rational r1);

// default print op prints the rational in rational form
::std::ostream& operator<<(::std::ostream& s, const rational r);

// prints the rational in decimal form
::std::ostream& operator<<(::std::ostream& s, const rational::printdecimal& r);

// inlines

#include "config.h"

inline rational::XR& rational::XR::operator=(const rational r)
{
  x_nm = r.numerator();
  x_dv = r.divisor();
  return *this;
}

inline rational::rational()
  : nm(0), dv(1) { }

inline rational::rational(long n)
  : nm(n), dv(1) { }

inline rational::rational(long n, long d)
  : nm(d < 0 ? -1*n : n), dv(d < 0 ? -1*d : d) { }

inline rational::rational(const rational& r)
  : nm(r.nm), dv(r.dv) { }

inline rational::rational(const rational::XR& x)
  : nm(x.x_nm), dv(x.x_dv) { };

inline long rational::numerator() const { return nm; }

inline long rational::divisor() const { return dv; }

inline long rational::sign() const
{
  return (nm < 0 ? -1 : (nm > 0 ? 1 : 0));
}

inline bool rational::zero() const
{
  return nm == 0;
}

inline bool rational::finite() const
{
  return dv != 0;
}

inline bool rational::infinite() const
{
  return dv == 0;
}

inline bool rational::integral() const
{
  return dv == 1;
}

inline rational rational::reduce(rational r)
{
  if (r.infinite()) return infinity(r.sign());
  if (r.sign() == 0) return rational(0,1);
  long c = gcd(r.nm, r.dv);
  return rational(r.nm / c, r.dv / c);
}

inline void rational::split(const rational r, long& s, long& k, long& m)
{
  s = r.sign();
  k = ((s * r.nm) / r.dv);
  m = (s * r.nm) - (k * r.dv);
}

inline rational rational::invert(const rational r)
{
  return rational(r.dv, r.nm);
}

inline rational rational::infinity(const long s)
{
  return rational((s < 0 ? -1 : (s > 0 ? 1 : 0)), 0);
}

inline rational rational::infinity(const rational r)
{
  return rational(r.sign(),0);
}

inline rational rational::min(const rational r0, const rational r1)
{
  if (r0.infinite()) {
    if (r0.sign() < 0) return r0;
    else return r1;
  }
  else if (r1.infinite()) {
    if (r1.sign() < 0) return r1;
    else return r0;
  }
  else if (r0 < r1) return r0;
  else return r1;
}

inline rational rational::max(const rational r0, const rational r1)
{
  if (r0.infinite()) {
    if (r0.sign() > 0) return r0;
    else return r1;
  }
  else if (r1.infinite()) {
    if (r1.sign() > 0) return r1;
    else return r0;
  }
  else if (r0 < r1) return r1;
  else return r0;
}

inline rational rational::reduce() const
{
  return reduce(*this);
}

inline void rational::split(long& s, long& k, long& m) const
{
  split(*this, s, k, m);
}

inline rational rational::invert() const
{
  return invert(*this);
}

inline rational rational::floor() const
{
  return floor_to(*this, 1);
}

inline rational rational::floor_to(long d) const
{
  return floor_to(*this, d);
}

inline rational rational::ceil() const
{
  return ceil_to(*this, 1);
}

inline rational rational::ceil_to(long d) const
{
  return ceil_to(*this, d);
}

inline rational rational::frac() const
{
  return frac(*this);
}

inline rational rational::round(long d_max) const
{
  return round(*this, d_max);
}

inline rational rational::round() const
{
  return round(*this, SAFE_RATIONAL_PRECISION);
}

inline rational rational::operator=(const rational r)
{
  nm = r.nm;
  dv = r.dv;
  return *this;
}

inline rational rational::operator=(long n)
{
  nm = n;
  dv = 1;
  return *this;
}

inline rational rational::operator+=(const rational r)
{
  return *this = (*this + r);
}

inline rational rational::operator-=(const rational r)
{
  return *this = (*this - r);
}

inline rational rational::operator*=(const rational r)
{
  return *this = (*this * r);
}

inline rational rational::operator/=(const rational r)
{
  return *this = (*this / r);
}

inline rational rational::operator+=(long n)
{
  return *this = (*this + n);
}

inline rational rational::operator-=(long n)
{
  return *this = (*this - n);
}

inline rational rational::operator*=(long n)
{
  return *this = (*this * n);
}

inline rational rational::operator/=(long n)
{
  return *this = (*this / n);
}

inline double rational::decimal() const { return nm/(double)dv; };

inline bool operator==(const rational r0, const rational r1)
{
  if (r0.infinite() && r1.infinite()) return r0.sign() == r1.sign();
  else return ((r0.numerator() == r1.numerator()) &&
	       (r0.divisor() == r1.divisor()));
}

inline bool operator==(const rational r0, long n1)
{
  return ((r0.numerator() == n1) && (r0.divisor() == 1));
}

inline bool operator==(long n0, const rational r1)
{
  return ((r1.numerator() == n0) && (r1.divisor() == 1));
}

inline bool operator!=(const rational r0, const rational r1)
{
  return !(r0 == r1);
}

inline bool operator!=(const rational r0, long n1)
{
  return !(r0 == n1);
}

inline bool operator!=(long n0, const rational r1)
{
  return !(n0 == r1);
}

inline bool operator<=(const rational r0, const rational r1)
{
  if (r0 == r1) return true;
  return (r0 < r1);
}

inline bool operator>=(const rational r0, const rational r1)
{
  if (r0 == r1) return true;
  return (r0 > r1);
}

// unary minus
inline rational operator-(const rational r0)
{
  return rational(-r0.numerator(), r0.divisor());
}

inline rational operator/(const rational r0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "rational(" << r0 << ") / rational(" << r1 << ")"
	      << ::std::endl;
#endif
  return (r0 * r1.invert());
}

inline rational operator+(long n0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "long(" << n0 << ") + rational(" << r1 << ")"
	    << ::std::endl;
#endif
  return (r1 + n0);
}

inline rational operator-(long n0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "long(" << n0 << ") - rational(" << r1 << ")"
	    << ::std::endl;
#endif
  return (-r1 + n0);
}

inline rational operator*(long n0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "long(" << n0 << ") * rational(" << r1 << ")"
	    << ::std::endl;
#endif
  return (r1 * n0);
}

inline rational operator/(long n0, const rational r1)
{
#ifdef TRACE_PRINT_RATIONAL_ARITHMETIC
  ::std::cerr << "long(" << n0 << ") / rational(" << r1 << ")"
	    << ::std::endl;
#endif
  return (r1.invert() * n0);
}

inline ::std::ostream& operator<<(::std::ostream& s, const rational r)
{
  if (r.infinite()) {
    if (r.sign() < 0) return s << "-INF";
    else return s << "INF";
  }
  else if (r.integral()) {
    return s << r.numerator();
  }
  else {
    return s << r.numerator() << '/' << r.divisor();
  }
}

inline ::std::ostream& operator<<
(::std::ostream& s, const rational::printdecimal& d)
{
  if (d.r.infinite()) {
    if (d.r.sign() < 0) return s << "-INF";
    else return s << "INF";
  }
  else {
    rational ipart = d.r.floor();
    rational dpart = d.r.sign() * d.r.frac();
    s << ipart;
    if (!dpart.zero()) {
      s << '.';
      unsigned int n = 0;
      while (!dpart.zero() && (n < RATIONAL_PRINT_MAX_DECIMALS)) {
	dpart = safemul(dpart, rational(10,1));
	rational digit = dpart.floor();
	s << digit;
	dpart = dpart.frac();
	n += 1;
      }
    }
    return s;
  }
}

END_HSPS_NAMESPACE

#endif
