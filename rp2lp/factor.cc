
#include "factor.h"
#include <stdlib.h>

BEGIN_HSPS_NAMESPACE

void factors(index_type n, index_vec& f)
{
  f.set_length(0);
  fmultiply(n, f);
}

index_type number(const index_vec& f)
{
  index_type n = 1;
  for (index_type k = 0; k < f.length(); k++)
    n = (n * f[k]);
  return n;
}

void fmultiply(index_type n, index_vec& f)
{
  bool prime = false;
  while (!prime) {
    prime = true;
    for (index_type k = 2; (k <= (n/k)) && prime; k++) {
      if ((n % k) == 0) {
	f.append(k);
	n = n/k;
	prime = false;
      }
    }
  }
  if ((n > 1) || (f.length() == 0)) f.append(n);
}

bool fdivide(index_type n, index_vec& f)
{
  index_vec fn(0, 0);
  fmultiply(n, fn);
  for (index_type k = 0; k < fn.length(); k++) {
    index_type p = f.first(fn[k]);
    if (p == no_such_index)
      return false;
    f.remove(p);
  }
  return true;
}

void ffactorial(index_type n, index_vec& f)
{
  f.assign_value(1, 1);
  for (index_type k = 2; k <= n; k++)
    fmultiply(k, f);
}

void fcombinations(index_type n, index_type k, index_vec& f)
{
  f.assign_value(1, 1);
  for (index_type i = n; i > (n - k); i--) {
    fmultiply(i, f);
  }
  for (index_type i = 2; i <= k; i++) {
    bool ok = fdivide(i, f);
    if (!ok) {
      std::cerr << "error in fcombinations: division failed!" << std::endl;
      std::cerr << "n = " << n << ", k = " << k << ", i = " << i
		<< ", f = " << f << std::endl;
      exit(255);
    }
  }
}

END_HSPS_NAMESPACE
