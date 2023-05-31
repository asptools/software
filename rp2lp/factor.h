
#ifndef FACTOR_H
#define FACTOR_H

#include "config.h"
#include "index_type.h"

BEGIN_HSPS_NAMESPACE

void factors(index_type n, index_vec& f);
index_type number(const index_vec& f);

void fmultiply(index_type n, index_vec& f);
bool fdivide(index_type n, index_vec& f);

void ffactorial(index_type n, index_vec& f);
void fcombinations(index_type n, index_type k, index_vec& f);

END_HSPS_NAMESPACE

#endif
