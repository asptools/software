
#include "numeric_type.h"

BEGIN_HSPS_NAMESPACE

cost_vec_util::decreasing_cost_order cost_vec_util::decreasing;
cost_vec_util::increasing_cost_order cost_vec_util::increasing;

NTYPE cost_vec_util::max(const cost_vec& v)
{
  index_type i = v.arg_max();
  return (i == no_such_index ? NEG_INF : v[i]);
}

NTYPE cost_vec_util::min(const cost_vec& v)
{
  index_type i = v.arg_min();
  return (i == no_such_index ? POS_INF : v[i]);
}

NTYPE cost_vec_util::sum(const cost_vec& v)
{
  NTYPE s = 0;
  for (index_type k = 0; k < v.size(); k++)
    s += v[k];
  return s;
}

NTYPE cost_vec_util::max(const cost_vec& v, const index_set& s)
{
  index_type i = v.arg_max(s);
  return (i == no_such_index ? NEG_INF : v[i]);
}

NTYPE cost_vec_util::max(const cost_vec& v, const bool_vec& s)
{
  index_type i = v.arg_max(s);
  return (i == no_such_index ? NEG_INF : v[i]);
}

NTYPE cost_vec_util::min(const cost_vec& v, const index_set& s)
{
  index_type i = v.arg_min(s);
  return (i == no_such_index ? POS_INF : v[i]);
}

NTYPE cost_vec_util::min(const cost_vec& v, const bool_vec& s)
{
  index_type i = v.arg_min();
  return (i == no_such_index ? POS_INF : v[i]);
}

NTYPE cost_vec_util::sum(const cost_vec& v, const index_set& s)
{
  NTYPE sm = 0;
  for (index_type k = 0; k < s.size(); k++) {
    assert(s[k] < v.size());
    sm += v[s[k]];
  }
  return sm;
}

NTYPE cost_vec_util::sum(const cost_vec& v, const bool_vec& s)
{
  NTYPE sm = 0;
  for (index_type k = 0; k < v.size(); k++)
    if (s[k]) sm += v[k];
  return sm;
}

void cost_matrix::row_max(cost_vec& m)
{
  m.set_length(rows());
  for (index_type i = 0; i < rows(); i++)
    m[i] = cost_vec_util::max((*this)[i]);
}

void cost_matrix::row_min(cost_vec& m)
{
  m.set_length(rows());
  for (index_type i = 0; i < rows(); i++)
    m[i] = cost_vec_util::min((*this)[i]);
}

void cost_matrix::column_max(cost_vec& m)
{
  if (rows() == 0) {
    m.set_length(0);
    return;
  }
  m.assign_copy((*this)[0]);
  for (index_type i = 1; i < rows(); i++)
    for (index_type j = 0; j < (*this)[i].size(); j++)
      m[j] = MAX((*this)[i][j], m[j]);
}

void cost_matrix::column_min(cost_vec& m)
{
  if (rows() == 0) {
    m.set_length(0);
    return;
  }
  m.assign_copy((*this)[0]);
  for (index_type i = 1; i < rows(); i++)
    for (index_type j = 0; j < (*this)[i].size(); j++)
      m[j] = MIN((*this)[i][j], m[j]);
}

void cost_matrix::transitive_closure()
{
  assert(rows() == columns());
  for (index_type k = 0; k < rows(); k++)
    for (index_type i = 0; i < rows(); i++)
      for (index_type j = 0; j < rows(); j++)
	(*this)[i][j] = MIN((*this)[i][j], (*this)[i][k] + (*this)[k][j]);
}

index_cost_vec_util::decreasing_cost_order
  index_cost_vec_util::decreasing_cost;
index_cost_vec_util::increasing_cost_order
  index_cost_vec_util::increasing_cost;

END_HSPS_NAMESPACE
