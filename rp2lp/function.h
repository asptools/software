#ifndef FUNCTION_OBJECT_H
#define FUNCTION_OBJECT_H

BEGIN_HSPS_NAMESPACE

template<class D, class R> class MonadicFunction {
 public:
  virtual ~MonadicFunction();
  virtual R operator()(const D& a) const = 0;
};

template<class D, class R>
MonadicFunction<D,R>::~MonadicFunction()
{
  // nada
}

template<class D, class R> class DyadicFunction {
 public:
  virtual ~DyadicFunction();
  virtual R operator()(const D& a, const D& b) const = 0;
};

template<class D, class R>
DyadicFunction<D,R>::~DyadicFunction()
{
  // nada
}

template<class T> class compare_as_less {
 public:
  bool operator()(const T& a, const T& b) const;
};

template<class T> 
bool compare_as_less<T>::operator()(const T& a, const T& b) const
{
  return (a.compare(b) < 0);
}

template<class T> class compare_as_less_ptr {
 public:
  bool operator()(const T& a, const T& b) const;
};

template<class T> 
bool compare_as_less_ptr<T>::operator()(const T& a, const T& b) const
{
  return (a->compare(*b) < 0);
}

template<class T> class compare_as_equal {
 public:
  bool operator()(const T& a, const T& b) const;
};

template<class T> 
bool compare_as_equal<T>::operator()(const T& a, const T& b) const
{
  return (a.compare(b) == 0);
}


END_HSPS_NAMESPACE

#endif
