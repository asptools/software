
#include "ptr_table.h"

BEGIN_HSPS_NAMESPACE

ptr_table::~ptr_table() {
  if (_next) delete _next;
  if (_left) delete _left;
  if (_right) delete _right;
}

void ptr_table::key_sequence(key_vec& k)
{
  if (_prev) {
    _prev->key_sequence(k);
  }
  else {
    k.set_length(0);
  }
  k.append(_key);
}

ptr_table::key_vec* ptr_table::key_sequence()
{
  key_vec* k = new key_vec((void*)0, 0);
  key_sequence(*k);
  return k;
}

ptr_table* ptr_table::root()
{
  ptr_table* p = this;
  while (p->_prev) {
    p = p->_prev;
  }
  while (p->_up) {
    p = p->_up;
  }
  return p;
}

ptr_table* ptr_table::enum_key_first()
{
  ptr_table* p = this;
  bool leaf = false;
  while (!leaf) {
    if (p->_left) p = p->_left;
    else if (p->_right) p = p->_right;
    else leaf = true;
  }
  return p;
}

ptr_table* ptr_table::enum_key_next()
{
  if (_up) {
    if (_up->_right == this) {
      return _up;
    }
    else if (_up->_right) {
      return _up->_right->enum_key_first();
    }
  }
  return 0;
}

ptr_table* ptr_table::insert(void* k) {
  if (_key == 0) {
    _key = k;
    _count += 1;
    return this;
  }
  if (_key == k) return this;
  else if (k < _key) {
    if (!_left) {
      _left = new ptr_table(k, this, _prev);
      _count += 1;
      return _left;
    }
    ptr_table* r = _left->insert(k);
    _count = _left->_count + (_right ? _right->_count : 0) + 1;
    return r;
  }
  else {
    if (!_right) {
      _right = new ptr_table(k, this, _prev);
      _count += 1;
      return _right;
    }
    ptr_table* r = _right->insert(k);
    _count = (_left ? _left->_count : 0) + _right->_count + 1;
    return r;
  }
}

ptr_table* ptr_table::insert_next(void* k) {
  if (!_next) {
    _next = _next = new ptr_table(k, 0, this);
    return _next;
  }
  return _next->insert(k);
}

ptr_table* ptr_table::insert(void** k, index_type l) {
  ptr_table* r = this;
  for (index_type i = 0; i < l; i++) {
    r = r->insert_next(k[i]);
  }
  return r;
}

ptr_table* ptr_table::insert(key_vec& k) {
  ptr_table* r = this;
  for (index_type i = 0; i < k.length(); i++) {
    r = r->insert_next(k[i]);
  }
  return r;
}

ptr_table* ptr_table::find(void* k) {
  if (_key == 0) return 0;
  if (_key == k) return this;
  if (k < _key) {
    if (!_left) return 0;
    return _left->find(k);
  }
  else {
    if (!_right) return 0;
    return _right->find(k);
  }
}

ptr_table* ptr_table::find_next(void* k) {
  if (!_next) return 0;
  return _next->find(k);
}

ptr_table* ptr_table::find(void** k, index_type l) {
  ptr_table* r = this;
  for (index_type i = 0; (i < l) && r; i++) {
    r = r->find_next(k[i]);
  }
  return r;
}

ptr_table* ptr_table::find(key_vec& k) {
  ptr_table* r = this;
  for (index_type i = 1; (i < k.length()) && r; i++) {
    r = r->find_next(k[i]);
  }
  return r;
}

index_type ptr_table::count_values() {
  return ((_left ? _left->count_values() : 0) +
	  (_right ? _right->count_values() : 0) +
	  (_next ? _next->count_values() : 0) +
	  (val != 0 ? 1 : 0));
}

index_type ptr_table::count_keys() {
  return ((_left ? _left->count_keys() : 0) +
	  (_right ? _right->count_keys() : 0) +
	  (_key != 0 ? 1 : 0));
}

void ptr_table::dump_values(value_vec& vec) {
  if (val != 0) vec.append(val);
  if (_left) _left->dump_values(vec);
  if (_right) _right->dump_values(vec);
  if (_next) _next->dump_values(vec);
}

void ptr_table::dump(cell_vec& vec) {
  vec.append(this);
  if (_left) _left->dump(vec);
  if (_right) _right->dump(vec);
  if (_next) _next->dump(vec);
}

void ptr_table::dump_keys(key_vec& vec) {
  if (_key != 0) vec.append(_key);
  if (_left) _left->dump_keys(vec);
  if (_right) _right->dump_keys(vec);
}

ptr_table::value_vec* ptr_table::values() {
  value_vec* vec = new value_vec((void*)0,0);
  dump_values(*vec);
  return vec;
}

ptr_table::key_vec* ptr_table::keys() {
  key_vec* vec = new key_vec((void*)0,0);
  dump_keys(*vec);
  return vec;
}

END_HSPS_NAMESPACE
