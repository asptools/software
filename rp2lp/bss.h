#ifndef BSS_H
#define BSS_H

#include "config.h"
#include "search_base.h"
#include "nodeset.h"

BEGIN_HSPS_NAMESPACE

struct signature {
  NTYPE f;
  NTYPE h;
  Node* parent;
  Transition* op;

  signature(NTYPE _f) : f(_f), h(_f), parent(0), op(0) { };
  signature(NTYPE _f, NTYPE _h) : f(_f), h(_h), parent(0), op(0) { };
  signature(NTYPE _f, NTYPE _h, Node* p, Transition* o)
    : f(_f), h(_h), parent(p), op(o) { };
  signature(const Node& n)
    : f(n.val), h(n.est), parent(n.bp_pre), op(n.bp_trans) { };

  bool operator<(const signature& s) const {
    if (f < s.f) return true;
    else if (f > s.f) return false;
    else {
      if (h < s.h) return true;
      else if (h > s.h) return false;
      else {
	// a null parent pointer is > any non-null
	if (parent != NULL && s.parent == NULL) return true;
	if (parent == NULL && s.parent != NULL) return false;
	if (parent != NULL && s.parent != NULL) {
	  int d1 = parent->state->compare(*(s.parent->state));
	  if (d1 < 0) return true;
	  if (d1 > 0) return false;
	  // parent states are equal, so we compare ops: again,
	  // null is considered > any non-null pointer.
	  if (op != NULL && s.op == NULL) return true;
	  if (op == NULL && s.op != NULL) return true;
	  if (op != NULL && s.op != NULL) {
	    int d2 = op->compare(*s.op);
	    if (d2 < 0) return true;
	    return false;
	  }
	  else return false; // both ops are NULL
	}
	else return false; // both parents are NULL
      }
      assert(false);
    }
  };

  bool operator==(const signature& s) const {
    if (f != s.f) return false;
    if (h != s.h) return false;
    if (parent == NULL && s.parent != NULL) return false;
    if (parent != NULL && s.parent == NULL) return false;
    if (parent == NULL) return true; // both parents are NULL
    // compare parents:
    int d1 = parent->state->compare(*(s.parent->state));
    if (d1 != 0) return false;
    if (op == NULL && s.op != NULL) return false;
    if (op != NULL && s.op == NULL) return false;
    if (op == NULL) return true; // both ops are NULL
    // compare ops:
    int d2 = op->compare(*s.op);
    if (d2 != 0) return false;
    return true;
  };

  signature clone() {
    return signature(f, h, parent, op != NULL ? op->copy() : NULL);
  };

}; // end class signature

class BeamStackSearch : public SearchAlgorithm {
 protected:

  struct Layer {
    HashNodeSet* nodes;
    node_vec open;
    index_type next;
    signature L;
    signature U;

    Layer() : nodes(NULL), next(0), L(ZERO), U(POS_INF) { };
    bool empty() const { return next >= open.size(); };
    NTYPE bound() const { return U.f; };
    signature min_open() const {
      if (open.size() == 0)
	return signature(ZERO);
      else
	return signature(*(open[0]));
    };
    signature max_open() const {
      if (open.size() == 0)
	return signature(POS_INF);
      else
	return signature(*(open[open.size() - 1]));
    };
  };

 protected:
  void make_new_layer(index_type l);
  void delete_layer(index_type l);
  void delete_all_layers();
  void reset_layer(index_type l);
  index_type insert_pos(const signature& nss,
			node_vec& open,
			index_type last);

 private:
  lvector<Layer> layers;
  index_type beam_width;
  index_type current; // current layer
  Node* current_node;

 public:
  BeamStackSearch(Statistics& s, SearchResult& r, index_type w);
  virtual ~BeamStackSearch();

  virtual NTYPE start(State& s, NTYPE b);
  virtual NTYPE start(State& s);
  NTYPE main2(State& s);

  NTYPE main(State& s);
  virtual NTYPE new_state(State& s, NTYPE bound);

  virtual NTYPE cost() const;
};

inline std::ostream& operator<<(std::ostream& to, const signature& sig) {
  to << "[" << sig.f << "|" << sig.h << "|";
  if (sig.parent != NULL)
    to << (void*)sig.parent;
  else
    to << "-";
  to << "|";
  if (sig.op != NULL)
    sig.op->write(to);
  else
    to << "-";
  to << "]";
  return to;
};

END_HSPS_NAMESPACE

#endif
