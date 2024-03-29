#ifndef COST_TABLE_H
#define COST_TABLE_H

#include "config.h"
#include "heuristic.h"

#ifdef BUILD_WITH_BOOSTING
#include "search_base.h"
#include "hashtable.h"
#endif

BEGIN_HSPS_NAMESPACE

void sort2(index_type& i, index_type& j);
void sort3(index_type& i, index_type& j, index_type& l);

class PairIndexFunction : index_vec {
 public:
  PairIndexFunction(index_type n);
  ~PairIndexFunction();

  index_type operator()(index_type i, index_type j) const;
  index_type max_in() const;
  index_type max_out() const;
};

class CostNode {
 public:
  struct Value {
    NTYPE     val;
    bool      opt;
    CostNode* next;
    bool      mark;
    count_type work;

    Value() : val(0), opt(false), next(0), mark(false), work(0) { };
    Value(NTYPE v, bool o) : val(v), opt(o), next(0), mark(false), work(0) { };

    Value& operator=(const Value& v)
    { val = v.val; opt = v.opt; return *this; };
    Value& operator=(NTYPE v) { val = v; return *this; };
    bool operator==(const Value& v) { return (val == v.val); };
    bool operator<(const Value& v)  { return (val < v.val); };
    bool operator<=(const Value& v) { return (val <= v.val); };
    bool operator>(const Value& v)  { return (val > v.val); };
    bool operator>=(const Value& v) { return (val >= v.val); };
    bool operator==(const NTYPE v) { return (val == v); };
    bool operator<(const NTYPE v)  { return (val < v); };
    bool operator<=(const NTYPE v) { return (val <= v); };
    bool operator>(const NTYPE v)  { return (val > v); };
    bool operator>=(const NTYPE v) { return (val >= v); };

    void clear() {
      val = 0;
      opt = false;
      next = 0;
      mark = false;
      work = 0;
    };
    void clear(NTYPE v_default) {
      val = v_default;
      opt = false;
      next = 0;
      mark = false;
      work = 0;
    };
  };

 private:
  index_type _size;
  index_type _depth;
  index_type _first;
  NTYPE      _default;
  NTYPE      _max;
  Value*     _store;
  CostNode*  _prev;
  CostNode(CostNode* p, index_type f); /* not a copy constructor! */

  NTYPE eval(const index_set& s, index_type i);
  NTYPE eval(const bool_vec& s, index_type i);
  NTYPE incremental_eval(const bool_vec& s, index_type i, index_type i_new);
  NTYPE eval_to_bound(const index_set& s, index_type i, NTYPE bound);
  NTYPE eval_to_bound(const bool_vec& s, index_type i, NTYPE bound);
  void  find_all(const index_set& s, index_type i, const index_set& c,
		 index_set_vec& a);

  void set_max(NTYPE v);
  void write(std::ostream& s, index_set q) const;
  void write_pddl(std::ostream& s, Instance& ins, index_set q) const;

 public:
  static bool eval_set_mark;
#ifdef EVAL_EXTRA_STATS
  static count_type eval_rec_count;
#endif

  CostNode(index_type size);
  CostNode(index_type size, NTYPE v_default);
  ~CostNode();

  index_type size() const { return _size; };
  index_type depth() const { return _depth; };
  index_type first() const { return _first; };
  index_type last() const { return _size - 1; };

  Value&    val(index_type p);
  const Value& val(index_type p) const;
  NTYPE     max() const { return _max; };
  NTYPE     max_finite() const;
  CostNode& next(index_type p);
  Value*    val_p(index_type p);
  CostNode* next_p(index_type p);

  Value*    find(const index_set& s);
  Value*    find(const bool_vec& s);
  void      find_all(const index_set& s, index_set_vec& a);
  Value&    insert(const index_set& s);
  Value&    insert(const bool_vec& s);

  void store(index_type p, NTYPE v, bool opt);
  void store(index_type p, NTYPE v);

  NTYPE eval(const index_set& s);
  NTYPE eval(const bool_vec& s);
  extended_cost extended_eval(const index_set& s);
  extended_cost extended_eval(const bool_vec& s);
  NTYPE incremental_eval(const index_set& s, index_type i_new);
  NTYPE incremental_eval(const bool_vec& s, index_type i_new);
  NTYPE eval_to_bound(const index_set& s, NTYPE bound);
  NTYPE eval_to_bound(const bool_vec& s, NTYPE bound);

  void store(const index_set& s, NTYPE v, bool opt);
  void store(const bool_vec& s, NTYPE v, bool opt);
  void store(const index_set& s, NTYPE v);
  void store(const bool_vec& s, NTYPE v);

  NTYPE eval_min(const index_set& s);

  void compute_max();
  void clear_marks();
  void clear();

  void write(std::ostream& s) const;
  void write_pddl(std::ostream& s, Instance& ins) const;
};

class CostTable : public Heuristic, public CostNode {
 public:

  struct Entry {
    index_set set;
    CostNode::Value& val;
    Entry*    prev;
    Entry*    next;

    Entry(index_set s, CostNode::Value& v)
      : set(s), val(v), prev(0), next(0) { };

    Entry& operator=(const Entry& e);

    bool operator==(const Entry& e) { return set == e.set; };
    bool operator<(const Entry& e) { return val < e.val; };
    bool operator<=(const Entry& e) { return val <= e.val; };
    bool operator>(const Entry& e) { return val > e.val; };
    bool operator>=(const Entry& e) { return val >= e.val; };

    Entry* first();
    Entry* last();
    Entry* move_down();
    Entry* move_up();
    void unlink();
    void place_before(Entry* p);
    void place_after(Entry* p);

    void delete_list();
    count_type list_length();
  };

 protected:
  Stopwatch& stats;
  NTYPE* pre_cost;
  NTYPE* per_cost;

  // value iteration update, returns true iff table was changed
  bool update(index_type i, NTYPE v, bool_vec& f);
  bool update(index_type i, NTYPE v);
  bool update(index_type i, index_type j, NTYPE v);
  bool update(index_type i, index_type j, index_type l, NTYPE v);

  Entry* entries(CostNode* n, index_set s, Entry* l,
		 const index_set* filter, bool ign_zero, bool ign_opt,
		 bool ign_inf, bool req_mark, bool sort);
  count_type count_entries(CostNode* n, index_set s, bool ign_zero,
			   bool ign_opt, bool ign_inf);
  void fill(index_set& s, index_type d);

  void compute_hm_graph_min
    (const ACF& cost, const index_set& s, bool strict, index_set_graph& g);
  void compute_hm_graph_max
    (const ACF& cost, const index_set& s, bool strict, index_set_graph& g);

  void extend_goal_set(const index_set& sgset,
		       const AnyACF& costs,
		       bool_vec& ext_goal_set);
  void find_cut(const bool_vec& ext_goal_set,
		const AnyACF& costs,
		bool_vec& cut);

  friend class ILB; // so ILB can use extend_goal_set and find_cut

 public:
  CostTable(const Instance& i, Stopwatch& s);
  CostTable(const Instance& i, Stopwatch& s, NTYPE v_default);
  ~CostTable();

  virtual NTYPE eval(index_type i);
  virtual NTYPE eval(index_type i, index_type j);
  virtual NTYPE eval(index_type i, index_type j, index_type k);
  virtual NTYPE eval(const index_set& s)
    { return CostNode::eval(s); };
  virtual NTYPE eval(const bool_vec& s)
    { return CostNode::eval(s); };

  virtual NTYPE incremental_eval(const index_set& s, index_type i_new)
    { return CostNode::incremental_eval(s, i_new); };
  virtual NTYPE incremental_eval(const bool_vec& s, index_type i_new)
    { return CostNode::incremental_eval(s, i_new); };
  virtual NTYPE eval_to_bound(const index_set& s, NTYPE bound)
    { return CostNode::eval_to_bound(s, bound); };
  virtual NTYPE eval_to_bound(const bool_vec& s, NTYPE bound)
    { return CostNode::eval_to_bound(s, bound); };

  virtual void store(const index_set& s, NTYPE v, bool opt)
    { CostNode::store(s, v, opt); };
  virtual void store(const bool_vec& s, NTYPE v, bool opt)
    { CostNode::store(s, v, opt); };
  virtual void store(const index_set& s, NTYPE v)
    { CostNode::store(s, v); };
  virtual void store(const bool_vec& s, NTYPE v)
    { CostNode::store(s, v); };

  void store_list(Entry* l);
  static Entry* insert_entry(Entry* e, Entry* l);
  static Entry* prepend_entry(Entry* e, Entry* l);
  static Entry* remove_entry(Entry* e, Entry* l);
  static Entry* find_entry(index_set& s, Entry* l);
  static Entry* copy_list(Entry* l);
  count_type    write_list(std::ostream& s, Entry* l);

  bool   interesting(const index_set& s);

  Entry* atoms(); // all atoms
  Entry* pairs(); // all pairs
  Entry* entries(); // all entries
  // entries satisfying some filtering conditions
  Entry* entries(bool ign_zero, bool ign_opt, bool ign_inf, bool req_mark);
  Entry* entries(const index_set& s);
  // entries satisfying filtering conditions, unsorted
  Entry* unsorted_entries(bool ign_zero, bool ign_opt, bool ign_inf,
			  bool req_mark);
  // some predefined combinations of filtering conditions
  Entry* boostable_entries();
  Entry* boostable_entries(const index_set& s);
  Entry* boost_cd_entries();
  Entry* marked_entries(bool ign_nb);
  // count entries satisfying some filtering conditions
  count_type count_entries(bool ign_zero, bool ign_opt, bool ign_inf);

  NTYPE action_pre_cost(index_type i);
  NTYPE action_per_cost(index_type i);

  bool compatible(index_type a, index_type b, bool opt_resources) const;

  void compute_H1(const ACF& cost);
  void compute_H1(const ACF& cost, const bool_vec& init,
		  const bool_vec* aa = 0);
  void compute_H1max(const ACF& cost, const bool_vec& init,
		     const bool_vec* aa = 0);
  void compute_H2(const ACF& cost);
  void compute_H2(const ACF& cost, const bool_vec& init,
		  const bool_vec* aa = 0);
  void compute_H3(const ACF& cost,
		  const bool_vec* aa = 0);
  void compute_H2C(const ACF& cost, bool opt_resources);
  // void compute_H2C(const ACF& cost) { compute_H2C(cost, false); };
  void verify_H2C(const ACF& cost, bool opt_resources, bool opt_verbose);
  void compute_action_cost();
  void fill(index_type to_depth);
  void clear();

  NTYPE compute_lmcut(const ACF& cost, const bool_vec& s,
		      const index_set& g, const bool_vec* a = 0);
  NTYPE compute_lmcut(const ACF& cost, const index_set& s,
		      const index_set& g, const bool_vec* a = 0);

  // what are these for??
  // void closure(const ACF& cost);
  // void closure(const ACF& cost, const index_set& nd);

  void compute_hm_graph(const ACF& cost, const index_set& s,
			bool strict, index_set_graph& g);
  void compute_complete_hm_graph(const ACF& cost, index_type m,
				 const index_set& s, bool strict,
				 index_set_graph& g);

  bool simple_cas(const index_set& g, const ACF& cost, index_type m,
		  index_set& cas);
  bool simple_cas_max(const index_set& g, NTYPE c, const ACF& cost,
		      index_type m, const index_set& cas_in,
		      index_set_vec& stk, index_set& cas_out);
  bool simple_cas_min(const index_set& g, NTYPE c, const ACF& cost,
		      index_type m, const index_set& cas_in,
		      index_set_vec& stk, index_set& cas_out);

  bool critical_tree(const index_set& g,
		     const ACF& cost,
		     index_type m,
		     index_set_vec& cta);

  bool critical_tree_max(const index_set& g, NTYPE c,
			 const ACF& cost, index_type m,
			 index_set_vec& ps,
			 bool_vec& ra,
			 index_set_vec& cta);
  bool critical_tree_min(const index_set& g, NTYPE c,
			 const ACF& cost, index_type m,
			 index_set_vec& ps,
			 bool_vec& ra,
			 index_set_vec& cta);
};

#ifdef BUILD_WITH_BOOSTING
class CostTableBoost : public CostTable {
 protected:
  Statistics& bstats;

  // internal method used by boost_cd
  Entry* create_new_entries(const index_set& set,
			    const index_set& conflict_set,
			    const index_set& cd_ignore,
			    CostTable::Entry* list);

 public:
  static bool strong_conflict_detection;
  static bool ultra_weak_conflict_detection;
  static bool boost_new_entries;
  static count_type n_entries_boosted;
  static count_type n_entries_solved;
  static count_type n_entries_discarded;
  static count_type n_entries_created;
  static count_type n_boost_searches;
  static count_type n_boost_searches_with_cd;

  CostTableBoost(Instance& i, Statistics& s) : CostTable(i, s), bstats(s) { };
  ~CostTableBoost() { };

  Entry* boost_cd(Entry* list, StateFactory& search_space,
		  HashTable* solved_tab,
		  index_type cd_size_limit, const index_set& cd_atoms_ignore,
		  NTYPE cost_limit, const index_set& cost_limit_set,
		  count_type wps_limit, bool scale_wps);
  Entry* boost(Entry* list, StateFactory& search_space, HashTable* solved_tab,
	       NTYPE cost_limit, const index_set& cost_limit_set,
	       count_type wps_limit, bool scale_wps);
  void quick_boost_pass(Entry* list, StateFactory& search_space,
			NTYPE cost_limit, const index_set& cost_limit_set);
};
#endif

class LinearScanH2Eval : public CostTable {
  index_type  n_pairs;
  index_pair* pairs;
  NTYPE*      values;
  bool        complete;
 public:
  LinearScanH2Eval(const Instance& i, Stopwatch& s);
  ~LinearScanH2Eval();

  void compile_finite();
  virtual NTYPE eval(const bool_vec& s);
};

class ForwardH1 : public Heuristic {
  index_set goals;
  const ACF& cost;
  CostTable* table;

 public:
  ForwardH1(const Instance& i, const index_set& g, const ACF& c, Stopwatch& s);
  virtual ~ForwardH1();

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
};

class ForwardH2 : public Heuristic { // aargh!
  index_set goals;
  const ACF& cost;
  CostTable* table;

 public:
  ForwardH2(const Instance& i, const index_set& g, const ACF& c, Stopwatch& s);
  virtual ~ForwardH2();

  virtual NTYPE eval(const index_set& s);
  virtual NTYPE eval(const bool_vec& s);
};


// inlines

inline void sort2(index_type& i, index_type& j)
{
  if (i > j) {
    index_type t = i;
    i = j;
    j = t;
  }
}

inline void sort3(index_type& i, index_type& j, index_type& l)
{
  sort2(i, l);
  sort2(i, j);
  sort2(j, l);
}

inline CostNode::Value& CostNode::val(index_type p)
{
#ifdef CHECK_TABLE_INDEX
  assert((p >= _first) && (p < _size));
#endif
  return _store[p];
}

inline const CostNode::Value& CostNode::val(index_type p) const
{
#ifdef CHECK_TABLE_INDEX
  assert((p >= _first) && (p < _size));
#endif
  return _store[p];
}

inline CostNode::Value* CostNode::val_p(index_type p)
{
#ifdef CHECK_TABLE_INDEX
  assert((p >= _first) && (p < _size));
#endif
  return &(_store[p]);
}

inline CostNode& CostNode::next(index_type p)
{
#ifdef CHECK_TABLE_INDEX
  assert((p >= _first) && (p < _size));
#endif
  if (!_store[p].next) _store[p].next = new CostNode(this, p + 1);
  return *(_store[p].next);
}

inline CostNode* CostNode::next_p(index_type p)
{
#ifdef CHECK_TABLE_INDEX
  assert((p >= _first) && (p < _size));
#endif
  if (!_store[p].next) _store[p].next = new CostNode(this, p + 1);
  return _store[p].next;
}

inline CostTable::Entry& CostTable::Entry::operator=(const Entry& e)
{
  set = e.set;
  val = e.val;
}

inline std::ostream& operator<<(std::ostream& s, const CostNode::Value& v)
{
  s << PRINT_NTYPE(v.val);
  if (v.opt) s << '*';
  return s;
}

inline std::ostream& operator<<(std::ostream& s, const CostNode& t)
{
  t.write(s);
  return s;
}

inline std::ostream& operator<<(std::ostream& s, const CostTable::Entry& e)
{
  if (e.set.length() > 0) {
    s << "{" << e.set[0];
    for (index_type k = 1; k < e.set.length(); k++) s << ',' << e.set[k];
    s << "}:" << e.val;
  }
  else {
    s << "{}:" << e.val;
  }
  return s;
}

END_HSPS_NAMESPACE

#endif
