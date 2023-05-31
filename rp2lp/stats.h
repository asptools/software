#ifndef STATS_H
#define STATS_H

// To use stats.h/.cc separate from the hsp* code:
// 1. comment out the following line:
#include "config.h"

// 2. uncomment _exactly_one_ of the following options:
// #define RSS_FROM_RUSAGE_MAX
// #define RSS_FROM_PSINFO
// #define RSS_FROM_PROCFS_STAT
// #define RSS_FROM_MALLINFO

// 3. uncomment/leave commented out the following as you wish:
// #define RESOURCE_CHECK_INTERVAL 1
// #define PERIODIC_RESOURCE_CHECK
// #define HI_FREQ_STACK_CHECK
// #define EXIT_ON_SECOND_INTERRUPT

#include <sys/time.h>
#include <vector>
#include <iostream>

BEGIN_HSPS_NAMESPACE

// util constant/functions for dealing with timevals
extern const timeval TV_ZERO;
// t = t1 - t2
void timeval_diff(const timeval& t1, const timeval& t2, timeval& t);
// t += t1
void timeval_add(const timeval& t1, timeval& t);


class Stopwatch {
  Stopwatch* top_level;
  Stopwatch* parent;
  std::vector<Stopwatch*> child_vec;

  // accumulated time
  timeval start_t;
  timeval stop_t;
  timeval total_t;

  // limits
  long          time_limit;
  unsigned long memory_limit;
  unsigned long stack_limit;

  // resource state at last read
  static timeval       last_time;
  static unsigned long peak_mem;
  static unsigned long peak_size;
  static unsigned long peak_stack;
  static unsigned long init_stack;

  // signal handlers and related flags
  static volatile bool interrupt_signal_trapped;
  static volatile bool interrupt_signal_raised;
  static volatile bool alarm_signal_trapped;
  static volatile bool alarm_signal_raised;

  static void alarm_handler(int sig);
  static void interrupt_handler(int sig);
  static void check_stack();

  // internal util methods to manage signals/flags/breaks
  void trap_interrupt();
  void untrap_interrupt();
  void set_recursive_breaks();
  void set_monitors(bool need_to_update_res);
  void set_monitors_actual();

  void start1();
  void stop1();
  void stopA();

  void register_child(Stopwatch* child);
  void deregister_child(Stopwatch* child);

  // state of flags, monitors and breaks
  unsigned long flag_state;
  unsigned long break_mask;
  unsigned long rec_break_mask;
  long min_time_limit;
  unsigned long run_level;

 protected:
  void set_flags(unsigned long flags);
  void clear_flags(unsigned long flags);
  void set_breaks(unsigned long flags);
  void clear_breaks(unsigned long flags);

  void* parent_as(long id);
  virtual void* downcast(long id);
  void update_parent_flag();

  // override to insert action to be performed on start/stop/reset
  virtual void on_start() { };
  virtual void on_stop() { };
  virtual void on_reset() { };

  // note: overrides MUST call these methods!
  virtual void check_for_update();
  virtual void update_flags(unsigned long changed);

 public:
  static const long FLAG_PARENT = 1;
  static const long FLAG_INTERRUPT = 2;
  static const long FLAG_TIME = 4;
  static const long FLAG_MEMORY = 8;
  static const long FLAG_STACK = 16;
  static const long FLAG_ERROR = 32;

  static const long STOPWATCH_FLAGS = (FLAG_PARENT | FLAG_INTERRUPT
				       | FLAG_TIME | FLAG_MEMORY
				       | FLAG_STACK | FLAG_ERROR);
  static const long OTHER_FLAGS = ~STOPWATCH_FLAGS;

  static double seconds();

  Stopwatch();
  Stopwatch(Stopwatch* p);
  ~Stopwatch();

  unsigned long flags() const
  {
    return flag_state;
  };

  bool flag_is_set(unsigned long fs) const
  {
    return (flag_state & fs);
  };

  bool enabled_breaks(unsigned long flags) const
  {
    return (break_mask & flags);
  };

  bool recursively_enabled_breaks(unsigned long flags) const
  {
    return (rec_break_mask & flags);
  };

  unsigned long running() const
  {
    return run_level;
  };

  unsigned int depth() const;

  void enable_interrupt();
  void disable_interrupt();
  void enable_time_out(long t);
  void disable_time_out();
  void enable_memory_limit(unsigned long l);
  void disable_memory_limit();
  void enable_stack_limit(unsigned long l);
  void disable_stack_limit();

  bool break_signal_raised();
  void error();

  void start();
  void stop();
  void reset();

  void add(const Stopwatch& sw);

  double time() const;
  double total_time() const;
  static unsigned long peak_memory();
  static unsigned long peak_total_size();
  static unsigned long peak_stack_size();

  virtual void print(::std::ostream& s) const;
  void print_hierarchy(::std::ostream& s) const;
};

inline ::std::ostream& operator<<(::std::ostream& s, Stopwatch& t) {
  t.print(s);
  return s;
}

END_HSPS_NAMESPACE

#endif
