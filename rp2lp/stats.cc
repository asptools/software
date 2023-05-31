
#include "stats.h"
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <signal.h>
#include <alloca.h>
#ifdef RSS_FROM_PROCFS_PSINFO
#include <fcntl.h>
#include <procfs.h>
#include <sstream>
#endif
#ifdef RSS_FROM_PROCFS_STAT
#include <sstream>
#include <fstream>
#endif
#ifdef RSS_FROM_MALLINFO
#include <malloc.h>
#endif
#include <limits.h>
#include <math.h>

#ifndef RESOURCE_CHECK_INTERVAL
#define RESOURCE_CHECK_INTERVAL 1
#endif

BEGIN_HSPS_NAMESPACE

const timeval TV_ZERO = {0, 0};

// t = t1 - t2
void timeval_diff(const timeval& t1, const timeval& t2, timeval& t)
{
  t.tv_sec = t1.tv_sec - t2.tv_sec;
  t.tv_usec = t1.tv_usec - t2.tv_usec;
  while (t.tv_usec < 0) {
    t.tv_sec -= 1;
    t.tv_usec += 1000000;
  }
  while (t.tv_usec > 1000000) {
    t.tv_sec += 1;
    t.tv_usec -= 1000000;
  }
}

// t += t1
void timeval_add(const timeval& t1, timeval& t)
{
  t.tv_sec += t1.tv_sec;
  t.tv_usec += t1.tv_usec;
  while (t.tv_usec > 1000000) {
    t.tv_sec += 1;
    t.tv_usec -= 1000000;
  }
}

timeval       Stopwatch::last_time;
unsigned long Stopwatch::peak_mem = 0;
unsigned long Stopwatch::peak_size = 0;
unsigned long Stopwatch::peak_stack = 0;
unsigned long Stopwatch::init_stack = 0;

volatile bool Stopwatch::interrupt_signal_trapped = false;
volatile bool Stopwatch::interrupt_signal_raised = false;
volatile bool Stopwatch::alarm_signal_trapped = false;
volatile bool Stopwatch::alarm_signal_raised = false;

double Stopwatch::seconds()
{
  rusage u;
  getrusage(RUSAGE_SELF, &u);
  last_time = u.ru_utime;
#ifdef RSS_FROM_PROCFS_PSINFO
  ::std::ostringstream fn;
  fn << "/proc/" << getpid() << "/psinfo";
  int fd = open(fn.str().c_str(), O_RDONLY);
  if (fd < 0) {
    ::std::cerr << "error " << errno << " opening " << fn.str() << ::std::endl;
    exit(255);
  }
  psinfo_t psinfo_data;
  int res = read(fd, &psinfo_data, sizeof(psinfo_t));
  if (res < 0) {
    ::std::cerr << "error " << errno << " reading " << fn.str() << ::std::endl;
    exit(255);
  }
  if (res < sizeof(psinfo_t)) {
    ::std::cerr << "error: only " << res << " bytes of " << sizeof(psinfo_t)
	      << " read from " << fn.str() << ::std::endl;
    exit(255);
  }
  close(fd);
  if (psinfo_data.pr_rssize > Stopwatch::peak_mem) {
    Stopwatch::peak_mem = psinfo_data.pr_rssize;
  }
  if (psinfo_data.pr_size > Stopwatch::peak_size) {
    Stopwatch::peak_size = psinfo_data.pr_size;
  }
#endif
#ifdef RSS_FROM_PROCFS_STAT
  ::std::ostringstream fn;
  fn << "/proc/" << getpid() << "/stat";
  std::ifstream ps_in(fn.str().c_str());
  std::string item;
  for (unsigned int n = 0; n < 23; n++) {
    ps_in >> item;
    // std::cerr << "item #" << n << " = \"" << item << "\"" << std::endl;
  }
  ps_in.close();
  long vm = atoi(item.c_str()) / 1024;
  if (vm > Stopwatch::peak_mem) {
    Stopwatch::peak_mem = vm;
  }
  if (vm > Stopwatch::peak_size) {
    Stopwatch::peak_size = vm;
  }
#endif
#ifdef RSS_FROM_MALLINFO
  struct mallinfo mi = mallinfo();
  unsigned long s = (mi.arena / 1024);
  if (s > Stopwatch::peak_mem) {
    Stopwatch::peak_mem = s;
  }
#endif
#ifdef RSS_FROM_RUSAGE_MAX
  if (u.ru_maxrss > Stopwatch::peak_mem) {
    Stopwatch::peak_mem = u.ru_maxrss;
  }
#endif
  check_stack();
  return (last_time.tv_sec + (last_time.tv_usec / 1000000.0));
}

void Stopwatch::alarm_handler(int sig)
{
  alarm_signal_raised = true;
}

void Stopwatch::interrupt_handler(int sig)
{
#ifdef EXIT_ON_SECOND_INTERRUPT
  if (interrupt_signal_raised) exit(255);
#endif
  interrupt_signal_raised = true;
#ifdef ABORT_ON_INTERRUPT
  abort();
#endif
#ifdef EXIT_ON_SECOND_INTERRUPT
  signal(SIGINT, Stopwatch::interrupt_handler);
#endif
}

void Stopwatch::check_stack()
{
  unsigned long current_stack = (unsigned long)alloca(4);
  if (init_stack == 0) {
    init_stack = current_stack;
  }
  else {
    if (current_stack > init_stack) {
      unsigned long used_stack = ((current_stack - init_stack) / 1024);
      if (used_stack > peak_stack) {
	peak_stack = used_stack;
      }
      init_stack = current_stack;
    }
    else {
      unsigned long used_stack = ((init_stack - current_stack) / 1024);
      if (used_stack > peak_stack) {
	peak_stack = used_stack;
      }
    }
  }
}

Stopwatch::Stopwatch()
  : top_level(0),
    parent(0),
    start_t(TV_ZERO),
    stop_t(TV_ZERO),
    total_t(TV_ZERO),
    time_limit(0),
    memory_limit(0),
    stack_limit(0),
    flag_state(0),
    break_mask(FLAG_PARENT | FLAG_ERROR),
    rec_break_mask(0),
    min_time_limit(0),
    run_level(0)
{
  // done
}

Stopwatch::Stopwatch(Stopwatch* p)
  : top_level(0),
    parent(p),
    start_t(TV_ZERO),
    stop_t(TV_ZERO),
    total_t(TV_ZERO),
    time_limit(0),
    memory_limit(0),
    stack_limit(0),
    flag_state(0),
    break_mask(FLAG_PARENT | FLAG_ERROR),
    rec_break_mask(0),
    min_time_limit(0),
    run_level(0)
{
  if (parent->top_level)
    top_level = parent->top_level;
  else
    top_level = parent;
  parent->register_child(this);
#ifdef TRACE_PRINT_STOPWATCH
  std::cerr << "Stopwatch " << this << " created as child of "
	    << parent << std::endl;
#endif
}

Stopwatch::~Stopwatch()
{
  // really should have a nice clean-up for this case,
  // but I can't be bothered to implement it right now
  // assert(child_vec.size() == 0);
  if (parent) {
    parent->deregister_child(this);
  }
#ifdef TRACE_PRINT_STOPWATCH
    std::cerr << "Stopwatch " << this << " destroyed" << std::endl;
#endif
}

void Stopwatch::trap_interrupt()
{
  signal(SIGINT, Stopwatch::interrupt_handler);
  interrupt_signal_trapped = true;
  interrupt_signal_raised = false;
}

void Stopwatch::untrap_interrupt()
{
  signal(SIGINT, SIG_DFL);
  interrupt_signal_trapped = false;
  interrupt_signal_raised = false;
}

void Stopwatch::check_for_update()
{
  unsigned long changed = 0;
  if (interrupt_signal_trapped && interrupt_signal_raised) {
    changed = (changed | FLAG_INTERRUPT);
  }
  if (alarm_signal_trapped && alarm_signal_raised) {
    seconds();
    // std::cerr << "RESOURCE CHECK: ";
    // print(std::cerr);
    // std::cerr << std::endl;
    changed = (changed | FLAG_TIME | FLAG_MEMORY | FLAG_STACK);
  }
#ifdef HI_FREQ_STACK_CHECK
  else {
    check_stack();
    changed = (changed | FLAG_STACK);
  }
#endif
  if (changed) {
    update_flags(changed);
  }
  if (alarm_signal_trapped && alarm_signal_raised) {
    alarm_signal_trapped = false;
    set_recursive_breaks(); // update recursive breaks and min_time_limit
    set_monitors_actual(); // set new alarm if required
  }
}

void Stopwatch::update_flags(unsigned long changed)
{
  // update parent flag to reflect value passed from above
  if (changed & FLAG_PARENT)
    flag_state = (flag_state | FLAG_PARENT);
  else
    flag_state = (flag_state & ~FLAG_PARENT);
  // check if error flag is set
  if (changed & FLAG_ERROR)
    flag_state = (flag_state | FLAG_ERROR);
  // check and set other flags we care about
  if ((run_level > 0) && enabled_breaks(changed)) {
    if (changed & FLAG_INTERRUPT) {
      if (enabled_breaks(FLAG_INTERRUPT) &&
	  ((flag_state & FLAG_INTERRUPT) == 0))
	std::cerr << "interrupt: " << *this << std::endl;
      flag_state = (flag_state | FLAG_INTERRUPT);
    }
    if (enabled_breaks(FLAG_TIME)) {
      long t = (last_time.tv_sec - start_t.tv_sec) + total_t.tv_sec;
      if (t >= time_limit) {
	if ((flag_state & FLAG_TIME) == 0)
	  std::cerr << "time-out " << total_time() << std::endl;
	flag_state = (flag_state | FLAG_TIME);
      }
    }
    if (enabled_breaks(FLAG_MEMORY))
      if (peak_memory() > memory_limit) {
	if ((flag_state & FLAG_MEMORY) == 0)
	  std::cerr << "memory limit exceeded: " << peak_memory() << std::endl;
	flag_state = (flag_state | FLAG_MEMORY);
      }
    if (enabled_breaks(FLAG_STACK))
      if (peak_stack_size() > stack_limit) {
	if ((flag_state & FLAG_STACK) == 0)
	  std::cerr << "stack limit exceeded: "
		    << peak_stack_size() << std::endl;
	flag_state = (flag_state | FLAG_STACK);
      }
  }
  // if we're breaking, pass change with parent flag set down
  if (flag_state & break_mask)
    changed = (changed | FLAG_PARENT);
  if (recursively_enabled_breaks(changed)) {
    for (unsigned int k = 0; k < child_vec.size(); k++)
      child_vec[k]->update_flags(changed);
  }
}

void Stopwatch::error()
{
  if (top_level) {
    top_level->error();
  }
  else {
    update_flags(FLAG_ERROR);
  }
}

void* Stopwatch::parent_as(long id)
{
  return (parent ? parent->downcast(id) : 0);
}

void* Stopwatch::downcast(long id)
{
  if (id == 0)
    return this;
  else
    return 0;
}

unsigned int Stopwatch::depth() const
{
  if (parent)
    return (parent->depth() + 1);
  else
    return 0;
}

void Stopwatch::update_parent_flag()
{
  if (parent) {
    if (parent->flag_state & parent->break_mask)
      flag_state = (flag_state | FLAG_PARENT);
  }
}

void Stopwatch::set_recursive_breaks()
{
  min_time_limit = LONG_MAX - 1;
  if (run_level > 0) {
    rec_break_mask = break_mask;
    if (enabled_breaks(FLAG_TIME)) {
      min_time_limit =
	time_limit - ((last_time.tv_sec - start_t.tv_sec) + total_t.tv_sec);
      if (min_time_limit <= 0) {
	flag_state = (flag_state | FLAG_TIME);
      }
      if (enabled_breaks(FLAG_MEMORY | FLAG_STACK))
	if (min_time_limit > RESOURCE_CHECK_INTERVAL)
	  min_time_limit = RESOURCE_CHECK_INTERVAL;
    }
    else if (enabled_breaks(FLAG_MEMORY | FLAG_STACK)) {
      min_time_limit = RESOURCE_CHECK_INTERVAL;
    }
  }
  else {
    rec_break_mask = (FLAG_PARENT | FLAG_ERROR);
  }
  for (unsigned int k = 0; k < child_vec.size(); k++) {
    child_vec[k]->set_recursive_breaks();
    rec_break_mask = (rec_break_mask | child_vec[k]->rec_break_mask);
//    std::cerr << "sw@" << this << ": child " << k << " = sw@" << child_vec[k]
//	      << " has " << child_vec[k]->rec_break_mask
// 	      << " / " << child_vec[k]->min_time_limit << std::endl;
    if (child_vec[k]->rec_break_mask & FLAG_TIME) {
      if (child_vec[k]->min_time_limit < min_time_limit)
	min_time_limit = child_vec[k]->min_time_limit;
    }
  }
}

void Stopwatch::set_monitors(bool need_to_update_res)
{
  if (top_level) {
    top_level->set_monitors(need_to_update_res);
  }
  else { // we are the top level
    if (need_to_update_res) { // check/update resource state
      seconds();
    }
    unsigned long changed = FLAG_TIME | FLAG_MEMORY | FLAG_STACK;
    if (interrupt_signal_trapped && interrupt_signal_raised) {
      changed = (changed | FLAG_INTERRUPT);
    }
    update_flags(changed);
    set_recursive_breaks(); // update recursive breaks and min_time_limit
    set_monitors_actual();
  }
}

void Stopwatch::set_monitors_actual()
{
//   std::cerr << "set_monitors_actual: rec_break_mask = " << rec_break_mask
// 	    << ", min_time_limit = " << min_time_limit << std::endl;
  if (rec_break_mask & FLAG_INTERRUPT) {
    if (!interrupt_signal_trapped)
      trap_interrupt();
  }
  else {
    if (interrupt_signal_trapped)
      untrap_interrupt();
  }
  if (rec_break_mask & (FLAG_TIME | FLAG_MEMORY | FLAG_STACK)) {
    if (min_time_limit > 0) {
      alarm_signal_trapped = true;
      alarm_signal_raised = false;
      signal(SIGALRM, Stopwatch::alarm_handler);
      alarm(min_time_limit);
    }
  }
#ifdef PERIODIC_RESOURCE_CHECK
  else if (!alarm_signal_trapped) {
    alarm_signal_trapped = true;
    alarm_signal_raised = false;
    signal(SIGALRM, Stopwatch::alarm_handler);
    alarm(RESOURCE_CHECK_INTERVAL);
  }
#endif
}

void Stopwatch::enable_interrupt()
{
  set_breaks(FLAG_INTERRUPT);
  if (running()) set_monitors(true);
}

void Stopwatch::disable_interrupt()
{
  clear_breaks(FLAG_INTERRUPT);
  if (running()) set_monitors(true);
}

void Stopwatch::enable_time_out(long t)
{
  time_limit = t;
  flag_state = (flag_state & ~FLAG_TIME);
  set_breaks(FLAG_TIME);
  if (running()) set_monitors(true);
}

void Stopwatch::disable_time_out()
{
  clear_breaks(FLAG_TIME);
  if (running()) set_monitors(true);
}

void Stopwatch::enable_memory_limit(unsigned long l)
{
  memory_limit = l;
  flag_state = (flag_state & ~FLAG_MEMORY);
  set_breaks(FLAG_MEMORY);
  if (running()) set_monitors(true);
}

void Stopwatch::disable_memory_limit()
{
  clear_breaks(FLAG_MEMORY);
  if (running()) set_monitors(true);
}

void Stopwatch::enable_stack_limit(unsigned long l)
{
  stack_limit = l;
  flag_state = (flag_state & ~FLAG_STACK);
  set_breaks(FLAG_STACK);
  if (running()) set_monitors(true);
}

void Stopwatch::disable_stack_limit()
{
  clear_breaks(FLAG_STACK);
  if (running()) set_monitors(true);
}

bool Stopwatch::break_signal_raised()
{
  if (top_level) {
    top_level->check_for_update();
  }
  else { // we are the top level
    check_for_update();
  }
  return (flag_state & break_mask);
}

void Stopwatch::start1()
{
  run_level += 1;
  if (run_level == 1) { // transit from stopped to running
#ifdef TRACE_PRINT_STOPWATCH
    std::cerr << "start1 " << this << " setting start_t" << std::endl;
#endif
    start_t = last_time;
    on_start();
  }
  if (parent)
    parent->start1();
}

void Stopwatch::start()
{
#ifdef TRACE_PRINT_STOPWATCH
  std::cerr << "enter start: " << this
	    << ", run_level = " << run_level
	    << std::endl;
#endif
  seconds(); // update clock
  start1();  // recursive start
  if (run_level == 1) // if transit from stopped, set monitors
    set_monitors(false);
#ifdef TRACE_PRINT_STOPWATCH
    std::cerr << "exit start: " << this << ", run_level = " << run_level
	      << ", start_t = "
	      << (start_t.tv_sec + (start_t.tv_usec / 1000000.0))
	      << std::endl;
#endif
}

void Stopwatch::stop1()
{
  if (run_level == 0) return;
  run_level -= 1;
  if (run_level == 0) { // transit from running to stopped
    timeval_diff(last_time, start_t, stop_t);
    timeval_add(stop_t, total_t);
    on_stop();
  }
  if (parent)
    parent->stop1();
}

void Stopwatch::stopA()
{
  for (unsigned int k = 0; k < child_vec.size(); k++) {
    child_vec[k]->stopA();
  }
  if (run_level == 0) return;
  run_level = 0;
  timeval_diff(last_time, start_t, stop_t);
  timeval_add(stop_t, total_t);
  on_stop();
}

void Stopwatch::stop()
{
#ifdef TRACE_PRINT_STOPWATCH
  std::cerr << "enter stop: " << this << ", run_level = " << run_level
	    << std::endl;
#endif
  if (run_level == 0) return;
  seconds();
  stop1();
  if (run_level == 0) {
    for (unsigned int k = 0; k < child_vec.size(); k++) {
      child_vec[k]->stopA();
    }
    set_monitors(false);
  }
#ifdef TRACE_PRINT_STOPWATCH
  std::cerr << "exit stop: " << this << ", run_level = " << run_level
	    << std::endl;
#endif
}

void Stopwatch::register_child(Stopwatch* child)
{
  child_vec.push_back(child);
}

void Stopwatch::deregister_child(Stopwatch* child)
{
  std::vector<Stopwatch*>::iterator p = child_vec.begin();
  while (p != child_vec.end()) {
    if (*p == child) {
      child_vec.erase(p);
      p = child_vec.end();
    }
    else {
      p++;
    }
  }
#ifdef TRACE_PRINT_STOPWATCH
  std::cerr << "child " << child << " dereged from " << this << std::endl;
#endif
}

void Stopwatch::reset()
{
  total_t.tv_sec = 0;
  total_t.tv_usec = 0;
  stop_t.tv_sec = 0;
  stop_t.tv_usec = 0;
  flag_state = 0;
  on_reset();
  if (run_level > 0) {
    seconds();
    start_t = last_time;
    set_monitors(false);
  }
}

void Stopwatch::add(const Stopwatch& sw)
{
  timeval_add(sw.total_t, total_t);
}

void Stopwatch::set_flags(unsigned long flags)
{
  flag_state = (flag_state | (flags & OTHER_FLAGS));
}

void Stopwatch::clear_flags(unsigned long flags)
{
  flag_state = (flag_state & (~flags | ~OTHER_FLAGS));
}

void Stopwatch::set_breaks(unsigned long flags)
{
  break_mask = (break_mask | flags | FLAG_PARENT | FLAG_ERROR);
}

void Stopwatch::clear_breaks(unsigned long flags)
{
  break_mask = ((break_mask & ~flags) | FLAG_PARENT | FLAG_ERROR);
}

double Stopwatch::time() const
{
  if (run_level > 0) {
    seconds();
    timeval current_t;
    timeval_diff(last_time, start_t, current_t);
    return (current_t.tv_sec + (current_t.tv_usec / 1000000.0));
  }
  else {
    return (stop_t.tv_sec + (stop_t.tv_usec / 1000000.0));
  }
}

double Stopwatch::total_time() const
{
  if (run_level > 0)
    return (total_t.tv_sec + (total_t.tv_usec / 1000000.0)) + time();
  else
    return (total_t.tv_sec + (total_t.tv_usec / 1000000.0));
}

unsigned long Stopwatch::peak_memory()
{
  return peak_mem;
}

unsigned long Stopwatch::peak_total_size()
{
#if (defined(RSS_FROM_PROCFS_PSINFO) || defined(RSS_FROM_PROCFS_STAT))
  return peak_size;
#else
  return peak_mem + peak_stack;
#endif
}

unsigned long Stopwatch::peak_stack_size()
{
  return peak_stack;
}

void Stopwatch::print_hierarchy(::std::ostream& s) const
{
  if (parent) {
    parent->print_hierarchy(s);
  }
  s << "sw@" << this << ": ";
  print(s);
  s << std::endl;
}

void Stopwatch::print(::std::ostream& s) const
{
  s << total_time() << " sec., " << peak_memory() << "k heap, "
    << peak_stack_size() << "k stack, flags: " << flags()
    << " (" << depth() << "/" << break_mask << ")";
}

END_HSPS_NAMESPACE
