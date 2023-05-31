
include makedefs


CLASPDIR = /home/qgmafe/TUNI/IJCAI23/clasp
CPLEXDIR = /home/qgmafe/TUNI/IJCAI23/CPLEX-ACADEMIC/cplex
CONCERTDIR = /home/qgmafe/TUNI/IJCAI23/CPLEX-ACADEMIC/concert

#C_INCLUDE_PATH += $(CPLEX_PATH)/cplex/include
#C_INCLUDE_PATH += $(CPLEX_PATH)/concert/include
#CPLUS_INCLUDE_PATH += $(CPLEX_PATH)/cplex/include
#CPLUS_INCLUDE_PATH += $(CPLEX_PATH)/concert/include
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CLASPINCDIR   = $(CLASPDIR)/
POTINCDIR   = $(CLASPDIR)/libpotassco
CCOPT = -m64 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG
COPT  = -m64 -fPIC
CLASP_PATH = /home/qgmafe/TUNI/IJCAI23/clasp
CPLEX_PATH = /home/qgmafe/TUNI/IJCAI23/CPLEX-ACADEMIC
CPLEX_ARCH = x86-64_linux
CLASP_LIB =  $(CLASP_PATH)/lib
CPLEX_LIB =  $(CPLEX_PATH)/cplex/lib/$(CPLEX_ARCH)/static_pic
CPLEX_LIB2 = $(CPLEX_PATH)/concert/lib/$(CPLEX_ARCH)/static_pic
LLIBS  = -ldl -lm -I$(CONCERTINCDIR) -I$(CPLEXINCDIR) -I$(CLASPINCDIR) -I$(POTINCDIR) #-lmalloc
CPLEX_CFLAGS = -DIL_STD -I$(CONCERTINCDIR) -I$(CPLEXINCDIR)
CLASP_CFLAGS = -DIL_STD -I$(CLASPINCDIR) -I$(POTINCDIR)
CPLEX_LFLAGS = -L$(CPLEX_LIB) -L$(CPLEX_LIB2) -I$(CONCERTINCDIR) -I$(CPLEXINCDIR)
CLASP_LFLAGS = -L$(CLASP_LIB) -I$(CLASPINCDIR) -I$(POTINCDIR)
CPLEX_LLIBS  = -lilocplex -lcplex -lconcert -lpthread -ldl
CLASP_LLIBS  = -lclasp -lpotassco -ldl



ifeq ($(buildwithglpk), 1)
  CFLAGS += $(GLPK_CFLAGS)
  CFLAGS += -DHAVE_GLPK
  LFLAGS += $(GLPK_LFLAGS)
  LLIBS += $(GLPK_LLIBS)
endif

ifeq ($(buildwithscip), 1)
  CFLAGS += $(SCIP_CFLAGS)
  CFLAGS += -DHAVE_SCIP
  LFLAGS += $(SCIP_LFLAGS)
  LLIBS += $(SCIP_LLIBS)
endif

ifeq ($(buildwithcplex), 1)
  CFLAGS += $(CLASP_CFLAGS)
  CFLAGS += $(CPLEX_CFLAGS)
  CFLAGS += -DHAVE_CPLEX
  LFLAGS += $(CPLEX_LFLAGS)
  LFLAGS += $(CLASP_LFLAGS)
  LLIBS += $(CPLEX_LLIBS)
  LLIBS += $(CLASP_LLIBS)
  
endif

ifeq ($(buildwithgecode), 1)
  CFLAGS += $(GECODE_CFLAGS)
  CFLAGS += -DHAVE_GECODE
  LFLAGS += $(GECODE_LFLAGS)
  LLIBS += $(GECODE_LLIBS)
endif

.SUFFIXES : .cc $(SUFFIXES)

.cc.o:
	$(CCC) $(CFLAGS) -c $*.cc

.c.o:
	$(CC) $(CFLAGS) -c $*.c

%.profile.o: %.cc
	$(CCC) $(CFLAGS) $(PROFILE) -c $*.cc -o $*.profile.o

%.shared.o: %.cc
	$(CCC) $(CFLAGS) $(SHARED) -c $*.cc -o $*.shared.o

%.d: %.cc
	$(CCDEP) $(CFLAGS) $< > $@

LIBUTILSRCS = rational.cc index_type.cc factor.cc numeric_type.cc \
 graph.cc hypergraph.cc enumerators.cc rng.cc stats.cc

LIBPLANSRCS = name.cc problem.cc search.cc search_base.cc heuristic.cc \
 exec.cc plans.cc mutex.cc preprocess.cc cost_table.cc htn.cc soft.cc \
 reduce.cc additive.cc lmcut.cc atomset.cc resource.cc seq_reg.cc \
 para_reg.cc temporal.cc forward.cc hashtable.cc ida.cc idao.cc minimax.cc \
 itest.cc nodeset.cc bfs.cc bb.cc bfhs.cc bss.cc collector.cc random.cc \
 explore.cc sas.cc sas_heuristic.cc pdb.cc pdb_construction.cc simplify.cc \
 tcn.cc pop.cc mdd.cc path.cc ilb.cc pcc.cc

LIBPDDLGENSRCS = grammar.cc scanner.cc
LIBPDDLNGSRCS = ptr_table.cc char_map.cc string_table.cc base.cc parser.cc \
 logic.cc ext_state.cc

SUPSRCS = graph_grammar.cc cda_hplus_interface.cc
GENSRCS = pddl2.y pddl2.lex

PRGSRCS = pddlcat.cc hsp0.cc hsp_f.cc hsp_a.cc hsp_b.cc hsp_c.cc d_hsp.cc \
 rpgc.cc hsp_p.cc
EMPRGSRCS = hsp_e.cc hsp_h.cc test_apx.cc strips.cc \
 convert_cost_to_preference.cc create_nb_problems.cc test_nbb.cc \
 test_ppc.cc test_sd.cc analyze_deadlines.cc random_sas.cc clockwork.cc \
 roverMSP_to_TSP.cc findprime.cc rnd.cc pddl2pdl.cc pddl2ctl.cc \
 pddl2ltl.cc pddl2qbf.cc schedule-with-til.cc test_ilb.cc \
 test_cegar.cc make-state-space-graph.cc
NMPRGSRCS = check_ra.cc sample.cc ggp.cc pddbg.cc sas_us_plan.cc sbsp.cc \
 test_reach.cc test_misc.cc

EXTRASRCS = locals/config.h function.h chopshop.c
DOCS =
UTILS = locals/fake-flex++ locals/fake-bison++ locals/scanner.*.h \
 locals/scanner.*.cc locals/grammar.*.h locals/grammar.*.cc \
 plan-to-xhtml.xsl mkbuildscript.sh ipc14-run-hpp.sh ipc14-run-hpp-ce.sh \
 ipc14-readme.txt

LIBTARGETS = libhspsutil.a
NORMALTARGETS =
EXTRATARGETS =
ifeq ($(buildlibplan), 1)
  LIBTARGETS += libhspsplan.a
endif
ifeq ($(buildlibpddl), 1)
  LIBTARGETS += libhspspddl.a
  NORMALTARGETS += $(PRGSRCS:.cc=$(EXEEXT))
  EXTRATARGETS += $(EMPRGSRCS:.cc=$(EXEEXT))
endif

all: $(LIBTARGETS) $(NORMALTARGETS)
libs: $(LIBTARGETS)
extra: $(EXTRATARGETS)

grammar.h: pddl2.y base.h
	$(BISON) -v -d -o grammar.cc pddl2.y
grammar.cc: pddl2.y base.h
	$(BISON) -v -d -o grammar.cc pddl2.y
scanner.cc: pddl2.lex grammar.h
	$(FLEX) -i -hscanner.h -oscanner.cc pddl2.lex
scanner.h: pddl2.lex grammar.h
	$(FLEX) -i -hscanner.h -oscanner.cc pddl2.lex

locals/grammar.linux.h: grammar.h
	cp grammar.h locals/grammar.linux.h
locals/grammar.linux.cc: grammar.cc
	cp grammar.cc locals/grammar.linux.cc
locals/scanner.linux.h: scanner.h
	cp scanner.h locals/scanner.linux.h
locals/scanner.linux.cc: scanner.cc
	cp scanner.cc locals/scanner.linux.cc

update: locals/grammar.linux.h locals/grammar.linux.cc \
 locals/scanner.linux.h locals/scanner.linux.cc

libhspsutil.a: $(LIBUTILSRCS:.cc=.o)
	ar cr libhspsutil.a $(LIBUTILSRCS:.cc=.o)
	ranlib libhspsutil.a

libhspsplan.a: $(LIBPLANSRCS:.cc=.o)
	ar cr libhspsplan.a $(LIBPLANSRCS:.cc=.o)
	ranlib libhspsplan.a

libhspspddl.a: $(LIBPDDLNGSRCS:.cc=.o) $(LIBPDDLGENSRCS:.cc=.o)
	ar cr libhspspddl.a $(LIBPDDLNGSRCS:.cc=.o) $(LIBPDDLGENSRCS:.cc=.o)
	ranlib libhspspddl.a

libcdahplus: $(LIBUTILSRCS:.cc=.shared.o) $(LIBPLANSRCS:.cc=.shared.o) $(LIBPDDLNGSRCS:.cc=.shared.o) $(LIBPDDLGENSRCS:.cc=.shared.o) cda_hplus_interface.shared.o
	gcc -shared -Wl,-soname,libcdahplus.so.1 -o libcdahplus.so.1.0.1 $(LIBUTILSRCS:.cc=.shared.o) $(LIBPLANSRCS:.cc=.shared.o) $(LIBPDDLNGSRCS:.cc=.shared.o) $(LIBPDDLGENSRCS:.cc=.shared.o) cda_hplus_interface.shared.o

PROFILEOBJS=$(LIBUTILSRCS:.cc=.profile.o) $(LIBPLANSRCS:.cc=.profile.o) $(LIBPDDLNGSRCS:.cc=.profile.o) $(LIBPDDLGENSRCS:.cc=.profile.o)

pddlcat$(EXEEXT): pddlcat.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) pddlcat.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o pddlcat$(EXEEXT)

hsp0$(EXEEXT): hsp0.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) hsp0.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o hsp0$(EXEEXT)

hsp_a$(EXEEXT): hsp_a.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) hsp_a.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o hsp_a$(EXEEXT)

hsp_b$(EXEEXT): hsp_b.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) hsp_b.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o hsp_b$(EXEEXT)

hsp_c$(EXEEXT): hsp_c.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) hsp_c.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o hsp_c$(EXEEXT)

d_hsp$(EXEEXT): d_hsp.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) d_hsp.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o d_hsp$(EXEEXT)

hsp_f$(EXEEXT): hsp_f.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) hsp_f.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o hsp_f$(EXEEXT)

hsp_h$(EXEEXT): hsp_h.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) hsp_h.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o hsp_h$(EXEEXT)

hsp_e$(EXEEXT): hsp_e.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) hsp_e.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o hsp_e$(EXEEXT)

hsp_g$(EXEEXT): hsp_g.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) hsp_g.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o hsp_g$(EXEEXT)

hsp_p$(EXEEXT): hsp_p.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) hsp_p.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o hsp_p$(EXEEXT)

rpgc$(EXEEXT): rpgc.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) rpgc.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o rpgc$(EXEEXT)

chopshop$(EXEEXT): chopshop.o
	$(LINK) $(LFLAGS) chopshop.o $(LLIBS) -o chopshop$(EXEEXT)

create_nb_problems$(EXEEXT): create_nb_problems.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) create_nb_problems.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o create_nb_problems$(EXEEXT)

test_nbb$(EXEEXT): test_nbb.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_nbb.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_nbb$(EXEEXT)

test_ppc$(EXEEXT): test_ppc.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_ppc.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_ppc$(EXEEXT)

test_ppc.profile$(EXEEXT): test_ppc.profile.o $(PROFILEOBJS)
	$(LINK) $(LFLAGS) $(PROFILE) test_ppc.profile.o $(PROFILEOBJS) $(LLIBS) $(PROFILELIBS) -o test_ppc.profile$(EXEEXT)

test_sd$(EXEEXT): test_sd.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_sd.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_sd$(EXEEXT)

analyze_deadlines$(EXEEXT): analyze_deadlines.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) analyze_deadlines.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o analyze_deadlines$(EXEEXT)

convert_cost_to_preference$(EXEEXT): convert_cost_to_preference.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) convert_cost_to_preference.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o convert_cost_to_preference$(EXEEXT)

random_sas$(EXEEXT): random_sas.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) random_sas.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o random_sas$(EXEEXT)

clockwork$(EXEEXT): clockwork.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) clockwork.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o clockwork$(EXEEXT)

roverMSP_to_TSP$(EXEEXT): roverMSP_to_TSP.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) roverMSP_to_TSP.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o roverMSP_to_TSP$(EXEEXT)

sas_us_plan$(EXEEXT): sas_us_plan.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) sas_us_plan.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o sas_us_plan$(EXEEXT)

pddl2pdl$(EXEEXT): pddl2pdl.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) pddl2pdl.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o pddl2pdl$(EXEEXT)

pddl2ctl$(EXEEXT): pddl2ctl.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) pddl2ctl.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o pddl2ctl$(EXEEXT)

pddl2ltl$(EXEEXT): pddl2ltl.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) pddl2ltl.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o pddl2ltl$(EXEEXT)

pddl2qbf$(EXEEXT): pddl2qbf.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) pddl2qbf.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o pddl2qbf$(EXEEXT)

schedule-with-til$(EXEEXT): schedule-with-til.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) schedule-with-til.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o schedule-with-til$(EXEEXT)

sbsp$(EXEEXT): sbsp.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) sbsp.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o sbsp$(EXEEXT)

make-state-space-graph$(EXEEXT): make-state-space-graph.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) make-state-space-graph.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o make-state-space-graph$(EXEEXT)

pddbg$(EXEEXT): pddbg.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) pddbg.o -lhspspddl -lhspsplan -lhspsutil -lreadline -lboost_regex $(LLIBS) -o pddbg$(EXEEXT)

# create sample_current_config.cc by copying if necessary
sample_current_config.cc: sample.cc
	cp -f sample.cc sample_current_config.cc

sample$(EXEEXT): sample_current_config.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) sample_current_config.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o sample$(EXEEXT)

##
## misc. test programs
##

check_ra$(EXEEXT): check_ra.o libhspsutil.a
	$(LINK) $(LFLAGS) check_ra.o -lhspsutil $(LLIBS) -o check_ra$(EXEEXT)

# test_misc$(EXEEXT): test_misc.o libhspsutil.a libhspsplan.a libhspspddl.a
# 	$(LINK) $(LFLAGS) test_misc.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_misc$(EXEEXT)
test_misc$(EXEEXT): test_misc.o libhspsutil.a
	$(LINK) $(LFLAGS) test_misc.o -lhspsutil $(LLIBS) -o test_misc$(EXEEXT)

test_prf$(EXEEXT): test_prf.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_prf.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_prf$(EXEEXT)

test_bdo$(EXEEXT): test_bdo.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_bdo.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_bdo$(EXEEXT)

test_apx$(EXEEXT): test_apx.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_apx.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_apx$(EXEEXT)

strips$(EXEEXT): strips.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) strips.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o strips$(EXEEXT)

findprime$(EXEEXT): findprime.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) findprime.o -lhspsutil $(LLIBS) -o findprime$(EXEEXT)

rnd$(EXEEXT): rnd.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) rnd.o -lhspsutil $(LLIBS) -o rnd$(EXEEXT)

ggp$(EXEEXT): ggp.o graph_grammar.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) ggp.o graph_grammar.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o ggp$(EXEEXT)

test_ilb$(EXEEXT): test_ilb.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_ilb.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_ilb$(EXEEXT)

test_cbd$(EXEEXT): test_cbd.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_cbd.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_cbd$(EXEEXT)

test_ilb.profile$(EXEEXT): test_ilb.profile.o $(PROFILEOBJS)
	$(LINK) $(LFLAGS) $(PROFILE) test_ilb.profile.o  $(PROFILEOBJS) $(LLIBS) $(PROFILELIBS) -o test_ilb.profile$(EXEEXT)

test_cegar$(EXEEXT): test_cegar.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_cegar.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_cegar$(EXEEXT)

bdpi$(EXEEXT): bdpi.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) bdpi.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o bdpi$(EXEEXT)

## reach test

REACH_OBJS = \
$(REACH_DIR)/reach_api.o $(REACH_DIR)/reach.o $(REACH_DIR)/initial.o \
$(REACH_DIR)/hs.o $(REACH_DIR)/list_utils.o $(REACH_DIR)/stats.o

test_reach$(EXEEXT): test_reach.o libhspsutil.a libhspsplan.a libhspspddl.a
	$(LINK) $(LFLAGS) test_reach.o $(REACH_OBJS) -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o test_reach$(EXEEXT)


##
## admin targets
##

clean:
	rm -f *.d *.o *.a $(PRGSRCS:.cc=$(EXEEXT)) $(EMPRGSRCS:.cc=$(EXEEXT)) $(NMPRGSRCS:.cc=$(EXEEXT)) chopshop$(EXEEXT)

distclean:
	make clean
	make -C amse21 clean
	make -C SPP clean
	make -C MOST clean
	make -C pddl3sim clean
	make -C hspIAO clean


remake: clean
	rm -f grammar.* scanner.*


ALLLIBSRCS = $(LIBUTILSRCS) $(LIBPLANSRCS) $(LIBPDDLNGSRCS) $(GENSRCS) $(SUPSRCS)

arc:
	rm -f hsps.tar hsps.tar.gz
	tar cvf hsps.tar $(ALLLIBSRCS:.cc=.h) $(ALLLIBSRCS) $(PRGSRCS) $(EMPRGSRCS) \
 $(NMPRGSRCS) $(EXTRASRCS) makefile
	gzip hsps.tar

unarc: hsps.tar.gz
	gzip -d hsps.tar.gz
	tar xvf hsps.tar
	rm -f hsps.tar

dist:
	-rm hsps.tar.gz
	tar cvf hsps.tar $(ALLLIBSRCS:.cc=.h) $(ALLLIBSRCS) $(PRGSRCS) $(EMPRGSRCS) \
 $(NMPRGSRCS) $(EXTRASRCS) $(DOCS) locals/makedefs $(UTILS) \
 makefile
	gzip hsps.tar

libutildist:
	-rm libhspsutil.tar.gz
	tar cvf libhspsutil.tar $(LIBUTILSRCS:.cc=.h) $(LIBUTILSRCS) \
          locals/config.h locals/makedefs makefile readme-libhspsutil
	gzip libhspsutil.tar


##
## build IPC6 dist
##

%.cc_esc: %.cc
	$(SED) -f escape.sed  $< > $@

%.h_esc: %.h
	$(SED) -f escape.sed  $< > $@

%.ccx: %.cc_esc
	$(CC) $(CFLAGS) -DIPC6 -E -P -x c++ $< | $(SED) -f unescape.sed > $@

escape: $(LIBSRCS:.cc=.h_esc) $(LIBSRCS:.cc=.cc_esc) $(LIBGENSRCS:.cc=.h_esc) $(LIBGENSRCS:.cc=.cc_esc) $(PRGSRCS:.cc=.cc_esc) function.h_esc config.h_esc

IPC6COMMONFILES = rational.cc index_type.cc name.cc ptr_table.cc char_map.cc \
 string_table.cc numeric_type.cc graph.cc enumerators.cc rng.cc problem.cc \
 soft.cc search.cc heuristic.cc stats.cc exec.cc plans.cc preprocess.cc \
 cost_table.cc additive.cc atomset.cc resource.cc seq_reg.cc forward.cc \
 search_base.cc hashtable.cc ida.cc nodeset.cc bfs.cc base.cc parser.cc \
 grammar.cc scanner.cc

ipc6-files: escape $(IPC6COMMONFILES:.cc=.ccx) hsp0.ccx hsp_f.ccx hsp_p.ccx

ipc6: escape ipc6-files
	-rm -rf seq-opt-hsps0
	-mkdir seq-opt-hsps0
	-cp *.ccx seq-opt-hsps0/
	./mkbuildscript.sh hsp0.ccx $(IPC6COMMONFILES:.cc=.ccx) > seq-opt-hsps0/build
	chmod 744 seq-opt-hsps0/build
	-cp ipc6-run-hsps0.sh seq-opt-hsps0/plan
	-rm -rf seq-opt-hspsf
	-mkdir seq-opt-hspsf
	-cp *.ccx seq-opt-hspsf/
	./mkbuildscript.sh hsp_f.ccx $(IPC6COMMONFILES:.cc=.ccx) > seq-opt-hspsf/build
	chmod 744 seq-opt-hspsf/build
	-cp ipc6-run-hspsf.sh seq-opt-hspsf/plan
	-rm -rf netben-opt-hspsp
	-mkdir netben-opt-hspsp
	-cp *.ccx netben-opt-hspsp/
	./mkbuildscript.sh hsp_p.ccx $(IPC6COMMONFILES:.cc=.ccx) > netben-opt-hspsp/build
	chmod 744 netben-opt-hspsp/build
	-cp ipc6-run-hspsp.sh netben-opt-hspsp/plan

ipc6-arc: ipc6
	tar cvf hsps-ipc6-all.tar seq-opt-hsps0/* seq-opt-hspsf/* netben-opt-hspsp/*
	gzip -9 hsps-ipc6-all.tar


##
## build IPC-14 dist
##

IPC14LIBFILES = ilb.cc problem.cc index_type.cc numeric_type.cc rational.cc \
 rng.cc name.cc graph.cc stats.cc preprocess.cc mutex.cc heuristic.cc \
 parser.cc scanner.cc string_table.cc char_map.cc grammar.cc base.cc \
 ptr_table.cc search.cc cost_table.cc lmcut.cc plans.cc exec.cc \
 search_base.cc enumerators.cc

IPC14FILES = test_ilb.cc config.h

ipc14:
	-rm -rf seq-opt-hpp
	-mkdir seq-opt-hpp
	-cp $(IPC14LIBFILES) $(IPC14LIBFILES:.cc=.h) $(IPC14FILES) seq-opt-hpp
	env CCC="$(CCC)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" LLIBS="$(LLIBS)" ./mkbuildscript.sh test_ilb.cc $(IPC14LIBFILES) > seq-opt-hpp/build
	chmod 744 seq-opt-hpp/build
	-cp ipc14-run-hpp.sh seq-opt-hpp/plan
	chmod 744 seq-opt-hpp/plan
	-cp ipc14-readme.txt seq-opt-hpp/README
	-rm -rf seq-opt-hpp-ce
	-mkdir seq-opt-hpp-ce
	-cp $(IPC14LIBFILES) $(IPC14LIBFILES:.cc=.h) $(IPC14FILES) seq-opt-hpp-ce
	env CCC="$(CCC)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" LLIBS="$(LLIBS)" ./mkbuildscript.sh test_ilb.cc $(IPC14LIBFILES) > seq-opt-hpp-ce/build
	chmod 744 seq-opt-hpp-ce/build
	-cp ipc14-run-hpp-ce.sh seq-opt-hpp-ce/plan
	chmod 744 seq-opt-hpp-ce/plan
	-cp ipc14-readme.txt seq-opt-hpp-ce/README

uipc16:
	-rm -rf uipc-hpp
	-mkdir uipc-hpp
	-cp $(IPC14LIBFILES) $(IPC14LIBFILES:.cc=.h) $(IPC14FILES) uipc-hpp
	env CCC="$(CCC)" CFLAGS="$(CFLAGS)" LFLAGS="$(LFLAGS)" LLIBS="$(LLIBS)" ./mkbuildscript.sh test_ilb.cc $(IPC14LIBFILES) > uipc-hpp/build
	chmod 744 uipc-hpp/build
	-cp uipc-run-hpp.sh uipc-hpp/plan
	chmod 744 uipc-hpp/plan

##
## experimental stuff
##

lgrammar.h: lisp.y lisp.h
	$(BISON) -v -d -o lgrammar.cc lisp.y
lgrammar.cc: lisp.y lisp.h
	$(BISON) -v -d -o lgrammar.cc lisp.y
lscanner.h: lisp.lex lgrammar.h
	$(FLEX) -i -hlscanner.h -olscanner.cc lisp.lex
lscanner.cc: lisp.lex lgrammar.h
	$(FLEX) -i -hlscanner.h -olscanner.cc lisp.lex

lscanner.o: lscanner.cc lscanner.h lisp.h input.h
lgrammar.o: lgrammar.cc lgrammar.h lisp.h input.h
lisp.o: lisp.cc lisp.h input.h
input.o: input.cc input.h
new_test.o: new_test.cc lscanner.h lgrammar.h lisp.h input.h
new_test: new_test.o input.o lisp.o lgrammar.o lscanner.o
	$(LINK) $(LFLAGS) new_test.o lscanner.o lgrammar.o lisp.o input.o -lhspspddl -lhspsplan -lhspsutil $(LLIBS) -o new_test$(EXEEXT)


ALLDEPFILES = $(LIBUTILSRCS:%.cc=%.d)
ifeq ($(buildlibplan), 1)
  ALLDEPFILES += $(LIBPLANSRCS:%.cc=%.d)
endif
ifeq ($(buildlibpddl), 1)
  ALLDEPFILES += $(EMPRGSRCS:%.cc=%.d)
  ALLDEPFILES += $(NMPRGSRCS:%.cc=%.d)
  ALLDEPFILES += $(PRGSRCS:%.cc=%.d)
  ALLDEPFILES += $(LIBPDDLNGSRCS:%.cc=%.d)
  ALLDEPFILES += $(LIBPDDLGENSRCS:%.cc=%.d)
endif

dep: $(ALLDEPFILES)


##
## include dependencies
##

-include $(LIBUTILSRCS:%.cc=%.d)

ifeq ($(buildlibplan), 1)
-include $(LIBPLANSRCS:%.cc=%.d)
endif

ifeq ($(buildlibpddl), 1)
-include $(EMPRGSRCS:%.cc=%.d)
-include $(NMPRGSRCS:%.cc=%.d)
-include $(PRGSRCS:%.cc=%.d)
-include $(LIBPDDLNGSRCS:%.cc=%.d)
-include $(LIBPDDLGENSRCS:%.cc=%.d)
endif
