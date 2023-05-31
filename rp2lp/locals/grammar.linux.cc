#define YY_PDDL_Parser_h_included

/*  A Bison++ parser, made from pddl2.y  */

 /* with Bison++ version bison++ Version 1.21-8, adapted from GNU bison by coetmeur@icdc.fr
  */


#line 1 "/usr/local/lib/bison.cc"
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Bob Corbett and Richard Stallman

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 1, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

/* HEADER SECTION */
#if defined( _MSDOS ) || defined(MSDOS) || defined(__MSDOS__) 
#define __MSDOS_AND_ALIKE
#endif
#if defined(_WINDOWS) && defined(_MSC_VER)
#define __HAVE_NO_ALLOCA
#define __MSDOS_AND_ALIKE
#endif

#ifndef alloca
#if defined( __GNUC__)
#define alloca __builtin_alloca

#elif (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc)  || defined (__sgi)
#include <alloca.h>

#elif defined (__MSDOS_AND_ALIKE)
#include <malloc.h>
#ifndef __TURBOC__
/* MS C runtime lib */
#define alloca _alloca
#endif

#elif defined(_AIX)
#include <malloc.h>
#pragma alloca

#elif defined(__hpux)
#ifdef __cplusplus
extern "C" {
void *alloca (unsigned int);
};
#else /* not __cplusplus */
void *alloca ();
#endif /* not __cplusplus */

#endif /* not _AIX  not MSDOS, or __TURBOC__ or _AIX, not sparc.  */
#endif /* alloca not defined.  */
#ifdef c_plusplus
#ifndef __cplusplus
#define __cplusplus
#endif
#endif
#ifdef __cplusplus
#ifndef YY_USE_CLASS
#define YY_USE_CLASS
#endif
#else
#ifndef __STDC__
#define const
#endif
#endif
#include <stdio.h>
#define YYBISON 1  

/* #line 73 "/usr/local/lib/bison.cc" */
#line 85 "grammar.cc"
#define YY_PDDL_Parser_ERROR  log_error
#define YY_PDDL_Parser_ERROR_BODY  = 0
#define YY_PDDL_Parser_ERROR_VERBOSE  1
#define YY_PDDL_Parser_LEX  next_token
#define YY_PDDL_Parser_LEX_BODY  = 0
#define YY_PDDL_Parser_DEBUG  1
#define YY_PDDL_Parser_INHERIT  : public HSPS::PDDL_Base
#define YY_PDDL_Parser_CONSTRUCTOR_PARAM  HSPS::StringTable& t
#define YY_PDDL_Parser_CONSTRUCTOR_INIT  : HSPS::PDDL_Base(t), error_flag(false), \
  current_param(0, 0), stored_n_param(0, 0), current_constant_decl(0, 0), \
  current_atom(0), current_atom_stack(0, 0), current_context(0), \
  stored_context(0, 0), current_item(0), current_goal(0, 0), \
  current_preference_name(0), current_entry(0), current_plan_file(0), \
  n_plans_in_current_file(0)
#define YY_PDDL_Parser_MEMBERS  \
public: \
virtual std::ostream& at_position(std::ostream& s) = 0; \
virtual const char*   current_file() = 0; \
bool                  error_flag; \
private: \
HSPS::PDDL_Base::variable_vec current_param; \
HSPS::index_vec               stored_n_param; \
HSPS::PDDL_Base::TypeSet      current_type_set; \
HSPS::PDDL_Base::symbol_vec   current_constant_decl; \
HSPS::index_type              last_n_functions; \
HSPS::PDDL_Base::AtomBase*    current_atom; \
HSPS::PDDL_Base::atom_base_vec current_atom_stack; \
HSPS::PDDL_Base::Context*     current_context; \
HSPS::lvector<Context*> stored_context; \
HSPS::PDDL_Base::DKEL_Item*   current_item; \
HSPS::lvector<ConjunctiveGoal*> current_goal; \
HSPS::PDDL_Base::Symbol*      current_preference_name; \
HSPS::PDDL_Base::HTableEntry* current_entry; \
const char* current_plan_file; \
HSPS::index_type n_plans_in_current_file;
#line 39 "pddl2.y"

#include <stdlib.h>
#include "base.h"
#include <sstream>

#line 45 "pddl2.y"
typedef union {
  HSPS::StringTable::Cell*         sym;
  HSPS::PDDL_Base::Expression*     exp;
  HSPS::PDDL_Base::ListExpression* lst;
  HSPS::PDDL_Base::Relation*       rel;
  HSPS::PDDL_Base::Goal*           goal;
  HSPS::PDDL_Base::Formula*        ff;
  HSPS::PDDL_Base::mode_keyword    tkw;
  HSPS::PDDL_Base::relation_type   rkw;
  HSPS::PDDL_Base::set_constraint_keyword sckw;
  NNTYPE rval;
  int    ival;
  char*  sval;
} yy_PDDL_Parser_stype;
#define YY_PDDL_Parser_STYPE yy_PDDL_Parser_stype

#line 73 "/usr/local/lib/bison.cc"
/* %{ and %header{ and %union, during decl */
#define YY_PDDL_Parser_BISON 1
#ifndef YY_PDDL_Parser_COMPATIBILITY
#ifndef YY_USE_CLASS
#define  YY_PDDL_Parser_COMPATIBILITY 1
#else
#define  YY_PDDL_Parser_COMPATIBILITY 0
#endif
#endif

#if YY_PDDL_Parser_COMPATIBILITY != 0
/* backward compatibility */
#ifdef YYLTYPE
#ifndef YY_PDDL_Parser_LTYPE
#define YY_PDDL_Parser_LTYPE YYLTYPE
#endif
#endif
#ifdef YYSTYPE
#ifndef YY_PDDL_Parser_STYPE 
#define YY_PDDL_Parser_STYPE YYSTYPE
#endif
#endif
#ifdef YYDEBUG
#ifndef YY_PDDL_Parser_DEBUG
#define  YY_PDDL_Parser_DEBUG YYDEBUG
#endif
#endif
#ifdef YY_PDDL_Parser_STYPE
#ifndef yystype
#define yystype YY_PDDL_Parser_STYPE
#endif
#endif
/* use goto to be compatible */
#ifndef YY_PDDL_Parser_USE_GOTO
#define YY_PDDL_Parser_USE_GOTO 1
#endif
#endif

/* use no goto to be clean in C++ */
#ifndef YY_PDDL_Parser_USE_GOTO
#define YY_PDDL_Parser_USE_GOTO 0
#endif

#ifndef YY_PDDL_Parser_PURE

/* #line 117 "/usr/local/lib/bison.cc" */
#line 191 "grammar.cc"

#line 117 "/usr/local/lib/bison.cc"
/*  YY_PDDL_Parser_PURE */
#endif

/* section apres lecture def, avant lecture grammaire S2 */

/* #line 121 "/usr/local/lib/bison.cc" */
#line 200 "grammar.cc"

#line 121 "/usr/local/lib/bison.cc"
/* prefix */
#ifndef YY_PDDL_Parser_DEBUG

/* #line 123 "/usr/local/lib/bison.cc" */
#line 207 "grammar.cc"

#line 123 "/usr/local/lib/bison.cc"
/* YY_PDDL_Parser_DEBUG */
#endif


#ifndef YY_PDDL_Parser_LSP_NEEDED

/* #line 128 "/usr/local/lib/bison.cc" */
#line 217 "grammar.cc"

#line 128 "/usr/local/lib/bison.cc"
 /* YY_PDDL_Parser_LSP_NEEDED*/
#endif



/* DEFAULT LTYPE*/
#ifdef YY_PDDL_Parser_LSP_NEEDED
#ifndef YY_PDDL_Parser_LTYPE
typedef
  struct yyltype
    {
      int timestamp;
      int first_line;
      int first_column;
      int last_line;
      int last_column;
      char *text;
   }
  yyltype;

#define YY_PDDL_Parser_LTYPE yyltype
#endif
#endif
/* DEFAULT STYPE*/
      /* We used to use `unsigned long' as YY_PDDL_Parser_STYPE on MSDOS,
	 but it seems better to be consistent.
	 Most programs should declare their own type anyway.  */

#ifndef YY_PDDL_Parser_STYPE
#define YY_PDDL_Parser_STYPE int
#endif
/* DEFAULT MISCELANEOUS */
#ifndef YY_PDDL_Parser_PARSE
#define YY_PDDL_Parser_PARSE yyparse
#endif
#ifndef YY_PDDL_Parser_LEX
#define YY_PDDL_Parser_LEX yylex
#endif
#ifndef YY_PDDL_Parser_LVAL
#define YY_PDDL_Parser_LVAL yylval
#endif
#ifndef YY_PDDL_Parser_LLOC
#define YY_PDDL_Parser_LLOC yylloc
#endif
#ifndef YY_PDDL_Parser_CHAR
#define YY_PDDL_Parser_CHAR yychar
#endif
#ifndef YY_PDDL_Parser_NERRS
#define YY_PDDL_Parser_NERRS yynerrs
#endif
#ifndef YY_PDDL_Parser_DEBUG_FLAG
#define YY_PDDL_Parser_DEBUG_FLAG yydebug
#endif
#ifndef YY_PDDL_Parser_ERROR
#define YY_PDDL_Parser_ERROR yyerror
#endif
#ifndef YY_PDDL_Parser_PARSE_PARAM
#ifndef __STDC__
#ifndef __cplusplus
#ifndef YY_USE_CLASS
#define YY_PDDL_Parser_PARSE_PARAM
#ifndef YY_PDDL_Parser_PARSE_PARAM_DEF
#define YY_PDDL_Parser_PARSE_PARAM_DEF
#endif
#endif
#endif
#endif
#ifndef YY_PDDL_Parser_PARSE_PARAM
#define YY_PDDL_Parser_PARSE_PARAM void
#endif
#endif
#if YY_PDDL_Parser_COMPATIBILITY != 0
/* backward compatibility */
#ifdef YY_PDDL_Parser_LTYPE
#ifndef YYLTYPE
#define YYLTYPE YY_PDDL_Parser_LTYPE
#else
/* WARNING obsolete !!! user defined YYLTYPE not reported into generated header */
#endif
#endif
#ifndef YYSTYPE
#define YYSTYPE YY_PDDL_Parser_STYPE
#else
/* WARNING obsolete !!! user defined YYSTYPE not reported into generated header */
#endif
#ifdef YY_PDDL_Parser_PURE
#ifndef YYPURE
#define YYPURE YY_PDDL_Parser_PURE
#endif
#endif
#ifdef YY_PDDL_Parser_DEBUG
#ifndef YYDEBUG
#define YYDEBUG YY_PDDL_Parser_DEBUG 
#endif
#endif
#ifndef YY_PDDL_Parser_ERROR_VERBOSE
#ifdef YYERROR_VERBOSE
#define YY_PDDL_Parser_ERROR_VERBOSE YYERROR_VERBOSE
#endif
#endif
#ifndef YY_PDDL_Parser_LSP_NEEDED
#ifdef YYLSP_NEEDED
#define YY_PDDL_Parser_LSP_NEEDED YYLSP_NEEDED
#endif
#endif
#endif
#ifndef YY_USE_CLASS
/* TOKEN C */

/* #line 236 "/usr/local/lib/bison.cc" */
#line 330 "grammar.cc"
#define	TK_OPEN	258
#define	TK_CLOSE	259
#define	TK_OPEN_SQ	260
#define	TK_CLOSE_SQ	261
#define	TK_GREATER	262
#define	TK_GREATEQ	263
#define	TK_LESS	264
#define	TK_LESSEQ	265
#define	TK_COLON	266
#define	TK_HASHT	267
#define	TK_EQ	268
#define	TK_HYPHEN	269
#define	TK_PLUS	270
#define	TK_MUL	271
#define	TK_DIV	272
#define	TK_UMINUS	273
#define	TK_NEW_SYMBOL	274
#define	TK_OBJ_SYMBOL	275
#define	TK_TYPE_SYMBOL	276
#define	TK_PRED_SYMBOL	277
#define	TK_OBJFUN_SYMBOL	278
#define	TK_FUN_SYMBOL	279
#define	TK_VAR_SYMBOL	280
#define	TK_ACTION_SYMBOL	281
#define	TK_MISC_SYMBOL	282
#define	TK_KEYWORD	283
#define	TK_NEW_VAR_SYMBOL	284
#define	TK_PREFERENCE_SYMBOL	285
#define	TK_SET_SYMBOL	286
#define	TK_FLOAT	287
#define	TK_INT	288
#define	TK_STRING	289
#define	KW_REQS	290
#define	KW_CONSTANTS	291
#define	KW_PREDS	292
#define	KW_FUNS	293
#define	KW_TYPES	294
#define	KW_DEFINE	295
#define	KW_DOMAIN	296
#define	KW_ACTION	297
#define	KW_PROCESS	298
#define	KW_EVENT	299
#define	KW_ARGS	300
#define	KW_PRE	301
#define	KW_COND	302
#define	KW_AT_START	303
#define	KW_AT_END	304
#define	KW_OVER_ALL	305
#define	KW_EFFECT	306
#define	KW_INVARIANT	307
#define	KW_DURATION	308
#define	KW_AND	309
#define	KW_OR	310
#define	KW_EXISTS	311
#define	KW_FORALL	312
#define	KW_IMPLY	313
#define	KW_NOT	314
#define	KW_WHEN	315
#define	KW_EITHER	316
#define	KW_PROBLEM	317
#define	KW_FORDOMAIN	318
#define	KW_OBJECTS	319
#define	KW_INIT	320
#define	KW_GOAL	321
#define	KW_LENGTH	322
#define	KW_SERIAL	323
#define	KW_PARALLEL	324
#define	KW_METRIC	325
#define	KW_MINIMIZE	326
#define	KW_MAXIMIZE	327
#define	KW_DURATION_VAR	328
#define	KW_TOTAL_TIME	329
#define	KW_INCREASE	330
#define	KW_DECREASE	331
#define	KW_SCALE_UP	332
#define	KW_SCALE_DOWN	333
#define	KW_ASSIGN	334
#define	KW_TAG	335
#define	KW_NAME	336
#define	KW_VARS	337
#define	KW_SET_CONSTRAINT	338
#define	KW_SETOF	339
#define	KW_AT_LEAST_N	340
#define	KW_AT_MOST_N	341
#define	KW_EXACTLY_N	342
#define	KW_CONTEXT	343
#define	KW_FORMULA	344
#define	KW_IRRELEVANT	345
#define	KW_PLAN	346
#define	KW_HEURISTIC	347
#define	KW_OPT	348
#define	KW_INF	349
#define	KW_FACT	350
#define	KW_SET	351
#define	KW_EXPANSION	352
#define	KW_TASKS	353
#define	KW_PREFERENCE	354
#define	KW_VIOLATED	355
#define	KW_WITHIN	356
#define	KW_ASSOC	357
#define	KW_CONSTRAINTS	358
#define	KW_ALWAYS	359
#define	KW_SOMETIME	360
#define	KW_AT_MOST_ONCE	361
#define	KW_SOMETIME_BEFORE	362
#define	KW_SOMETIME_AFTER	363
#define	KW_ALWAYS_WITHIN	364
#define	KW_IFF	365
#define	KW_FALSE	366
#define	KW_TRUE	367
#define	KW_NUMBER	368
#define	KW_UNDEFINED	369
#define	KW_MINOP	370
#define	KW_MAXOP	371
#define	KW_DERIVED	372


#line 236 "/usr/local/lib/bison.cc"
 /* #defines tokens */
#else
/* CLASS */
#ifndef YY_PDDL_Parser_CLASS
#define YY_PDDL_Parser_CLASS PDDL_Parser
#endif
#ifndef YY_PDDL_Parser_INHERIT
#define YY_PDDL_Parser_INHERIT
#endif
#ifndef YY_PDDL_Parser_MEMBERS
#define YY_PDDL_Parser_MEMBERS 
#endif
#ifndef YY_PDDL_Parser_LEX_BODY
#define YY_PDDL_Parser_LEX_BODY  
#endif
#ifndef YY_PDDL_Parser_ERROR_BODY
#define YY_PDDL_Parser_ERROR_BODY  
#endif
#ifndef YY_PDDL_Parser_CONSTRUCTOR_PARAM
#define YY_PDDL_Parser_CONSTRUCTOR_PARAM
#endif
#ifndef YY_PDDL_Parser_CONSTRUCTOR_CODE
#define YY_PDDL_Parser_CONSTRUCTOR_CODE
#endif
#ifndef YY_PDDL_Parser_CONSTRUCTOR_INIT
#define YY_PDDL_Parser_CONSTRUCTOR_INIT
#endif
/* choose between enum and const */
#ifndef YY_PDDL_Parser_USE_CONST_TOKEN
#define YY_PDDL_Parser_USE_CONST_TOKEN 0
/* yes enum is more compatible with flex,  */
/* so by default we use it */ 
#endif
#if YY_PDDL_Parser_USE_CONST_TOKEN != 0
#ifndef YY_PDDL_Parser_ENUM_TOKEN
#define YY_PDDL_Parser_ENUM_TOKEN yy_PDDL_Parser_enum_token
#endif
#endif

class YY_PDDL_Parser_CLASS YY_PDDL_Parser_INHERIT
{
public: 
#if YY_PDDL_Parser_USE_CONST_TOKEN != 0
/* static const int token ... */

/* #line 280 "/usr/local/lib/bison.cc" */
#line 495 "grammar.cc"
static const int TK_OPEN;
static const int TK_CLOSE;
static const int TK_OPEN_SQ;
static const int TK_CLOSE_SQ;
static const int TK_GREATER;
static const int TK_GREATEQ;
static const int TK_LESS;
static const int TK_LESSEQ;
static const int TK_COLON;
static const int TK_HASHT;
static const int TK_EQ;
static const int TK_HYPHEN;
static const int TK_PLUS;
static const int TK_MUL;
static const int TK_DIV;
static const int TK_UMINUS;
static const int TK_NEW_SYMBOL;
static const int TK_OBJ_SYMBOL;
static const int TK_TYPE_SYMBOL;
static const int TK_PRED_SYMBOL;
static const int TK_OBJFUN_SYMBOL;
static const int TK_FUN_SYMBOL;
static const int TK_VAR_SYMBOL;
static const int TK_ACTION_SYMBOL;
static const int TK_MISC_SYMBOL;
static const int TK_KEYWORD;
static const int TK_NEW_VAR_SYMBOL;
static const int TK_PREFERENCE_SYMBOL;
static const int TK_SET_SYMBOL;
static const int TK_FLOAT;
static const int TK_INT;
static const int TK_STRING;
static const int KW_REQS;
static const int KW_CONSTANTS;
static const int KW_PREDS;
static const int KW_FUNS;
static const int KW_TYPES;
static const int KW_DEFINE;
static const int KW_DOMAIN;
static const int KW_ACTION;
static const int KW_PROCESS;
static const int KW_EVENT;
static const int KW_ARGS;
static const int KW_PRE;
static const int KW_COND;
static const int KW_AT_START;
static const int KW_AT_END;
static const int KW_OVER_ALL;
static const int KW_EFFECT;
static const int KW_INVARIANT;
static const int KW_DURATION;
static const int KW_AND;
static const int KW_OR;
static const int KW_EXISTS;
static const int KW_FORALL;
static const int KW_IMPLY;
static const int KW_NOT;
static const int KW_WHEN;
static const int KW_EITHER;
static const int KW_PROBLEM;
static const int KW_FORDOMAIN;
static const int KW_OBJECTS;
static const int KW_INIT;
static const int KW_GOAL;
static const int KW_LENGTH;
static const int KW_SERIAL;
static const int KW_PARALLEL;
static const int KW_METRIC;
static const int KW_MINIMIZE;
static const int KW_MAXIMIZE;
static const int KW_DURATION_VAR;
static const int KW_TOTAL_TIME;
static const int KW_INCREASE;
static const int KW_DECREASE;
static const int KW_SCALE_UP;
static const int KW_SCALE_DOWN;
static const int KW_ASSIGN;
static const int KW_TAG;
static const int KW_NAME;
static const int KW_VARS;
static const int KW_SET_CONSTRAINT;
static const int KW_SETOF;
static const int KW_AT_LEAST_N;
static const int KW_AT_MOST_N;
static const int KW_EXACTLY_N;
static const int KW_CONTEXT;
static const int KW_FORMULA;
static const int KW_IRRELEVANT;
static const int KW_PLAN;
static const int KW_HEURISTIC;
static const int KW_OPT;
static const int KW_INF;
static const int KW_FACT;
static const int KW_SET;
static const int KW_EXPANSION;
static const int KW_TASKS;
static const int KW_PREFERENCE;
static const int KW_VIOLATED;
static const int KW_WITHIN;
static const int KW_ASSOC;
static const int KW_CONSTRAINTS;
static const int KW_ALWAYS;
static const int KW_SOMETIME;
static const int KW_AT_MOST_ONCE;
static const int KW_SOMETIME_BEFORE;
static const int KW_SOMETIME_AFTER;
static const int KW_ALWAYS_WITHIN;
static const int KW_IFF;
static const int KW_FALSE;
static const int KW_TRUE;
static const int KW_NUMBER;
static const int KW_UNDEFINED;
static const int KW_MINOP;
static const int KW_MAXOP;
static const int KW_DERIVED;


#line 280 "/usr/local/lib/bison.cc"
 /* decl const */
#else
enum YY_PDDL_Parser_ENUM_TOKEN { YY_PDDL_Parser_NULL_TOKEN=0

/* #line 283 "/usr/local/lib/bison.cc" */
#line 619 "grammar.cc"
	,TK_OPEN=258
	,TK_CLOSE=259
	,TK_OPEN_SQ=260
	,TK_CLOSE_SQ=261
	,TK_GREATER=262
	,TK_GREATEQ=263
	,TK_LESS=264
	,TK_LESSEQ=265
	,TK_COLON=266
	,TK_HASHT=267
	,TK_EQ=268
	,TK_HYPHEN=269
	,TK_PLUS=270
	,TK_MUL=271
	,TK_DIV=272
	,TK_UMINUS=273
	,TK_NEW_SYMBOL=274
	,TK_OBJ_SYMBOL=275
	,TK_TYPE_SYMBOL=276
	,TK_PRED_SYMBOL=277
	,TK_OBJFUN_SYMBOL=278
	,TK_FUN_SYMBOL=279
	,TK_VAR_SYMBOL=280
	,TK_ACTION_SYMBOL=281
	,TK_MISC_SYMBOL=282
	,TK_KEYWORD=283
	,TK_NEW_VAR_SYMBOL=284
	,TK_PREFERENCE_SYMBOL=285
	,TK_SET_SYMBOL=286
	,TK_FLOAT=287
	,TK_INT=288
	,TK_STRING=289
	,KW_REQS=290
	,KW_CONSTANTS=291
	,KW_PREDS=292
	,KW_FUNS=293
	,KW_TYPES=294
	,KW_DEFINE=295
	,KW_DOMAIN=296
	,KW_ACTION=297
	,KW_PROCESS=298
	,KW_EVENT=299
	,KW_ARGS=300
	,KW_PRE=301
	,KW_COND=302
	,KW_AT_START=303
	,KW_AT_END=304
	,KW_OVER_ALL=305
	,KW_EFFECT=306
	,KW_INVARIANT=307
	,KW_DURATION=308
	,KW_AND=309
	,KW_OR=310
	,KW_EXISTS=311
	,KW_FORALL=312
	,KW_IMPLY=313
	,KW_NOT=314
	,KW_WHEN=315
	,KW_EITHER=316
	,KW_PROBLEM=317
	,KW_FORDOMAIN=318
	,KW_OBJECTS=319
	,KW_INIT=320
	,KW_GOAL=321
	,KW_LENGTH=322
	,KW_SERIAL=323
	,KW_PARALLEL=324
	,KW_METRIC=325
	,KW_MINIMIZE=326
	,KW_MAXIMIZE=327
	,KW_DURATION_VAR=328
	,KW_TOTAL_TIME=329
	,KW_INCREASE=330
	,KW_DECREASE=331
	,KW_SCALE_UP=332
	,KW_SCALE_DOWN=333
	,KW_ASSIGN=334
	,KW_TAG=335
	,KW_NAME=336
	,KW_VARS=337
	,KW_SET_CONSTRAINT=338
	,KW_SETOF=339
	,KW_AT_LEAST_N=340
	,KW_AT_MOST_N=341
	,KW_EXACTLY_N=342
	,KW_CONTEXT=343
	,KW_FORMULA=344
	,KW_IRRELEVANT=345
	,KW_PLAN=346
	,KW_HEURISTIC=347
	,KW_OPT=348
	,KW_INF=349
	,KW_FACT=350
	,KW_SET=351
	,KW_EXPANSION=352
	,KW_TASKS=353
	,KW_PREFERENCE=354
	,KW_VIOLATED=355
	,KW_WITHIN=356
	,KW_ASSOC=357
	,KW_CONSTRAINTS=358
	,KW_ALWAYS=359
	,KW_SOMETIME=360
	,KW_AT_MOST_ONCE=361
	,KW_SOMETIME_BEFORE=362
	,KW_SOMETIME_AFTER=363
	,KW_ALWAYS_WITHIN=364
	,KW_IFF=365
	,KW_FALSE=366
	,KW_TRUE=367
	,KW_NUMBER=368
	,KW_UNDEFINED=369
	,KW_MINOP=370
	,KW_MAXOP=371
	,KW_DERIVED=372


#line 283 "/usr/local/lib/bison.cc"
 /* enum token */
     }; /* end of enum declaration */
#endif
public:
 int YY_PDDL_Parser_PARSE (YY_PDDL_Parser_PARSE_PARAM);
 virtual void YY_PDDL_Parser_ERROR(char *msg) YY_PDDL_Parser_ERROR_BODY;
#ifdef YY_PDDL_Parser_PURE
#ifdef YY_PDDL_Parser_LSP_NEEDED
 virtual int  YY_PDDL_Parser_LEX (YY_PDDL_Parser_STYPE *YY_PDDL_Parser_LVAL,YY_PDDL_Parser_LTYPE *YY_PDDL_Parser_LLOC) YY_PDDL_Parser_LEX_BODY;
#else
 virtual int  YY_PDDL_Parser_LEX (YY_PDDL_Parser_STYPE *YY_PDDL_Parser_LVAL) YY_PDDL_Parser_LEX_BODY;
#endif
#else
 virtual int YY_PDDL_Parser_LEX() YY_PDDL_Parser_LEX_BODY;
 YY_PDDL_Parser_STYPE YY_PDDL_Parser_LVAL;
#ifdef YY_PDDL_Parser_LSP_NEEDED
 YY_PDDL_Parser_LTYPE YY_PDDL_Parser_LLOC;
#endif
 int   YY_PDDL_Parser_NERRS;
 int    YY_PDDL_Parser_CHAR;
#endif
#if YY_PDDL_Parser_DEBUG != 0
 int YY_PDDL_Parser_DEBUG_FLAG;   /*  nonzero means print parse trace     */
#endif
public:
 YY_PDDL_Parser_CLASS(YY_PDDL_Parser_CONSTRUCTOR_PARAM);
public:
 YY_PDDL_Parser_MEMBERS 
};
/* other declare folow */
#if YY_PDDL_Parser_USE_CONST_TOKEN != 0

/* #line 314 "/usr/local/lib/bison.cc" */
#line 771 "grammar.cc"
const int YY_PDDL_Parser_CLASS::TK_OPEN=258;
const int YY_PDDL_Parser_CLASS::TK_CLOSE=259;
const int YY_PDDL_Parser_CLASS::TK_OPEN_SQ=260;
const int YY_PDDL_Parser_CLASS::TK_CLOSE_SQ=261;
const int YY_PDDL_Parser_CLASS::TK_GREATER=262;
const int YY_PDDL_Parser_CLASS::TK_GREATEQ=263;
const int YY_PDDL_Parser_CLASS::TK_LESS=264;
const int YY_PDDL_Parser_CLASS::TK_LESSEQ=265;
const int YY_PDDL_Parser_CLASS::TK_COLON=266;
const int YY_PDDL_Parser_CLASS::TK_HASHT=267;
const int YY_PDDL_Parser_CLASS::TK_EQ=268;
const int YY_PDDL_Parser_CLASS::TK_HYPHEN=269;
const int YY_PDDL_Parser_CLASS::TK_PLUS=270;
const int YY_PDDL_Parser_CLASS::TK_MUL=271;
const int YY_PDDL_Parser_CLASS::TK_DIV=272;
const int YY_PDDL_Parser_CLASS::TK_UMINUS=273;
const int YY_PDDL_Parser_CLASS::TK_NEW_SYMBOL=274;
const int YY_PDDL_Parser_CLASS::TK_OBJ_SYMBOL=275;
const int YY_PDDL_Parser_CLASS::TK_TYPE_SYMBOL=276;
const int YY_PDDL_Parser_CLASS::TK_PRED_SYMBOL=277;
const int YY_PDDL_Parser_CLASS::TK_OBJFUN_SYMBOL=278;
const int YY_PDDL_Parser_CLASS::TK_FUN_SYMBOL=279;
const int YY_PDDL_Parser_CLASS::TK_VAR_SYMBOL=280;
const int YY_PDDL_Parser_CLASS::TK_ACTION_SYMBOL=281;
const int YY_PDDL_Parser_CLASS::TK_MISC_SYMBOL=282;
const int YY_PDDL_Parser_CLASS::TK_KEYWORD=283;
const int YY_PDDL_Parser_CLASS::TK_NEW_VAR_SYMBOL=284;
const int YY_PDDL_Parser_CLASS::TK_PREFERENCE_SYMBOL=285;
const int YY_PDDL_Parser_CLASS::TK_SET_SYMBOL=286;
const int YY_PDDL_Parser_CLASS::TK_FLOAT=287;
const int YY_PDDL_Parser_CLASS::TK_INT=288;
const int YY_PDDL_Parser_CLASS::TK_STRING=289;
const int YY_PDDL_Parser_CLASS::KW_REQS=290;
const int YY_PDDL_Parser_CLASS::KW_CONSTANTS=291;
const int YY_PDDL_Parser_CLASS::KW_PREDS=292;
const int YY_PDDL_Parser_CLASS::KW_FUNS=293;
const int YY_PDDL_Parser_CLASS::KW_TYPES=294;
const int YY_PDDL_Parser_CLASS::KW_DEFINE=295;
const int YY_PDDL_Parser_CLASS::KW_DOMAIN=296;
const int YY_PDDL_Parser_CLASS::KW_ACTION=297;
const int YY_PDDL_Parser_CLASS::KW_PROCESS=298;
const int YY_PDDL_Parser_CLASS::KW_EVENT=299;
const int YY_PDDL_Parser_CLASS::KW_ARGS=300;
const int YY_PDDL_Parser_CLASS::KW_PRE=301;
const int YY_PDDL_Parser_CLASS::KW_COND=302;
const int YY_PDDL_Parser_CLASS::KW_AT_START=303;
const int YY_PDDL_Parser_CLASS::KW_AT_END=304;
const int YY_PDDL_Parser_CLASS::KW_OVER_ALL=305;
const int YY_PDDL_Parser_CLASS::KW_EFFECT=306;
const int YY_PDDL_Parser_CLASS::KW_INVARIANT=307;
const int YY_PDDL_Parser_CLASS::KW_DURATION=308;
const int YY_PDDL_Parser_CLASS::KW_AND=309;
const int YY_PDDL_Parser_CLASS::KW_OR=310;
const int YY_PDDL_Parser_CLASS::KW_EXISTS=311;
const int YY_PDDL_Parser_CLASS::KW_FORALL=312;
const int YY_PDDL_Parser_CLASS::KW_IMPLY=313;
const int YY_PDDL_Parser_CLASS::KW_NOT=314;
const int YY_PDDL_Parser_CLASS::KW_WHEN=315;
const int YY_PDDL_Parser_CLASS::KW_EITHER=316;
const int YY_PDDL_Parser_CLASS::KW_PROBLEM=317;
const int YY_PDDL_Parser_CLASS::KW_FORDOMAIN=318;
const int YY_PDDL_Parser_CLASS::KW_OBJECTS=319;
const int YY_PDDL_Parser_CLASS::KW_INIT=320;
const int YY_PDDL_Parser_CLASS::KW_GOAL=321;
const int YY_PDDL_Parser_CLASS::KW_LENGTH=322;
const int YY_PDDL_Parser_CLASS::KW_SERIAL=323;
const int YY_PDDL_Parser_CLASS::KW_PARALLEL=324;
const int YY_PDDL_Parser_CLASS::KW_METRIC=325;
const int YY_PDDL_Parser_CLASS::KW_MINIMIZE=326;
const int YY_PDDL_Parser_CLASS::KW_MAXIMIZE=327;
const int YY_PDDL_Parser_CLASS::KW_DURATION_VAR=328;
const int YY_PDDL_Parser_CLASS::KW_TOTAL_TIME=329;
const int YY_PDDL_Parser_CLASS::KW_INCREASE=330;
const int YY_PDDL_Parser_CLASS::KW_DECREASE=331;
const int YY_PDDL_Parser_CLASS::KW_SCALE_UP=332;
const int YY_PDDL_Parser_CLASS::KW_SCALE_DOWN=333;
const int YY_PDDL_Parser_CLASS::KW_ASSIGN=334;
const int YY_PDDL_Parser_CLASS::KW_TAG=335;
const int YY_PDDL_Parser_CLASS::KW_NAME=336;
const int YY_PDDL_Parser_CLASS::KW_VARS=337;
const int YY_PDDL_Parser_CLASS::KW_SET_CONSTRAINT=338;
const int YY_PDDL_Parser_CLASS::KW_SETOF=339;
const int YY_PDDL_Parser_CLASS::KW_AT_LEAST_N=340;
const int YY_PDDL_Parser_CLASS::KW_AT_MOST_N=341;
const int YY_PDDL_Parser_CLASS::KW_EXACTLY_N=342;
const int YY_PDDL_Parser_CLASS::KW_CONTEXT=343;
const int YY_PDDL_Parser_CLASS::KW_FORMULA=344;
const int YY_PDDL_Parser_CLASS::KW_IRRELEVANT=345;
const int YY_PDDL_Parser_CLASS::KW_PLAN=346;
const int YY_PDDL_Parser_CLASS::KW_HEURISTIC=347;
const int YY_PDDL_Parser_CLASS::KW_OPT=348;
const int YY_PDDL_Parser_CLASS::KW_INF=349;
const int YY_PDDL_Parser_CLASS::KW_FACT=350;
const int YY_PDDL_Parser_CLASS::KW_SET=351;
const int YY_PDDL_Parser_CLASS::KW_EXPANSION=352;
const int YY_PDDL_Parser_CLASS::KW_TASKS=353;
const int YY_PDDL_Parser_CLASS::KW_PREFERENCE=354;
const int YY_PDDL_Parser_CLASS::KW_VIOLATED=355;
const int YY_PDDL_Parser_CLASS::KW_WITHIN=356;
const int YY_PDDL_Parser_CLASS::KW_ASSOC=357;
const int YY_PDDL_Parser_CLASS::KW_CONSTRAINTS=358;
const int YY_PDDL_Parser_CLASS::KW_ALWAYS=359;
const int YY_PDDL_Parser_CLASS::KW_SOMETIME=360;
const int YY_PDDL_Parser_CLASS::KW_AT_MOST_ONCE=361;
const int YY_PDDL_Parser_CLASS::KW_SOMETIME_BEFORE=362;
const int YY_PDDL_Parser_CLASS::KW_SOMETIME_AFTER=363;
const int YY_PDDL_Parser_CLASS::KW_ALWAYS_WITHIN=364;
const int YY_PDDL_Parser_CLASS::KW_IFF=365;
const int YY_PDDL_Parser_CLASS::KW_FALSE=366;
const int YY_PDDL_Parser_CLASS::KW_TRUE=367;
const int YY_PDDL_Parser_CLASS::KW_NUMBER=368;
const int YY_PDDL_Parser_CLASS::KW_UNDEFINED=369;
const int YY_PDDL_Parser_CLASS::KW_MINOP=370;
const int YY_PDDL_Parser_CLASS::KW_MAXOP=371;
const int YY_PDDL_Parser_CLASS::KW_DERIVED=372;


#line 314 "/usr/local/lib/bison.cc"
 /* const YY_PDDL_Parser_CLASS::token */
#endif
/*apres const  */
YY_PDDL_Parser_CLASS::YY_PDDL_Parser_CLASS(YY_PDDL_Parser_CONSTRUCTOR_PARAM) YY_PDDL_Parser_CONSTRUCTOR_INIT
{
#if YY_PDDL_Parser_DEBUG != 0
YY_PDDL_Parser_DEBUG_FLAG=0;
#endif
YY_PDDL_Parser_CONSTRUCTOR_CODE;
};
#endif

/* #line 325 "/usr/local/lib/bison.cc" */
#line 903 "grammar.cc"


#define	YYFINAL		1256
#define	YYFLAG		-32768
#define	YYNTBASE	118

#define YYTRANSLATE(x) ((unsigned)(x) <= 372 ? yytranslate[x] : 390)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     1,     2,     3,     4,     5,
     6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
    16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
    26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
    36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
    46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
    56,    57,    58,    59,    60,    61,    62,    63,    64,    65,
    66,    67,    68,    69,    70,    71,    72,    73,    74,    75,
    76,    77,    78,    79,    80,    81,    82,    83,    84,    85,
    86,    87,    88,    89,    90,    91,    92,    93,    94,    95,
    96,    97,    98,    99,   100,   101,   102,   103,   104,   105,
   106,   107,   108,   109,   110,   111,   112,   113,   114,   115,
   116,   117
};

#if YY_PDDL_Parser_DEBUG != 0
static const short yyprhs[] = {     0,
     0,     3,     6,     9,    12,    15,    16,    22,    25,    28,
    31,    34,    37,    40,    43,    46,    47,    52,    54,    56,
    58,    60,    62,    64,    66,    68,    70,    72,    74,    79,
    82,    85,    88,    89,    94,    97,    98,    99,   105,   110,
   111,   120,   123,   124,   127,   130,   132,   134,   137,   139,
   140,   146,   147,   153,   154,   160,   161,   170,   172,   173,
   176,   178,   179,   185,   186,   192,   197,   202,   205,   206,
   209,   212,   213,   214,   220,   221,   227,   232,   235,   236,
   239,   241,   244,   246,   248,   250,   252,   254,   255,   262,
   263,   270,   274,   278,   282,   286,   290,   293,   297,   298,
   299,   305,   306,   312,   314,   316,   321,   324,   327,   329,
   331,   333,   335,   337,   339,   341,   343,   345,   347,   348,
   354,   360,   361,   370,   379,   380,   389,   398,   399,   411,
   423,   426,   427,   430,   431,   433,   435,   437,   438,   444,
   446,   448,   454,   463,   465,   467,   469,   471,   473,   478,
   484,   490,   496,   502,   508,   514,   518,   520,   522,   526,
   531,   533,   535,   538,   540,   543,   548,   550,   553,   556,
   557,   558,   559,   571,   572,   573,   583,   584,   591,   592,
   593,   608,   609,   610,   623,   624,   634,   640,   642,   647,
   649,   652,   654,   659,   661,   662,   668,   669,   678,   684,
   693,   694,   700,   701,   702,   703,   713,   722,   727,   729,
   732,   734,   735,   741,   747,   748,   757,   766,   768,   773,
   777,   780,   781,   782,   783,   793,   794,   801,   802,   803,
   816,   817,   826,   828,   830,   836,   838,   843,   845,   847,
   852,   855,   857,   859,   861,   862,   868,   869,   879,   880,
   889,   891,   893,   894,   903,   904,   916,   917,   927,   928,
   938,   939,   949,   950,   963,   964,   977,   978,   991,   993,
   998,  1000,  1003,  1005,  1007,  1009,  1015,  1017,  1023,  1025,
  1027,  1033,  1035,  1037,  1038,  1049,  1052,  1053,  1054,  1060,
  1061,  1062,  1073,  1078,  1079,  1088,  1091,  1092,  1095,  1098,
  1101,  1104,  1106,  1108,  1110,  1112,  1118,  1123,  1129,  1132,
  1135,  1138,  1141,  1144,  1147,  1150,  1153,  1156,  1157,  1162,
  1165,  1168,  1171,  1172,  1173,  1179,  1180,  1190,  1191,  1204,
  1205,  1215,  1216,  1226,  1232,  1237,  1245,  1250,  1258,  1261,
  1262,  1264,  1270,  1271,  1277,  1278,  1291,  1293,  1294,  1300,
  1301,  1307,  1313,  1314,  1315,  1326,  1327,  1336,  1337,  1346,
  1349,  1350,  1353,  1354,  1355,  1361,  1362,  1371,  1377,  1386,
  1388,  1390,  1395,  1400,  1405,  1411,  1417,  1423,  1430,  1436,
  1438,  1440,  1442,  1448,  1454,  1456,  1460,  1462,  1464,  1465,
  1471,  1472,  1481,  1485,  1487,  1489,  1490,  1496,  1502,  1507,
  1512,  1517,  1523,  1529,  1530,  1539,  1540,  1549,  1552,  1554,
  1555,  1561,  1563,  1565,  1566,  1574,  1575,  1583,  1584,  1590,
  1591,  1600,  1604,  1607,  1610,  1613,  1616,  1617,  1620,  1623,
  1626,  1627,  1633,  1636,  1642,  1647,  1649,  1650,  1654,  1661,
  1665,  1672,  1676,  1683,  1684,  1687,  1693,  1694,  1697,  1703,
  1706,  1712,  1713,  1716,  1718,  1720,  1722,  1724,  1726,  1728,
  1729,  1735,  1736,  1745,  1746,  1755,  1761,  1762,  1771,  1772,
  1784,  1785,  1797,  1798,  1810,  1819,  1820,  1829,  1830,  1842,
  1843,  1855,  1856,  1871,  1876,  1881,  1883,  1885,  1887,  1890,
  1893,  1894,  1895,  1901,  1902,  1911,  1912,  1913,  1914,  1930,
  1931,  1937,  1940,  1943,  1946,  1947,  1950,  1951,  1959,  1961,
  1963,  1966,  1968,  1969,  1978,  1982,  1983,  1986,  1988,  1989,
  1995,  2000,  2003,  2004,  2005,  2011,  2014,  2016,  2019,  2021,
  2023,  2025,  2026,  2032,  2033,  2042,  2043,  2052,  2055,  2058,
  2059,  2062,  2065,  2066,  2067,  2073,  2075,  2076,  2084,  2085,
  2091,  2092
};

static const short yyrhs[] = {   118,
   119,     0,   118,   260,     0,   118,   354,     0,   118,   360,
     0,   118,   368,     0,     0,     3,    40,   121,   120,     4,
     0,   120,   124,     0,   120,   143,     0,   120,   147,     0,
   120,   126,     0,   120,   134,     0,   120,   152,     0,   120,
   379,     0,   120,   273,     0,     0,     3,    41,   122,     4,
     0,    19,     0,    20,     0,    21,     0,    22,     0,    24,
     0,    25,     0,    26,     0,    27,     0,    30,     0,    19,
     0,    26,     0,     3,    35,   125,     4,     0,   125,    39,
     0,   125,   103,     0,   125,    28,     0,     0,     3,    37,
   127,     4,     0,   127,   128,     0,     0,     0,     3,    19,
   129,   130,     4,     0,   130,   132,    14,    21,     0,     0,
   130,   132,    14,     3,    61,   131,   133,     4,     0,   130,
   132,     0,     0,   132,    29,     0,   132,    25,     0,    29,
     0,    25,     0,   133,    21,     0,    21,     0,     0,     3,
    38,   135,   136,     4,     0,     0,   140,    14,   113,   137,
   136,     0,     0,   140,    14,    21,   138,   136,     0,     0,
   140,    14,     3,    61,   133,     4,   139,   136,     0,   140,
     0,     0,   141,   140,     0,   141,     0,     0,     3,    19,
   142,   130,     4,     0,     0,     3,    39,   144,   145,     4,
     0,   145,   146,    14,    21,     0,   145,   146,    14,    19,
     0,   145,   146,     0,     0,   146,    21,     0,   146,    19,
     0,     0,     0,     3,    36,   148,   150,     4,     0,     0,
     3,    64,   149,   150,     4,     0,   150,   151,    14,    21,
     0,   150,   151,     0,     0,   151,    19,     0,    19,     0,
   151,    20,     0,    20,     0,   153,     0,   254,     0,   296,
     0,   305,     0,     0,     3,    42,   123,   154,   155,     4,
     0,     0,   155,    45,     3,   156,   130,     4,     0,   155,
    96,   157,     0,   155,    51,   212,     0,   155,    46,   160,
     0,   155,    47,   160,     0,   155,    53,   241,     0,   155,
   249,     0,   155,   102,    34,     0,     0,     0,     3,    19,
   158,   171,     4,     0,     0,     3,    31,   159,   171,     4,
     0,   161,     0,   163,     0,     3,    54,   162,     4,     0,
     3,     4,     0,   162,   163,     0,   163,     0,   165,     0,
   168,     0,   184,     0,   202,     0,   204,     0,   177,     0,
    48,     0,    49,     0,    50,     0,     0,     3,    22,   166,
   172,     4,     0,     3,    13,   176,   176,     4,     0,     0,
     3,   164,     3,    22,   167,   172,     4,     4,     0,     3,
   164,     3,    13,   176,   176,     4,     4,     0,     0,     3,
    59,     3,    22,   169,   172,     4,     4,     0,     3,    59,
     3,    13,   176,   176,     4,     4,     0,     0,     3,   164,
     3,    59,     3,    22,   170,   172,     4,     4,     4,     0,
     3,   164,     3,    59,     3,    13,   176,   176,     4,     4,
     4,     0,   171,   173,     0,     0,   172,   176,     0,     0,
    25,     0,    20,     0,    19,     0,     0,     3,    23,   175,
   172,     4,     0,   173,     0,   174,     0,     3,   178,   179,
   179,     4,     0,     3,   164,     3,   178,   179,   179,     4,
     4,     0,     7,     0,     8,     0,     9,     0,    10,     0,
    13,     0,     3,    14,   179,     4,     0,     3,    15,   179,
   180,     4,     0,     3,    14,   179,   179,     4,     0,     3,
    16,   179,   181,     4,     0,     3,    17,   179,   179,     4,
     0,     3,   115,   179,   179,     4,     0,     3,   116,   179,
   179,     4,     0,   179,    17,   179,     0,    33,     0,    32,
     0,     3,    74,     4,     0,     3,   100,    30,     4,     0,
   182,     0,   179,     0,   179,   180,     0,   179,     0,   179,
   181,     0,     3,    24,   183,     4,     0,    24,     0,    25,
   183,     0,    20,   183,     0,     0,     0,     0,     3,    84,
   185,   322,   325,     3,    22,   186,   172,     4,     4,     0,
     0,     0,     3,    57,     3,   187,   130,     4,   188,   195,
     4,     0,     0,     3,    58,   189,   198,   196,     4,     0,
     0,     0,     3,   164,     3,    84,   190,   322,   325,     3,
    22,   191,   172,     4,     4,     4,     0,     0,     0,     3,
   164,     3,    57,     3,   192,   130,     4,   193,   195,     4,
     4,     0,     0,     3,   164,     3,    58,   194,   198,   196,
     4,     4,     0,     3,    58,   198,   196,     4,     0,   196,
     0,     3,    54,   197,     4,     0,   199,     0,   199,   197,
     0,   199,     0,     3,    54,   327,     4,     0,   328,     0,
     0,     3,    22,   200,   172,     4,     0,     0,     3,    59,
     3,    22,   201,   172,     4,     4,     0,     3,    13,   176,
   176,     4,     0,     3,    59,     3,    13,   176,   176,     4,
     4,     0,     0,     3,    55,   203,   208,     4,     0,     0,
     0,     0,     3,    56,     3,   205,   130,     4,   206,   207,
     4,     0,     3,    54,   327,     3,    55,   208,     4,     4,
     0,     3,    55,   208,     4,     0,   209,     0,   209,   208,
     0,   209,     0,     0,     3,    22,   210,   172,     4,     0,
     3,    13,   176,   176,     4,     0,     0,     3,    59,     3,
    22,   211,   172,     4,     4,     0,     3,    59,     3,    13,
   176,   176,     4,     4,     0,   214,     0,     3,    54,   213,
     4,     0,     3,    54,     4,     0,   213,   214,     0,     0,
     0,     0,     3,    57,     3,   215,   130,     4,   216,   221,
     4,     0,     0,     3,    60,   217,   223,   222,     4,     0,
     0,     0,     3,   164,     3,    57,     3,   218,   130,     4,
   219,   221,     4,     4,     0,     0,     3,   164,     3,    60,
   220,   223,   222,     4,     0,   225,     0,   234,     0,     3,
    60,   223,   222,     4,     0,   222,     0,     3,    54,   224,
     4,     0,   225,     0,   328,     0,     3,    54,   327,     4,
     0,   225,   224,     0,   225,     0,   226,     0,   231,     0,
     0,     3,    22,   227,   172,     4,     0,     0,     3,    79,
     3,    23,   228,   172,     4,   230,     4,     0,     0,     3,
   164,     3,    22,   229,   172,     4,     4,     0,   176,     0,
   114,     0,     0,     3,    59,     3,    22,   232,   172,     4,
     4,     0,     0,     3,   164,     3,    59,     3,    22,   233,
   172,     4,     4,     4,     0,     0,     3,    75,     3,    24,
   235,   172,     4,   179,     4,     0,     0,     3,    76,     3,
    24,   236,   172,     4,   179,     4,     0,     0,     3,    79,
     3,    24,   237,   172,     4,   179,     4,     0,     0,     3,
   164,     3,    75,     3,    24,   238,   172,     4,   179,     4,
     4,     0,     0,     3,   164,     3,    76,     3,    24,   239,
   172,     4,   179,     4,     4,     0,     0,     3,   164,     3,
    79,     3,    24,   240,   172,     4,   179,     4,     4,     0,
   243,     0,     3,    54,   242,     4,     0,   243,     0,   243,
   242,     0,   244,     0,   245,     0,   247,     0,     3,    13,
    73,   179,     4,     0,   179,     0,     3,   246,    73,   179,
     4,     0,     9,     0,    10,     0,     3,   248,    73,   179,
     4,     0,     7,     0,     8,     0,     0,     3,    97,   250,
   323,   324,    98,     3,   251,     4,     4,     0,   252,   251,
     0,     0,     0,     3,    26,   253,   172,     4,     0,     0,
     0,     3,   117,     3,    22,   255,   257,     4,   256,   300,
     4,     0,   257,   259,    14,    21,     0,     0,   257,   259,
    14,     3,    61,   258,   133,     4,     0,   257,   259,     0,
     0,   259,    29,     0,   259,    25,     0,   259,    20,     0,
   259,    19,     0,    29,     0,    25,     0,    20,     0,    19,
     0,     3,    40,   261,   262,     4,     0,     3,    62,   122,
     4,     0,   262,     3,    63,   122,     4,     0,   262,   124,
     0,   262,   147,     0,   262,   263,     0,   262,   273,     0,
   262,   291,     0,   262,   294,     0,   262,   312,     0,   262,
   305,     0,   262,   379,     0,     0,     3,    65,   264,     4,
     0,   264,   271,     0,   264,   269,     0,   264,   265,     0,
     0,     0,     3,    22,   266,   172,     4,     0,     0,     3,
   122,   295,     3,    22,   267,   172,     4,     4,     0,     0,
     3,   122,   295,     3,    59,     3,    22,   268,   172,     4,
     4,     4,     0,     0,     3,    13,     3,    23,   270,   172,
     4,    20,     4,     0,     0,     3,    13,     3,    24,   272,
   172,     4,   295,     4,     0,     3,    13,    24,   295,     4,
     0,     3,    66,   275,     4,     0,     3,    66,     3,    54,
   274,     4,     4,     0,     3,   103,   275,     4,     0,     3,
   103,     3,    54,   274,     4,     4,     0,   274,   275,     0,
     0,   290,     0,     3,    99,    19,   278,     4,     0,     0,
     3,    55,   276,   285,     4,     0,     0,     3,    57,     3,
   277,   130,     4,     3,    99,    19,   278,     4,     4,     0,
   290,     0,     0,     3,    54,   279,   285,     4,     0,     0,
     3,    55,   280,   285,     4,     0,     3,    58,   287,   287,
     4,     0,     0,     0,     3,    58,     3,    54,   281,   286,
     4,   282,   287,     4,     0,     0,     3,    57,     3,   283,
   130,     4,   278,     4,     0,     0,     3,    56,     3,   284,
   130,     4,   278,     4,     0,   290,   285,     0,     0,   287,
   286,     0,     0,     0,     3,    22,   288,   172,     4,     0,
     0,     3,    59,     3,    22,   289,   172,     4,     4,     0,
     3,    13,   176,   176,     4,     0,     3,    59,     3,    13,
   176,   176,     4,     4,     0,   287,     0,   177,     0,     3,
   104,   278,     4,     0,     3,   105,   278,     4,     0,     3,
   106,   278,     4,     0,     3,   107,   278,   278,     4,     0,
     3,   108,   278,   278,     4,     0,     3,   101,   295,   278,
     4,     0,     3,   109,   295,   278,   278,     4,     0,     3,
    70,   292,   293,     4,     0,    71,     0,    72,     0,   179,
     0,     3,    67,    68,    33,     4,     0,     3,    67,    69,
    33,     4,     0,    33,     0,    33,    17,    33,     0,    32,
     0,    94,     0,     0,     3,    52,   297,   316,   298,     0,
     0,    83,   299,     3,   345,    33,   346,     4,     4,     0,
    89,   300,     4,     0,   111,     0,   112,     0,     0,     3,
    22,   301,   171,     4,     0,     3,    13,   173,   173,     4,
     0,     3,    59,   300,     4,     0,     3,    54,   304,     4,
     0,     3,    55,   304,     4,     0,     3,    58,   300,   300,
     4,     0,     3,   110,   300,   300,     4,     0,     0,     3,
    57,     3,   302,   130,     4,   300,     4,     0,     0,     3,
    56,     3,   303,   130,     4,   300,     4,     0,   304,   300,
     0,   300,     0,     0,     3,    90,   306,   316,   307,     0,
   308,     0,   310,     0,     0,    42,     3,    26,   309,   172,
     4,     4,     0,     0,    95,     3,    22,   311,   172,     4,
     4,     0,     0,     3,    52,   313,   316,   314,     0,     0,
    83,   315,     3,   345,    33,   346,     4,     4,     0,    89,
   300,     4,     0,   316,   317,     0,   316,   318,     0,   316,
   319,     0,   316,   321,     0,     0,    80,   122,     0,    81,
    19,     0,    81,    27,     0,     0,    82,     3,   320,   130,
     4,     0,    88,   328,     0,    88,     3,    54,   327,     4,
     0,    82,     3,   130,     4,     0,   322,     0,     0,    88,
   328,   326,     0,    88,     3,    54,   327,     4,   326,     0,
    46,   328,   325,     0,    46,     3,    54,   327,     4,   325,
     0,    47,   328,   325,     0,    47,     3,    54,   327,     4,
   325,     0,     0,    88,   328,     0,    88,     3,    54,   327,
     4,     0,     0,    46,   328,     0,    46,     3,    54,   327,
     4,     0,    47,   328,     0,    47,     3,    54,   327,     4,
     0,     0,   327,   328,     0,   328,     0,   329,     0,   333,
     0,   338,     0,   341,     0,   344,     0,     0,     3,    22,
   330,   172,     4,     0,     0,     3,   164,     3,    22,   331,
   172,     4,     4,     0,     0,     3,    65,     3,    22,   332,
   172,     4,     4,     0,     3,    13,   176,   176,     4,     0,
     0,     3,    59,     3,    22,   334,   172,     4,     4,     0,
     0,     3,   164,     3,    59,     3,    22,   335,   172,     4,
     4,     4,     0,     0,     3,    65,     3,    59,     3,    22,
   336,   172,     4,     4,     4,     0,     0,     3,    59,     3,
    65,     3,    22,   337,   172,     4,     4,     4,     0,     3,
    59,     3,    13,   176,   176,     4,     4,     0,     0,     3,
    66,     3,    22,   339,   172,     4,     4,     0,     0,     3,
    66,     3,    59,     3,    22,   340,   172,     4,     4,     4,
     0,     0,     3,    59,     3,    66,     3,    22,   342,   172,
     4,     4,     4,     0,     0,     3,    59,     3,    66,     3,
    59,     3,    22,   343,   172,     4,     4,     4,     4,     0,
     3,    21,    25,     4,     0,     3,    21,   174,     4,     0,
    85,     0,    86,     0,    87,     0,   346,   347,     0,   346,
   350,     0,     0,     0,     3,    22,   348,   172,     4,     0,
     0,     3,    59,     3,    22,   349,   172,     4,     4,     0,
     0,     0,     0,     3,    84,    82,     3,   351,   130,     4,
   352,   325,     3,    22,   353,   172,     4,     4,     0,     0,
     3,    91,   355,   356,     4,     0,   357,   356,     0,    93,
   356,     0,   358,   356,     0,     0,    81,   122,     0,     0,
   295,    11,     3,    26,   359,   172,     4,     0,   361,     0,
   365,     0,   362,   361,     0,   362,     0,     0,   295,    11,
     3,    26,   363,   172,     4,   364,     0,     5,   295,     6,
     0,     0,   366,   365,     0,   366,     0,     0,     3,    26,
   367,   172,     4,     0,     3,    92,   369,     4,     0,   369,
   370,     0,     0,     0,     3,   371,   373,     4,   372,     0,
    93,   295,     0,   295,     0,   373,   374,     0,   374,     0,
   375,     0,   377,     0,     0,     3,    22,   376,   172,     4,
     0,     0,     3,    59,     3,    22,   378,   172,     4,     4,
     0,     0,     3,    96,   380,   381,   323,   325,   382,     4,
     0,    81,    19,     0,    81,    27,     0,     0,   382,   385,
     0,   382,   383,     0,     0,     0,     3,   122,   384,   172,
     4,     0,   122,     0,     0,     3,    84,   386,   322,   325,
   387,     4,     0,     0,     3,   122,   388,   172,     4,     0,
     0,     3,    59,     3,   122,   389,   172,     4,     4,     0
};

#endif

#if YY_PDDL_Parser_DEBUG != 0
static const short yyrline[] = { 0,
   105,   107,   108,   109,   110,   111,   114,   118,   120,   121,
   122,   123,   124,   125,   126,   127,   130,   138,   140,   141,
   142,   143,   144,   145,   146,   147,   150,   152,   167,   171,
   173,   174,   175,   180,   184,   186,   189,   194,   204,   209,
   213,   217,   221,   224,   237,   244,   256,   265,   270,   277,
   282,   285,   291,   291,   321,   321,   351,   351,   355,   358,
   360,   363,   368,   380,   385,   388,   407,   426,   434,   437,
   443,   449,   454,   459,   460,   464,   467,   488,   499,   502,
   512,   521,   528,   539,   541,   542,   543,   548,   553,   562,
   567,   571,   572,   573,   574,   575,   576,   577,   582,   585,
   594,   598,   603,   609,   611,   612,   615,   619,   621,   624,
   626,   631,   632,   633,   634,   640,   645,   649,   655,   660,
   666,   674,   678,   684,   695,   700,   706,   714,   718,   724,
   779,   786,   789,   796,   799,   807,   811,   819,   825,   845,
   850,   855,   862,   871,   876,   880,   884,   888,   894,   899,
   903,   907,   911,   915,   919,   923,   927,   931,   935,   939,
   943,   949,   954,   960,   965,   971,   976,   982,   987,   991,
   997,  1003,  1007,  1022,  1027,  1036,  1046,  1050,  1056,  1061,
  1065,  1081,  1086,  1095,  1107,  1111,  1120,  1122,  1125,  1127,
  1130,  1132,  1135,  1137,  1140,  1145,  1152,  1156,  1163,  1172,
  1209,  1214,  1222,  1223,  1228,  1237,  1249,  1251,  1252,  1255,
  1257,  1260,  1265,  1271,  1279,  1283,  1289,  1299,  1301,  1302,
  1312,  1314,  1317,  1323,  1331,  1342,  1346,  1351,  1356,  1364,
  1377,  1381,  1388,  1389,  1392,  1394,  1397,  1399,  1402,  1404,
  1407,  1409,  1412,  1414,  1417,  1423,  1433,  1437,  1464,  1469,
  1481,  1486,  1492,  1498,  1508,  1513,  1525,  1530,  1537,  1541,
  1548,  1552,  1559,  1563,  1570,  1574,  1581,  1585,  1594,  1596,
  1599,  1601,  1604,  1606,  1607,  1610,  1616,  1623,  1630,  1632,
  1635,  1642,  1644,  1647,  1657,  1680,  1682,  1685,  1690,  1701,
  1708,  1717,  1728,  1733,  1737,  1741,  1745,  1748,  1763,  1768,
  1773,  1780,  1794,  1799,  1804,  1816,  1820,  1828,  1830,  1831,
  1832,  1833,  1834,  1835,  1836,  1837,  1838,  1839,  1842,  1846,
  1848,  1849,  1850,  1853,  1858,  1869,  1873,  1886,  1890,  1905,
  1910,  1922,  1927,  1940,  1956,  1958,  1959,  1960,  1963,  1965,
  1968,  1973,  1979,  1983,  1989,  1993,  2012,  2017,  2021,  2027,
  2031,  2037,  2053,  2057,  2075,  2084,  2088,  2104,  2108,  2126,
  2134,  2137,  2145,  2148,  2153,  2160,  2164,  2171,  2179,  2188,
  2193,  2197,  2201,  2205,  2209,  2213,  2217,  2221,  2257,  2261,
  2271,  2282,  2297,  2302,  2315,  2317,  2318,  2319,  2324,  2330,
  2333,  2338,  2345,  2353,  2358,  2362,  2366,  2371,  2375,  2379,
  2384,  2389,  2393,  2397,  2401,  2416,  2420,  2437,  2442,  2449,
  2456,  2459,  2461,  2464,  2469,  2481,  2486,  2498,  2505,  2508,
  2513,  2520,  2528,  2530,  2531,  2532,  2533,  2536,  2543,  2550,
  2557,  2562,  2568,  2570,  2573,  2583,  2585,  2588,  2590,  2591,
  2592,  2593,  2594,  2595,  2598,  2600,  2601,  2604,  2606,  2607,
  2608,  2609,  2612,  2614,  2617,  2619,  2620,  2621,  2622,  2625,
  2630,  2637,  2641,  2647,  2651,  2657,  2667,  2672,  2679,  2683,
  2689,  2693,  2699,  2703,  2709,  2719,  2724,  2730,  2734,  2742,
  2747,  2753,  2757,  2765,  2772,  2779,  2784,  2788,  2794,  2796,
  2797,  2799,  2804,  2809,  2813,  2820,  2826,  2834,  2838,  2854,
  2864,  2877,  2879,  2884,  2885,  2888,  2896,  2901,  2913,  2919,
  2926,  2928,  2931,  2936,  2967,  2969,  2972,  2974,  2977,  2982,
  3014,  3018,  3020,  3023,  3028,  3035,  3041,  3047,  3049,  3052,
  3054,  3057,  3063,  3072,  3078,  3087,  3096,  3103,  3109,  3113,
  3116,  3118,  3119,  3122,  3132,  3139,  3154,  3161,  3175,  3185,
  3193,  3202
};

static const char * const yytname[] = {   "$","error","$illegal.","TK_OPEN",
"TK_CLOSE","TK_OPEN_SQ","TK_CLOSE_SQ","TK_GREATER","TK_GREATEQ","TK_LESS","TK_LESSEQ",
"TK_COLON","TK_HASHT","TK_EQ","TK_HYPHEN","TK_PLUS","TK_MUL","TK_DIV","TK_UMINUS",
"TK_NEW_SYMBOL","TK_OBJ_SYMBOL","TK_TYPE_SYMBOL","TK_PRED_SYMBOL","TK_OBJFUN_SYMBOL",
"TK_FUN_SYMBOL","TK_VAR_SYMBOL","TK_ACTION_SYMBOL","TK_MISC_SYMBOL","TK_KEYWORD",
"TK_NEW_VAR_SYMBOL","TK_PREFERENCE_SYMBOL","TK_SET_SYMBOL","TK_FLOAT","TK_INT",
"TK_STRING","KW_REQS","KW_CONSTANTS","KW_PREDS","KW_FUNS","KW_TYPES","KW_DEFINE",
"KW_DOMAIN","KW_ACTION","KW_PROCESS","KW_EVENT","KW_ARGS","KW_PRE","KW_COND",
"KW_AT_START","KW_AT_END","KW_OVER_ALL","KW_EFFECT","KW_INVARIANT","KW_DURATION",
"KW_AND","KW_OR","KW_EXISTS","KW_FORALL","KW_IMPLY","KW_NOT","KW_WHEN","KW_EITHER",
"KW_PROBLEM","KW_FORDOMAIN","KW_OBJECTS","KW_INIT","KW_GOAL","KW_LENGTH","KW_SERIAL",
"KW_PARALLEL","KW_METRIC","KW_MINIMIZE","KW_MAXIMIZE","KW_DURATION_VAR","KW_TOTAL_TIME",
"KW_INCREASE","KW_DECREASE","KW_SCALE_UP","KW_SCALE_DOWN","KW_ASSIGN","KW_TAG",
"KW_NAME","KW_VARS","KW_SET_CONSTRAINT","KW_SETOF","KW_AT_LEAST_N","KW_AT_MOST_N",
"KW_EXACTLY_N","KW_CONTEXT","KW_FORMULA","KW_IRRELEVANT","KW_PLAN","KW_HEURISTIC",
"KW_OPT","KW_INF","KW_FACT","KW_SET","KW_EXPANSION","KW_TASKS","KW_PREFERENCE",
"KW_VIOLATED","KW_WITHIN","KW_ASSOC","KW_CONSTRAINTS","KW_ALWAYS","KW_SOMETIME",
"KW_AT_MOST_ONCE","KW_SOMETIME_BEFORE","KW_SOMETIME_AFTER","KW_ALWAYS_WITHIN",
"KW_IFF","KW_FALSE","KW_TRUE","KW_NUMBER","KW_UNDEFINED","KW_MINOP","KW_MAXOP",
"KW_DERIVED","pddl_declarations","pddl_domain","domain_elements","domain_name",
"any_symbol","action_symbol","domain_requires","require_list","domain_predicates",
"predicate_list","predicate_decl","@1","typed_param_list","@2","typed_param_sym_list",
"non_empty_type_name_list","domain_functions","@3","function_list","@4","@5",
"@6","function_decl_list","function_decl","@7","domain_types","@8","typed_type_list",
"primitive_type_list","domain_constants","@9","@10","typed_constant_list","ne_constant_sym_list",
"domain_structure","action_declaration","@11","action_elements","@12","action_set_name",
"@13","@14","action_condition","empty_action_condition","action_condition_list",
"single_action_condition","timing_keyword","positive_atom_condition","@15","@16",
"negative_atom_condition","@17","@18","flat_atom_argument_list","atom_argument_list",
"flat_atom_argument","functional_atom_argument","@19","atom_argument","numeric_condition",
"relation_keyword","d_expression","d_sum","d_product","d_function","d_argument_list",
"set_condition","@20","@21","@22","@23","@24","@25","@26","@27","@28","@29",
"universal_condition_body","one_or_more_condition_atoms","quantified_condition_atom_list",
"one_or_more_context_atoms","quantified_condition_atom","@30","@31","disjunctive_condition",
"@32","disjunctive_set_condition","@33","@34","existential_condition_body","disjunction_list",
"disjunction_atom","@35","@36","action_effect","action_effect_list","single_action_effect",
"@37","@38","@39","@40","@41","@42","quantified_effect_body","one_or_more_atomic_effects",
"effect_conditions","atomic_effect_list","atomic_effect","positive_atom_effect",
"@43","@44","@45","object_assignment_value","negative_atom_effect","@46","@47",
"numeric_effect","@48","@49","@50","@51","@52","@53","action_duration","action_duration_list",
"action_duration_exp","action_exact_duration","action_min_duration","less_or_lesseq",
"action_max_duration","greater_or_greatereq","action_expansion","@54","task_list",
"task","@55","axiom_declaration","@56","@57","atom_schema_typed_arg_list","@58",
"atom_schema_arg_list","pddl_problem","problem_name","problem_elements","initial_state",
"init_elements","init_atom","@59","@60","@61","init_object_function","@62","init_function",
"@63","goal_spec","goal_spec_list","single_goal_spec","@64","@65","new_goal",
"@66","@67","@68","@69","@70","@71","new_goal_list","new_atomic_goal_list","new_atomic_goal",
"@72","@73","new_single_goal","metric_spec","metric_keyword","metric_expression",
"length_spec","numeric_value","domain_invariant","@74","domain_invariant_body",
"@75","fol_formula","@76","@77","@78","fol_formula_list","irrelevant_item","@79",
"irrelevant_item_content","irrelevant_action","@80","irrelevant_atom","@81",
"problem_invariant","@82","problem_invariant_body","@83","dkel_element_list",
"dkel_tag_spec","dkel_name_spec","dkel_vars_spec","@84","dkel_context_spec",
"req_vars_spec","opt_vars_spec","opt_context_and_precondition_spec","opt_context_spec",
"opt_precondition_spec","context_list","context_atom","positive_context_atom",
"@85","@86","@87","negative_context_atom","@88","@89","@90","@91","positive_context_goal_atom",
"@92","@93","negative_context_goal_atom","@94","@95","type_constraint_atom",
"set_constraint_type","invariant_set","invariant_atom","@96","@97","invariant_set_of_atoms",
"@98","@99","@100","pddl_plan","@101","plan_elements","plan_name","plan_step",
"@102","ipc_plan","ipc_plan_step_list","ipc_plan_step","@103","opt_step_duration",
"ipc_plan_step_seq","simple_plan_step","@104","heuristic_table","table_entry_list",
"table_entry","@105","entry_value","ne_entry_atom_list","entry_atom","pos_entry_atom",
"@106","neg_entry_atom","@107","reference_set","@108","opt_set_name","reference_list",
"reference","@109","simple_reference_set","@110","simple_reference_set_body",
"@111","@112",""
};
#endif

static const short yyr1[] = {     0,
   118,   118,   118,   118,   118,   118,   119,   120,   120,   120,
   120,   120,   120,   120,   120,   120,   121,   122,   122,   122,
   122,   122,   122,   122,   122,   122,   123,   123,   124,   125,
   125,   125,   125,   126,   127,   127,   129,   128,   130,   131,
   130,   130,   130,   132,   132,   132,   132,   133,   133,   135,
   134,   137,   136,   138,   136,   139,   136,   136,   136,   140,
   140,   142,   141,   144,   143,   145,   145,   145,   145,   146,
   146,   146,   148,   147,   149,   147,   150,   150,   150,   151,
   151,   151,   151,   152,   152,   152,   152,   154,   153,   156,
   155,   155,   155,   155,   155,   155,   155,   155,   155,   158,
   157,   159,   157,   160,   160,   160,   161,   162,   162,   163,
   163,   163,   163,   163,   163,   164,   164,   164,   166,   165,
   165,   167,   165,   165,   169,   168,   168,   170,   168,   168,
   171,   171,   172,   172,   173,   173,   173,   175,   174,   176,
   176,   177,   177,   178,   178,   178,   178,   178,   179,   179,
   179,   179,   179,   179,   179,   179,   179,   179,   179,   179,
   179,   180,   180,   181,   181,   182,   182,   183,   183,   183,
   185,   186,   184,   187,   188,   184,   189,   184,   190,   191,
   184,   192,   193,   184,   194,   184,   195,   195,   196,   196,
   197,   197,   198,   198,   200,   199,   201,   199,   199,   199,
   203,   202,   204,   205,   206,   204,   207,   207,   207,   208,
   208,   210,   209,   209,   211,   209,   209,   212,   212,   212,
   213,   213,   215,   216,   214,   217,   214,   218,   219,   214,
   220,   214,   214,   214,   221,   221,   222,   222,   223,   223,
   224,   224,   225,   225,   227,   226,   228,   226,   229,   226,
   230,   230,   232,   231,   233,   231,   235,   234,   236,   234,
   237,   234,   238,   234,   239,   234,   240,   234,   241,   241,
   242,   242,   243,   243,   243,   244,   244,   245,   246,   246,
   247,   248,   248,   250,   249,   251,   251,   253,   252,   255,
   256,   254,   257,   258,   257,   257,   257,   259,   259,   259,
   259,   259,   259,   259,   259,   260,   261,   262,   262,   262,
   262,   262,   262,   262,   262,   262,   262,   262,   263,   264,
   264,   264,   264,   266,   265,   267,   265,   268,   265,   270,
   269,   272,   271,   271,   273,   273,   273,   273,   274,   274,
   275,   275,   276,   275,   277,   275,   278,   279,   278,   280,
   278,   278,   281,   282,   278,   283,   278,   284,   278,   285,
   285,   286,   286,   288,   287,   289,   287,   287,   287,   290,
   290,   290,   290,   290,   290,   290,   290,   290,   291,   292,
   292,   293,   294,   294,   295,   295,   295,   295,   297,   296,
   299,   298,   298,   300,   300,   301,   300,   300,   300,   300,
   300,   300,   300,   302,   300,   303,   300,   304,   304,   306,
   305,   307,   307,   309,   308,   311,   310,   313,   312,   315,
   314,   314,   316,   316,   316,   316,   316,   317,   318,   318,
   320,   319,   321,   321,   322,   323,   323,   324,   324,   324,
   324,   324,   324,   324,   325,   325,   325,   326,   326,   326,
   326,   326,   327,   327,   328,   328,   328,   328,   328,   330,
   329,   331,   329,   332,   329,   329,   334,   333,   335,   333,
   336,   333,   337,   333,   333,   339,   338,   340,   338,   342,
   341,   343,   341,   344,   344,   345,   345,   345,   346,   346,
   346,   348,   347,   349,   347,   351,   352,   353,   350,   355,
   354,   356,   356,   356,   356,   357,   359,   358,   360,   360,
   361,   361,   363,   362,   364,   364,   365,   365,   367,   366,
   368,   369,   369,   371,   370,   372,   372,   373,   373,   374,
   374,   376,   375,   378,   377,   380,   379,   381,   381,   381,
   382,   382,   382,   384,   383,   383,   386,   385,   388,   387,
   389,   387
};

static const short yyr2[] = {     0,
     2,     2,     2,     2,     2,     0,     5,     2,     2,     2,
     2,     2,     2,     2,     2,     0,     4,     1,     1,     1,
     1,     1,     1,     1,     1,     1,     1,     1,     4,     2,
     2,     2,     0,     4,     2,     0,     0,     5,     4,     0,
     8,     2,     0,     2,     2,     1,     1,     2,     1,     0,
     5,     0,     5,     0,     5,     0,     8,     1,     0,     2,
     1,     0,     5,     0,     5,     4,     4,     2,     0,     2,
     2,     0,     0,     5,     0,     5,     4,     2,     0,     2,
     1,     2,     1,     1,     1,     1,     1,     0,     6,     0,
     6,     3,     3,     3,     3,     3,     2,     3,     0,     0,
     5,     0,     5,     1,     1,     4,     2,     2,     1,     1,
     1,     1,     1,     1,     1,     1,     1,     1,     0,     5,
     5,     0,     8,     8,     0,     8,     8,     0,    11,    11,
     2,     0,     2,     0,     1,     1,     1,     0,     5,     1,
     1,     5,     8,     1,     1,     1,     1,     1,     4,     5,
     5,     5,     5,     5,     5,     3,     1,     1,     3,     4,
     1,     1,     2,     1,     2,     4,     1,     2,     2,     0,
     0,     0,    11,     0,     0,     9,     0,     6,     0,     0,
    14,     0,     0,    12,     0,     9,     5,     1,     4,     1,
     2,     1,     4,     1,     0,     5,     0,     8,     5,     8,
     0,     5,     0,     0,     0,     9,     8,     4,     1,     2,
     1,     0,     5,     5,     0,     8,     8,     1,     4,     3,
     2,     0,     0,     0,     9,     0,     6,     0,     0,    12,
     0,     8,     1,     1,     5,     1,     4,     1,     1,     4,
     2,     1,     1,     1,     0,     5,     0,     9,     0,     8,
     1,     1,     0,     8,     0,    11,     0,     9,     0,     9,
     0,     9,     0,    12,     0,    12,     0,    12,     1,     4,
     1,     2,     1,     1,     1,     5,     1,     5,     1,     1,
     5,     1,     1,     0,    10,     2,     0,     0,     5,     0,
     0,    10,     4,     0,     8,     2,     0,     2,     2,     2,
     2,     1,     1,     1,     1,     5,     4,     5,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     0,     4,     2,
     2,     2,     0,     0,     5,     0,     9,     0,    12,     0,
     9,     0,     9,     5,     4,     7,     4,     7,     2,     0,
     1,     5,     0,     5,     0,    12,     1,     0,     5,     0,
     5,     5,     0,     0,    10,     0,     8,     0,     8,     2,
     0,     2,     0,     0,     5,     0,     8,     5,     8,     1,
     1,     4,     4,     4,     5,     5,     5,     6,     5,     1,
     1,     1,     5,     5,     1,     3,     1,     1,     0,     5,
     0,     8,     3,     1,     1,     0,     5,     5,     4,     4,
     4,     5,     5,     0,     8,     0,     8,     2,     1,     0,
     5,     1,     1,     0,     7,     0,     7,     0,     5,     0,
     8,     3,     2,     2,     2,     2,     0,     2,     2,     2,
     0,     5,     2,     5,     4,     1,     0,     3,     6,     3,
     6,     3,     6,     0,     2,     5,     0,     2,     5,     2,
     5,     0,     2,     1,     1,     1,     1,     1,     1,     0,
     5,     0,     8,     0,     8,     5,     0,     8,     0,    11,
     0,    11,     0,    11,     8,     0,     8,     0,    11,     0,
    11,     0,    14,     4,     4,     1,     1,     1,     2,     2,
     0,     0,     5,     0,     8,     0,     0,     0,    15,     0,
     5,     2,     2,     2,     0,     2,     0,     7,     1,     1,
     2,     1,     0,     8,     3,     0,     2,     1,     0,     5,
     4,     2,     0,     0,     5,     2,     1,     2,     1,     1,
     1,     0,     5,     0,     8,     0,     8,     2,     2,     0,
     2,     2,     0,     0,     5,     1,     0,     7,     0,     5,
     0,     8
};

static const short yydefact[] = {     6,
     0,     0,   387,   385,   388,     1,     2,     0,     3,     4,
   509,   512,   510,   518,     5,   519,     0,   500,   523,     0,
     0,   511,     0,   517,   134,     0,    16,   318,   505,     0,
   386,     0,     0,     0,     0,     0,     0,     0,   505,     0,
     0,   505,   505,   524,   521,   522,   513,     0,   520,   137,
   136,   135,   140,   141,   133,    18,    19,    20,    21,    22,
    23,    24,    25,    26,     0,     0,     0,     7,     8,    11,
    12,     9,    10,    13,    84,    85,    15,    86,    87,    14,
     0,   306,   309,   310,   311,   312,   313,   314,   316,   315,
   317,   506,   503,     0,   501,   502,   504,     0,   134,   138,
    17,   307,    33,    73,    36,    50,    64,     0,   389,    75,
     0,   410,   536,     0,     0,   418,     0,   323,     0,     0,
     0,     0,     0,   529,   530,   531,     0,   134,     0,    79,
     0,    59,    69,    27,    28,    88,   427,    79,     0,   371,
     0,   370,   341,   427,   540,     0,     0,     0,   427,     0,
     0,     0,     0,   380,   381,     0,   507,   532,     0,     0,
   528,   516,     0,    29,    32,    30,    31,     0,     0,    34,
    35,     0,     0,    58,    61,    72,    99,     0,     0,   144,
   145,   146,   147,   148,   364,   116,   117,   118,   340,   343,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,   335,     0,     0,   437,   340,   337,   290,     0,
   308,     0,   319,   322,   321,   320,     0,     0,     0,   167,
   158,   157,   382,   161,     0,   134,   134,     0,     0,   527,
   525,     0,   514,   139,    74,    81,    83,    78,    37,    62,
    51,     0,    60,    65,    68,     0,     0,     0,     0,   391,
     0,     0,   390,   423,   424,   425,   426,    76,     0,   134,
     0,   361,   345,     0,     0,     0,     0,     0,   347,     0,
     0,     0,     0,     0,     0,     0,     0,     0,   411,   412,
   413,   538,   539,     0,   436,   447,     0,   297,   420,     0,
   419,     0,   324,     0,   383,   384,     0,     0,     0,     0,
   170,     0,     0,     0,     0,     0,   379,     0,     0,   534,
   526,     0,     0,    80,    82,    43,    43,     0,    54,    52,
     0,    71,    70,     0,    89,     0,   203,   203,     0,     0,
     0,     0,    97,   428,   429,   430,   431,     0,     0,   433,
   455,   456,   457,   458,   459,     0,   394,   395,     0,     0,
     0,     0,     0,   339,     0,     0,   361,    43,     0,   366,
     0,     0,   348,   350,     0,     0,     0,   372,   373,   374,
     0,     0,     0,   148,     0,     0,     0,     0,    43,     0,
   543,     0,     0,     0,     0,     0,     0,   134,     0,     0,
     0,     0,     0,   170,   170,     0,   159,     0,     0,     0,
   156,   508,   533,   134,   515,    77,     0,     0,     0,    59,
    59,    67,    66,   284,    90,     0,    94,   104,   105,   110,
   111,   115,   112,   113,   114,    95,     0,    93,   218,   233,
   243,   244,   234,     0,   277,    96,   269,   273,   274,   275,
     0,    92,    98,    43,     0,     0,     0,   460,     0,     0,
     0,     0,     0,     0,   396,     0,     0,     0,     0,     0,
     0,     0,   393,   368,   365,   336,   344,   360,     0,     0,
   134,   342,   377,   361,   361,   358,   356,     0,     0,   375,
   376,     0,     0,   142,   414,   416,     0,     0,   445,     0,
   338,   291,   305,   304,   303,   302,   296,     0,   422,   330,
   332,     0,     0,     0,   149,     0,   162,     0,   164,     0,
     0,   169,   168,   166,   160,     0,     0,     0,    38,    47,
    46,    42,    63,    49,     0,    55,    53,   437,    43,   107,
   148,   119,   203,   201,     0,     0,   177,     0,   171,     0,
   245,   222,     0,     0,   226,     0,     0,     0,     0,   282,
   283,   279,   280,     0,     0,     0,     0,   100,   102,     0,
   486,   487,   488,     0,     0,     0,     0,   134,     0,     0,
   454,     0,     0,     0,     0,     0,   132,   409,     0,     0,
   406,   404,     0,     0,     0,     0,     0,     0,     0,     0,
    43,    43,     0,   353,     0,     0,   378,     0,   134,   134,
   435,     0,     0,   537,   546,   542,   541,     0,     0,   301,
   300,   299,   298,     0,   134,   134,   334,   325,   326,     0,
   151,   163,   150,   165,   152,   153,   154,   155,     0,     0,
    45,    44,    56,    48,   444,     0,     0,   134,     0,     0,
   109,     0,   204,   174,     0,     0,     0,     0,   134,   220,
     0,   223,     0,     0,     0,     0,     0,     0,     0,     0,
     0,   271,     0,     0,   132,   132,   432,   491,     0,   484,
   485,     0,   434,   453,     0,   467,     0,     0,   464,     0,
   476,     0,   462,     0,     0,     0,   400,   408,   401,    43,
    43,     0,   399,     0,     0,     0,     0,   349,   351,     0,
     0,   363,   352,     0,     0,     0,     0,   547,   544,     0,
     0,   293,   491,     0,     0,   134,     0,   535,     0,    39,
    59,     0,     0,     0,     0,    91,     0,     0,   106,   108,
     0,     0,   211,    43,    43,     0,     0,   194,     0,   125,
   447,   148,   122,     0,   185,     0,   179,     0,     0,   219,
   221,    43,   253,     0,     0,   239,   257,   259,   247,   261,
   249,     0,     0,   231,     0,     0,     0,     0,   270,   272,
     0,     0,     0,     0,     0,   466,   461,     0,   134,     0,
     0,   134,     0,   134,     0,   134,     0,   398,   397,   131,
     0,     0,   402,   403,     0,   369,   367,     0,     0,     0,
   363,   143,     0,     0,   446,     0,   134,   292,   294,     0,
     0,     0,     0,   328,    40,    57,     0,   447,     0,   447,
     0,   452,     0,   121,   120,     0,   212,     0,   202,   210,
     0,     0,     0,     0,     0,   190,     0,   134,     0,     0,
   134,   182,     0,     0,     0,   246,     0,   134,     0,     0,
     0,   238,   134,   134,   134,   134,   134,   228,     0,     0,
     0,     0,     0,   276,   278,   281,   101,   103,     0,     0,
   489,   490,     0,     0,   473,   480,     0,     0,   471,     0,
   478,     0,   469,     0,     0,     0,     0,     0,   354,   362,
   415,   417,   447,     0,     0,     0,     0,     0,     0,   134,
     0,     0,   440,     0,   442,     0,     0,     0,   438,   287,
     0,   134,     0,   205,   175,     0,     0,   195,     0,     0,
   178,     0,     0,     0,     0,     0,    43,     0,     0,   128,
   447,   224,     0,     0,     0,     0,     0,   227,     0,     0,
     0,     0,     0,    43,   255,     0,   263,   265,   267,   492,
     0,     0,   392,     0,     0,   134,   134,     0,     0,   134,
     0,   134,     0,   134,     0,     0,     0,   359,   357,     0,
     0,   545,     0,   421,   331,   333,   327,     0,     0,     0,
     0,     0,     0,   448,     0,   450,     0,     0,   287,     0,
     0,     0,   215,     0,     0,   193,     0,   134,     0,     0,
   192,     0,     0,     0,   172,     0,     0,     0,     0,     0,
   134,     0,     0,     0,   240,     0,     0,   242,     0,     0,
     0,     0,     0,     0,     0,     0,   134,     0,   134,   134,
   134,   134,     0,     0,   475,   468,     0,     0,   482,   465,
     0,   477,     0,   463,     0,   407,   405,     0,     0,     0,
     0,   295,     0,    41,   447,   447,   452,     0,     0,   288,
     0,   286,   214,   213,     0,   134,     0,     0,   209,     0,
     0,   188,     0,     0,   189,   191,     0,   197,   127,   126,
   134,   124,   123,   183,     0,     0,     0,     0,     0,     0,
   236,   254,   237,   241,     0,     0,   252,   251,     0,     0,
   250,   229,     0,   232,     0,     0,     0,     0,   494,   496,
     0,     0,   134,     0,     0,     0,   346,   355,     0,   549,
   548,     0,   441,   443,   439,     0,     0,   134,   285,     0,
     0,     0,     0,   206,     0,   176,   199,   196,     0,   134,
     0,     0,   186,     0,     0,   180,     0,   225,   258,   260,
   248,   262,     0,     0,     0,     0,     0,   493,   134,    43,
     0,     0,     0,     0,     0,     0,     0,   134,   329,   449,
   451,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,   134,     0,     0,     0,     0,     0,     0,
     0,     0,   474,   481,     0,   472,   479,   470,   551,     0,
   289,   217,   216,     0,   208,     0,     0,     0,   173,     0,
   130,   129,     0,     0,     0,   256,     0,     0,     0,     0,
   497,     0,   134,   550,     0,   187,   200,   198,   184,     0,
   235,   230,   264,   266,   268,   495,   447,     0,     0,     0,
     0,     0,   483,     0,     0,   181,     0,   552,   207,   498,
   134,     0,     0,   499,     0,     0
};

static const short yydefgoto[] = {     1,
     6,    36,    27,    65,   136,    69,   129,    70,   131,   171,
   316,   407,   901,   522,   525,    71,   132,   173,   411,   410,
   721,   174,   175,   317,    72,   133,   176,   245,    73,   130,
   138,   168,   238,    74,    75,   177,   246,   529,   442,   665,
   666,   417,   418,   640,   419,   453,   420,   638,   841,   421,
   838,  1011,   686,    33,    53,    54,   128,    55,   140,   202,
   435,   508,   510,   224,   396,   423,   647,  1081,   735,   995,
   645,   845,  1184,   927,  1142,   843,  1071,  1072,  1000,   737,
   836,   998,  1140,   424,   642,   425,   734,   994,  1068,   732,
   733,   912,  1066,   428,   651,   429,   752,  1013,   654,   944,
  1153,   860,  1090,  1091,   755,  1017,   852,   431,   649,   855,
   857,  1099,   432,   848,  1027,   433,   853,   854,   856,  1029,
  1030,  1031,   436,   661,   662,   438,   439,   556,   440,   557,
   333,   528,   988,   989,  1128,    76,   288,   608,   383,   895,
   497,     7,    28,    37,    85,   151,   214,   388,   716,   900,
   215,   615,   216,   616,    77,   261,   354,   262,   358,   268,
   474,   475,   702,   970,   592,   591,   356,   800,   142,   260,
   471,   269,    87,   156,   225,    88,    40,    78,   137,   253,
   338,   578,   577,   691,   690,   579,    79,   144,   279,   280,
   599,   281,   600,    90,   149,   291,   384,   178,   254,   255,
   256,   444,   257,   285,   286,   725,   381,   909,   570,   571,
   341,   568,   786,   782,   342,   779,   964,   960,   956,   343,
   784,   962,   344,   957,  1113,   345,   564,   775,   871,  1032,
  1159,   872,  1160,  1237,  1251,     9,    29,    41,    42,    43,
   226,    10,    11,    12,    99,   233,    13,    14,    25,    15,
    30,    46,    98,   231,   123,   124,   125,   227,   126,   404,
    80,   145,   206,   490,   606,   807,   607,   806,  1051,  1168,
  1223
};

static const short yypact[] = {-32768,
    66,    16,-32768,    28,-32768,-32768,-32768,   193,-32768,-32768,
-32768,    73,-32768,    62,-32768,-32768,    89,-32768,-32768,   267,
   189,-32768,   264,-32768,-32768,     5,-32768,-32768,    81,   223,
-32768,   297,   612,  1173,  1173,   425,   481,  1173,    81,   317,
   345,    81,    81,-32768,-32768,-32768,-32768,   346,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,   387,   405,   260,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
   385,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,   421,-32768,-32768,-32768,   455,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,   518,-32768,-32768,
   471,-32768,-32768,   508,   549,-32768,  1173,-32768,   438,   620,
   530,   112,   602,-32768,-32768,-32768,   935,-32768,    36,-32768,
   701,   561,-32768,-32768,-32768,-32768,-32768,-32768,   486,-32768,
   562,-32768,-32768,-32768,   489,   777,   580,   588,-32768,   622,
   758,   596,   603,-32768,-32768,   390,-32768,-32768,   642,    92,
-32768,   671,   960,-32768,-32768,-32768,-32768,   268,   646,-32768,
-32768,   659,   677,   679,   561,   698,-32768,  1006,   434,-32768,
-32768,-32768,-32768,   359,-32768,-32768,-32768,-32768,-32768,-32768,
   697,   706,   704,    73,   710,   710,   710,   710,   710,    73,
   722,   390,-32768,   473,   303,   649,-32768,-32768,-32768,  1026,
-32768,  1497,-32768,-32768,-32768,-32768,   731,   734,    95,-32768,
-32768,-32768,   728,-32768,   744,-32768,-32768,   773,    73,-32768,
-32768,    73,-32768,-32768,-32768,-32768,-32768,   648,-32768,-32768,
-32768,    40,-32768,-32768,   682,   577,  1173,   503,   747,-32768,
   750,    44,-32768,-32768,-32768,-32768,-32768,-32768,   359,-32768,
   760,   806,-32768,   120,   710,   710,   470,   762,-32768,   808,
   812,   710,   710,   710,   798,   727,   820,   827,-32768,-32768,
-32768,-32768,-32768,   838,-32768,   759,   815,-32768,-32768,    44,
-32768,   113,    97,    73,-32768,-32768,   390,   390,   390,   390,
   143,   851,   829,   390,   390,   390,-32768,   988,  1029,-32768,
-32768,   856,   843,-32768,-32768,-32768,-32768,   826,-32768,-32768,
   257,-32768,-32768,   772,-32768,   876,   888,   888,   907,   499,
   918,   855,-32768,-32768,-32768,-32768,-32768,   921,   892,-32768,
-32768,-32768,-32768,-32768,-32768,   176,-32768,-32768,   925,   930,
  1055,   358,   932,-32768,   550,   941,   806,-32768,   359,-32768,
   948,   952,-32768,-32768,   958,   963,   972,-32768,-32768,-32768,
   955,   982,   710,-32768,   390,   197,   971,   978,-32768,  1025,
-32768,  1016,   623,  1034,  1043,   757,    73,-32768,  1050,   386,
   727,   727,   727,   143,   143,  1059,-32768,  1064,   727,   727,
-32768,-32768,-32768,-32768,-32768,-32768,   332,   578,  1049,   561,
   561,-32768,-32768,-32768,-32768,   662,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,   844,-32768,-32768,-32768,
-32768,-32768,-32768,    64,   728,-32768,-32768,-32768,-32768,-32768,
   168,-32768,-32768,-32768,   462,   359,    34,-32768,  1076,  1095,
  1101,  1115,  1119,   763,-32768,    44,    44,  1124,  1125,    44,
    44,    44,-32768,-32768,-32768,-32768,-32768,-32768,   579,   359,
-32768,-32768,-32768,   806,   806,-32768,-32768,   195,  1126,-32768,
-32768,  1109,   727,-32768,-32768,-32768,   657,   922,-32768,  1186,
-32768,-32768,-32768,-32768,-32768,-32768,  1128,   462,-32768,-32768,
-32768,  1137,  1198,   114,-32768,   277,   727,  1139,   727,  1141,
   315,-32768,-32768,-32768,-32768,   327,   383,  1211,-32768,-32768,
-32768,   241,-32768,-32768,   111,-32768,-32768,   649,-32768,-32768,
   359,-32768,  1143,-32768,  1147,  1148,-32768,  1151,-32768,  1163,
-32768,  1160,  1164,  1165,-32768,  1167,  1168,  1170,  1171,-32768,
-32768,-32768,-32768,  1104,   564,  1111,  1123,-32768,-32768,   670,
-32768,-32768,-32768,  1176,   359,  1200,  1215,-32768,  1113,   868,
-32768,    31,   159,   166,   181,   763,-32768,-32768,    47,    54,
-32768,-32768,    44,  1222,    44,  1224,  1230,  1218,  1231,  1235,
-32768,-32768,   359,-32768,   183,  1238,-32768,   394,-32768,-32768,
-32768,  1076,   818,-32768,-32768,-32768,-32768,    44,   100,-32768,
-32768,-32768,-32768,  1214,-32768,-32768,-32768,-32768,-32768,  1250,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,  1253,   240,
-32768,-32768,-32768,-32768,    38,   718,   359,-32768,   719,   871,
-32768,  1257,-32768,-32768,  1258,   312,   649,   631,-32768,-32768,
   903,-32768,  1242,  1273,  1270,  1274,   872,   411,   390,   243,
  1279,   564,   390,   390,-32768,-32768,-32768,-32768,  1286,-32768,
-32768,  1221,-32768,-32768,   359,-32768,  1294,  1298,-32768,  1310,
-32768,  1317,-32768,  1324,  1327,   964,-32768,-32768,-32768,-32768,
-32768,  1330,-32768,  1331,  1239,  1346,  1353,-32768,-32768,   769,
   775,  1126,-32768,  1360,  1225,  1229,   905,-32768,-32768,  1364,
  1311,-32768,-32768,  1248,  1252,-32768,  1349,-32768,  1314,-32768,
   561,  1384,  1391,  1398,  1307,-32768,  1404,  1255,-32768,-32768,
   252,  1405,  1257,-32768,-32768,   956,  1409,-32768,   359,-32768,
   759,   359,-32768,  1421,-32768,  1428,-32768,  1259,   326,-32768,
-32768,-32768,-32768,   977,  1435,-32768,-32768,-32768,-32768,-32768,
-32768,  1439,  1442,-32768,  1443,  1446,  1458,   423,-32768,-32768,
   430,   472,  1015,  1077,   923,-32768,-32768,   359,-32768,  1453,
   226,-32768,  1457,-32768,  1460,-32768,  1461,-32768,-32768,-32768,
   785,   792,-32768,-32768,  1449,-32768,-32768,   710,   710,  1482,
  1126,-32768,  1494,  1501,-32768,   649,-32768,-32768,-32768,   944,
  1491,    73,  1262,-32768,-32768,-32768,   996,   759,  1017,   759,
  1051,   927,  1509,-32768,-32768,   359,-32768,  1511,-32768,-32768,
   824,   825,  1076,   279,  1516,-32768,   359,-32768,  1512,   359,
-32768,-32768,  1258,   351,   649,-32768,   831,-32768,  1076,   743,
  1521,-32768,-32768,-32768,-32768,-32768,-32768,-32768,  1504,  1273,
  1505,  1506,  1507,-32768,-32768,-32768,-32768,-32768,    17,  1524,
-32768,-32768,  1528,  1266,-32768,-32768,  1530,  1285,-32768,  1289,
-32768,  1292,-32768,    44,    44,   710,  1531,  1532,-32768,-32768,
-32768,-32768,   759,  1296,  1049,  1533,  1534,  1535,  1536,-32768,
  1049,  1076,-32768,  1076,-32768,  1076,  1538,  1539,-32768,  1540,
   359,-32768,   469,-32768,-32768,   946,   359,-32768,  1541,  1542,
-32768,  1543,  1299,  1526,  1545,  1303,-32768,  1409,   359,-32768,
   759,-32768,  1322,   990,  1547,  1548,  1549,-32768,  1326,  1329,
  1333,  1336,  1340,-32768,-32768,  1435,-32768,-32768,-32768,-32768,
  1550,  1452,-32768,  1551,  1552,-32768,-32768,  1537,  1553,-32768,
  1554,-32768,  1556,-32768,  1557,  1558,  1559,-32768,-32768,  1126,
  1561,-32768,   136,-32768,-32768,-32768,-32768,  1359,   172,   992,
  1008,  1020,  1071,-32768,  1090,-32768,  1520,  1562,  1540,  1563,
  1363,   359,-32768,  1565,  1566,-32768,   359,-32768,   262,  1567,
  1541,   520,  1568,  1569,-32768,  1570,  1571,   848,  1572,   359,
-32768,  1574,  1575,  1576,-32768,    41,  1577,  1547,  1560,   251,
   390,   390,    35,   390,  1578,   886,-32768,  1580,-32768,-32768,
-32768,-32768,  1564,  1582,-32768,-32768,  1366,  1370,-32768,-32768,
  1373,-32768,  1377,-32768,  1396,-32768,-32768,  1583,  1584,  1161,
  1585,-32768,  1586,-32768,   759,   759,   927,  1076,  1076,-32768,
  1587,-32768,-32768,-32768,   359,-32768,   249,  1588,-32768,   255,
  1589,-32768,  1590,  1400,-32768,-32768,   359,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,  1591,  1592,  1403,  1579,   382,  1593,
-32768,-32768,-32768,-32768,   534,   569,-32768,-32768,  1594,   585,
-32768,-32768,  1407,-32768,  1410,  1414,  1433,  1437,-32768,-32768,
  1595,  1596,-32768,  1598,  1599,  1600,-32768,-32768,  1602,-32768,
-32768,  1603,-32768,-32768,-32768,  1048,  1053,-32768,-32768,  1604,
  1440,  1076,  1257,-32768,  1258,-32768,-32768,-32768,   359,-32768,
  1444,  1566,-32768,  1605,  1606,-32768,  1273,-32768,-32768,-32768,
-32768,-32768,  1575,  1607,   390,   390,   390,-32768,-32768,-32768,
  1608,  1609,  1447,  1610,  1611,  1612,  1173,-32768,-32768,-32768,
-32768,  1451,  1613,  1614,  1616,  1617,  1409,  1618,  1470,  1619,
  1620,  1621,  1622,-32768,  1435,  1623,  1624,   608,   617,   660,
  1474,   893,-32768,-32768,  1625,-32768,-32768,-32768,-32768,  1477,
-32768,-32768,-32768,  1110,-32768,  1626,  1627,  1628,-32768,  1629,
-32768,-32768,  1481,  1630,  1631,-32768,  1632,  1633,  1634,  1635,
-32768,  1636,-32768,-32768,  1257,-32768,-32768,-32768,-32768,  1637,
-32768,-32768,-32768,-32768,-32768,-32768,   759,  1638,  1484,  1639,
  1640,  1642,-32768,  1643,  1644,-32768,  1641,-32768,-32768,-32768,
-32768,  1488,  1645,-32768,  1646,-32768
};

static const short yypgoto[] = {-32768,
-32768,-32768,-32768,   -34,-32768,  1517,-32768,-32768,-32768,-32768,
-32768,  -309,-32768,-32768,  -833,-32768,-32768,  -384,-32768,-32768,
-32768,  1390,-32768,-32768,-32768,-32768,-32768,-32768,  1615,-32768,
-32768,  1432,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,  1251,-32768,-32768,  -499,  -136,-32768,-32768,-32768,-32768,
-32768,-32768,   412,   -99,  -422,  1159,-32768,  -177,  -305,  -254,
  -154,  1144,  1145,-32768,   696,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,   478,  -731,   652,  -828,
  -884,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,  -721,
   656,-32768,-32768,-32768,-32768,  1004,-32768,-32768,-32768,-32768,
-32768,-32768,   504,  -746,  -832,   638,  -324,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,   997,  1328,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,   672,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,  1647,  1455,   281,-32768,-32768,  -179,
-32768,-32768,-32768,-32768,-32768,-32768,  -327,   859,  -356,-32768,
-32768,   -78,-32768,-32768,-32768,-32768,    12,-32768,-32768,-32768,
-32768,  -238,-32768,-32768,-32768,  1207,  1648,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,    58,-32768,-32768,
-32768,-32768,-32768,  -616,  1138,-32768,  -716,   610,  -586,  -210,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,  1172,   959,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,   694,-32768,-32768,
-32768,-32768,  1653,-32768,-32768,-32768,  1654,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,  1546,-32768,-32768,-32768,-32768,
  1649,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768
};


#define	YYLAST		1686


static const short yytable[] = {   127,
    66,   223,   201,    92,   430,   835,   259,   408,   851,   201,
   479,   830,     8,   349,   928,   707,   270,   271,   272,   273,
   375,   422,   422,     8,   839,   526,   527,   946,   163,   468,
   741,   576,   143,   641,  1001,   143,    48,    48,   950,   164,
   340,    16,   318,   675,    20,    34,   346,   276,   469,   346,
   687,   385,   676,    50,    51,    17,   346,   689,   566,    52,
   319,   973,   541,   165,    23,  1255,    35,   979,     2,   487,
   550,   551,   552,   553,   166,   951,   554,   297,   298,   299,
   300,   350,   150,   722,   723,   361,   362,   301,   186,   187,
   188,    26,   371,   372,   373,   677,   678,     3,     4,   544,
   952,   903,   711,   905,     3,     4,    18,    19,   297,   298,
   299,   300,     3,     4,   633,   386,  1001,   555,   301,   936,
   712,   376,   596,     3,     4,   724,   308,   309,   -21,   -21,
   201,   634,   359,   158,   560,   619,   387,   302,   167,  1052,
   730,   360,   390,   391,   392,   393,   589,   590,  1097,   399,
   400,   401,   320,   685,   347,   348,   634,   347,   348,     5,
   351,    38,   394,   303,   347,   348,     5,   395,   302,   489,
   159,   230,   620,    39,     5,  1054,   971,   294,   304,   305,
   679,   470,   143,   357,   229,     5,   558,   681,   454,   893,
   -21,    32,   634,   482,   303,   593,  1009,   455,   559,  1028,
   484,   204,   683,    21,   185,   266,   210,   593,   143,   304,
   305,   274,   334,   306,  1012,   201,   185,   680,   201,   636,
   483,   583,   584,   585,   682,    44,    45,   422,   931,   456,
   457,   458,   459,   460,   461,   506,   507,   509,   511,   684,
   311,   192,   719,   312,   516,   517,   916,   876,   594,   550,
   551,   552,   553,   192,   630,   554,   297,   298,   299,   300,
   720,   826,   934,   790,   826,   631,   301,   917,   565,   632,
   827,   235,   761,   827,   917,   412,   918,   413,   357,   540,
   621,   700,   701,   918,   877,   462,   236,   237,   503,    16,
   549,   917,   587,   306,   103,   104,   105,   106,   107,    31,
   918,   108,  1132,  1133,   518,   389,  1177,   828,   919,   763,
   828,   109,  1135,   920,  1185,   980,   302,   981,   626,   982,
   920,   282,    47,   110,   739,   111,   430,    94,   598,   283,
   627,   306,   919,   740,   422,   519,   816,   920,  1123,  1124,
   688,   688,   303,   306,   692,   801,   694,   541,    95,   112,
   790,   790,   507,   637,   509,   113,   520,   304,   305,   674,
   521,    48,   114,   929,   180,   181,   182,   183,   100,   710,
   184,   588,   930,   186,   187,   188,   115,    50,    51,   185,
   791,   792,   543,    52,   544,   545,   628,   669,   219,   505,
   101,   141,   219,   375,   147,   357,   357,   704,   502,   306,
   546,   547,   306,   541,   548,   186,   187,   188,   102,   220,
   306,  1176,   190,   220,   191,   259,   192,   221,   222,   103,
   104,   221,   222,   121,   831,   832,   864,    67,    68,   186,
   187,   188,   761,   865,   738,   935,   116,   258,  1214,   306,
   544,  1147,   847,   756,   801,  1206,   306,   117,   110,   118,
   111,   119,   236,   237,   120,   605,   193,   122,   194,   727,
   936,   195,   196,   197,   198,   199,   200,   762,   672,   763,
   764,  1126,  1127,   139,   112,   866,   180,   181,   182,   183,
   113,   992,   184,    81,    82,   765,   766,   114,   306,   767,
   993,   185,   180,   181,   182,   183,   674,   778,   184,   705,
   706,   434,   540,  1240,   768,   152,   153,   185,   771,   772,
   146,   818,   820,   822,   277,   714,   715,   186,   187,   188,
  1242,   335,   220,   363,   364,   365,   366,   367,   192,   336,
   221,   222,  1077,   186,   187,   188,   134,  1149,   728,   189,
   190,  1078,   191,   135,   192,  1175,   561,   562,   563,   748,
   306,   148,   247,   248,   249,   157,   180,   181,   182,   183,
   251,   837,   184,   172,   840,   203,   660,   278,   709,   205,
   194,   185,  1150,   195,   196,   197,   198,   199,   200,   324,
   325,   523,   586,   208,   193,   306,   194,   220,  1152,   195,
   196,   197,   198,   199,   200,   221,   222,   186,   187,   188,
   873,   306,   520,   520,   122,   160,   521,   521,   192,   209,
  1018,  1217,   549,  1049,    48,    49,   813,  1008,   887,   888,
  1218,   326,   327,   328,   306,   211,   492,   329,   217,   330,
    50,    51,   738,   306,  1026,   218,    52,   180,   181,   182,
   183,   493,   494,   742,   228,   965,   966,   495,   911,   756,
   194,   496,   743,   195,   196,   197,   198,   199,   200,   922,
   601,   313,   925,  1219,   239,   530,   314,   315,   180,   181,
   182,   183,   331,   667,   531,   232,   306,   240,   332,   874,
   241,   520,   878,   532,   880,   521,   882,   744,   745,   746,
   154,   155,   242,  1018,   520,   321,   984,   986,   521,   263,
   322,   244,   323,   169,   170,   674,   967,   894,   264,   186,
   187,   188,   267,   937,   747,   533,   534,   535,   536,   537,
   538,   726,   265,   674,   275,   180,   181,   182,   183,   219,
   284,   531,    93,   990,   295,    96,    97,   296,   923,   997,
   532,   926,   520,   306,   306,   539,   521,   307,   933,   337,
   220,  1010,   339,   939,   940,   941,   942,   943,   221,   222,
   212,   213,   352,   353,   541,   368,   186,   187,   188,   674,
   674,   674,   798,   534,   535,   536,   537,   538,   799,   500,
   501,    50,    51,   180,   181,   182,   183,    52,   884,   184,
   186,   187,   188,   520,   310,   885,   935,   521,   185,   520,
   978,   544,   539,   521,   180,   181,   182,   183,   355,   520,
   374,   369,   991,   521,  1065,   370,   520,   352,   382,  1073,
   521,   936,   377,   898,   186,   187,   188,   914,   915,   378,
   207,   190,  1086,   191,   932,   192,    56,    57,    58,    59,
   379,    60,    61,    62,    63,  1098,   380,    64,   520,   520,
  1192,  1084,   521,   521,   397,   520,  1037,  1038,   398,   521,
  1041,   405,  1043,   406,  1045,   541,  1095,  1096,   414,  1100,
   569,   673,   520,   639,   729,   193,   521,   194,   415,   937,
   195,   196,   197,   198,   199,   200,   409,  1130,   443,  1102,
   416,   186,   187,   188,   759,   760,  1221,   542,  1074,  1139,
   543,   708,   544,   545,   446,   749,   750,   569,   805,   427,
   520,  1087,   447,   448,   521,   674,   674,   520,   546,   547,
   441,   521,   548,   445,   738,   869,   870,  1103,   463,  1105,
  1106,  1107,  1108,   464,   446,   466,   756,    48,   162,   186,
   187,   188,   447,   448,   467,   449,   869,   896,   569,   996,
   450,   472,   937,    50,    51,   473,   451,   452,   480,    52,
   476,  1178,    48,   234,   674,   477,  1131,   789,   446,   186,
   187,   188,   907,   908,   478,   602,   447,   448,    50,    51,
   450,  1141,    50,    51,    52,   481,   451,   452,    52,   446,
    48,   402,   569,  1015,   569,  1055,   485,   447,   448,   486,
  1188,  1189,  1190,   186,   187,   188,    50,    51,   446,   833,
   569,  1056,    52,  1163,   450,  1120,   447,   448,   867,   491,
   451,   452,   569,  1057,   186,   187,   188,   488,  1172,   446,
   849,    48,   403,    50,    51,   450,   498,   447,   448,    52,
  1179,   451,   452,   186,   187,   188,   499,    50,    51,   902,
   569,  1170,   504,    52,   450,   569,  1171,    48,   465,  1191,
   451,   452,   514,   446,   186,   187,   188,   515,  1200,   524,
   904,   447,   448,    50,    51,   450,   773,   774,   569,    52,
   868,   451,   452,   446,  1213,   247,   248,   249,   250,   512,
   513,   447,   448,   251,   252,    50,    51,   572,   186,   187,
   188,    52,   446,   573,   906,   247,   248,   249,   289,   450,
   447,   448,   597,   251,   290,   451,   452,   574,   186,   187,
   188,   575,   446,  1239,  1058,   446,   581,   582,   595,   450,
   447,   448,  1199,   447,   448,   451,   452,   186,   187,   188,
   617,   609,   623,  1059,   625,   639,   610,   611,   450,   643,
   644,  1252,   612,   646,   451,   452,   613,   186,   187,   188,
   186,   187,   188,   650,  1225,   648,   652,   653,   450,   655,
   656,   450,   657,   658,   451,   452,   659,   451,   452,    56,
    57,    58,    59,   663,    60,    61,    62,    63,   603,   604,
    64,    56,    57,    58,    59,   664,    60,    61,    62,    63,
    48,   618,    64,   670,    56,    57,    58,    59,   668,    60,
    61,    62,    63,    48,   629,    64,    50,    51,   671,  1119,
    48,   697,    52,    48,   777,   693,   695,    48,   803,    50,
    51,    48,   804,   696,   698,    52,    50,    51,   699,    50,
    51,   703,    52,    50,    51,    52,   713,    50,    51,    52,
    48,   811,   717,    52,    48,   812,   718,    48,   825,   731,
   736,    48,   846,   753,    48,   899,    50,    51,    48,   955,
    50,    51,    52,    50,    51,   754,    52,    50,    51,    52,
    50,    51,   769,    52,    50,    51,    52,    48,   959,   776,
    52,    48,   961,   757,    48,   963,   780,   758,    48,   972,
   781,    48,  1004,    50,    51,    48,  1007,    50,    51,    52,
    50,    51,   783,    52,    50,    51,    52,    50,    51,   785,
    52,    50,    51,    52,    48,  1014,   787,    52,    48,  1021,
   788,    48,  1022,   793,   794,    48,  1023,   795,    48,  1024,
    50,    51,    48,  1025,    50,    51,    52,    50,    51,   796,
    52,    50,    51,    52,    50,    51,   797,    52,    50,    51,
    52,    48,  1053,   802,    52,    48,  1064,   808,    48,  1111,
   814,   809,    48,  1112,   815,    48,  1114,    50,    51,    48,
  1115,    50,    51,    52,    50,    51,   817,    52,    50,    51,
    52,    50,    51,   819,    52,    50,    51,    52,    48,  1116,
   821,    52,    48,  1138,   823,    48,  1145,   824,   829,    48,
  1154,   834,    48,  1155,    50,    51,    48,  1156,    50,    51,
    52,    50,    51,   842,    52,    50,    51,    52,    50,    51,
   844,    52,    50,    51,    52,    48,  1157,   850,    52,    48,
  1158,   858,    48,  1174,   859,   861,    48,  1180,   862,    48,
  1195,    50,    51,    48,  1201,    50,    51,    52,    50,    51,
   863,    52,    50,    51,    52,    50,    51,   886,    52,    50,
    51,    52,    48,  1208,   875,    52,    48,  1220,   879,    48,
  1224,   881,   883,    48,  1230,   889,    48,  1244,    50,    51,
    48,  1253,    50,    51,    52,    50,    51,   891,    52,    50,
    51,    52,    50,    51,   892,    52,    50,    51,    52,   292,
   897,   910,    52,   913,   924,    56,    57,    58,   293,   921,
    60,    61,    62,    63,   938,   945,    64,   953,   947,   948,
   949,   954,   958,  1034,   968,   969,   974,   975,   976,   977,
   983,   985,   987,   999,  1002,  1060,  1003,  1005,  1006,  1016,
  1019,  1020,  1033,    83,  1035,  1036,  1040,  1042,  1039,  1044,
  1046,  1047,  1048,  1050,   243,  1061,  1063,  1067,  1070,   179,
  1075,  1079,  1080,  1082,  1083,  1085,  1088,  1089,   426,  1092,
  1093,  1101,   759,  1104,  1110,  1109,  1117,  1118,  1121,  1122,
  1129,  1134,  1136,  1137,  1143,  1144,  1148,  1151,  1161,  1162,
  1146,  1164,  1165,  1166,  1167,   567,  1169,  1173,  1182,  1183,
  1187,  1193,  1194,  1196,  1197,  1198,  1202,  1203,  1204,  1181,
  1205,  1207,  1209,  1210,  1211,  1212,  1215,  1216,  1222,  1226,
  1227,  1228,  1229,  1231,  1232,  1233,  1234,  1235,  1236,  1238,
  1241,  1243,  1245,  1246,  1247,  1256,  1248,  1249,  1254,  1069,
   622,    84,  1076,   624,   751,  1094,  1186,   437,   770,   890,
  1062,   287,  1250,   580,    22,   635,  1125,    24,   161,   614,
     0,   810,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,    86,    89,    91
};

static const short yycheck[] = {    99,
    35,   156,   139,    38,   329,   737,   184,   317,   755,   146,
   367,   733,     1,   252,   843,   602,   196,   197,   198,   199,
   275,   327,   328,    12,   741,   410,   411,   860,   128,   357,
   647,   454,   111,   533,   919,   114,     3,     3,    22,     4,
   251,    26,     3,    13,    17,    41,     3,   202,   358,     3,
     4,   290,    22,    19,    20,    40,     3,     4,    25,    25,
    21,   895,    22,    28,     3,     0,    62,   901,     3,   379,
     7,     8,     9,    10,    39,    59,    13,    14,    15,    16,
    17,   259,   117,    46,    47,   265,   266,    24,    48,    49,
    50,     3,   272,   273,   274,    65,    66,    32,    33,    59,
    84,   818,     3,   820,    32,    33,    91,    92,    14,    15,
    16,    17,    32,    33,     4,     3,  1001,    54,    24,    79,
    21,   276,   479,    32,    33,    88,   226,   227,    32,    33,
   267,    21,    13,    22,   444,    22,    24,    74,   103,     4,
   640,    22,   297,   298,   299,   300,   474,   475,   114,   304,
   305,   306,   113,   576,   111,   112,    21,   111,   112,    94,
   260,    81,    20,   100,   111,   112,    94,    25,    74,   380,
    59,   160,    59,    93,    94,     4,   893,   212,   115,   116,
    22,   359,   261,   262,    93,    94,    19,    22,    13,   806,
    94,     3,    21,   373,   100,    13,   928,    22,    31,   946,
     4,   144,    22,    11,    22,   194,   149,    13,   287,   115,
   116,   200,   247,    17,   931,   352,    22,    59,   355,   529,
   375,   460,   461,   462,    59,     3,     4,   533,   845,    54,
    55,    56,    57,    58,    59,   390,   391,   392,   393,    59,
   229,    59,     3,   232,   399,   400,   833,    22,    54,     7,
     8,     9,    10,    59,    14,    13,    14,    15,    16,    17,
    21,    13,   849,   686,    13,    25,    24,    13,   446,    29,
    22,     4,    22,    22,    13,    19,    22,    21,   357,   416,
     4,   591,   592,    22,    59,   110,    19,    20,   388,    26,
   427,    13,   470,    17,    35,    36,    37,    38,    39,    33,
    22,    42,    54,    55,   404,   294,  1135,    59,    54,    59,
    59,    52,    58,    59,  1147,   902,    74,   904,     4,   906,
    59,    19,    26,    64,    13,    66,   651,    11,   483,    27,
     4,    17,    54,    22,   640,     4,   721,    59,  1055,  1056,
   579,   580,   100,    17,   583,   702,   585,    22,     4,    90,
   773,   774,   507,   531,   509,    96,    25,   115,   116,   570,
    29,     3,   103,    13,     7,     8,     9,    10,    23,   608,
    13,   471,    22,    48,    49,    50,   117,    19,    20,    22,
   690,   691,    57,    25,    59,    60,     4,   565,     3,     4,
     4,   111,     3,   648,   114,   474,   475,     4,   387,    17,
    75,    76,    17,    22,    79,    48,    49,    50,     4,    24,
    17,  1133,    55,    24,    57,   593,    59,    32,    33,    35,
    36,    32,    33,     3,   734,   735,     4,     3,     4,    48,
    49,    50,    22,     4,   645,    54,    52,     4,  1185,    17,
    59,    60,   752,   654,   801,  1177,    17,    63,    64,    65,
    66,    67,    19,    20,    70,   490,    99,     3,   101,   637,
    79,   104,   105,   106,   107,   108,   109,    57,   568,    59,
    60,  1058,  1059,     3,    90,     4,     7,     8,     9,    10,
    96,    13,    13,     3,     4,    75,    76,   103,    17,    79,
    22,    22,     7,     8,     9,    10,   707,   675,    13,   599,
   600,     3,   639,  1225,   659,    68,    69,    22,   663,   664,
     3,   722,   723,   724,    42,   615,   616,    48,    49,    50,
  1237,    19,    24,    54,    55,    56,    57,    58,    59,    27,
    32,    33,    13,    48,    49,    50,    19,     4,   638,    54,
    55,    22,    57,    26,    59,  1132,    85,    86,    87,   649,
    17,     3,    80,    81,    82,    26,     7,     8,     9,    10,
    88,   739,    13,     3,   742,     4,     3,    95,   603,    81,
   101,    22,     4,   104,   105,   106,   107,   108,   109,     3,
     4,     4,     4,     4,    99,    17,   101,    24,     4,   104,
   105,   106,   107,   108,   109,    32,    33,    48,    49,    50,
   778,    17,    25,    25,     3,     4,    29,    29,    59,    22,
   935,     4,   749,   970,     3,     4,   716,   927,   798,   799,
     4,    45,    46,    47,    17,     4,     4,    51,    33,    53,
    19,    20,   843,    17,   944,    33,    25,     7,     8,     9,
    10,    19,    20,    13,     3,   884,   885,    25,   826,   860,
   101,    29,    22,   104,   105,   106,   107,   108,   109,   837,
     4,    14,   840,     4,    19,     4,    19,    20,     7,     8,
     9,    10,    96,     4,    13,     5,    17,    19,   102,   779,
     4,    25,   782,    22,   784,    29,   786,    57,    58,    59,
    71,    72,    14,  1018,    25,    14,   907,   908,    29,     3,
    19,     4,    21,     3,     4,   916,   886,   807,     3,    48,
    49,    50,     3,   850,    84,    54,    55,    56,    57,    58,
    59,     4,    19,   934,     3,     7,     8,     9,    10,     3,
    82,    13,    39,   911,     4,    42,    43,     4,   838,   917,
    22,   841,    25,    17,    17,    84,    29,     4,   848,     3,
    24,   929,     3,   853,   854,   855,   856,   857,    32,    33,
     3,     4,     3,     4,    22,     4,    48,    49,    50,   980,
   981,   982,     4,    55,    56,    57,    58,    59,     4,    23,
    24,    19,    20,     7,     8,     9,    10,    25,     4,    13,
    48,    49,    50,    25,    22,     4,    54,    29,    22,    25,
   900,    59,    84,    29,     7,     8,     9,    10,     3,    25,
    13,     4,   912,    29,   992,     4,    25,     3,     4,   997,
    29,    79,     3,   812,    48,    49,    50,     4,     4,     3,
    54,    55,  1010,    57,     4,    59,    19,    20,    21,    22,
     3,    24,    25,    26,    27,  1023,    88,    30,    25,    25,
  1160,     4,    29,    29,     4,    25,   956,   957,    30,    29,
   960,     6,   962,    21,   964,    22,  1021,  1022,    97,  1024,
     3,     4,    25,     3,     4,    99,    29,   101,     3,  1016,
   104,   105,   106,   107,   108,   109,    61,  1065,    34,     4,
     3,    48,    49,    50,    23,    24,     4,    54,   998,  1077,
    57,    84,    59,    60,    13,     3,     4,     3,     4,     3,
    25,  1011,    21,    22,    29,  1126,  1127,    25,    75,    76,
     3,    29,    79,     3,  1135,     3,     4,  1027,     4,  1029,
  1030,  1031,  1032,     4,    13,     4,  1147,     3,     4,    48,
    49,    50,    21,    22,     4,    54,     3,     4,     3,     4,
    59,     4,  1089,    19,    20,     4,    65,    66,     4,    25,
     3,  1139,     3,     4,  1175,     3,  1066,     4,    13,    48,
    49,    50,    46,    47,     3,    54,    21,    22,    19,    20,
    59,  1081,    19,    20,    25,     4,    65,    66,    25,    13,
     3,     4,     3,     4,     3,     4,    26,    21,    22,    22,
  1155,  1156,  1157,    48,    49,    50,    19,    20,    13,    54,
     3,     4,    25,  1113,    59,  1050,    21,    22,     4,     4,
    65,    66,     3,     4,    48,    49,    50,     3,  1128,    13,
    54,     3,     4,    19,    20,    59,     3,    21,    22,    25,
  1140,    65,    66,    48,    49,    50,     4,    19,    20,    54,
     3,     4,     3,    25,    59,     3,     4,     3,     4,  1159,
    65,    66,     4,    13,    48,    49,    50,     4,  1168,    21,
    54,    21,    22,    19,    20,    59,   665,   666,     3,    25,
     4,    65,    66,    13,  1184,    80,    81,    82,    83,   394,
   395,    21,    22,    88,    89,    19,    20,     3,    48,    49,
    50,    25,    13,     3,    54,    80,    81,    82,    83,    59,
    21,    22,     4,    88,    89,    65,    66,     3,    48,    49,
    50,     3,    13,  1223,    54,    13,     3,     3,     3,    59,
    21,    22,  1167,    21,    22,    65,    66,    48,    49,    50,
     4,    14,     4,    54,     4,     3,    19,    20,    59,     3,
     3,  1251,    25,     3,    65,    66,    29,    48,    49,    50,
    48,    49,    50,     4,    55,     3,     3,     3,    59,     3,
     3,    59,     3,     3,    65,    66,    73,    65,    66,    19,
    20,    21,    22,    73,    24,    25,    26,    27,     3,     4,
    30,    19,    20,    21,    22,    73,    24,    25,    26,    27,
     3,     4,    30,     4,    19,    20,    21,    22,    33,    24,
    25,    26,    27,     3,     4,    30,    19,    20,     4,    59,
     3,     4,    25,     3,     4,     4,     3,     3,     4,    19,
    20,     3,     4,     4,     4,    25,    19,    20,     4,    19,
    20,     4,    25,    19,    20,    25,    33,    19,    20,    25,
     3,     4,     3,    25,     3,     4,     4,     3,     4,     3,
     3,     3,     4,    22,     3,     4,    19,    20,     3,     4,
    19,    20,    25,    19,    20,     3,    25,    19,    20,    25,
    19,    20,     4,    25,    19,    20,    25,     3,     4,     4,
    25,     3,     4,    24,     3,     4,     3,    24,     3,     4,
     3,     3,     4,    19,    20,     3,     4,    19,    20,    25,
    19,    20,     3,    25,    19,    20,    25,    19,    20,     3,
    25,    19,    20,    25,     3,     4,     3,    25,     3,     4,
     4,     3,     4,     4,     4,     3,     4,    99,     3,     4,
    19,    20,     3,     4,    19,    20,    25,    19,    20,     4,
    25,    19,    20,    25,    19,    20,     4,    25,    19,    20,
    25,     3,     4,     4,    25,     3,     4,     4,     3,     4,
    22,    61,     3,     4,    61,     3,     4,    19,    20,     3,
     4,    19,    20,    25,    19,    20,     3,    25,    19,    20,
    25,    19,    20,     3,    25,    19,    20,    25,     3,     4,
     3,    25,     3,     4,    98,     3,     4,     4,     4,     3,
     4,     3,     3,     4,    19,    20,     3,     4,    19,    20,
    25,    19,    20,     3,    25,    19,    20,    25,    19,    20,
     3,    25,    19,    20,    25,     3,     4,     3,    25,     3,
     4,     3,     3,     4,     3,     3,     3,     4,     3,     3,
     4,    19,    20,     3,     4,    19,    20,    25,    19,    20,
     3,    25,    19,    20,    25,    19,    20,    19,    25,    19,
    20,    25,     3,     4,    22,    25,     3,     4,    22,     3,
     4,    22,    22,     3,     4,     4,     3,     4,    19,    20,
     3,     4,    19,    20,    25,    19,    20,     4,    25,    19,
    20,    25,    19,    20,     4,    25,    19,    20,    25,    13,
    20,     3,    25,     3,     3,    19,    20,    21,    22,     4,
    24,    25,    26,    27,     4,    22,    30,     4,    24,    24,
    24,     4,     3,    82,     4,     4,     4,     4,     4,     4,
     3,     3,     3,     3,     3,    26,     4,    22,     4,     3,
     3,     3,     3,    37,     4,     4,     4,     4,    22,     4,
     4,     4,     4,     3,   175,     4,     4,     3,     3,   138,
     4,     4,     4,     4,     4,     4,     3,     3,   328,     4,
     4,     4,    23,     4,     3,    22,     4,     4,     4,     4,
     4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
    22,     4,     4,     4,     3,   447,     4,     4,     4,     4,
     4,     4,     4,     4,     4,     4,     4,     4,     3,  1142,
     4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
     4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
     4,     4,     4,     4,     3,     0,     4,     4,     4,   994,
   507,    37,  1001,   509,   651,  1018,  1153,   330,   662,   801,
   989,   207,    22,   457,    12,   528,  1057,    14,   123,   498,
    -1,   713,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    37,    37,    37
};

#line 325 "/usr/local/lib/bison.cc"
 /* fattrs + tables */

/* parser code folow  */


/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

/* Note: dollar marks section change
   the next  is replaced by the list of actions, each action
   as one case of the switch.  */ 

#if YY_PDDL_Parser_USE_GOTO != 0
/* 
 SUPRESSION OF GOTO : on some C++ compiler (sun c++)
  the goto is strictly forbidden if any constructor/destructor
  is used in the whole function (very stupid isn't it ?)
 so goto are to be replaced with a 'while/switch/case construct'
 here are the macro to keep some apparent compatibility
*/
#define YYGOTO(lb) {yy_gotostate=lb;continue;}
#define YYBEGINGOTO  enum yy_labels yy_gotostate=yygotostart; \
                     for(;;) switch(yy_gotostate) { case yygotostart: {
#define YYLABEL(lb) } case lb: {
#define YYENDGOTO } } 
#define YYBEGINDECLARELABEL enum yy_labels {yygotostart
#define YYDECLARELABEL(lb) ,lb
#define YYENDDECLARELABEL  };
#else
/* macro to keep goto */
#define YYGOTO(lb) goto lb
#define YYBEGINGOTO 
#define YYLABEL(lb) lb:
#define YYENDGOTO
#define YYBEGINDECLARELABEL 
#define YYDECLARELABEL(lb)
#define YYENDDECLARELABEL 
#endif
/* LABEL DECLARATION */
YYBEGINDECLARELABEL
  YYDECLARELABEL(yynewstate)
  YYDECLARELABEL(yybackup)
/* YYDECLARELABEL(yyresume) */
  YYDECLARELABEL(yydefault)
  YYDECLARELABEL(yyreduce)
  YYDECLARELABEL(yyerrlab)   /* here on detecting error */
  YYDECLARELABEL(yyerrlab1)   /* here on error raised explicitly by an action */
  YYDECLARELABEL(yyerrdefault)  /* current state does not do anything special for the error token. */
  YYDECLARELABEL(yyerrpop)   /* pop the current state because it cannot handle the error token */
  YYDECLARELABEL(yyerrhandle)  
YYENDDECLARELABEL
/* ALLOCA SIMULATION */
/* __HAVE_NO_ALLOCA */
#ifdef __HAVE_NO_ALLOCA
int __alloca_free_ptr(char *ptr,char *ref)
{if(ptr!=ref) free(ptr);
 return 0;}

#define __ALLOCA_alloca(size) malloc(size)
#define __ALLOCA_free(ptr,ref) __alloca_free_ptr((char *)ptr,(char *)ref)

#ifdef YY_PDDL_Parser_LSP_NEEDED
#define __ALLOCA_return(num) \
            return( __ALLOCA_free(yyss,yyssa)+\
		    __ALLOCA_free(yyvs,yyvsa)+\
		    __ALLOCA_free(yyls,yylsa)+\
		   (num))
#else
#define __ALLOCA_return(num) \
            return( __ALLOCA_free(yyss,yyssa)+\
		    __ALLOCA_free(yyvs,yyvsa)+\
		   (num))
#endif
#else
#define __ALLOCA_return(num) return(num)
#define __ALLOCA_alloca(size) alloca(size)
#define __ALLOCA_free(ptr,ref) 
#endif

/* ENDALLOCA SIMULATION */

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (YY_PDDL_Parser_CHAR = YYEMPTY)
#define YYEMPTY         -2
#define YYEOF           0
#define YYACCEPT        __ALLOCA_return(0)
#define YYABORT         __ALLOCA_return(1)
#define YYERROR         YYGOTO(yyerrlab1)
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL          YYGOTO(yyerrlab)
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(token, value) \
do                                                              \
  if (YY_PDDL_Parser_CHAR == YYEMPTY && yylen == 1)                               \
    { YY_PDDL_Parser_CHAR = (token), YY_PDDL_Parser_LVAL = (value);                 \
      yychar1 = YYTRANSLATE (YY_PDDL_Parser_CHAR);                                \
      YYPOPSTACK;                                               \
      YYGOTO(yybackup);                                            \
    }                                                           \
  else                                                          \
    { YY_PDDL_Parser_ERROR ("syntax error: cannot back up"); YYERROR; }   \
while (0)

#define YYTERROR        1
#define YYERRCODE       256

#ifndef YY_PDDL_Parser_PURE
/* UNPURE */
#define YYLEX           YY_PDDL_Parser_LEX()
#ifndef YY_USE_CLASS
/* If nonreentrant, and not class , generate the variables here */
int     YY_PDDL_Parser_CHAR;                      /*  the lookahead symbol        */
YY_PDDL_Parser_STYPE      YY_PDDL_Parser_LVAL;              /*  the semantic value of the */
				/*  lookahead symbol    */
int YY_PDDL_Parser_NERRS;                 /*  number of parse errors so far */
#ifdef YY_PDDL_Parser_LSP_NEEDED
YY_PDDL_Parser_LTYPE YY_PDDL_Parser_LLOC;   /*  location data for the lookahead     */
			/*  symbol                              */
#endif
#endif


#else
/* PURE */
#ifdef YY_PDDL_Parser_LSP_NEEDED
#define YYLEX           YY_PDDL_Parser_LEX(&YY_PDDL_Parser_LVAL, &YY_PDDL_Parser_LLOC)
#else
#define YYLEX           YY_PDDL_Parser_LEX(&YY_PDDL_Parser_LVAL)
#endif
#endif
#ifndef YY_USE_CLASS
#if YY_PDDL_Parser_DEBUG != 0
int YY_PDDL_Parser_DEBUG_FLAG;                    /*  nonzero means print parse trace     */
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif
#endif



/*  YYINITDEPTH indicates the initial size of the parser's stacks       */

#ifndef YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif


#if __GNUC__ > 1                /* GNU C and GNU C++ define this.  */
#define __yy_bcopy(FROM,TO,COUNT)       __builtin_memcpy(TO,FROM,COUNT)
#else                           /* not GNU C or C++ */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */

#ifdef __cplusplus
static void __yy_bcopy (char *from, char *to, int count)
#else
#ifdef __STDC__
static void __yy_bcopy (char *from, char *to, int count)
#else
static void __yy_bcopy (from, to, count)
     char *from;
     char *to;
     int count;
#endif
#endif
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}
#endif

int
#ifdef YY_USE_CLASS
 YY_PDDL_Parser_CLASS::
#endif
     YY_PDDL_Parser_PARSE(YY_PDDL_Parser_PARSE_PARAM)
#ifndef __STDC__
#ifndef __cplusplus
#ifndef YY_USE_CLASS
/* parameter definition without protypes */
YY_PDDL_Parser_PARSE_PARAM_DEF
#endif
#endif
#endif
{
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YY_PDDL_Parser_STYPE *yyvsp;
  int yyerrstatus;      /*  number of tokens to shift before error messages enabled */
  int yychar1=0;          /*  lookahead token as an internal (translated) token number */

  short yyssa[YYINITDEPTH];     /*  the state stack                     */
  YY_PDDL_Parser_STYPE yyvsa[YYINITDEPTH];        /*  the semantic value stack            */

  short *yyss = yyssa;          /*  refer to the stacks thru separate pointers */
  YY_PDDL_Parser_STYPE *yyvs = yyvsa;     /*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YY_PDDL_Parser_LSP_NEEDED
  YY_PDDL_Parser_LTYPE yylsa[YYINITDEPTH];        /*  the location stack                  */
  YY_PDDL_Parser_LTYPE *yyls = yylsa;
  YY_PDDL_Parser_LTYPE *yylsp;

#define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
#define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  int yystacksize = YYINITDEPTH;

#ifdef YY_PDDL_Parser_PURE
  int YY_PDDL_Parser_CHAR;
  YY_PDDL_Parser_STYPE YY_PDDL_Parser_LVAL;
  int YY_PDDL_Parser_NERRS;
#ifdef YY_PDDL_Parser_LSP_NEEDED
  YY_PDDL_Parser_LTYPE YY_PDDL_Parser_LLOC;
#endif
#endif

  YY_PDDL_Parser_STYPE yyval;             /*  the variable used to return         */
				/*  semantic values from the action     */
				/*  routines                            */

  int yylen;
/* start loop, in which YYGOTO may be used. */
YYBEGINGOTO

#if YY_PDDL_Parser_DEBUG != 0
  if (YY_PDDL_Parser_DEBUG_FLAG)
    fprintf(stderr, "Starting parse\n");
#endif
  yystate = 0;
  yyerrstatus = 0;
  YY_PDDL_Parser_NERRS = 0;
  YY_PDDL_Parser_CHAR = YYEMPTY;          /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YY_PDDL_Parser_LSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state, which is found in  yystate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
YYLABEL(yynewstate)

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      YY_PDDL_Parser_STYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;
#ifdef YY_PDDL_Parser_LSP_NEEDED
      YY_PDDL_Parser_LTYPE *yyls1 = yyls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = yyssp - yyss + 1;

#ifdef yyoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
#ifdef YY_PDDL_Parser_LSP_NEEDED
      /* This used to be a conditional around just the two extra args,
	 but that might be undefined if yyoverflow is a macro.  */
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yyls1, size * sizeof (*yylsp),
		 &yystacksize);
#else
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yystacksize);
#endif

      yyss = yyss1; yyvs = yyvs1;
#ifdef YY_PDDL_Parser_LSP_NEEDED
      yyls = yyls1;
#endif
#else /* no yyoverflow */
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	{
	  YY_PDDL_Parser_ERROR("parser stack overflow");
	  __ALLOCA_return(2);
	}
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;
      yyss = (short *) __ALLOCA_alloca (yystacksize * sizeof (*yyssp));
      __yy_bcopy ((char *)yyss1, (char *)yyss, size * sizeof (*yyssp));
      __ALLOCA_free(yyss1,yyssa);
      yyvs = (YY_PDDL_Parser_STYPE *) __ALLOCA_alloca (yystacksize * sizeof (*yyvsp));
      __yy_bcopy ((char *)yyvs1, (char *)yyvs, size * sizeof (*yyvsp));
      __ALLOCA_free(yyvs1,yyvsa);
#ifdef YY_PDDL_Parser_LSP_NEEDED
      yyls = (YY_PDDL_Parser_LTYPE *) __ALLOCA_alloca (yystacksize * sizeof (*yylsp));
      __yy_bcopy ((char *)yyls1, (char *)yyls, size * sizeof (*yylsp));
      __ALLOCA_free(yyls1,yylsa);
#endif
#endif /* no yyoverflow */

      yyssp = yyss + size - 1;
      yyvsp = yyvs + size - 1;
#ifdef YY_PDDL_Parser_LSP_NEEDED
      yylsp = yyls + size - 1;
#endif

#if YY_PDDL_Parser_DEBUG != 0
      if (YY_PDDL_Parser_DEBUG_FLAG)
	fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

#if YY_PDDL_Parser_DEBUG != 0
  if (YY_PDDL_Parser_DEBUG_FLAG)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

  YYGOTO(yybackup);
YYLABEL(yybackup)

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* YYLABEL(yyresume) */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    YYGOTO(yydefault);

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (YY_PDDL_Parser_CHAR == YYEMPTY)
    {
#if YY_PDDL_Parser_DEBUG != 0
      if (YY_PDDL_Parser_DEBUG_FLAG)
	fprintf(stderr, "Reading a token: ");
#endif
      YY_PDDL_Parser_CHAR = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (YY_PDDL_Parser_CHAR <= 0)           /* This means end of input. */
    {
      yychar1 = 0;
      YY_PDDL_Parser_CHAR = YYEOF;                /* Don't call YYLEX any more */

#if YY_PDDL_Parser_DEBUG != 0
      if (YY_PDDL_Parser_DEBUG_FLAG)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      yychar1 = YYTRANSLATE(YY_PDDL_Parser_CHAR);

#if YY_PDDL_Parser_DEBUG != 0
      if (YY_PDDL_Parser_DEBUG_FLAG)
	{
	  fprintf (stderr, "Next token is %d (%s", YY_PDDL_Parser_CHAR, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise meaning
	     of a token, for further debugging info.  */
#ifdef YYPRINT
	  YYPRINT (stderr, YY_PDDL_Parser_CHAR, YY_PDDL_Parser_LVAL);
#endif
	  fprintf (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    YYGOTO(yydefault);

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	YYGOTO(yyerrlab);
      yyn = -yyn;
      YYGOTO(yyreduce);
    }
  else if (yyn == 0)
    YYGOTO(yyerrlab);

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YY_PDDL_Parser_DEBUG != 0
  if (YY_PDDL_Parser_DEBUG_FLAG)
    fprintf(stderr, "Shifting token %d (%s), ", YY_PDDL_Parser_CHAR, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (YY_PDDL_Parser_CHAR != YYEOF)
    YY_PDDL_Parser_CHAR = YYEMPTY;

  *++yyvsp = YY_PDDL_Parser_LVAL;
#ifdef YY_PDDL_Parser_LSP_NEEDED
  *++yylsp = YY_PDDL_Parser_LLOC;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (yyerrstatus) yyerrstatus--;

  yystate = yyn;
  YYGOTO(yynewstate);

/* Do the default action for the current state.  */
YYLABEL(yydefault)

  yyn = yydefact[yystate];
  if (yyn == 0)
    YYGOTO(yyerrlab);

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
YYLABEL(yyreduce)
  yylen = yyr2[yyn];
  if (yylen > 0)
    yyval = yyvsp[1-yylen]; /* implement default value of the action */

#if YY_PDDL_Parser_DEBUG != 0
  if (YY_PDDL_Parser_DEBUG_FLAG)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
	fprintf (stderr, "%s ", yytname[yyrhs[i]]);
      fprintf (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif


/* #line 811 "/usr/local/lib/bison.cc" */
#line 2631 "grammar.cc"

  switch (yyn) {

case 17:
#line 132 "pddl2.y"
{
  domain_name = yyvsp[-1].sym->text;
  if (current_file()) domain_file = strdup(current_file());
;
    break;}
case 18:
#line 139 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 19:
#line 140 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 20:
#line 141 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 21:
#line 142 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 22:
#line 143 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 23:
#line 144 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 24:
#line 145 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 25:
#line 146 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 26:
#line 147 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 27:
#line 151 "pddl2.y"
{ yyval.sym = yyvsp[0].sym; ;
    break;}
case 28:
#line 153 "pddl2.y"
{
  if (warning_level > 0) {
    std::cerr << "warning: redeclaration of action " << yyvsp[0].sym->text << std::endl;
    if (!best_effort) exit(1);
    yyval.sym = tab.gensym(yyvsp[0].sym->text);
  }
  else {
    yyval.sym = yyvsp[0].sym;
  }
;
    break;}
case 37:
#line 191 "pddl2.y"
{
  current_param.clear()
;
    break;}
case 38:
#line 195 "pddl2.y"
{
  PredicateSymbol* p = new PredicateSymbol(yyvsp[-3].sym->text);
  p->param = current_param;
  dom_predicates.append(p);
  clear_context(current_param);
  yyvsp[-3].sym->val = p;
;
    break;}
case 39:
#line 206 "pddl2.y"
{
  set_variable_type(current_param, (TypeSymbol*)yyvsp[0].sym->val);
;
    break;}
case 40:
#line 210 "pddl2.y"
{
  current_type_set.clear();
;
    break;}
case 41:
#line 214 "pddl2.y"
{
  set_variable_type(current_param, current_type_set);
;
    break;}
case 42:
#line 218 "pddl2.y"
{
  set_variable_type(current_param, dom_top_type);
;
    break;}
case 44:
#line 226 "pddl2.y"
{
  yyvsp[0].sym->val = new VariableSymbol(yyvsp[0].sym->text);
  current_param.append((VariableSymbol*)yyvsp[0].sym->val);
  if (trace_print_context) {
    std::cerr << "variable ";
    current_param[current_param.length() - 1]->print(std::cerr);
    std::cerr << " added to context (now "
	      << current_param.length() << " variables)"
	      << std::endl;
  }
;
    break;}
case 45:
#line 238 "pddl2.y"
{
  std::cerr << "error: variable ";
  ((VariableSymbol*)yyvsp[0].sym->val)->print(std::cerr);
  std::cerr << " redeclared" << std::endl;
  exit(255);
;
    break;}
case 46:
#line 245 "pddl2.y"
{
  yyvsp[0].sym->val = new VariableSymbol(yyvsp[0].sym->text);
  current_param.append((VariableSymbol*)yyvsp[0].sym->val);
  if (trace_print_context) {
    std::cerr << "variable ";
    current_param[current_param.length() - 1]->print(std::cerr);
    std::cerr << " added to context (now "
	      << current_param.length() << " variables)"
	      << std::endl;
  }
;
    break;}
case 47:
#line 257 "pddl2.y"
{
  std::cerr << "error: variable ";
  ((VariableSymbol*)yyvsp[0].sym->val)->print(std::cerr);
  std::cerr << " redeclared" << std::endl;
  exit(255);
;
    break;}
case 48:
#line 267 "pddl2.y"
{
  current_type_set.append((TypeSymbol*)yyvsp[0].sym->val);
;
    break;}
case 49:
#line 271 "pddl2.y"
{
  current_type_set.append((TypeSymbol*)yyvsp[0].sym->val);
;
    break;}
case 50:
#line 279 "pddl2.y"
{
  last_n_functions = dom_functions.length();
;
    break;}
case 52:
#line 287 "pddl2.y"
{
  last_n_functions = dom_functions.length();
;
    break;}
case 54:
#line 292 "pddl2.y"
{
  TypeSymbol* t = (TypeSymbol*)yyvsp[0].sym->val;
  for (HSPS::index_type k = last_n_functions; k < dom_functions.length(); k++){
    if (write_info) {
      std::cerr << "info: converting ";
      dom_functions[k]->print(std::cerr);
      std::cerr << " to object function with type " << t->print_name
		<< std::endl;
    }
    HSPS::PDDL_Base::ObjectFunctionSymbol* f =
      new HSPS::PDDL_Base::ObjectFunctionSymbol(dom_functions[k]->print_name);
    f->param = dom_functions[k]->param;
    f->sym_types.assign_value(t, 1);
    dom_object_functions.append(f);
    HSPS::StringTable::Cell* c =
      (HSPS::StringTable::Cell*)tab.find(f->print_name);
    if (c == 0) {
      std::cerr << "very bad error: function "
		<< dom_functions[k]->print_name
		<< " declared but not found in string table!"
		<< std::endl;
      exit(255);
    }
    c->val = f;
  }
  dom_functions.set_length(last_n_functions);
  // last_n_functions = dom_functions.length();
;
    break;}
case 56:
#line 322 "pddl2.y"
{
  for (HSPS::index_type k = last_n_functions; k < dom_functions.length(); k++){
    if (write_info) {
      std::cerr << "info: converting ";
      dom_functions[k]->print(std::cerr);
      std::cerr << " to object function with type ";
      current_type_set.write_type(std::cerr);
      std::cerr << std::endl;
    }
    HSPS::PDDL_Base::ObjectFunctionSymbol* f =
      new HSPS::PDDL_Base::ObjectFunctionSymbol(dom_functions[k]->print_name);
    f->param = dom_functions[k]->param;
    f->sym_types.assign_copy(current_type_set);
    dom_object_functions.append(f);
    HSPS::StringTable::Cell* c =
      (HSPS::StringTable::Cell*)tab.find(f->print_name);
    if (c == 0) {
      std::cerr << "very bad error: function "
		<< dom_functions[k]->print_name
		<< " declared but not found in string table!"
		<< std::endl;
      exit(255);
    }
    c->val = f;
  }
  dom_functions.set_length(last_n_functions);
  // last_n_functions = dom_functions.length();
;
    break;}
case 58:
#line 352 "pddl2.y"
{
  last_n_functions = dom_functions.length();
;
    break;}
case 62:
#line 365 "pddl2.y"
{
  current_param.clear();
;
    break;}
case 63:
#line 369 "pddl2.y"
{
  FunctionSymbol* f = new FunctionSymbol(yyvsp[-3].sym->text);
  f->param = current_param;
  dom_functions.append(f);
  clear_context(current_param);
  yyvsp[-3].sym->val = f;
;
    break;}
case 64:
#line 382 "pddl2.y"
{
  current_type_set.clear();
;
    break;}
case 66:
#line 390 "pddl2.y"
{
  // set_type_type(dom_types, (TypeSymbol*)$4->val);
  for (HSPS::index_type k = 0; k < current_type_set.length(); k++) {
    if (current_type_set[k] == (TypeSymbol*)yyvsp[0].sym->val) {
      log_error("cyclic type declaration");
    }
    else {
      //current_type_set[k]->sym_types.assign_value((TypeSymbol*)$4->val, 1);
      if ((current_type_set[k]->sym_types.length() == 1) &&
	  current_type_set[k]->sym_types[0] == dom_top_type)
	current_type_set[k]->sym_types[0] = ((TypeSymbol*)yyvsp[0].sym->val);
      else 
	current_type_set[k]->sym_types.append((TypeSymbol*)yyvsp[0].sym->val);
    }
  }
  current_type_set.clear();
;
    break;}
case 67:
#line 408 "pddl2.y"
{
  yyvsp[0].sym->val = new TypeSymbol(yyvsp[0].sym->text);
  ((TypeSymbol*)yyvsp[0].sym->val)->sym_types.assign_value(dom_top_type, 1);
  // if (warning_level > 0)
  // std::cerr << "warning: assuming " << $4->text << " - object" << std::endl;
  // ((TypeSymbol*)$4->val)->sym_types.assign_value(dom_top_type, 1);
  // set_type_type(dom_types, (TypeSymbol*)$4->val);
  dom_types.append((TypeSymbol*)yyvsp[0].sym->val);
  for (HSPS::index_type k = 0; k < current_type_set.length(); k++) {
    //current_type_set[k]->sym_types.assign_value((TypeSymbol*)$4->val, 1);
    if ((current_type_set[k]->sym_types.length() == 1) &&
	current_type_set[k]->sym_types[0] == dom_top_type)
      current_type_set[k]->sym_types[0] = ((TypeSymbol*)yyvsp[0].sym->val);
    else 
      current_type_set[k]->sym_types.append((TypeSymbol*)yyvsp[0].sym->val);
  }
  current_type_set.clear();
;
    break;}
case 68:
#line 427 "pddl2.y"
{
  // set_type_type(dom_types, dom_top_type);
  for (HSPS::index_type k = 0; k < current_type_set.length(); k++)
    if (current_type_set[k] != dom_top_type)
      current_type_set[k]->sym_types.assign_value(dom_top_type, 1);
  current_type_set.clear();
;
    break;}
case 70:
#line 439 "pddl2.y"
{
  /* the type is already (implicitly) declared */
  current_type_set.append((TypeSymbol*)yyvsp[0].sym->val);
;
    break;}
case 71:
#line 444 "pddl2.y"
{
  yyvsp[0].sym->val = new TypeSymbol(yyvsp[0].sym->text);
  dom_types.append((TypeSymbol*)yyvsp[0].sym->val);
  current_type_set.append((TypeSymbol*)yyvsp[0].sym->val);
;
    break;}
case 73:
#line 456 "pddl2.y"
{
  current_constant_decl.clear();
;
    break;}
case 75:
#line 461 "pddl2.y"
{
  current_constant_decl.clear();
;
    break;}
case 77:
#line 469 "pddl2.y"
{
  TypeSymbol* t = (TypeSymbol*)yyvsp[0].sym->val;
  for (HSPS::index_type k = 0; k < current_constant_decl.size(); k++) {
    HSPS::bool_vec d(false, current_constant_decl[k]->sym_types.length());
    bool is_d = false;
    for (HSPS::index_type i = 0; i < current_constant_decl[k]->sym_types.length(); i++)
      if (current_constant_decl[k]->sym_types[i]->subtype_or_equal(t))
	is_d = true;
      else if (t->subtype_or_equal(current_constant_decl[k]->sym_types[i]))
	d[i] = true;
    if (!is_d) {
      current_constant_decl[k]->sym_types.remove(d);
      current_constant_decl[k]->sym_types.append(t);
      t->add_element(current_constant_decl[k]);
    }
  }
  current_constant_decl.clear();
  // set_constant_type(dom_constants, (TypeSymbol*)$4->val);
;
    break;}
case 78:
#line 489 "pddl2.y"
{
  for (HSPS::index_type k = 0; k < current_constant_decl.size(); k++) {
    if (current_constant_decl[k]->sym_types.empty()) {
      current_constant_decl[k]->sym_types.append(dom_top_type);
      dom_top_type->add_element(current_constant_decl[k]);
    }
  }
  current_constant_decl.clear();
  // set_constant_type(dom_constants, dom_top_type);
;
    break;}
case 80:
#line 504 "pddl2.y"
{
  yyvsp[0].sym->val = new Symbol(yyvsp[0].sym->text);
  if (problem_name) {
    ((Symbol*)yyvsp[0].sym->val)->defined_in_problem = true;
  }
  dom_constants.append((Symbol*)yyvsp[0].sym->val);
  current_constant_decl.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 81:
#line 513 "pddl2.y"
{
  yyvsp[0].sym->val = new Symbol(yyvsp[0].sym->text);
  if (problem_name) {
    ((Symbol*)yyvsp[0].sym->val)->defined_in_problem = true;
  }
  dom_constants.append((Symbol*)yyvsp[0].sym->val);
  current_constant_decl.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 82:
#line 522 "pddl2.y"
{
  if (warning_level > 0) {
    std::cerr << "warning: redeclaration of constant " << yyvsp[0].sym->text << std::endl;
  }
  current_constant_decl.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 83:
#line 529 "pddl2.y"
{
  if (warning_level > 0) {
    std::cerr << "warning: redeclaration of constant " << yyvsp[0].sym->text << std::endl;
  }
  current_constant_decl.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 88:
#line 550 "pddl2.y"
{
  dom_actions.append(new ActionSymbol(yyvsp[0].sym->text));
;
    break;}
case 89:
#line 554 "pddl2.y"
{
  // post-processing should be done on all actions after the complete
  // domain and problem have been read (calling PDDL_Base::post_process())
  clear_context(current_param);
  yyvsp[-3].sym->val = dom_actions[dom_actions.length() - 1];
;
    break;}
case 90:
#line 564 "pddl2.y"
{
  current_param.clear();
;
    break;}
case 91:
#line 568 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->param = current_param;
;
    break;}
case 98:
#line 578 "pddl2.y"
{
  //std::cerr << "read assoc string: [" << $3 << "]" << std::endl;
  dom_actions[dom_actions.length() - 1]->assoc = yyvsp[0].sval;
;
    break;}
case 100:
#line 587 "pddl2.y"
{
  SetSymbol* ssym = new SetSymbol(yyvsp[0].sym->text);
  yyvsp[0].sym->val = ssym;
  partitions.append(ssym);
  SetName* s = new SetName(ssym);
  current_atom = s;
;
    break;}
case 101:
#line 595 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->part = (SetName*)current_atom;
;
    break;}
case 102:
#line 599 "pddl2.y"
{
  SetName* s = new SetName((SetSymbol*)yyvsp[0].sym->val);
  current_atom = s;
;
    break;}
case 103:
#line 604 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->part = (SetName*)current_atom;
;
    break;}
case 115:
#line 635 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->num_pre.append(yyvsp[0].rel);
;
    break;}
case 116:
#line 642 "pddl2.y"
{
  yyval.tkw = PDDL_Base::md_start;
;
    break;}
case 117:
#line 646 "pddl2.y"
{
  yyval.tkw = PDDL_Base::md_end;
;
    break;}
case 118:
#line 650 "pddl2.y"
{
  yyval.tkw = PDDL_Base::md_all;
;
    break;}
case 119:
#line 657 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 120:
#line 661 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  dom_actions[dom_actions.length() - 1]->pos_pre.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->pos_pre = true;
;
    break;}
case 121:
#line 667 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-2].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-1].sym->val;
  dom_actions[dom_actions.length() - 1]->pos_pre.append(eq_atom);
;
    break;}
case 122:
#line 675 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, yyvsp[-2].tkw);
;
    break;}
case 123:
#line 679 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  dom_actions[dom_actions.length() - 1]->pos_pre.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->pos_pre = true;
;
    break;}
case 124:
#line 686 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred, yyvsp[-6].tkw);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-3].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-2].sym->val;
  dom_actions[dom_actions.length() - 1]->pos_pre.append(eq_atom);
;
    break;}
case 125:
#line 697 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 126:
#line 701 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  dom_actions[dom_actions.length() - 1]->neg_pre.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 127:
#line 707 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-3].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-2].sym->val;
  dom_actions[dom_actions.length() - 1]->neg_pre.append(eq_atom);
;
    break;}
case 128:
#line 715 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, yyvsp[-4].tkw);
;
    break;}
case 129:
#line 719 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  dom_actions[dom_actions.length() - 1]->neg_pre.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 130:
#line 726 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred, yyvsp[-9].tkw);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-4].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-3].sym->val;
  dom_actions[dom_actions.length() - 1]->neg_pre.append(eq_atom);
;
    break;}
case 131:
#line 781 "pddl2.y"
{
  if (current_atom != 0) {
    current_atom->param.append((Symbol*)yyvsp[0].sym->val);
  }
;
    break;}
case 133:
#line 791 "pddl2.y"
{
  if (current_atom != 0) {
    current_atom->param.append((Symbol*)yyvsp[0].sym->val);
  }
;
    break;}
case 135:
#line 801 "pddl2.y"
{
  if (yyvsp[0].sym->val == 0) {
    log_error("undeclared variable in atom argument");
  }
  yyval.sym = yyvsp[0].sym;
;
    break;}
case 136:
#line 808 "pddl2.y"
{
  yyval.sym = yyvsp[0].sym;
;
    break;}
case 137:
#line 812 "pddl2.y"
{
  yyvsp[0].sym->val = new Symbol(yyvsp[0].sym->text);
  dom_constants.append((Symbol*)yyvsp[0].sym->val);
  yyval.sym = yyvsp[0].sym;
;
    break;}
case 138:
#line 821 "pddl2.y"
{
  current_atom_stack.append(current_atom);
  current_atom = new FTerm((ObjectFunctionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 139:
#line 826 "pddl2.y"
{
  ObjectFunctionSymbol* f = (ObjectFunctionSymbol*)yyvsp[-3].sym->val;
  VariableSymbol* v =
    (VariableSymbol*)gensym(sym_variable, "?omsk", f->sym_types);
  v->binding = (FTerm*)current_atom;
  assert(current_atom_stack.length() > 0);
  current_atom = current_atom_stack[current_atom_stack.length() - 1];
  current_atom_stack.dec_length();
  HSPS::StringTable::Cell* c =
    (HSPS::StringTable::Cell*)tab.find(v->print_name);
  if (c == 0) {
    std::cerr << "very bad error: omsk symbol " << v->print_name
	      << " generated but not found in string table!"
	      << std::endl;
    exit(255);
  }
  yyval.sym = c;
;
    break;}
case 140:
#line 847 "pddl2.y"
{
  yyval.sym = yyvsp[0].sym;
;
    break;}
case 141:
#line 851 "pddl2.y"
{
  yyval.sym = yyvsp[0].sym;
;
    break;}
case 142:
#line 857 "pddl2.y"
{
  yyval.rel = new Relation(yyvsp[-3].rkw, yyvsp[-2].exp, yyvsp[-1].exp);
  yyvsp[-2].exp->mark_functions_in_condition();
  yyvsp[-1].exp->mark_functions_in_condition();
;
    break;}
case 143:
#line 864 "pddl2.y"
{
  yyval.rel = new Relation(yyvsp[-4].rkw, yyvsp[-6].tkw, yyvsp[-3].exp, yyvsp[-2].exp);
  yyvsp[-3].exp->mark_functions_in_condition();
  yyvsp[-2].exp->mark_functions_in_condition();
;
    break;}
case 144:
#line 873 "pddl2.y"
{
  yyval.rkw = rel_greater;
;
    break;}
case 145:
#line 877 "pddl2.y"
{
  yyval.rkw = rel_greater_equal;
;
    break;}
case 146:
#line 881 "pddl2.y"
{
  yyval.rkw = rel_less;
;
    break;}
case 147:
#line 885 "pddl2.y"
{
  yyval.rkw = rel_less_equal;
;
    break;}
case 148:
#line 889 "pddl2.y"
{
  yyval.rkw = rel_equal;
;
    break;}
case 149:
#line 896 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_sub, new ConstantExpression(0), yyvsp[-1].exp);
;
    break;}
case 150:
#line 900 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_add, yyvsp[-2].exp, yyvsp[-1].exp);
;
    break;}
case 151:
#line 904 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_sub, yyvsp[-2].exp, yyvsp[-1].exp);
;
    break;}
case 152:
#line 908 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_mul, yyvsp[-2].exp, yyvsp[-1].exp);
;
    break;}
case 153:
#line 912 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_div, yyvsp[-2].exp, yyvsp[-1].exp);
;
    break;}
case 154:
#line 916 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_min, yyvsp[-2].exp, yyvsp[-1].exp);
;
    break;}
case 155:
#line 920 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_max, yyvsp[-2].exp, yyvsp[-1].exp);
;
    break;}
case 156:
#line 924 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_div, yyvsp[-2].exp, yyvsp[0].exp);
;
    break;}
case 157:
#line 928 "pddl2.y"
{
  yyval.exp = new ConstantExpression(yyvsp[0].ival);
;
    break;}
case 158:
#line 932 "pddl2.y"
{
  yyval.exp = new ConstantExpression(NN_TO_N(yyvsp[0].rval));
;
    break;}
case 159:
#line 936 "pddl2.y"
{
  yyval.exp = new TimeExpression();
;
    break;}
case 160:
#line 940 "pddl2.y"
{
  yyval.exp = new PreferenceExpression((Symbol*)yyvsp[-1].sym->val);
;
    break;}
case 161:
#line 944 "pddl2.y"
{
  yyval.exp = yyvsp[0].exp;
;
    break;}
case 162:
#line 951 "pddl2.y"
{
  yyval.exp = yyvsp[0].exp;
;
    break;}
case 163:
#line 955 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_add, yyvsp[-1].exp, yyvsp[0].exp);
;
    break;}
case 164:
#line 962 "pddl2.y"
{
  yyval.exp = yyvsp[0].exp;
;
    break;}
case 165:
#line 966 "pddl2.y"
{
  yyval.exp = new BinaryExpression(exp_mul, yyvsp[-1].exp, yyvsp[0].exp);
;
    break;}
case 166:
#line 973 "pddl2.y"
{
  yyval.exp = new FunctionExpression((FunctionSymbol*)yyvsp[-2].sym->val, yyvsp[-1].lst);
;
    break;}
case 167:
#line 977 "pddl2.y"
{
  yyval.exp = new FunctionExpression((FunctionSymbol*)yyvsp[0].sym->val, 0);
;
    break;}
case 168:
#line 984 "pddl2.y"
{
  yyval.lst = new ListExpression((VariableSymbol*)yyvsp[-1].sym->val, yyvsp[0].lst);
;
    break;}
case 169:
#line 988 "pddl2.y"
{
  yyval.lst = new ListExpression((Symbol*)yyvsp[-1].sym->val, yyvsp[0].lst);
;
    break;}
case 170:
#line 992 "pddl2.y"
{
  yyval.lst = 0;
;
    break;}
case 171:
#line 999 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  current_context = new SetOf();
;
    break;}
case 172:
#line 1004 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 173:
#line 1008 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  SetOf* s = (SetOf*)current_context;
  s->pos_atoms.append((Atom*)current_atom);
  dom_actions[dom_actions.length() - 1]->set_pre.append(s);
  ((Atom*)current_atom)->pred->pos_pre = true;
  assert(stored_n_param.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = 0;
;
    break;}
case 174:
#line 1023 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  current_context = new SetOf();
;
    break;}
case 175:
#line 1028 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    s->param.append(current_param[k]);
  }
  dom_actions[dom_actions.length() - 1]->set_pre.append(s);
;
    break;}
case 176:
#line 1037 "pddl2.y"
{
  assert(stored_n_param.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = 0;
;
    break;}
case 177:
#line 1047 "pddl2.y"
{
  current_context = new SetOf();
;
    break;}
case 178:
#line 1051 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  dom_actions[dom_actions.length() - 1]->set_pre.append(s);
  current_context = 0;
;
    break;}
case 179:
#line 1057 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  current_context = new SetOf();
;
    break;}
case 180:
#line 1062 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 181:
#line 1066 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  SetOf* s = (SetOf*)current_context;
  s->pos_atoms.append((Atom*)current_atom);
  s->set_mode(yyvsp[-12].tkw);
  dom_actions[dom_actions.length() - 1]->set_pre.append(s);
  ((Atom*)current_atom)->pred->pos_pre = true;
  assert(stored_n_param.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = 0;
;
    break;}
case 182:
#line 1082 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  current_context = new SetOf();
;
    break;}
case 183:
#line 1087 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    s->param.append(current_param[k]);
  }
  dom_actions[dom_actions.length() - 1]->set_pre.append(s);
;
    break;}
case 184:
#line 1096 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  s->set_mode(yyvsp[-10].tkw);
  assert(stored_n_param.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = 0;
;
    break;}
case 185:
#line 1108 "pddl2.y"
{
  current_context = new SetOf();
;
    break;}
case 186:
#line 1112 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  s->set_mode(yyvsp[-7].tkw);
  dom_actions[dom_actions.length() - 1]->set_pre.append(s);
  current_context = 0;
;
    break;}
case 195:
#line 1142 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 196:
#line 1146 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  SetOf* s = (SetOf*)current_context;
  s->pos_atoms.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->pos_pre = true;
;
    break;}
case 197:
#line 1153 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 198:
#line 1157 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  SetOf* s = (SetOf*)current_context;
  s->neg_atoms.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 199:
#line 1164 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-2].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-1].sym->val;
  SetOf* s = (SetOf*)current_context;
  s->pos_atoms.append(eq_atom);
;
    break;}
case 200:
#line 1173 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-3].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-2].sym->val;
  SetOf* s = (SetOf*)current_context;
  s->neg_atoms.append(eq_atom);
;
    break;}
case 201:
#line 1211 "pddl2.y"
{
  current_context = new SetOf();
;
    break;}
case 202:
#line 1215 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  dom_actions[dom_actions.length() - 1]->dis_pre.append(s);
  current_context = 0;
;
    break;}
case 204:
#line 1224 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  current_context = new SetOf();
;
    break;}
case 205:
#line 1229 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    s->param.append(current_param[k]);
  }
  dom_actions[dom_actions.length() - 1]->dis_pre.append(s);
;
    break;}
case 206:
#line 1238 "pddl2.y"
{
  assert(stored_n_param.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = 0;
;
    break;}
case 212:
#line 1262 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 213:
#line 1266 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  ((SetOf*)current_context)->pos_atoms.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->pos_pre = true;
;
    break;}
case 214:
#line 1272 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-2].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-1].sym->val;
  ((SetOf*)current_context)->pos_atoms.append(eq_atom);
;
    break;}
case 215:
#line 1280 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 216:
#line 1284 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  ((SetOf*)current_context)->neg_atoms.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 217:
#line 1290 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-3].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-2].sym->val;
  ((SetOf*)current_context)->neg_atoms.append(eq_atom);
;
    break;}
case 223:
#line 1319 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  current_context = new SetOf();
;
    break;}
case 224:
#line 1324 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    s->param.append(current_param[k]);
  }
;
    break;}
case 225:
#line 1332 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->set_eff.append((SetOf*)current_context);
  assert(stored_n_param.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = 0;
;
    break;}
case 226:
#line 1343 "pddl2.y"
{
  current_context = new SetOf();
;
    break;}
case 227:
#line 1347 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->set_eff.append((SetOf*)current_context);
  current_context = 0;
;
    break;}
case 228:
#line 1352 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  current_context = new SetOf();
;
    break;}
case 229:
#line 1357 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    s->param.append(current_param[k]);
  }
;
    break;}
case 230:
#line 1365 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  s->set_mode(yyvsp[-10].tkw);
  dom_actions[dom_actions.length() - 1]->set_eff.append(s);
  assert(stored_n_param.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = 0;
;
    break;}
case 231:
#line 1378 "pddl2.y"
{
  current_context = new SetOf();
;
    break;}
case 232:
#line 1382 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  s->set_mode(yyvsp[-6].tkw);
  dom_actions[dom_actions.length() - 1]->set_eff.append(s);
  current_context = 0;
;
    break;}
case 245:
#line 1419 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
  ((PredicateSymbol*)yyvsp[0].sym->val)->modded = true;
;
    break;}
case 246:
#line 1424 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  if (current_context != 0) {
    ((SetOf*)current_context)->pos_atoms.append((Atom*)current_atom);
  }
  else {
    dom_actions[dom_actions.length() - 1]->adds.append((Atom*)current_atom);
  }
;
    break;}
case 247:
#line 1434 "pddl2.y"
{
  current_atom = new FTerm((ObjectFunctionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 248:
#line 1438 "pddl2.y"
{
  FTerm* ft = (FTerm*)current_atom;
  ft->fun->modded = true;
  VariableSymbol* v =
    (VariableSymbol*)gensym(sym_variable,"?omsk",ft->fun->sym_types);
  v->binding = ft;
  Atom* a = new Atom(dom_assign_pred);
  a->param.set_length(2);
  a->param[0] = v;
  a->param[1] = (Symbol*)yyvsp[-1].sym->val;
  if (current_context != 0) {
    if ((warning_level > 0) || !best_effort) {
      std::cerr << "warning: object function assignment ";
      a->print(std::cerr);
      std::cerr << " in quantified/conditional effect ignored"
		<< std::endl;
      std::cerr << "context: ";
      current_context->print(std::cerr);
      std::cerr << std::endl;
    }
    if (!best_effort) exit(1);
  }
  else {
    dom_actions[dom_actions.length() - 1]->adds.append(a);
  }
;
    break;}
case 249:
#line 1465 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, yyvsp[-2].tkw);
  ((PredicateSymbol*)yyvsp[0].sym->val)->modded = true;
;
    break;}
case 250:
#line 1470 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  if (current_context != 0) {
    ((SetOf*)current_context)->pos_atoms.append((Atom*)current_atom);
  }
  else {
    dom_actions[dom_actions.length() - 1]->adds.append((Atom*)current_atom);
  }
;
    break;}
case 251:
#line 1483 "pddl2.y"
{
  yyval.sym = yyvsp[0].sym;
;
    break;}
case 252:
#line 1487 "pddl2.y"
{
  yyval.sym = (HSPS::StringTable::Cell*)tab.find("undefined");
  assert(yyval.sym != 0);
;
    break;}
case 253:
#line 1494 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
  ((PredicateSymbol*)yyvsp[0].sym->val)->modded = true;
;
    break;}
case 254:
#line 1499 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  if (current_context != 0) {
    ((SetOf*)current_context)->neg_atoms.append((Atom*)current_atom);
  }
  else {
    dom_actions[dom_actions.length() - 1]->dels.append((Atom*)current_atom);
  }
;
    break;}
case 255:
#line 1509 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, yyvsp[-4].tkw);
  ((PredicateSymbol*)yyvsp[0].sym->val)->modded = true;
;
    break;}
case 256:
#line 1514 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  if (current_context != 0) {
    ((SetOf*)current_context)->neg_atoms.append((Atom*)current_atom);
  }
  else {
    dom_actions[dom_actions.length() - 1]->dels.append((Atom*)current_atom);
  }
;
    break;}
case 257:
#line 1527 "pddl2.y"
{
  current_atom = new FChangeAtom((FunctionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 258:
#line 1531 "pddl2.y"
{
  ((FChangeAtom*)current_atom)->val = yyvsp[-1].exp;
  ((FChangeAtom*)current_atom)->fun->modified = true;
  if (!yyvsp[-1].exp->is_integral()) ((FChangeAtom*)current_atom)->fun->integral = false;
  dom_actions[dom_actions.length() - 1]->incs.append((FChangeAtom*)current_atom);
;
    break;}
case 259:
#line 1538 "pddl2.y"
{
  current_atom = new FChangeAtom((FunctionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 260:
#line 1542 "pddl2.y"
{
  ((FChangeAtom*)current_atom)->val = yyvsp[-1].exp;
  ((FChangeAtom*)current_atom)->fun->modified = true;
  if (!yyvsp[-1].exp->is_integral()) ((FChangeAtom*)current_atom)->fun->integral = false;
  dom_actions[dom_actions.length() - 1]->decs.append((FChangeAtom*)current_atom);
;
    break;}
case 261:
#line 1549 "pddl2.y"
{
  current_atom = new FChangeAtom((FunctionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 262:
#line 1553 "pddl2.y"
{
  ((FChangeAtom*)current_atom)->val = yyvsp[-1].exp;
  ((FChangeAtom*)current_atom)->fun->modified = true;
  if (!yyvsp[-1].exp->is_integral()) ((FChangeAtom*)current_atom)->fun->integral = false;
  dom_actions[dom_actions.length() - 1]->fass.append((FChangeAtom*)current_atom);
;
    break;}
case 263:
#line 1560 "pddl2.y"
{
  current_atom = new FChangeAtom((FunctionSymbol*)yyvsp[0].sym->val, yyvsp[-4].tkw);
;
    break;}
case 264:
#line 1564 "pddl2.y"
{
  ((FChangeAtom*)current_atom)->val = yyvsp[-2].exp;
  ((FChangeAtom*)current_atom)->fun->modified = true;
  if (!yyvsp[-2].exp->is_integral()) ((FChangeAtom*)current_atom)->fun->integral = false;
  dom_actions[dom_actions.length() - 1]->incs.append((FChangeAtom*)current_atom);
;
    break;}
case 265:
#line 1571 "pddl2.y"
{
  current_atom = new FChangeAtom((FunctionSymbol*)yyvsp[0].sym->val, yyvsp[-4].tkw);
;
    break;}
case 266:
#line 1575 "pddl2.y"
{
  ((FChangeAtom*)current_atom)->val = yyvsp[-2].exp;
  ((FChangeAtom*)current_atom)->fun->modified = true;
  if (!yyvsp[-2].exp->is_integral()) ((FChangeAtom*)current_atom)->fun->integral = false;
  dom_actions[dom_actions.length() - 1]->decs.append((FChangeAtom*)current_atom);
;
    break;}
case 267:
#line 1582 "pddl2.y"
{
  current_atom = new FChangeAtom((FunctionSymbol*)yyvsp[0].sym->val, yyvsp[-4].tkw);
;
    break;}
case 268:
#line 1586 "pddl2.y"
{
  ((FChangeAtom*)current_atom)->val = yyvsp[-2].exp;
  ((FChangeAtom*)current_atom)->fun->modified = true;
  if (!yyvsp[-2].exp->is_integral()) ((FChangeAtom*)current_atom)->fun->integral = false;
  dom_actions[dom_actions.length() - 1]->fass.append((FChangeAtom*)current_atom);
;
    break;}
case 276:
#line 1612 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->dmin = yyvsp[-1].exp;
  dom_actions[dom_actions.length() - 1]->dmax = yyvsp[-1].exp;
;
    break;}
case 277:
#line 1617 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->dmin = yyvsp[0].exp;
  dom_actions[dom_actions.length() - 1]->dmax = yyvsp[0].exp;
;
    break;}
case 278:
#line 1625 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->dmax = yyvsp[-1].exp;
;
    break;}
case 281:
#line 1637 "pddl2.y"
{
  dom_actions[dom_actions.length() - 1]->dmin = yyvsp[-1].exp;
;
    break;}
case 284:
#line 1649 "pddl2.y"
{
  current_context = new SequentialTaskNet();
  stored_n_param.append(current_param.length());
  if (trace_print_context) {
    std::cerr << "pushed context (" << current_param.length() << " variables)"
	      << std::endl;
  }
;
    break;}
case 285:
#line 1659 "pddl2.y"
{
  SequentialTaskNet* n = (SequentialTaskNet*)current_context;
  n->abs_act = dom_actions[dom_actions.length() - 1];
  dom_actions[dom_actions.length() - 1]->exps.append(n);
  if (trace_print_context) {
    std::cerr << "poping context from "
	      << current_param.length()
	      << " to "
	      << stored_n_param[stored_n_param.length() - 1]
	      << " variables..." << std::endl;
  }
  assert(stored_n_param.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = 0;
;
    break;}
case 288:
#line 1687 "pddl2.y"
{
  current_atom = new Reference((ActionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 289:
#line 1691 "pddl2.y"
{
  Reference* ref = (Reference*)current_atom;
  // ActionSymbol* act = (ActionSymbol*)$2->val;
  // act->refs[act->n_refs] = ref;
  // act->n_refs += 1;
  SequentialTaskNet* task_net = (SequentialTaskNet*)current_context;
  task_net->tasks.append(ref);
;
    break;}
case 290:
#line 1703 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
  ((PredicateSymbol*)yyvsp[0].sym->val)->modded = true;
;
    break;}
case 291:
#line 1709 "pddl2.y"
{
  Axiom* new_axiom = new Axiom((Atom*)current_atom, 0);
  assert(stored_n_param.length() > 0);
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++)
    new_axiom->param.append(current_param[k]);
  dom_axioms.append(new_axiom);
;
    break;}
case 292:
#line 1718 "pddl2.y"
{
  Axiom* new_axiom = dom_axioms[dom_axioms.length() - 1];
  new_axiom->body = yyvsp[-1].ff;
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
;
    break;}
case 293:
#line 1730 "pddl2.y"
{
  set_variable_type(current_param, (TypeSymbol*)yyvsp[0].sym->val);
;
    break;}
case 294:
#line 1734 "pddl2.y"
{
  current_type_set.clear();
;
    break;}
case 295:
#line 1738 "pddl2.y"
{
  set_variable_type(current_param, current_type_set);
;
    break;}
case 296:
#line 1742 "pddl2.y"
{
  set_variable_type(current_param, dom_top_type);
;
    break;}
case 298:
#line 1750 "pddl2.y"
{
  yyvsp[0].sym->val = new VariableSymbol(yyvsp[0].sym->text);
  current_param.append((VariableSymbol*)yyvsp[0].sym->val);
  if (trace_print_context) {
    std::cerr << "variable ";
    current_param[current_param.length() - 1]->print(std::cerr);
    std::cerr << " added to context (now "
	      << current_param.length() << " variables)"
	      << std::endl;
  }
  assert(current_atom != 0);
  current_atom->param.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 299:
#line 1764 "pddl2.y"
{
  assert(current_atom != 0);
  current_atom->param.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 300:
#line 1769 "pddl2.y"
{
  assert(current_atom != 0);
  current_atom->param.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 301:
#line 1774 "pddl2.y"
{
  yyvsp[0].sym->val = new Symbol(yyvsp[0].sym->text);
  dom_constants.append((Symbol*)yyvsp[0].sym->val);
  assert(current_atom != 0);
  current_atom->param.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 302:
#line 1781 "pddl2.y"
{
  yyvsp[0].sym->val = new VariableSymbol(yyvsp[0].sym->text);
  current_param.append((VariableSymbol*)yyvsp[0].sym->val);
  if (trace_print_context) {
    std::cerr << "variable ";
    current_param[current_param.length() - 1]->print(std::cerr);
    std::cerr << " added to context (now "
	      << current_param.length() << " variables)"
	      << std::endl;
  }
  assert(current_atom != 0);
  current_atom->param.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 303:
#line 1795 "pddl2.y"
{
  assert(current_atom != 0);
  current_atom->param.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 304:
#line 1800 "pddl2.y"
{
  assert(current_atom != 0);
  current_atom->param.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 305:
#line 1805 "pddl2.y"
{
  yyvsp[0].sym->val = new Symbol(yyvsp[0].sym->text);
  dom_constants.append((Symbol*)yyvsp[0].sym->val);
  assert(current_atom != 0);
  current_atom->param.append((Symbol*)yyvsp[0].sym->val);
;
    break;}
case 307:
#line 1822 "pddl2.y"
{
  problem_name = yyvsp[-1].sym->text;
  if (current_file()) problem_file = strdup(current_file());
;
    break;}
case 324:
#line 1855 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_init);
;
    break;}
case 325:
#line 1859 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  PredicateSymbol* p = (PredicateSymbol*)yyvsp[-3].sym->val;
  if (p->param.length() != current_atom->param.length()) {
    log_error("wrong number of arguments for predicate in (:init ...");
  }
  ((Atom*)current_atom)->insert(p->init);
  current_atom->at_time = 0;
  dom_init.append((Atom*)current_atom);
;
    break;}
case 326:
#line 1870 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_init);
;
    break;}
case 327:
#line 1874 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  PredicateSymbol* p = (PredicateSymbol*)yyvsp[-4].sym->val;
  if (p->param.length() != current_atom->param.length()) {
    log_error("wrong number of arguments for predicate in (:init ...");
  }
  ((Atom*)current_atom)->insert(p->init);
  current_atom->at_time = yyvsp[-6].rval;
  p->added = true;
  p->modded = true;
  dom_pos_til.append((Atom*)current_atom);
;
    break;}
case 328:
#line 1887 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_init);
;
    break;}
case 329:
#line 1891 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  PredicateSymbol* p = (PredicateSymbol*)yyvsp[-5].sym->val;
  if (p->param.length() != current_atom->param.length()) {
    log_error("wrong number of arguments for predicate in (:init ...");
  }
  ((Atom*)current_atom)->insert(p->init);
  current_atom->at_time = yyvsp[-9].rval;
  p->deleted = true;
  p->modded = true;
  dom_neg_til.append((Atom*)current_atom);
;
    break;}
case 330:
#line 1907 "pddl2.y"
{
  current_atom = new OInitAtom((ObjectFunctionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 331:
#line 1911 "pddl2.y"
{
  ObjectFunctionSymbol* f = (ObjectFunctionSymbol*)yyvsp[-5].sym->val;
  if (f->param.length() != current_atom->param.length()) {
    log_error("wrong number of arguments for object function in (:init ...");
  }
  ((OInitAtom*)current_atom)->val = (Symbol*)yyvsp[-1].sym->val;
  current_atom->at_time = 0;
  ((OInitAtom*)current_atom)->insert(f->init);
  dom_obj_init.append((OInitAtom*)current_atom);
;
    break;}
case 332:
#line 1924 "pddl2.y"
{
  current_atom = new FInitAtom((FunctionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 333:
#line 1928 "pddl2.y"
{
  FunctionSymbol* f = (FunctionSymbol*)yyvsp[-5].sym->val;
  if (f->param.length() != current_atom->param.length()) {
    log_error("wrong number of arguments for function in (:init ...");
  }
  ((FInitAtom*)current_atom)->val = NN_TO_N(yyvsp[-1].rval);
  if (!INTEGRAL(((FInitAtom*)current_atom)->val))
    ((FInitAtom*)current_atom)->fun->integral = false;
  ((FInitAtom*)current_atom)->insert(f->init);
  current_atom->at_time = 0;
  dom_fun_init.append((FInitAtom*)current_atom);
;
    break;}
case 334:
#line 1941 "pddl2.y"
{
  FunctionSymbol* f = (FunctionSymbol*)yyvsp[-2].sym->val;
  current_atom = new FInitAtom((FunctionSymbol*)yyvsp[-2].sym->val);
  if (f->param.length() != 0) {
    log_error("wrong number of arguments for function in (:init ...");
  }
  ((FInitAtom*)current_atom)->val = NN_TO_N(yyvsp[-1].rval);
  if (!INTEGRAL(((FInitAtom*)current_atom)->val))
    ((FInitAtom*)current_atom)->fun->integral = false;
  current_atom->at_time = 0;
  ((FInitAtom*)current_atom)->insert(f->init);
  dom_fun_init.append((FInitAtom*)current_atom);
;
    break;}
case 341:
#line 1970 "pddl2.y"
{
  dom_goals.append(yyvsp[0].goal);
;
    break;}
case 342:
#line 1974 "pddl2.y"
{
  Symbol* name = new Symbol(sym_preference, yyvsp[-2].sym->text);
  yyvsp[-2].sym->val = name;
  dom_preferences.append(new Preference(name, yyvsp[-1].goal));
;
    break;}
case 343:
#line 1980 "pddl2.y"
{
  current_goal.append(new DisjunctiveGoal());
;
    break;}
case 344:
#line 1984 "pddl2.y"
{
  assert(current_goal.length() > 0);
  dom_goals.append(current_goal[current_goal.length() - 1]);
  current_goal.dec_length();
;
    break;}
case 345:
#line 1990 "pddl2.y"
{
  stored_n_param.append(current_param.length());
;
    break;}
case 346:
#line 1994 "pddl2.y"
{
  Symbol* name = new Symbol(sym_preference, yyvsp[-3].sym->text);
  yyvsp[-3].sym->val = name;
  Preference* qp = new Preference(name, yyvsp[-2].goal);
  assert(stored_n_param.length() > 0);
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    qp->param.append(current_param[k]);
  }
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  dom_preferences.append(qp);
;
    break;}
case 347:
#line 2014 "pddl2.y"
{
  yyval.goal = yyvsp[0].goal;
;
    break;}
case 348:
#line 2018 "pddl2.y"
{
  current_goal.append(new ConjunctiveGoal());
;
    break;}
case 349:
#line 2022 "pddl2.y"
{
  assert(current_goal.length() > 0);
  yyval.goal = current_goal[current_goal.length() - 1];
  current_goal.dec_length();
;
    break;}
case 350:
#line 2028 "pddl2.y"
{
  current_goal.append(new DisjunctiveGoal());
;
    break;}
case 351:
#line 2032 "pddl2.y"
{
  assert(current_goal.length() > 0);
  yyval.goal = current_goal[current_goal.length() - 1];
  current_goal.dec_length();
;
    break;}
case 352:
#line 2038 "pddl2.y"
{
  assert((yyvsp[-2].goal->g_class == goal_pos_atom) || (yyvsp[-2].goal->g_class == goal_neg_atom));
  DisjunctiveGoal* dg = new DisjunctiveGoal();
  if (yyvsp[-2].goal->g_class == goal_neg_atom) {
    yyvsp[-2].goal->g_class = goal_pos_atom;
    ((AtomicGoal*)yyvsp[-2].goal)->atom->pred->pos_pre = true;
  }
  else {
    yyvsp[-2].goal->g_class = goal_neg_atom;
    ((AtomicGoal*)yyvsp[-2].goal)->atom->pred->neg_pre = true;
  }
  dg->goals.append(yyvsp[-2].goal);
  dg->goals.append(yyvsp[-1].goal);
  yyval.goal = dg;
;
    break;}
case 353:
#line 2054 "pddl2.y"
{
  current_goal.append(new DisjunctiveGoal());
;
    break;}
case 354:
#line 2058 "pddl2.y"
{
  assert(current_goal.length() > 0);
  DisjunctiveGoal* dg =
    (DisjunctiveGoal*)current_goal[current_goal.length() - 1];
  for (HSPS::index_type k = 0; k < dg->goals.size(); k++) {
    assert((dg->goals[k]->g_class == goal_pos_atom) ||
	   (dg->goals[k]->g_class == goal_neg_atom));
    if (dg->goals[k]->g_class == goal_neg_atom) {
      dg->goals[k]->g_class = goal_pos_atom;
      ((AtomicGoal*)dg->goals[k])->atom->pred->pos_pre = true;
    }
    else {
      dg->goals[k]->g_class = goal_neg_atom;
      ((AtomicGoal*)dg->goals[k])->atom->pred->neg_pre = true;
    }
  }
;
    break;}
case 355:
#line 2076 "pddl2.y"
{
  assert(current_goal.length() > 0);
  DisjunctiveGoal* dg =
    (DisjunctiveGoal*)current_goal[current_goal.length() - 1];
  dg->goals.append(yyvsp[-1].goal);
  yyval.goal = current_goal[current_goal.length() - 1];
  current_goal.dec_length();
;
    break;}
case 356:
#line 2085 "pddl2.y"
{
  stored_n_param.append(current_param.length());
;
    break;}
case 357:
#line 2089 "pddl2.y"
{
  QuantifiedGoal* qg = new QuantifiedGoal(goal_forall);
  qg->goal = yyvsp[-1].goal;
  assert(stored_n_param.length() > 0);
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    qg->param.append(current_param[k]);
  }
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  yyval.goal = qg;
;
    break;}
case 358:
#line 2105 "pddl2.y"
{
  stored_n_param.append(current_param.length());
;
    break;}
case 359:
#line 2109 "pddl2.y"
{
  QuantifiedGoal* qg = new QuantifiedGoal(goal_exists);
  qg->goal = yyvsp[-1].goal;
  assert(stored_n_param.length() > 0);
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    qg->param.append(current_param[k]);
  }
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  yyval.goal = qg;
;
    break;}
case 360:
#line 2128 "pddl2.y"
{
  assert(current_goal.length() > 0);
  ConjunctiveGoal* cg = current_goal[current_goal.length() - 1];
  assert(cg != 0);
  cg->goals.append(yyvsp[-1].goal);
;
    break;}
case 362:
#line 2139 "pddl2.y"
{
  assert(current_goal.length() > 0);
  ConjunctiveGoal* cg = current_goal[current_goal.length() - 1];
  assert(cg != 0);
  cg->goals.append(yyvsp[-1].goal);
;
    break;}
case 364:
#line 2150 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 365:
#line 2154 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  yyval.goal = new AtomicGoal((Atom*)current_atom, false);
  ((Atom*)current_atom)->pred->pos_pre = true;
  ((Atom*)current_atom)->insert(((Atom*)current_atom)->pred->pos_goal);
;
    break;}
case 366:
#line 2161 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 367:
#line 2165 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  yyval.goal = new AtomicGoal((Atom*)current_atom, true);
  ((Atom*)current_atom)->pred->neg_pre = true;
  ((Atom*)current_atom)->insert(((Atom*)current_atom)->pred->neg_goal);
;
    break;}
case 368:
#line 2172 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-2].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-1].sym->val;
  yyval.goal = new AtomicGoal(eq_atom, false);
;
    break;}
case 369:
#line 2180 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-3].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-2].sym->val;
  yyval.goal = new AtomicGoal(eq_atom, true);
;
    break;}
case 370:
#line 2190 "pddl2.y"
{
  yyval.goal = yyvsp[0].goal;
;
    break;}
case 371:
#line 2194 "pddl2.y"
{
  yyval.goal = new NumericGoal(yyvsp[0].rel);
;
    break;}
case 372:
#line 2198 "pddl2.y"
{
  yyval.goal = new SimpleSequenceGoal(goal_always, yyvsp[-1].goal);
;
    break;}
case 373:
#line 2202 "pddl2.y"
{
  yyval.goal = new SimpleSequenceGoal(goal_sometime, yyvsp[-1].goal);
;
    break;}
case 374:
#line 2206 "pddl2.y"
{
  yyval.goal = new SimpleSequenceGoal(goal_at_most_once, yyvsp[-1].goal);
;
    break;}
case 375:
#line 2210 "pddl2.y"
{
  yyval.goal = new TriggeredSequenceGoal(goal_sometime_before, yyvsp[-2].goal, yyvsp[-1].goal);
;
    break;}
case 376:
#line 2214 "pddl2.y"
{
  yyval.goal = new TriggeredSequenceGoal(goal_sometime_after, yyvsp[-2].goal, yyvsp[-1].goal);
;
    break;}
case 377:
#line 2218 "pddl2.y"
{
  yyval.goal = new DeadlineGoal(yyvsp[-2].rval, yyvsp[-1].goal);
;
    break;}
case 378:
#line 2222 "pddl2.y"
{
  yyval.goal = new TriggeredDeadlineGoal(yyvsp[-2].goal, yyvsp[-3].rval, yyvsp[-1].goal);
;
    break;}
case 380:
#line 2263 "pddl2.y"
{
  if (metric_type != metric_none) {
    if (warning_level > 0) {
      std::cerr << "warning: multiple :metric expressions - overwriting previous definition" << std::endl;
    }
  }
  metric_type = metric_minimize;
;
    break;}
case 381:
#line 2272 "pddl2.y"
{
  if (metric_type != metric_none) {
    if (warning_level > 0) {
      std::cerr << "warning: multiple :metric expressions - overwriting previous definition" << std::endl;
    }
  }
  metric_type = metric_maximize;
;
    break;}
case 382:
#line 2284 "pddl2.y"
{
  if (yyvsp[0].exp->exp_class == exp_time) {
    metric = 0;
    metric_type = metric_makespan;
    yyval.exp = 0;
  }
  else {
    metric = yyvsp[0].exp;
    yyval.exp = yyvsp[0].exp;
  }
;
    break;}
case 383:
#line 2299 "pddl2.y"
{
  serial_length = I_TO_N(yyvsp[-1].ival);
;
    break;}
case 384:
#line 2303 "pddl2.y"
{
  parallel_length = I_TO_N(yyvsp[-1].ival);
;
    break;}
case 385:
#line 2316 "pddl2.y"
{ yyval.rval = N_TO_NN(yyvsp[0].ival); ;
    break;}
case 386:
#line 2317 "pddl2.y"
{ yyval.rval = N_TO_NN(R_TO_N(yyvsp[-2].ival,yyvsp[0].ival)); ;
    break;}
case 387:
#line 2318 "pddl2.y"
{ yyval.rval = yyvsp[0].rval; ;
    break;}
case 388:
#line 2319 "pddl2.y"
{ yyval.rval = POS_INF; ;
    break;}
case 389:
#line 2326 "pddl2.y"
{
  current_item = new DKEL_Item(":invariant");
  current_context = current_item;
;
    break;}
case 391:
#line 2335 "pddl2.y"
{
  dom_sc_invariants.append(new SetConstraint(current_item));
;
    break;}
case 392:
#line 2339 "pddl2.y"
{
  dom_sc_invariants[dom_sc_invariants.length() - 1]->sc_type = yyvsp[-4].sckw;
  dom_sc_invariants[dom_sc_invariants.length() - 1]->sc_count = yyvsp[-3].ival;
  clear_context(current_param);
  current_context = 0;
;
    break;}
case 393:
#line 2346 "pddl2.y"
{
  dom_f_invariants.append(new InvariantFormula(current_item, yyvsp[-1].ff));
  clear_context(current_param);
  current_context = 0;
;
    break;}
case 394:
#line 2355 "pddl2.y"
{
  yyval.ff = new Formula(fc_false);
;
    break;}
case 395:
#line 2359 "pddl2.y"
{
  yyval.ff = new Formula(fc_true);
;
    break;}
case 396:
#line 2363 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 397:
#line 2367 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  yyval.ff = new AFormula((Atom*)current_atom);
;
    break;}
case 398:
#line 2372 "pddl2.y"
{
  yyval.ff = new EqFormula((Symbol*)yyvsp[-2].sym->val, (Symbol*)yyvsp[-1].sym->val);
;
    break;}
case 399:
#line 2376 "pddl2.y"
{
  yyval.ff = new NFormula(yyvsp[-1].ff);
;
    break;}
case 400:
#line 2380 "pddl2.y"
{
  yyval.ff = yyvsp[-1].ff;
  yyval.ff->fc = fc_conjunction;
;
    break;}
case 401:
#line 2385 "pddl2.y"
{
  yyval.ff = yyvsp[-1].ff;
  yyval.ff->fc = fc_disjunction;
;
    break;}
case 402:
#line 2390 "pddl2.y"
{
  yyval.ff = new BFormula(fc_implication, yyvsp[-2].ff, yyvsp[-1].ff);
;
    break;}
case 403:
#line 2394 "pddl2.y"
{
  yyval.ff = new BFormula(fc_equivalence, yyvsp[-2].ff, yyvsp[-1].ff);
;
    break;}
case 404:
#line 2398 "pddl2.y"
{
  stored_n_param.append(current_param.length());
;
    break;}
case 405:
#line 2402 "pddl2.y"
{
  QFormula* qf = new QFormula(fc_universal, yyvsp[-1].ff);
  assert(stored_n_param.length() > 0);
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    qf->vars.append(current_param[k]);
  }
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  yyval.ff = qf;
;
    break;}
case 406:
#line 2417 "pddl2.y"
{
  stored_n_param.append(current_param.length());
;
    break;}
case 407:
#line 2421 "pddl2.y"
{
  QFormula* qf = new QFormula(fc_existential, yyvsp[-1].ff);
  assert(stored_n_param.length() > 0);
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    qf->vars.append(current_param[k]);
  }
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  yyval.ff = qf;
;
    break;}
case 408:
#line 2439 "pddl2.y"
{
  ((CFormula*)yyval.ff)->add(yyvsp[0].ff);
;
    break;}
case 409:
#line 2443 "pddl2.y"
{
  yyval.ff = new CFormula(fc_list);
  ((CFormula*)yyval.ff)->add(yyvsp[0].ff);
;
    break;}
case 410:
#line 2451 "pddl2.y"
{
  current_item = new IrrelevantItem();
  if (problem_name) current_item->defined_in_problem = true;
  current_context = current_item;
;
    break;}
case 414:
#line 2466 "pddl2.y"
{
  current_atom = new Reference((ActionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 415:
#line 2470 "pddl2.y"
{
  IrrelevantItem* item = (IrrelevantItem*)current_item;
  item->entity = (Reference*)current_atom;
  dom_irrelevant.append(item);
  ActionSymbol* act = (ActionSymbol*)yyvsp[-4].sym->val;
  act->irr_ins.append(item);
  clear_context(current_param);
  current_context = 0;
;
    break;}
case 416:
#line 2483 "pddl2.y"
{
  current_atom = new Reference((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 417:
#line 2487 "pddl2.y"
{
  IrrelevantItem* item = (IrrelevantItem*)current_item;
  item->entity = (Reference*)current_atom;
  dom_irrelevant.append(item);
  PredicateSymbol* pred = (PredicateSymbol*)yyvsp[-4].sym->val;
  pred->irr_ins.append(item);
  clear_context(current_param);
  current_context = 0;
;
    break;}
case 418:
#line 2500 "pddl2.y"
{
  current_item = new DKEL_Item(":invariant");
  current_item->defined_in_problem = true;
  current_context = current_item;
;
    break;}
case 420:
#line 2510 "pddl2.y"
{
  dom_sc_invariants.append(new SetConstraint(current_item));
;
    break;}
case 421:
#line 2514 "pddl2.y"
{
  dom_sc_invariants[dom_sc_invariants.length() - 1]->sc_type = yyvsp[-4].sckw;
  dom_sc_invariants[dom_sc_invariants.length() - 1]->sc_count = yyvsp[-3].ival;
  clear_context(current_param);
  current_context = 0;
;
    break;}
case 422:
#line 2521 "pddl2.y"
{
  dom_f_invariants.append(new InvariantFormula(current_item, yyvsp[-1].ff));
  clear_context(current_param);
  current_context = 0;
;
    break;}
case 428:
#line 2538 "pddl2.y"
{
  current_item->item_tags.insert(yyvsp[0].sym->text);
;
    break;}
case 429:
#line 2545 "pddl2.y"
{
  yyvsp[0].sym->val = new Symbol(sym_misc, yyvsp[0].sym->text);
  current_item->name = (Symbol*)yyvsp[0].sym->val;
  current_item->item_tags.insert(yyvsp[0].sym->text);
;
    break;}
case 430:
#line 2551 "pddl2.y"
{
  current_item->name = (Symbol*)yyvsp[0].sym->val;
  current_item->item_tags.insert(yyvsp[0].sym->text);
;
    break;}
case 431:
#line 2559 "pddl2.y"
{
  current_param.clear();
;
    break;}
case 432:
#line 2563 "pddl2.y"
{
  current_context->param = current_param;
;
    break;}
case 435:
#line 2575 "pddl2.y"
{
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    current_context->param.append(current_param[k]);
  }
;
    break;}
case 460:
#line 2627 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 461:
#line 2631 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->pos_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->pos_pre = true;
;
    break;}
case 462:
#line 2638 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 463:
#line 2642 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->pos_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->pos_pre = true;
;
    break;}
case 464:
#line 2648 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_init);
;
    break;}
case 465:
#line 2652 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->pos_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->pos_pre = true;
;
    break;}
case 466:
#line 2658 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-2].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-1].sym->val;
  current_context->pos_con.append(eq_atom);
;
    break;}
case 467:
#line 2669 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 468:
#line 2673 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->neg_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 469:
#line 2680 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 470:
#line 2684 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->neg_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 471:
#line 2690 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_init);
;
    break;}
case 472:
#line 2694 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->neg_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 473:
#line 2700 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_init);
;
    break;}
case 474:
#line 2704 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->neg_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 475:
#line 2710 "pddl2.y"
{
  Atom* eq_atom = new Atom(dom_eq_pred);
  eq_atom->param.set_length(2);
  eq_atom->param[0] = (Symbol*)yyvsp[-3].sym->val;
  eq_atom->param[1] = (Symbol*)yyvsp[-2].sym->val;
  current_context->neg_con.append(eq_atom);
;
    break;}
case 476:
#line 2721 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_pos_goal);
;
    break;}
case 477:
#line 2725 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->pos_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->pos_pre = true;
;
    break;}
case 478:
#line 2731 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_neg_goal);
;
    break;}
case 479:
#line 2735 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->pos_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 480:
#line 2744 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_pos_goal);
;
    break;}
case 481:
#line 2748 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->neg_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->pos_pre = true;
;
    break;}
case 482:
#line 2754 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val, md_neg_goal);
;
    break;}
case 483:
#line 2758 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_context->neg_con.append((Atom*)current_atom);
  ((Atom*)current_atom)->pred->neg_pre = true;
;
    break;}
case 484:
#line 2767 "pddl2.y"
{
  TypeConstraint* c =
    new TypeConstraint((VariableSymbol*)yyvsp[-1].sym->val, (TypeSymbol*)yyvsp[-2].sym->val);
  current_context->type_con.append(c);
;
    break;}
case 485:
#line 2773 "pddl2.y"
{
  TypeConstraint* c =
    new TypeConstraint((VariableSymbol*)yyvsp[-1].sym->val, (TypeSymbol*)yyvsp[-2].sym->val);
  current_context->type_con.append(c);
;
    break;}
case 486:
#line 2781 "pddl2.y"
{
  yyval.sckw = sc_at_least;
;
    break;}
case 487:
#line 2785 "pddl2.y"
{
  yyval.sckw = sc_at_most;
;
    break;}
case 488:
#line 2789 "pddl2.y"
{
  yyval.sckw = sc_exactly;
;
    break;}
case 492:
#line 2801 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 493:
#line 2805 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  dom_sc_invariants[dom_sc_invariants.length()-1]->pos_atoms.append((Atom*)current_atom);
;
    break;}
case 494:
#line 2810 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 495:
#line 2814 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  dom_sc_invariants[dom_sc_invariants.length()-1]->neg_atoms.append((Atom*)current_atom);
;
    break;}
case 496:
#line 2822 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  current_context = new SetOf();
;
    break;}
case 497:
#line 2827 "pddl2.y"
{
  SetOf* s = (SetOf*)current_context;
  for (HSPS::index_type k = stored_n_param[stored_n_param.length() - 1];
       k < current_param.length(); k++) {
    s->param.append(current_param[k]);
  }
;
    break;}
case 498:
#line 2835 "pddl2.y"
{
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 499:
#line 2839 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  SetOf* s = (SetOf*)current_context;
  s->pos_atoms.append((Atom*)current_atom);
  dom_sc_invariants[dom_sc_invariants.length()-1]->atom_sets.append(s);
  assert(stored_n_param.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = 0;
;
    break;}
case 500:
#line 2856 "pddl2.y"
{
  input_plans.append(new InputPlan());
  input_plans[input_plans.length() - 1]->src = current_file();
  if (current_plan_file != current_file()) {
    n_plans_in_current_file = 0;
    current_plan_file = current_file();
  }
;
    break;}
case 501:
#line 2865 "pddl2.y"
{
  if (input_plans[input_plans.length() - 1]->name == 0)
    if (current_plan_file) {
      std::ostringstream pn;
      pn << current_plan_file << ":" << n_plans_in_current_file;
      Symbol* plan_file_name = new Symbol(sym_misc, strdup(pn.str().c_str()));
      input_plans[input_plans.length() - 1]->name = plan_file_name;
    }
  n_plans_in_current_file += 1;
;
    break;}
case 503:
#line 2880 "pddl2.y"
{
  assert(input_plans.length() > 0);
  input_plans[input_plans.length() - 1]->is_opt = true;
;
    break;}
case 506:
#line 2890 "pddl2.y"
{
  assert(input_plans.length() > 0);
  input_plans[input_plans.length() - 1]->name = new Symbol(yyvsp[0].sym->text);
;
    break;}
case 507:
#line 2898 "pddl2.y"
{
  current_atom = new Reference((ActionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 508:
#line 2902 "pddl2.y"
{
  Reference* ref = (Reference*)current_atom;
  ActionSymbol* act = (ActionSymbol*)yyvsp[-3].sym->val;
  act->refs.append(ref);
  assert(input_plans.length() > 0);
  input_plans[input_plans.length() - 1]->
    steps.append(new InputPlanStep(ref, yyvsp[-6].rval));
  clear_context(current_param);
;
    break;}
case 509:
#line 2915 "pddl2.y"
{
  // input_plans.append(0);
  current_plan_file = 0;
;
    break;}
case 510:
#line 2920 "pddl2.y"
{
  // input_plans.append(0);
  current_plan_file = 0;
;
    break;}
case 513:
#line 2933 "pddl2.y"
{
  current_atom = new Reference((ActionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 514:
#line 2937 "pddl2.y"
{
  Reference* ref = (Reference*)current_atom;
  ActionSymbol* act = (ActionSymbol*)yyvsp[-4].sym->val;
  act->refs.append(ref);
  if (input_plans.length() == 0) {
    at_position(std::cerr) << "beginning of new plan" << std::endl;
    input_plans.append(new InputPlan());
    input_plans[input_plans.length() - 1]->src = current_file();
    if (current_file()) {
      Symbol* plan_file_name = new Symbol(sym_misc, current_file());
      input_plans[input_plans.length() - 1]->name = plan_file_name;
    }
    current_plan_file = current_file();
  }
  else if (current_file() != current_plan_file) {
    at_position(std::cerr) << "beginning of new plan (new file)" << std::endl;
    input_plans.append(new InputPlan());
    input_plans[input_plans.length() - 1]->src = current_file();
    if (current_file()) {
      Symbol* plan_file_name = new Symbol(sym_misc, current_file());
      input_plans[input_plans.length() - 1]->name = plan_file_name;
    }
    current_plan_file = current_file();
  }
  input_plans[input_plans.length() - 1]->
    steps.append(new InputPlanStep(ref, yyvsp[-7].rval));
  clear_context(current_param);
;
    break;}
case 519:
#line 2979 "pddl2.y"
{
  current_atom = new Reference((ActionSymbol*)yyvsp[0].sym->val);
;
    break;}
case 520:
#line 2983 "pddl2.y"
{
  Reference* ref = (Reference*)current_atom;
  ActionSymbol* act = (ActionSymbol*)yyvsp[-3].sym->val;
  act->refs.append(ref);
  if (input_plans.length() == 0) {
    at_position(std::cerr) << "beginning of new plan" << std::endl;
    input_plans.append(new InputPlan());
    input_plans[input_plans.length() - 1]->src = current_file();
    if (current_file()) {
      Symbol* plan_file_name = new Symbol(sym_misc, current_file());
      input_plans[input_plans.length() - 1]->name = plan_file_name;
    }
    current_plan_file = current_file();
  }
  else if (current_file() != current_plan_file) {
    at_position(std::cerr) << "beginning of new plan (new file)" << std::endl;
    input_plans.append(new InputPlan());
    input_plans[input_plans.length() - 1]->src = current_file();
    if (current_file()) {
      Symbol* plan_file_name = new Symbol(sym_misc, current_file());
      input_plans[input_plans.length() - 1]->name = plan_file_name;
    }
    current_plan_file = current_file();
  }
  HSPS::index_type n = input_plans[input_plans.length() - 1]->steps.length();
  input_plans[input_plans.length() - 1]->
    steps.append(new InputPlanStep(ref, n));
  clear_context(current_param);
;
    break;}
case 524:
#line 3025 "pddl2.y"
{
  current_entry = new HTableEntry();
;
    break;}
case 525:
#line 3029 "pddl2.y"
{
  h_table.append(current_entry);
  current_entry = 0;
;
    break;}
case 526:
#line 3037 "pddl2.y"
{
  current_entry->cost = yyvsp[0].rval;
  current_entry->opt = true;
;
    break;}
case 527:
#line 3042 "pddl2.y"
{
  current_entry->cost = yyvsp[0].rval;
;
    break;}
case 532:
#line 3059 "pddl2.y"
{
  assert(current_entry);
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 533:
#line 3064 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_entry->atoms.append((Atom*)current_atom);
  current_entry->neg.append(false);
  assert(current_entry->atoms.length() == current_entry->neg.length());
;
    break;}
case 534:
#line 3074 "pddl2.y"
{
  assert(current_entry);
  current_atom = new Atom((PredicateSymbol*)yyvsp[0].sym->val);
;
    break;}
case 535:
#line 3079 "pddl2.y"
{
  ((Atom*)current_atom)->check();
  current_entry->atoms.append((Atom*)current_atom);
  current_entry->neg.append(true);
  assert(current_entry->atoms.length() == current_entry->neg.length());
;
    break;}
case 536:
#line 3089 "pddl2.y"
{
  ReferenceSet* set = new ReferenceSet();
  current_context = set;
  stored_n_param.assign_value(0, 1);
  current_param.clear();
  input_sets.append(set);
;
    break;}
case 537:
#line 3097 "pddl2.y"
{
  clear_context(current_param);
  current_context = 0;
;
    break;}
case 538:
#line 3105 "pddl2.y"
{
  yyvsp[0].sym->val = new Symbol(sym_misc, yyvsp[0].sym->text);
  ((ReferenceSet*)current_context)->name = (Symbol*)yyvsp[0].sym->val;
;
    break;}
case 539:
#line 3110 "pddl2.y"
{
  ((ReferenceSet*)current_context)->name = (Symbol*)yyvsp[0].sym->val;
;
    break;}
case 544:
#line 3124 "pddl2.y"
{
  if (yyvsp[0].sym->val) {
    current_atom = new Reference((Symbol*)yyvsp[0].sym->val, false, true);
  }
  else {
    current_atom = new Reference(new Symbol(sym_misc, yyvsp[0].sym->text), false, true);
  }
;
    break;}
case 545:
#line 3133 "pddl2.y"
{
  assert(input_sets.length() > 0);
  assert(input_sets[input_sets.length() - 1] != 0);
  input_sets[input_sets.length() - 1]->
    add(new SimpleReferenceSet((Reference*)current_atom));
;
    break;}
case 546:
#line 3140 "pddl2.y"
{
  assert(input_sets.length() > 0);
  assert(input_sets[input_sets.length() - 1] != 0);
  if (yyvsp[0].sym->val) {
    input_sets[input_sets.length() - 1]->add
      (new SimpleReferenceSet(new Reference((Symbol*)yyvsp[0].sym->val, false, false)));
  }
  else {
    input_sets[input_sets.length() - 1]->add
      (new SimpleReferenceSet(new Reference(new Symbol(sym_misc, yyvsp[0].sym->text), false, false)));
  }
;
    break;}
case 547:
#line 3156 "pddl2.y"
{
  stored_n_param.append(current_param.length());
  stored_context.append(current_context);
  current_context = new SimpleReferenceSet(0);
;
    break;}
case 548:
#line 3162 "pddl2.y"
{
  assert(stored_n_param.length() > 0);
  assert(stored_context.length() > 0);
  clear_context(current_param,
		stored_n_param[stored_n_param.length() - 1],
		current_param.length());
  current_param.set_length(stored_n_param[stored_n_param.length() - 1]);
  stored_n_param.dec_length();
  current_context = stored_context[stored_context.length() - 1];
  stored_context.dec_length();
;
    break;}
case 549:
#line 3177 "pddl2.y"
{
  if (yyvsp[0].sym->val) {
    current_atom = new Reference((Symbol*)yyvsp[0].sym->val, false, true);
  }
  else {
    current_atom = new Reference(new Symbol(sym_misc, yyvsp[0].sym->text), false, true);
  }
;
    break;}
case 550:
#line 3186 "pddl2.y"
{
  SimpleReferenceSet* s = (SimpleReferenceSet*)current_context;
  s->ref = (Reference*)current_atom;
  assert(input_sets.length() > 0);
  assert(input_sets[input_sets.length() - 1] != 0);
  input_sets[input_sets.length() - 1]->add(s);
;
    break;}
case 551:
#line 3194 "pddl2.y"
{
  if (yyvsp[0].sym->val) {
    current_atom = new Reference((Symbol*)yyvsp[0].sym->val, true, true);
  }
  else {
    current_atom = new Reference(new Symbol(sym_misc, yyvsp[0].sym->text), true, true);
  }
;
    break;}
case 552:
#line 3203 "pddl2.y"
{
  SimpleReferenceSet* s = (SimpleReferenceSet*)current_context;
  s->ref = (Reference*)current_atom;
  assert(input_sets.length() > 0);
  assert(input_sets[input_sets.length() - 1] != 0);
  input_sets[input_sets.length() - 1]->add(s);
;
    break;}
}

#line 811 "/usr/local/lib/bison.cc"
   /* the action file gets copied in in place of this dollarsign  */
  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YY_PDDL_Parser_LSP_NEEDED
  yylsp -= yylen;
#endif

#if YY_PDDL_Parser_DEBUG != 0
  if (YY_PDDL_Parser_DEBUG_FLAG)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;

#ifdef YY_PDDL_Parser_LSP_NEEDED
  yylsp++;
  if (yylen == 0)
    {
      yylsp->first_line = YY_PDDL_Parser_LLOC.first_line;
      yylsp->first_column = YY_PDDL_Parser_LLOC.first_column;
      yylsp->last_line = (yylsp-1)->last_line;
      yylsp->last_column = (yylsp-1)->last_column;
      yylsp->text = 0;
    }
  else
    {
      yylsp->last_line = (yylsp+yylen-1)->last_line;
      yylsp->last_column = (yylsp+yylen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  YYGOTO(yynewstate);

YYLABEL(yyerrlab)   /* here on detecting error */

  if (! yyerrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++YY_PDDL_Parser_NERRS;

#ifdef YY_PDDL_Parser_ERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  /* Start X at -yyn if nec to avoid negative indexes in yycheck.  */
	  for (x = (yyn < 0 ? -yyn : 0);
	       x < (sizeof(yytname) / sizeof(char *)); x++)
	    if (yycheck[x + yyn] == x)
	      size += strlen(yytname[x]) + 15, count++;
	  msg = (char *) malloc(size + 15);
	  if (msg != 0)
	    {
	      strcpy(msg, "parse error");

	      if (count < 5)
		{
		  count = 0;
		  for (x = (yyn < 0 ? -yyn : 0);
		       x < (sizeof(yytname) / sizeof(char *)); x++)
		    if (yycheck[x + yyn] == x)
		      {
			strcat(msg, count == 0 ? ", expecting `" : " or `");
			strcat(msg, yytname[x]);
			strcat(msg, "'");
			count++;
		      }
		}
	      YY_PDDL_Parser_ERROR(msg);
	      free(msg);
	    }
	  else
	    YY_PDDL_Parser_ERROR ("parse error; also virtual memory exceeded");
	}
      else
#endif /* YY_PDDL_Parser_ERROR_VERBOSE */
	YY_PDDL_Parser_ERROR("parse error");
    }

  YYGOTO(yyerrlab1);
YYLABEL(yyerrlab1)   /* here on error raised explicitly by an action */

  if (yyerrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (YY_PDDL_Parser_CHAR == YYEOF)
	YYABORT;

#if YY_PDDL_Parser_DEBUG != 0
      if (YY_PDDL_Parser_DEBUG_FLAG)
	fprintf(stderr, "Discarding token %d (%s).\n", YY_PDDL_Parser_CHAR, yytname[yychar1]);
#endif

      YY_PDDL_Parser_CHAR = YYEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3;              /* Each real token shifted decrements this */

  YYGOTO(yyerrhandle);

YYLABEL(yyerrdefault)  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (yyn) YYGOTO(yydefault);
#endif

YYLABEL(yyerrpop)   /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss) YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YY_PDDL_Parser_LSP_NEEDED
  yylsp--;
#endif

#if YY_PDDL_Parser_DEBUG != 0
  if (YY_PDDL_Parser_DEBUG_FLAG)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

YYLABEL(yyerrhandle)

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    YYGOTO(yyerrdefault);

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    YYGOTO(yyerrdefault);

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	YYGOTO(yyerrpop);
      yyn = -yyn;
      YYGOTO(yyreduce);
    }
  else if (yyn == 0)
    YYGOTO(yyerrpop);

  if (yyn == YYFINAL)
    YYACCEPT;

#if YY_PDDL_Parser_DEBUG != 0
  if (YY_PDDL_Parser_DEBUG_FLAG)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++yyvsp = YY_PDDL_Parser_LVAL;
#ifdef YY_PDDL_Parser_LSP_NEEDED
  *++yylsp = YY_PDDL_Parser_LLOC;
#endif

  yystate = yyn;
  YYGOTO(yynewstate);
/* end loop, in which YYGOTO may be used. */
  YYENDGOTO
}

/* END */

/* #line 1010 "/usr/local/lib/bison.cc" */
#line 5875 "grammar.cc"
#line 3212 "pddl2.y"

