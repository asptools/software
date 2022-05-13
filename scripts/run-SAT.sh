#!/bin/bash

BIN=../bin
DIR=.
KISSAT=../bin #kissat executable address	
GRINGO=../bin #gringo executable address


$GRINGO/gringo --output smodels $* \
| $BIN/lpstrip \
| $BIN/lpcat \
| $BIN/lpshift \
| $BIN/lp2acyc \
| $BIN/lp2normal2 \
| $BIN/lp2sat -g \
| $BIN/graph2sat \
| $KISSAT/kissat -n --sat --relaxed 







