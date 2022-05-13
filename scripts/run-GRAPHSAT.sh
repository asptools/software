#!/bin/bash

BIN=../bin
DIR=.
GRAPHSAT=../bin #graphsat executable address	
GRINGO=../bin #gringo executable address


$GRINGO/gringo --output smodels $* \
| $BIN/lpstrip \
| $BIN/lpcat \
| $BIN/lpshift \
| $BIN/lp2acyc \
| $BIN/lp2normal2 \
| $BIN/lp2sat -g \
| $GRAPHSAT/graphsat









