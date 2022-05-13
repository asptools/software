#!/bin/bash

BIN=../bin
DIR=.
CLASP=../bin #clasp executable address	
GRINGO=../bin #gringo executable address


$GRINGO/gringo --output smodels $* \
| $CLASP/clasp










