#!/bin/bash

BIN=../bin
DIR=.
WASP=../bin #wasp executable address	
GRINGO=../bin #gringo executable address


$GRINGO/gringo --output smodels $* \
| $WASP/wasp --stats=0






