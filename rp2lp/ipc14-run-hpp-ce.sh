#!/bin/bash

whereami=`dirname $0`
${whereami}/test_ilb -ce -x -a -s "$1" "$2" > "$3"
