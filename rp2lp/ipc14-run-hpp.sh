#!/bin/bash

whereami=`dirname $0`
${whereami}/test_ilb -c 3 -x -a -s "$1" "$2" > "$3"
