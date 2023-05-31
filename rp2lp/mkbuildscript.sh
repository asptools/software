#!/bin/bash

#CC=g++
#CFLAGS="-O6 -m32 -x c++"

echo "set -x"

echo "CFLAGS=\"$CFLAGS\""
echo "LFLAGS=\"$LFLAGS\""
echo "LLIBS=\"$LLIBS\""

for file in $*; do
    echo "$CCC -c \$CFLAGS $file"
done

echo -n "$CCC \$LFLAGS -o `echo $1 | sed -e 's/\.ccx\?//'` ";

for file in $*; do
    echo -n " `echo $file | sed -e 's/\.ccx\?/.o/'` ";
done

echo " \$LLIBS"
