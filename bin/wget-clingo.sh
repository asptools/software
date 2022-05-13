#!/bin/bash

# Download gringo, clasp and clingo

if test -x gringo -a -x clasp -a -x clingo
then
  echo "gringo, clasp, and clingo are already available!"
  exit 0
fi

VER=5.3.0
PACK=clingo-$VER-linux-x86_64.tar.gz

if test ! -f $PACK
then
  wget https://github.com/potassco/clingo/releases/download/v$VER/$PACK
fi

if test -f $PACK
then
  # Extract

  tar xvfz $PACK \
    ${PACK%%.tar.gz}/gringo ${PACK%%.tar.gz}/clasp ${PACK%%.tar.gz}/clingo
  mv ${PACK%%.tar.gz}/gringo .
  mv ${PACK%%.tar.gz}/clasp .
  mv ${PACK%%.tar.gz}/clingo .
  
  # Cleanup

  rm -f $PACK
  rm -f ${PACK%%.tar.gz}/*
  rmdir ${PACK%%.tar.gz}
else
  echo "Downloading gringo/clasp/clingo (v$VER) failed!"
fi
