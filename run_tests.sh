#!/bin/sh

OWD=$PWD
TESTDIR=$PWD/test
TEST="casim_test"

cd $TESTDIR

timestamp=$(date +%FT%T)
LOG="${OWD}/${TEST}@${timestamp}.log"

python $TEST.py -v
