#!/bin/bash

# execute as ./run_brazilian.sh

dir="/uscms_data/d3/cihar29/Analysis/CMSSW_8_1_0/src/theta/utils2/2017/"
plotObs="false"

signals=( "gkk" "zp1" "zp10" "zp30" )
variables=( "mass" "st" "sumrmin" )

for sig in "${signals[@]}" ; do
  for var in "${variables[@]}" ; do

    folder="${sig}_${var}"
    cmd="root -l -b -q 'brazilian.c (\"$dir\", \"$folder\", $plotObs)'"
    eval $cmd

  done
done
