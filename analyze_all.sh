#!/bin/bash

# execute as ./analyze_all.sh

args=("$@")

if [ $# -eq 0 ] ; then

  mcfiles=( "TTbar" "lowDY" "highDY" "STtchannel" "SaTtchannel" "STschannel" "STtWchannel" "SaTtWchannel" "Wjet" )

  line1="ISMC       true"
  line4="channel    mm"
  line5="era        Summer16_25nsV5"
  line6="jet_type   AK4PFchs"

  for i in "${mcfiles[@]}"
  do

    line2="inName     "$i".root"
    line3="outName    "$i"_mm.root"
    echo -e "$line1\n$line2\n$line3\n$line4\n$line5\n$line6" > pars.txt

    analyze pars.txt mc_weights.txt

  done

fi
