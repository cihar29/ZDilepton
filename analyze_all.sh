#!/bin/bash

# execute as ./analyze_all.sh

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a channel (mm, ee, em)."
fi

if [ $# -eq 1 ] ; then

  channel=${args[0]}

  mcfiles=( "TTbar" "lowDY" "highDY" "STtchannel" "SaTtchannel" "STschannel" "STtWchannel" "SaTtWchannel" "Wjet"
            "zprime-M3000-W300" "gluon-M3000" )

  line1="ISMC       true"
  line4="channel    $channel"
  line5="eras       Summer16_23Sep2016V4_MC"
  line6="jet_type   AK4PFchs"

  for i in "${mcfiles[@]}"
  do

    line2="inName     "$i".root"
    line3="outName    "$i"_"$channel".root"
    echo -e "$line1\n$line2\n$line3\n$line4\n$line5\n$line6" > parsMC.txt

    analyze parsMC.txt mc_weights.txt

  done

fi
