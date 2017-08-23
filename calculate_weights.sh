#!/bin/bash

# execute as ./calculate_weights.sh lumi_value

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a lumi value."
fi

if [ $# -eq 1 ] ; then

  lumi=${args[0]}

  dataset=(
    "DYhigh"
    "DYlow"
    "STschannel"
    "STtWchannel"
    "STtchannel"
    "SaTtWchannel"
    "SaTtchannel"
    "TTbar"
    "WJets"
    "WW"
    "WZ"
    "ZZ"
    "gluon_M-1000"
    "gluon_M-1250"
    "gluon_M-1500"
    "gluon_M-2000"
    "gluon_M-2500"
    "gluon_M-3000"
    "gluon_M-3500"
    "gluon_M-4000"
    "gluon_M-500"
    "gluon_M-750"
    "zprime_M-1000_W-10"
    "zprime_M-1000_W-100"
    "zprime_M-1000_W-300"
    "zprime_M-1250_W-125"
    "zprime_M-1250_W-12p5"
    "zprime_M-1500_W-15"
    "zprime_M-1500_W-150"
    "zprime_M-3000_W-300"
  )

  xs=( 5765400 18610000 26380 35850 136020 35850 80950 831760 61526700 118700 47100 16500
       10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 )

  events=( 145407553 35291457 2989179 6952779 67240483 6933037 38810801 154857762 86731469 7981093 3995807 1988085
           98558 99995 99993 99995 99990 99740 99488 99113 99997 100000 103781 101055 79474 96841 102827 99684 111101 189120 )

  file="dataset lumi xs lumi*xs events weight"

  for ((i=0;i<${#dataset[@]};++i)); do

    lumi_xs=$(echo "$lumi*${xs[i]}" | bc)

    line="${dataset[i]} $lumi ${xs[i]} $lumi_xs ${events[i]} $(echo "scale=5; $lumi_xs/${events[i]}" | bc)"

    file="$file\n$line"

  done

  echo -e "$file" | column -t > mc_weights.txt

fi
