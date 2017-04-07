#!/bin/bash

# execute as ./calculate_weights.sh lumi_value

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a lumi value."
fi

if [ $# -eq 1 ] ; then

  lumi=${args[0]}

  dataset=( "ttbar" "lowDY" "highDY" "STtchannel" "SaTtchannel" "STschannel" "STtWchannel" "SaTtWchannel" "Wjet"
            "zprime-M3000-W300" "gluon-M3000" )
  xs=( 831760 18610000 5765400 136020 80950 26380 35850 35850 61526700
       415880 415880 )
  events=( 154989126 35291457 145802430 67240483 38810801 2989179 6952779 6933037 86410031
           189120 99740 )

  file="dataset lumi xs lumi*xs events weight"

  for ((i=0;i<${#dataset[@]};++i)); do

    lumi_xs=$(echo "$lumi*${xs[i]}" | bc)

    line="${dataset[i]} $lumi ${xs[i]} $lumi_xs ${events[i]} $(echo "scale=5; $lumi_xs/${events[i]}" | bc)"

    file="$file\n$line"

  done

  echo -e "$file" | column -t > mc_weights.txt

fi
