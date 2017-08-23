#!/bin/bash

# execute as nohup ./run_all.sh &

#START=$(date +%s.%N)

channels=( "mm" "em" "ee" )
uncerts=( "topPt_weight" "jec" "jer" "pdf" "q2" "btagSF" "mistagSF" "pileup" )
types=( "UP" "DOWN" )

for chan in "${channels[@]}" ; do
  echo "Processing" $chan
  sh ./analyze_all.sh $chan > "logs/logMC_${chan}.txt"

  for uncert in "${uncerts[@]}" ; do
    for type in "${types[@]}" ; do
      echo "Processing" $chan $uncert $type
      sh ./analyze_all.sh $chan $uncert $type > "logs/logMC_${chan}_${uncert}${type}.txt"
    done
  done
done

#END=$(date +%s.%N)
#DIFF=$(echo "$END - $START" | bc)
#echo $DIFF
