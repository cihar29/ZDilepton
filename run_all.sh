#!/bin/bash

# execute as nohup ./run_all.sh > & run.txt &

start=$(date +%s.%N)

channels=( "mm" "em" "ee" )
uncerts=( "topPtWeight" "jec" "jer" "pdf" "q2ttbar" "q2dy" "q2st" "q2signal" "btagSF" "mistagSF" "pileup" )
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

end=$(date +%s.%N)
echo "$(echo "$end - $start" | bc) seconds"
