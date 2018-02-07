#!/bin/bash

# execute as nohup ./run_all.sh > & run.txt &

start=$(date +%s.%N)

channels=( "mm" "em" "ee" )
uncerts=( "topPtWeight" "jec" "jer" "pdf" "q2ttbar" "q2dy" "q2st" "q2signal" "btagSF" "mistagSF" "pileup" )
types=( "UP" "DOWN" )

dir="/uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/root_trees"

mkdir logs

for chan in "${channels[@]}" ; do
  echo "Processing" $chan
  mkdir $chan
  sh ./analyze_all.sh $chan > "logs/logMC_${chan}.txt"

  for uncert in "${uncerts[@]}" ; do
    for type in "${types[@]}" ; do
      echo "Processing" $chan $uncert $type
      sh ./analyze_all.sh $chan $uncert $type > "logs/logMC_${chan}_${uncert}${type}.txt"
    done
  done

  dname="Muon"
  if [[ $chan = "ee" ]] ; then
    dname="Ele"
  fi
  lines=( "isMC         false"
          "setSUMDRCut  OFF"
          "inName       ${dir}/${dname}.root"
          "outName      ${chan}/${dname}_${chan}.root"
          "channel      ${chan}"
          "eras         Summer16_23Sep2016BCDV6_DATA Summer16_23Sep2016EFV6_DATA Summer16_23Sep2016GV6_DATA Summer16_23Sep2016HV6_DATA"
          "jet_type     AK4PFchs"
        )
  out=""

  for line in "${lines[@]}" ; do
    out="$out$line\n"
  done
  echo -e "$out" > parsData.txt
  analyze parsData.txt > "logs/logData_${chan}.txt"

done

end=$(date +%s.%N)
echo "$(echo "$end - $start" | bc) seconds"
