#!/bin/bash

# execute as nohup ./run_all.sh OFF > & run.txt &

args=("$@")

if [ $# -ne 1 ] ; then
  echo "Please provide a DR cut (OFF, ON, ONbt, ONnb, REVERSE)."
  exit
fi

dr=${args[0]}

channels=( "mm" "em" "ee" )
uncerts=( "jec" "jer" "btagSF" "mistagSF" "pileup" "topPtWeight" "pdf" "q2ttbar" "muTrigSys" "muIdSys" "eleTrigSys" "eleIdSys" )
types=( "UP" "DOWN" )

dir="/uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/root_trees"

mkdir logs
start=$(date +%s.%N)

for chan in "${channels[@]}" ; do
  echo "Processing" $chan
  mkdir $chan
  sh ./analyze_all.sh $chan $dr > "logs/logMC_${chan}.txt"

  for uncert in "${uncerts[@]}" ; do
    for type in "${types[@]}" ; do
      echo "Processing" $chan $uncert $type
      sh ./analyze_all.sh $chan $dr $uncert $type > "logs/logMC_${chan}_${uncert}${type}.txt"
    done
  done

  dname="Muon"
  if [[ $chan = "ee" ]] ; then
    dname="Ele"
  fi
  lines=( "isMC         false"
          "setSUMDRCut  $dr"
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
