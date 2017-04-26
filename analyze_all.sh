#!/bin/bash

# execute as ./analyze_all.sh

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a channel (mm, ee, em)."
fi

if [ $# -eq 1 ] ; then

  channel=${args[0]}

  mcfiles=( "TTbar" "lowDY" "highDY" "STtchannel" "SaTtchannel" "STschannel" "STtWchannel" "SaTtWchannel" "Wjet"
            "zprime-M3000-W300" "gluonkk-M3000" )

  for i in "${mcfiles[@]}"
  do

    lines=( "ISMC           true"
            "inName         ${i}.root"
            "outName        ${i}_${channel}.root"
            "channel        ${channel}"
            "eras           Summer16_23Sep2016V4_MC"
            "jet_type       AK4PFchs"
            "muTrigSfName   Trigger_EfficienciesAndSF_Period4.root"
            "muIdSfName     ID_EfficienciesAndSF_GH.root"
            "muTrackSfName  Tracking_EfficienciesAndSF_GH.root"
            "eRecoSfName    e_Reco_efficiency.root"
            "eIdSfName      e_MediumID_efficiency.root"
          )
    line=""

    for j in "${lines[@]}"
    do
      line="$line$j\n"
    done

    echo -e "$line" > parsMC.txt

    analyze parsMC.txt mc_weights.txt

  done

fi
