#!/bin/bash

# execute as ./analyze_all.sh

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a channel (mm, ee, em)."
fi

if [ $# -eq 1 ] ; then

  channel=${args[0]}

  pileupName="/uscms/home/cihar29/nobackup/Analysis/CMSSW_8_0_20/src/analysis/ZDilepton/mu_weights_SingleMuon.root"
  btagName="/uscms/home/cihar29/nobackup/Analysis/CMSSW_8_0_20/src/analysis/ZDilepton/btag_eff_SingleMuon.root"

  if [ "$channel" = "ee" ] ; then
    pileupName="/uscms/home/broozbah/nobackup/AnalysisZP/CMSSW_8_0_19/src/Analysis/ZDilepton/mu_weights_DoubleEle.root"
    btagName="/uscms/home/broozbah/nobackup/AnalysisZP/CMSSW_8_0_19/src/Analysis/ZDilepton/btag_eff_DE_ee.root"
  fi

  mcfiles=( "TTbar" "lowDY" "highDY" "STtchannel" "SaTtchannel" "STschannel" "STtWchannel" "SaTtWchannel" "Wjet"
            "zprime-M3000-W300" "gluonkk-M3000" )

  for i in "${mcfiles[@]}"
  do

    lines=( "isMC           true"
            "topPt_weight   NOMINAL" 
            "setDRCut       OFF"   
            "inName         /uscms_data/d3/broozbah/AnalysisZP/CMSSW_8_0_19/src/Analysis/ZDilepton/Chads_root/${i}.root"
            "outName        ./rootFiles_SM/Jun15/${i}_${channel}.root"
            "channel        ${channel}"
            "eras           Summer16_23Sep2016V4_MC"
            "res_era        Spring16_25nsV10_MC"
            "jet_type       AK4PFchs"
            "muTrigSfName   Trigger_EfficienciesAndSF_Period4.root"
            "muIdSfName     ID_EfficienciesAndSF_GH.root"
            "muTrackSfName  Tracking_EfficienciesAndSF_GH.root"
            "eRecoSfName    e_Reco_efficiency.root"
            "eIdSfName      e_MediumID_efficiency.root"
            "btagName       ${btagName}"
            "pileupName     ${pileupName}"
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

