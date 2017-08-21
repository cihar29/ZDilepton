#!/bin/bash

# execute as ./analyze_all.sh

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a channel (mm, ee, em)."
fi

if [ $# -eq 1 ] ; then

  channel=${args[0]}

  dir="/uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/root_trees/"
  mcfiles=(
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
    "gluon_M-3000"
    "zprime_M-1000_W-10"
    "zprime_M-1000_W-100"
    "zprime_M-1000_W-300"
    "zprime_M-1250_W-125"
    "zprime_M-1250_W-12p5"
    "zprime_M-1500_W-15"
    "zprime_M-1500_W-150"
    "zprime_M-3000_W-300"
  )

  for i in "${mcfiles[@]}"
  do

    lines=( "isMC           true"
            "topPt_weight   NOMINAL" 
            "jec            NOMINAL"
            "jer            NOMINAL"
            "pdf            NOMINAL"
            "q2             NOMINAL"
            "btagSF         NOMINAL"
            "mistagSF       NOMINAL"
            "pileup         NOMINAL"
            "setDRCut       OFF"   
            "inName         ${dir}${i}.root"
            "outName        ./${channel}/${i}_${channel}.root"
            "channel        ${channel}"
            "eras           Summer16_23Sep2016V4_MC"
            "res_era        Spring16_25nsV10_MC"
            "jet_type       AK4PFchs"
            "muTrigSfName   Trigger_EfficienciesAndSF_Period4.root"
            "muIdSfName     ID_EfficienciesAndSF_GH.root"
            "muTrackSfName  Tracking_EfficienciesAndSF_GH.root"
            "eRecoSfName    e_Reco_efficiency.root"
            "eIdSfName      e_MediumID_efficiency.root"
            "btagName       /uscms/home/broozbah/nobackup/AnalysisZP/CMSSW_8_0_19/src/Analysis/ZDilepton/btag_eff_default.root"
            "pileupName     /uscms/home/cihar29/nobackup/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/mu_weights.root"
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
