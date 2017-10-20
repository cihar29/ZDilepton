#!/bin/bash

# execute as ./analyze_all.sh

isValid()
{
  found=0
  check=$1
  array=("${!2}")
  for a in "${array[@]}" ; do
    if [[ $a = "$check" ]] ; then
        found=1
        break
    fi
  done
}

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a channel (mm, ee, em)."

else

  channel=${args[0]}

  uncert=""
  type=""
  underscore=""
  if [ $# -eq 3 ] ; then

    uncert=${args[1]}
    type=${args[2]}

    uncerts=( "topPtWeight" "jec" "jer" "pdf" "q2ttbar" "q2dy" "q2st" "q2signal" "btagSF" "mistagSF" "pileup" )
    types=( "UP" "DOWN" )

    isValid $uncert uncerts[@]
    if [ $found -eq 0 ] ; then
      echo "Invalid uncertainty"
      exit
    fi
    isValid $type types[@]
    if [ $found -eq 0 ] ; then
      echo "Invalid type"
      exit
    fi
    underscore="_"
  fi

  dir="/uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/root_trees/"
  mcfiles=(
    "DYhigh"
    "DYlow"
    "STschannel"
    "STtWchannel"
    "STtchannel"
    "SaTtWchannel"
    "SaTtchannel"
    "TTbar0-700"
    "TTbar700-1000"
    "TTbar1000-inf"
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
    "zprime_M-500_W-5"
    "zprime_M-500_W-50"
    "zprime_M-750_W-7p5"
    "zprime_M-750_W-75"
    "zprime_M-2000_W-20"
    "zprime_M-2000_W-200"
    "zprime_M-2000_W-600"
    "zprime_M-2500_W-25"
    "zprime_M-2500_W-250"
    "zprime_M-3000_W-30"
    "zprime_M-3500_W-35"
    "zprime_M-3500_W-350"
    "zprime_M-4000_W-40"
    "zprime_M-4000_W-400"
    "zprime_M-4000_W-1200"
    "zprime_M-4500_W-45"
    "zprime_M-4500_W-450"
    "zprime_M-5000_W-50"
    "zprime_M-5000_W-500"
    "zprime_M-5000_W-1500"
  )

  for file in "${mcfiles[@]}" ; do

    lines=( "isMC           true"
            "${uncert}      ${type}"
            "setDRCut       OFF"
            "inName         ${dir}${file}.root"
            "outName        ${channel}/${file}_${channel}${underscore}${uncert}${type}.root"
            "channel        ${channel}"
            "eras           Summer16_23Sep2016V4_MC"
            "res_era        Spring16_25nsV10_MC"
            "jet_type       AK4PFchs"
            "muTrigSfName   /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/Trigger_EfficienciesAndSF_Period4.root"
            "muIdSfName     /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/ID_EfficienciesAndSF_GH.root"
            "muTrackSfName  /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/Tracking_EfficienciesAndSF_GH.root"
            "eTrigSfName    /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/electronTrigSF.root"
            "eRecoSfName    /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/e_Reco_efficiency.root"
            "eIdSfName      /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/e_MediumID_efficiency.root"
            "btagName       /uscms/home/cihar29/nobackup/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/btag_eff.root"
            "pileupName     /uscms/home/cihar29/nobackup/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/mu_weights.root"
          )
    out=""

    for line in "${lines[@]}"
    do
      out="$out$line\n"
    done

    echo -e "$out" | column -t > parsMC.txt

    analyze parsMC.txt mc_weights.txt

  done

fi
