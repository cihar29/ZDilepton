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

    uncerts=( "topPt_weight" "jec" "jer" "pdf" "q2ttbar" "q2dy" "q2st" "q2signal" "btagSF" "mistagSF" "pileup" )
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

  for file in "${mcfiles[@]}" ; do

    lines=( "isMC           true"
            "${uncert}      ${type}"
            "setDRCut       OFF"   
            "inName         ${dir}${file}.root"
            "outName        ./root_Sep28/${channel}/${file}_${channel}${underscore}${uncert}${type}.root"
            "channel        ${channel}"
            "eras           Summer16_23Sep2016V4_MC"
            "res_era        Spring16_25nsV10_MC"
            "jet_type       AK4PFchs"
            "muTrigSfName   /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/Trigger_EfficienciesAndSF_Period4.root"
            "muIdSfName     /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/ID_EfficienciesAndSF_GH.root"
            "muTrackSfName  /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_26_patch2/src/analysis/ZDilepton/Tracking_EfficienciesAndSF_GH.root"
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
