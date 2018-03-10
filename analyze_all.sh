#!/bin/bash

# execute as ./analyze_all.sh mm OFF

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

if [ $# -lt 2 ] ; then
  echo "Please provide a channel (mm, ee, em) and a DR cut (OFF, ON, ONbt, ONnb, REVERSE)."
  exit
fi

channel=${args[0]}
dr=${args[1]}

channels=( "mm" "ee" "em" )
drs=( "OFF" "ON" "ONbt" "ONnb" "REVERSE" )

isValid $channel channels[@]
if [ $found -eq 0 ] ; then
  echo "Invalid channel"
  exit
fi
isValid $dr drs[@]
if [ $found -eq 0 ] ; then
  echo "Invalid DR cut"
  exit
fi

uncert=""
type=""
underscore=""
if [ $# -eq 4 ] ; then

  uncert=${args[2]}
  type=${args[3]}

  uncerts=( "jec" "jer" "btagSF" "mistagSF" "pileup" "topPtWeight" "pdf" "q2ttbar" "muTrigSys" "muIdSys" "eleTrigSys" "eleIdSys" )
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
  "gluon_M-4500"
  "gluon_M-5000"
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
  "zprime_M-3000_W-900"
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

: <<'COMMENT'

mcfiles=(
  "qcd30-50_Mu"
  "qcd50-80_Mu"
  "qcd80-120_Mu"
  "qcd120-170_Mu"
  "qcd170-300_Mu"
  "qcd300-470_Mu"
  "qcd470-600_Mu"
  "qcd600-800_Mu"
  "qcd800-1000_Mu"
  "qcd1000-inf_Mu"
  "qcd30-50_EM"
  "qcd50-80_EM"
  "qcd80-120_EM"
  "qcd120-170_EM"
  "qcd170-300_EM"
  "qcd300-inf_EM"
  "qcd30-80_bcToE"
  "qcd80-170_bcToE"
  "qcd170-250_bcToE"
  "qcd250-inf_bcToE"
)

COMMENT

for file in "${mcfiles[@]}" ; do

  if [[ ( "$uncert" == "topPtWeight" || "$uncert" == "pdf" || "$uncert" == "q2ttbar" ) && "$file" != "TTbar"* ]] ; then
    continue
  elif [[ ("$uncert" == "muTrigSys" || "$uncert" == "muIdSys") && "$channel" == "ee" ]] ; then
    continue
  elif [[ "$uncert" == "eleIdSys" && "$channel" == "mm" ]] ; then
    continue
  elif [[ "$uncert" == "eleTrigSys" && "$channel" != "ee" ]] ; then
    continue
  elif [[ "$uncert" == "topPtWeight" && "$type" == "DOWN" ]] ; then
    continue
  fi

  lines=( "isMC           true"
          "$uncert        $type"
          "setSUMDRCut    $dr"
          "inName         ${dir}${file}.root"
          "outName        ${channel}/${file}_${channel}${underscore}${uncert}${type}.root"
          "channel        ${channel}"
          "eras           Summer16_23Sep2016V4_MC"
          "res_era        Spring16_25nsV10_MC"
          "jet_type       AK4PFchs"
          "muTrigSfName   /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/Trigger_EfficienciesAndSF_BCDEF.root:0.549697 /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/Trigger_EfficienciesAndSF_GH.root:0.450303"
          "muIdSfName     /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/ID_EfficienciesAndSF_BCDEF.root:0.549697 /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/ID_EfficienciesAndSF_GH.root:0.450303"
          "muTrackSfName  /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/Tracking_EfficienciesAndSF_BCDEFGH.root"
          "eTrigSfName    /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/electronTrigSF.root"
          "eRecoSfName    /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/e_Reco_efficiency.root"
          "eIdSfName      /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/ele_id_effSF.root"
          "btagName       /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/btag_eff.root"
          "pileupName     /uscms_data/d3/cihar29/Analysis/CMSSW_8_0_30/src/analysis/ZDilepton/mu_weights.root"
        )
  out=""

  for line in "${lines[@]}" ; do
    out="$out$line\n"
  done

  #parFile="parsMC_${channel}.txt"
  echo -e "$out" | column -t > parsMC.txt

  analyze parsMC.txt mc_weights.txt

done
