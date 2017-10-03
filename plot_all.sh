#!/bin/bash

# execute as ./plot_all.sh channel cut#

args=("$@")

if [ $# -lt 2 ] ; then
  echo "Please provide a channel (mm, ee, em) and cut# (0, 1, ...)"
fi

if [ $# -eq 2 ] ; then

  channel=${args[0]}
  cut=${args[1]}
  dir=""

  if [[ $channel = "mm" ]] ; then
    dir="${dir}mm/"
    dataset="Muon"
  elif [[ $channel = "ee" ]] ; then
    dir="${dir}ee/"
    dataset="Ele"
  else
    dir="${dir}em/"
    dataset="Muon"
  fi

  subplot="ratio"
  if [[ "${subplot}" == "ratio" ]] ; then
    subymin="0"
  else
    subymin="-2"
  fi

  hists=( "dilepmass" "jet0btag" "jet0eta" "jet0pt" "jet1btag" "jet1eta" "jet1pt" "jethT" "lep0eta" "lep0pt" "lep1eta" "lep1pt" "metcorrpt"
          "sT" "masslljjm" "lep0perp" "lep1perp" "lepept" "lepmpt" "metpt" "nGoodEle" "nGoodMuon" "nGoodJet" "minjet0pt" "minjet1pt" "cleanjet0pt" "cleanjet1pt"
          "nbtag" "rl0cleanj" "rl1cleanj" "rl0l1" "rmin0" "rmin1" "sumrmin" "rbl" "masslmin0" "masslmin1" "sT_met" "dphi_jet0met" "dphi_jet1met" )

  xmins=( 0 0 -5 0 0 -5 0 0 -5 0 -5 0 0
          0 0 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 0 -6 -6)

  xmaxs=( 1000 1 5 2000 1 5 2000 3000 5 1000 5 600 1000
          5000 5000 1000 1000 1000 1000 1000 5 5 25 1000 1000 1000 1000
          5 5 5 5 5 5 10 5 2000 2000 5000 6 6)

  for ((i=0;i<${#hists[@]};++i)); do

    logy="false"
    if [[ "${hists[i]}" == *"pt"* || "${hists[i]}" == "sT" || "${hists[i]}" == "jethT" || "${hists[i]}" == "masslljjm" ]] ; then
      logy="true"
    fi

    rebin="2"
    if [[ "${hists[i]}" == *"met"* ]] ; then
      rebin="0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2500 5000"
    elif [[ "${hists[i]}" == "sT" || "${hists[i]}" == "jethT" || "${hists[i]}" == "masslljjm" ]] ; then
      rebin="0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1600 1800 2000 2500 3000 5000"
    elif [[ "${hists[i]}" == *"Good"* || "${hists[i]}" == "nbtag" ]] ; then
      rebin="1"
    fi

    lines=( "dataFileName   ${dir}${dataset}_${channel}.root"
            "dataName       Data"
            "plotData       false"
            "mcFileNames    ${dir}DYhigh_${channel}.root ${dir}DYlow_${channel}.root ${dir}STschannel_${channel}.root ${dir}STtWchannel_${channel}.root ${dir}STtchannel_${channel}.root ${dir}SaTtWchannel_${channel}.root ${dir}SaTtchannel_${channel}.root ${dir}TTbar_${channel}.root ${dir}WJets_${channel}.root ${dir}WW_${channel}.root ${dir}WZ_${channel}.root ${dir}ZZ_${channel}.root"
            "sigFileNames   ${dir}gluon_M-1000_${channel}.root ${dir}gluon_M-1250_${channel}.root ${dir}gluon_M-1500_${channel}.root ${dir}gluon_M-2000_${channel}.root ${dir}gluon_M-2500_${channel}.root ${dir}gluon_M-3000_${channel}.root ${dir}gluon_M-3500_${channel}.root ${dir}gluon_M-4000_${channel}.root ${dir}gluon_M-500_${channel}.root ${dir}gluon_M-750_${channel}.root ${dir}zprime_M-1000_W-10_${channel}.root ${dir}zprime_M-1000_W-100_${channel}.root ${dir}zprime_M-1000_W-300_${channel}.root ${dir}zprime_M-1250_W-125_${channel}.root ${dir}zprime_M-1250_W-12p5_${channel}.root ${dir}zprime_M-1500_W-15_${channel}.root ${dir}zprime_M-1500_W-150_${channel}.root ${dir}zprime_M-3000_W-300_${channel}.root"
            "sigScales      10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10"
            "leftText       CMS"
            "rightText      Run 2016 - 35.9 fb^{-1} (13 TeV)"
            "logx           false"
            "logy           ${logy}"
            "subplot        ${subplot}"
            "outName        ${channel}/${cut}/${hists[i]}"
            "hname          ${cut}_${hists[i]}"
            "xmin           ${xmins[i]}"
            "xmax           ${xmaxs[i]}"
            "ymin           0"
            "subymin        ${subymin}"
            "subymax        2"
            "systematics    topPt_weight jec jer btagSF mistagSF pileup pdf q2ttbar q2dy q2st q2signal"
            "sys_norm       lumi:0.025 sig_st:0.16 sig_db:0.15 mutrig:0.005 muid:0.01 muiso:0.01 eltrig:0.05 elid:0.01 eliso:0.01"
            "rebin          ${rebin}"
            "plotImpact     false"
            "theta          false" # zp1, zp10, zp30, or gkk
          )
    out=""

    for line in "${lines[@]}"
    do
      out="$out$line\n"
    done

    echo -e "$out" > plot_pars.txt

    plot plot_pars.txt

  done

fi
