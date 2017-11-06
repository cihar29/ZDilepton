#!/bin/bash

# execute as ./plot_all.sh channel cut#

args=("$@")

if [ $# -lt 2 ] ; then
  echo "Please provide a channel (mm, ee, em) and cut# (0, 1, ...)"
fi

if [ $# -eq 2 ] ; then

  channel=${args[0]}
  cut=${args[1]}
  dir="off/"

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

  fit=true

  hists=( "dilepmass" "jet0btag" "jet0eta" "jet0pt" "jet1btag" "jet1eta" "jet1pt" "jethT" "lep0eta" "lep0pt" "lep1eta" "lep1pt" "metcorrpt" "sT" "masslljjm"
          "lep0perp" "lep1perp" "lepept" "lepmpt" "metpt" "nGoodEle" "nGoodMuon" "nGoodJet" "minjet0pt" "minjet1pt" "cleanjet0pt" "cleanjet1pt" "nbtag"
          "rl0cleanj" "rl1cleanj" "rl0l1" "rmin0" "rmin1" "sumrmin" "rbl" "masslmin0" "masslmin1" "sT_met" "dphi_jet0met" "dphi_jet1met" "nPV" "lep0perp_in" "lep1perp_in" )

  xmins=( 0 0 -5 0 0 -5 0 0 -5 0 -5 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)

  xmaxs=( 1000 1 5 2000 1 5 2000 3000 5 1000 5 600 1000 5000 5000
          200 200 1000 1000 1000 5 5 25 1000 1000 1000 1000 5
          5 5 5 5 5 10 5 2000 2000 5000 3 3 80 200 200)

  for ((i=0;i<${#hists[@]};++i)); do

    if [ "$fit" = true ] ; then
      if [[ "${hists[i]}" != "jet"*"pt" && "${hists[i]}" != "lep0pt" && "${hists[i]}" != "lep1pt" && "${hists[i]}" != "sT_met" ]] ; then
        continue
      fi
      if [[ $cut != 3 && $cut != 6 ]] ; then
        continue
      fi
    fi

    logy="false"
    if [[ "${hists[i]}" == *"pt"* || "${hists[i]}" == "sT"* || "${hists[i]}" == "jethT" || "${hists[i]}" == "masslljjm" ]] ; then
      logy="true"
    fi

    rebin="2"
    if [[ "${hists[i]}" == "metcorrpt" || "${hists[i]}" == "metpt" ]] ; then
      rebin="0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2500 5000"
    elif [[ "${hists[i]}" == "sT" || "${hists[i]}" == "jethT" || "${hists[i]}" == "masslljjm" ]] ; then
      rebin="0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1600 1800 2000 2500 3000 5000"
    elif [[ "${hists[i]}" == *"Good"* || "${hists[i]}" == "nbtag" || "${hists[i]}" == *"perp"* ]] ; then
      rebin="1"
    elif [[ "${hists[i]}" == "sT_met" ]] ; then
      if [[ $channel = "mm" ]] ; then
        rebin="0 20 40 60 80 100 120 140 160 180 200 220 240 260 300 340 380 420 460 500 540 580 620 660 700 740 780 820 860 900 940 980 1020 1060 1100 1140 1180 1220 1260 1300 1340 1380 1420 1460 1500 1540 1580 1620 1660 1700 1740 1780 1820 1860 1900 1940 1980 2020 2060 2100 2140 2180 2220 2260 2340 2420 2500 2660 2820 5000"
      elif [[ $channel = "ee" ]] ; then
        rebin="0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 320 360 400 440 480 520 560 600 640 680 720 760 800 840 880 920 960 1000 1040 1080 1120 1160 1200 1240 1280 1320 1360 1400 1440 1480 1520 1560 1600 1640 1680 1720 1760 1800 1840 1880 1920 2000 2080 2240 2560 5000"
      else
        rebin="0 20 40 60 80 100 120 140 160 180 200 220 240 260 300 340 380 420 460 500 540 580 620 660 700 740 780 820 860 900 940 980 1020 1060 1100 1140 1180 1220 1260 1300 1340 1380 1420 1460 1500 1540 1580 1620 1660 1700 1740 1780 1820 1860 1900 1940 1980 2020 2060 2100 2140 2180 2220 2260 2340 2420 2500 2580 2740 5000"
      fi
    fi

    lines=( "dataFileName   ${dir}${dataset}_${channel}.root"
            "dataName       Data"
            "plotData       true"
            "mcFileNames    ${dir}TTbar0-700_${channel}.root ${dir}TTbar700-1000_${channel}.root ${dir}TTbar1000-inf_${channel}.root ${dir}DYhigh_${channel}.root ${dir}DYlow_${channel}.root ${dir}STschannel_${channel}.root ${dir}STtWchannel_${channel}.root ${dir}STtchannel_${channel}.root ${dir}SaTtWchannel_${channel}.root ${dir}SaTtchannel_${channel}.root ${dir}WJets_${channel}.root ${dir}WW_${channel}.root ${dir}WZ_${channel}.root ${dir}ZZ_${channel}.root"
            "sigFileNames   ${dir}gluon_M-4500_${channel}.root ${dir}gluon_M-5000_${channel}.root ${dir}gluon_M-1000_${channel}.root ${dir}gluon_M-1250_${channel}.root ${dir}gluon_M-1500_${channel}.root ${dir}gluon_M-2000_${channel}.root ${dir}gluon_M-2500_${channel}.root ${dir}gluon_M-3000_${channel}.root ${dir}gluon_M-3500_${channel}.root ${dir}gluon_M-4000_${channel}.root ${dir}gluon_M-500_${channel}.root ${dir}gluon_M-750_${channel}.root ${dir}zprime_M-1000_W-10_${channel}.root ${dir}zprime_M-1000_W-100_${channel}.root ${dir}zprime_M-1000_W-300_${channel}.root ${dir}zprime_M-1250_W-125_${channel}.root ${dir}zprime_M-1250_W-12p5_${channel}.root ${dir}zprime_M-1500_W-15_${channel}.root ${dir}zprime_M-1500_W-150_${channel}.root ${dir}zprime_M-3000_W-300_${channel}.root ${dir}zprime_M-500_W-5_${channel}.root ${dir}zprime_M-500_W-50_${channel}.root ${dir}zprime_M-750_W-7p5_${channel}.root ${dir}zprime_M-750_W-75_${channel}.root ${dir}zprime_M-2000_W-20_${channel}.root ${dir}zprime_M-2000_W-200_${channel}.root ${dir}zprime_M-2000_W-600_${channel}.root ${dir}zprime_M-2500_W-25_${channel}.root ${dir}zprime_M-2500_W-250_${channel}.root ${dir}zprime_M-3000_W-30_${channel}.root ${dir}zprime_M-3500_W-35_${channel}.root ${dir}zprime_M-3500_W-350_${channel}.root ${dir}zprime_M-4000_W-40_${channel}.root ${dir}zprime_M-4000_W-400_${channel}.root ${dir}zprime_M-4000_W-1200_${channel}.root ${dir}zprime_M-4500_W-45_${channel}.root ${dir}zprime_M-4500_W-450_${channel}.root ${dir}zprime_M-5000_W-50_${channel}.root ${dir}zprime_M-5000_W-500_${channel}.root ${dir}zprime_M-5000_W-1500_${channel}.root"
            "sigScales      10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10"
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
            "systematics    topPtWeight jec jer btagSF mistagSF pileup pdf q2ttbar q2dy q2st q2signal"
            "sys_norm       lumi:0.025 sig_st:0.16 sig_db:0.15 mutrig:0.005 muid:0.01 muiso:0.01 eltrig:0.01 elid:0.01 eliso:0.01"
            "rebin          ${rebin}"
            "plotImpact     false"
            "theta          false" # zp1, zp10, zp30, or gkk
            "fit            $fit"
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
