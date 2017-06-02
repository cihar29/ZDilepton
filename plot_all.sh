#!/bin/bash

# execute as ./plot_all.sh channel cut#

args=("$@")

if [ $# -eq 0 -o $# -eq 1 ] ; then
  echo "Please provide a channel (mm, ee, em) and cut# (0, 1, ...)"
fi

if [ $# -eq 2 ] ; then

  channel=${args[0]}
  cut=${args[1]}
  scale="1"

  hists=( "dilepmass" "jet0btag" "jet0eta" "jet0pt" "jet1btag" "jet1eta" "jet1pt" "jethT" "lep0eta" "lep0pt" "lep1eta" "lep1pt" "metcorrpt" "sT" "masslljjm" )
         # "lep0perp" "lep1perp" "lepept" "lepmpt" "metpt" "nGoodEle" "nGoodMuon" "nGoodJet" "minjet0pt" "minjet1pt"
         # "cleanjet0pt" "cleanjet1pt" "nbtag" "rl0cleanj" "rl1cleanj" "rl0l1" "rmin0" "rmin1" "rbl" "masslmin0" "masslmin1" )

  xmins=( 0 0 -5 0 0 -5 0 0 -5 0 -5 0 0 0 0
          0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 0)
  xmaxs=( 400 1 5 750 1 5 500 1000 5 500 5 300 300 1000 3000
          200 200 300 300 300 5 5 15 1000 1000 1000 1000
          5 5 5 5 5 5 5 500 500 )

  ymins=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 0 )

  subymins=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0 )
  subymaxs=( 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
             3 3 3 3 3 3 3 3 3 3
             3 3 3 3 3 3 3 3 3 3 3 )

  for ((i=0;i<${#hists[@]};++i)); do

    lines=( "dataFileName   Muon_${channel}.root"
            "dataName       Data"
            "mcFileNames    TTbar_${channel}.root highDY_${channel}.root lowDY_${channel}.root STschannel_${channel}.root STtchannel_${channel}.root SaTtchannel_${channel}.root STtWchannel_${channel}.root SaTtWchannel_${channel}.root Wjet_${channel}.root"
            "sigFileNames   zprime-M3000-W300_${channel}.root gluonkk-M3000_${channel}.root"
            "sigScales      1 1"
            "leftText       CMS"
            "rightText      Run 2016 - 35.9 fb^{-1} (13 TeV)"
            "logx           false"
            "subplot        ratio"
            "outName        ${hists[i]}"
            "hname          ${cut}_${hists[i]}"
            "xmin           ${xmins[i]}"
            "xmax           ${xmaxs[i]}"
            "ymin           ${ymins[i]}"
            "subymin        ${subymins[i]}"
            "subymax        ${subymaxs[i]}"
          )
    line=""

    for j in "${lines[@]}"
    do
      line="$line$j\n"
    done

    echo -e "$line" > plot_pars.txt

    plot plot_pars.txt

  done

fi
