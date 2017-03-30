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

  hists=( "dilepmass" "jet0btag" "jet0eta" "jet0pt" "jet1btag" "jet1eta" "jet1pt" "jethT" "lep0eta" "lep0perp" "lep0pt" "lep1eta" "lep1perp" "lep1pt"
          "jet0phi" "jet1phi" "lep0phi" "lep1phi" "lepept" "lepmpt" "metpt" "metcorrpt" "sT" "nEle" "nGoodEle" "nGoodJet" "nGoodMuon" "nJet" "nMuon"
          "nbtag" "rl0j" "rl1j" "rl0l1" "rmin0" "rmin1" "rbl" )

  xmins=( 0 0 -5 0 0 -5 0 0 -5 0 0 -5 0 0
          -4 -4 -4 -4 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 )
  xmaxs=( 400 1 5 750 1 5 500 1000 5 200 500 5 200 300
          4 4 4 4 300 300 300 300 1500 5 5 15 5 30 5 30 5
          5 5 5 5 1 1 5 )

  ymins=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 )

  subymins=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 )
  subymaxs=( 3 3 3 3 3 3 3 3 3 3 3 3 3 3 
             3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
             3 3 3 3 3 3 3 )

  for ((i=0;i<${#hists[@]};++i)); do

    lines=( "dataFileName   MuonGH_${channel}.root"
            "dataName       Data"
            "mcFileNames    TTbar_${channel}.root highDY_${channel}.root lowDY_${channel}.root STschannel_${channel}.root STtchannel_${channel}.root SaTtchannel_${channel}.root STtWchannel_${channel}.root SaTtWchannel_${channel}.root Wjet_${channel}.root"
            "mcScales       1 1 1 1 1 1 1 1 1"
            "leftText       CMS"
            "rightText      Run 2016GH - 16.1 fb^{-1} (13 TeV)"
            "logx           false"
            "subplot        ratio"
            "outName        MuonGH/${hists[i]}"
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
