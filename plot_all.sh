#!/bin/bash

# execute as ./plot_all.sh cut#

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a cut number (0, 1, or 2)."
fi

if [ $# -eq 1 ] ; then

  cut=${args[0]}
  channel="em"

  hists=( "dilepmass" "jet0btag" "jet0eta" "jet0pt" "jet1btag" "jet1eta" "jet1pt" "jethT" "lep0eta" "lep0perp" "lep0pt" "lep1eta" "lep1perp" "lep1pt" "metcorrpt"
          "metpt" "muonD0" "muonDz" "nEle" "nEleDiff" "nGoodEle" "nGoodJet" "nGoodMuon" "nJet" "nJetDiff" "nMuon" "nMuonDiff" "nbtag" "rl0j" "rl1j" "rl0l1"
          "rmin0" "rmin1" "sT" )

  xmins=( 0 0 -5 0 0 -5 0 0 -5 0 0 -5 0 0 0
          0 -2 -2 0 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 )
  ymins=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
          0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
          0 0 0 )

  subymins=( 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
             0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
             0 0 0 )
  subymaxs=( 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
             3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
             3 3 3 )

  if [ "$cut" = "0" ] ; then
    xmaxs=( 400 1 5 500 1 5 500 1000 5 200 500 5 200 500 300
            300 2 2 5 5 5 30 5 30 10 5 5 5 1 1 5
            1 1 1000 )
    ymaxs=( 100000 16000 10000 50000 16000 10000 50000 25000 10000 70000 25000 10000 150000 50000 30000
            30000 300000 300000 300000 300000 300000 30000 300000 30000 300000 300000 300000 250000 40000 40000 5000
            70000 70000 12000 )
  fi

  if [ "$cut" = "1" ] ; then
    xmaxs=( 300 1 5 400 1 5 400 800 5 100 400 5 100 300 300
            300 2 2 5 5 5 30 5 30 10 5 5 5 1 1 5
            1 1 1000 )
    ymaxs=( 20000 10000 4000 15000 10000 4000 10000 5000 4000 30000 10000 4000 40000 15000 10000
            10000 100000 100000 100000 100000 100000 10000 100000 10000 100000 100000 100000 100000 20000 20000 5000
            30000 30000 4000 )
  fi

  if [ "$cut" = "2" ] ; then
    xmaxs=( 300 1 5 400 1 5 400 800 5 100 400 5 50 200 300
            300 2 2 5 5 5 30 5 30 10 5 5 5 1 1 5
            1 1 1000 )
    ymaxs=( 4000 7000 2000 8000 7000 2000 5000 2000 2000 20000 4000 3000 40000 10000 5000
            5000 50000 50000 50000 50000 50000 10000 50000 10000 50000 50000 50000 50000 5000 5000 5000
            20000 20000 4000 )
  fi

  for ((i=0;i<${#hists[@]};++i)); do

    lines=( "dataFileName   MuonG_${channel}.root"
            "dataName       Data"
            "mcFileNames    TTbar_${channel}.root highDY_${channel}.root lowDY_${channel}.root STschannel_${channel}.root STtchannel_${channel}.root SaTtchannel_${channel}.root STtWchannel_${channel}.root SaTtWchannel_${channel}.root Wjet_${channel}.root"
            "mcScales       1 1 1 1 1 1 1 1 1"
            "leftText       CMS"
            "rightText      Run 2016G - 7.54 fb^{-1} (13 TeV)"
            "logx           false"
            "subplot        ratio"
            "outName        MuonGH/${hists[i]}"
            "hname          ${cut}_${hists[i]}"
            "xmin           ${xmins[i]}"
            "xmax           ${xmaxs[i]}"
            "ymin           ${ymins[i]}"
            "ymax           ${ymaxs[i]}"
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
