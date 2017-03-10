#!/bin/bash

# execute as ./plot_all.sh

args=("$@")

if [ $# -eq 0 ] ; then

  hists=( "nEleDiff" "nMuonDiff" "nJetDiff" "nEle" "nMuon" "nJet" "nGoodEle" "nGoodMuon" "nGoodJet" "lep0pt" "lep0eta" "lep1pt" "lep1eta" "dilepmass"
          "jet0pt" "jet0eta" "jet0btag" "jet1pt" "jet1eta" "jet1btag" "jethT" "metpt" "sT" )
  xmins=( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5, 0, -5, 0, 0, -5, 0, 0, -5, 0, 0, 0, 0 )
  xmaxs=( 10, 10, 20, 10, 10, 20, 10, 10, 20, 500, 5, 500, 5, 400, 500, 5, 1, 500, 5, 1, 500, 500, 500 )
  ymins=( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
  ymaxs=( 50000, 50000, 50000, 50000, 50000, 50000, 50000, 50000, 50000, 30000, 20000, 30000, 20000, 100000, 50000, 50000, 50000,
          50000, 50000, 50000, 50000, 50000, 50000 )
  subymins=(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 )
  subymaxs=(2, 2, 3, 2, 2, 3, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 )

  for ((i=0;i<${#hists[@]};++i)); do

    lines=( "dataFileName   MuonGH_mm.root"
            "dataName       Data"
            "mcFileNames    TTbar_mm.root highDY_mm.root lowDY_mm.root STschannel_mm.root STtchannel_mm.root SaTtchannel_mm.root STtWchannel_mm.root SaTtWchannel_mm.root Wjet_mm.root"
            "mcScales       1 1 1 1 1 1 1 1 1"
            "leftText       CMS"
            "rightText      Run 2016GH - 16.1 fb^{-1} (13 TeV)"
            "logx           false"
            "subplot        ratio"
            "outName        MuonGH/${hists[i]}"
            "hname          0_${hists[i]}"
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
