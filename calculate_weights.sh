#!/bin/bash

# execute as ./calculate_weights.sh lumi_value

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a lumi value."
fi

if [ $# -eq 1 ] ; then

  lumi=${args[0]}

  dataset=(
    "DYhigh"
    "DYlow"
    "STschannel"
    "STtWchannel"
    "STtchannel"
    "SaTtWchannel"
    "SaTtchannel"
    "TTbar_incl"
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
    "zprime_M-3000_W-900"
    "TTbar0-700_topPtWeightUP"    # no weights
    "TTbar700-1000_topPtWeightUP"
    "TTbar1000-inf_topPtWeightUP"
    "TTbar0-700"                  # topPtWeightNOM
    "TTbar700-1000"
    "TTbar1000-inf"
    "TTbar0-700_topPtWeightDN"
    "TTbar700-1000_topPtWeightDN"
    "TTbar1000-inf_topPtWeightDN"
    "TTbar0-700_pdfUP"
    "TTbar700-1000_pdfUP"
    "TTbar1000-inf_pdfUP"
    "TTbar0-700_pdfDN"
    "TTbar700-1000_pdfDN"
    "TTbar1000-inf_pdfDN"
    "TTbar0-700_q2UP"
    "TTbar700-1000_q2UP"
    "TTbar1000-inf_q2UP"
    "TTbar0-700_q2DN"
    "TTbar700-1000_q2DN"
    "TTbar1000-inf_q2DN"
  )

  xs=( 5765400 18610000 11360 35850 136020 35850 80950 831760 61526700 118700 47100 16500
       1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000
       1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000
       729955 80466.2 21338.6 734268 77540.2 19952.2 738224 74795.9 18740.1 733443 77969.8 20346.8 734904 77233.4 19622.3 731539 79256.9 20964.3 737582 75409.6 18768.5 )

  events=( 145407553 35291457 2989179 6952779 67240483 6933037 38810801 154947474 86731469 7981093 3995807 1988085
           98558 99995 99993 99995 99990 99740 99488 99113 99580 98404 99997 100000 103781 101055 79474 96841 102827 99684 111101 189120
           102459 101188 94011 100817 202525 104107 113996 100268 96063 91874 90960 100989 107851 102388 88030 100229 84470 111944 107118 91015 89110
           135982394 53567728 28536040 135982394 53567728 28536040 135982394 53567728 28536040 135982394 53567728 28536040 135982394 53567728 28536040 135982394 53567728 28536040 135982394 53567728 28536040 )

  file="dataset lumi xs lumi*xs events weight"

  for ((i=0;i<${#dataset[@]};++i)); do

    lumi_xs=$(echo "$lumi*${xs[i]}" | bc)

    line="${dataset[i]} $lumi ${xs[i]} $lumi_xs ${events[i]} $(echo "scale=5; $lumi_xs/${events[i]}" | bc)"

    file="$file\n$line"

  done

  echo -e "$file" | column -t > mc_weights.txt

fi
