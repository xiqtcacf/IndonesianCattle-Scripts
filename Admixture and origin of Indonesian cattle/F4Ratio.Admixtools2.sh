#!/bin/bash

batch=$1

export DIR=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/qpadmixture/using_Taurine/$batch
cd $DIR

inds=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/qpadmixture/using_Taurine/$batch/$batch

while read file ; do
  Rscript F4Ratio.Admixtools2.admixtureP.2refSource.Waterbuffol.R $file
  mv EUTaurine.weights.tsv $file.EUTaurine.weights.tsv
done < $inds
