#!/bin/bash

pop=$1 ### give a list of individuals for specific population e.g. Madura Pesisir

export DIR=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix_new/$pop
cd $DIR
##outgroup
Zebu="Nelore_BINE1,Nelore_BINE2,Bhagnari_23,Gir_BIGI3,Gir_BIGI4,Hariana_Har03,Kangayam_SAMN10131257,Lohani_18,Red_Sindhi_303,Sahiwal_Sha3b,Tharparkar_Thar1,N_22,N_23,N_3,N_45,N_6,N_1,N_11,N_14,N_16,N_17,N_24,N_38_copy,N_8,N_9"
###input bcf per chr
vcf_input=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix_new/bcf_perChr
###good sites, weights
Sites_filtering=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/Good_sites_Bostau9/autosome_rep_het_dep_map.bed
###two admixed pops
admixbcf_Banteng=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix_new/admixPop/Ref.Banteng.314inds.imputated.BosTau9.bcf
admixbcf_Gaur=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix_new/admixPop/Ref.Gaur.314inds.imputated.BosTau9.bcf

input=/home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix_new/$pop

###using like this: bash Using.10000.allChr.PerInd.hmmix.annotation.sh Jabres Banteng
while read sample ; do
  for i in {1..29}
  do
#######use weights, using modified jason
#    hmmix create_outgroup -ind=$Zebu -vcf=$vcf_input/Chr${i}.314inds.imputated.BosTau9.bcf -out=Chr${i}.Weights.outgroup.txt -weights=$Sites_filtering
#    hmmix mutation_rate -outgroup=Chr${i}.Weights.outgroup.txt -window_size=1000000 -out Chr${i}.Weights.mutationrate.bed -weights=$Sites_filtering
#    hmmix create_ingroup -ind=$sample -vcf=$vcf_input/Chr${i}.314inds.imputated.BosTau9.bcf -outgroup=Chr${i}.Weights.outgroup.txt -out=Chr${i}.obs.Weights -weights=$Sites_filtering
#    hmmix train -obs=Chr${i}.obs.Weights.$sample.txt -mutrates=Chr${i}.Weights.mutationrate.bed -param=10000.$pop.Cattle.Initialguesses.json -window_size=10000 -out=Chr${i}.Weights.trained.$sample.json -weights=$Sites_filtering
#    hmmix decode -obs=Chr${i}.obs.Weights.$sample.txt -mutrates=Chr${i}.Weights.mutationrate.bed -param=Chr${i}.Weights.trained.$sample.json -window_size=10000 -out=Chr${i}.Weights.$sample.decoded -weights=$Sites_filtering
    #### below is to explain introgressed snps from which population
     hmmix decode -obs=$input/Chr${i}.obs.Weights.$sample.txt -mutrates=$input/Chr${i}.Weights.mutationrate.bed -param=$input/Chr${i}.Weights.trained.$sample.json -window_size=10000 -out=Banteng.Chr${i}.$sample.decode.annotation -weights=$Sites_filtering -admixpop=$admixbcf_Banteng
     hmmix decode -obs=$input/Chr${i}.obs.Weights.$sample.txt -mutrates=$input/Chr${i}.Weights.mutationrate.bed -param=$input/Chr${i}.Weights.trained.$sample.json -window_size=10000 -out=Gaur.Chr${i}.$sample.decode.annotation -weights=$Sites_filtering -admixpop=$admixbcf_Gaur
  done
 for j in {2..29}
  do
  sed -i '1d' Banteng.Chr${j}.$sample.decode.annotation.diploid.txt
  sed -i '1d' Gaur.Chr${j}.$sample.decode.annotation.diploid.txt
  done
  cat Banteng.Chr1.$sample.decode.annotation.diploid.txt \
      Banteng.Chr2.$sample.decode.annotation.diploid.txt \
      Banteng.Chr3.$sample.decode.annotation.diploid.txt \
      Banteng.Chr4.$sample.decode.annotation.diploid.txt \
      Banteng.Chr5.$sample.decode.annotation.diploid.txt \
      Banteng.Chr6.$sample.decode.annotation.diploid.txt \
      Banteng.Chr7.$sample.decode.annotation.diploid.txt \
      Banteng.Chr8.$sample.decode.annotation.diploid.txt \
      Banteng.Chr9.$sample.decode.annotation.diploid.txt \
      Banteng.Chr10.$sample.decode.annotation.diploid.txt \
      Banteng.Chr11.$sample.decode.annotation.diploid.txt \
      Banteng.Chr12.$sample.decode.annotation.diploid.txt \
      Banteng.Chr13.$sample.decode.annotation.diploid.txt \
      Banteng.Chr14.$sample.decode.annotation.diploid.txt \
      Banteng.Chr15.$sample.decode.annotation.diploid.txt \
      Banteng.Chr16.$sample.decode.annotation.diploid.txt \
      Banteng.Chr17.$sample.decode.annotation.diploid.txt \
      Banteng.Chr18.$sample.decode.annotation.diploid.txt \
      Banteng.Chr19.$sample.decode.annotation.diploid.txt \
      Banteng.Chr20.$sample.decode.annotation.diploid.txt \
      Banteng.Chr21.$sample.decode.annotation.diploid.txt \
      Banteng.Chr22.$sample.decode.annotation.diploid.txt \
      Banteng.Chr23.$sample.decode.annotation.diploid.txt \
      Banteng.Chr24.$sample.decode.annotation.diploid.txt \
      Banteng.Chr25.$sample.decode.annotation.diploid.txt \
      Banteng.Chr26.$sample.decode.annotation.diploid.txt \
      Banteng.Chr27.$sample.decode.annotation.diploid.txt \
      Banteng.Chr28.$sample.decode.annotation.diploid.txt \
      Banteng.Chr29.$sample.decode.annotation.diploid.txt > Banteng.AllChr.$sample.decode.annotation.diploid.txt
  rm Banteng.Chr*.$sample.decode.annotation.diploid.txt
  mv Banteng.AllChr.$sample.decode.annotation.diploid.txt /home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix_new/$pop.annotation

  cat Gaur.Chr1.$sample.decode.annotation.diploid.txt \
      Gaur.Chr2.$sample.decode.annotation.diploid.txt \
      Gaur.Chr3.$sample.decode.annotation.diploid.txt \
      Gaur.Chr4.$sample.decode.annotation.diploid.txt \
      Gaur.Chr5.$sample.decode.annotation.diploid.txt \
      Gaur.Chr6.$sample.decode.annotation.diploid.txt \
      Gaur.Chr7.$sample.decode.annotation.diploid.txt \
      Gaur.Chr8.$sample.decode.annotation.diploid.txt \
      Gaur.Chr9.$sample.decode.annotation.diploid.txt \
      Gaur.Chr10.$sample.decode.annotation.diploid.txt \
      Gaur.Chr11.$sample.decode.annotation.diploid.txt \
      Gaur.Chr12.$sample.decode.annotation.diploid.txt \
      Gaur.Chr13.$sample.decode.annotation.diploid.txt \
      Gaur.Chr14.$sample.decode.annotation.diploid.txt \
      Gaur.Chr15.$sample.decode.annotation.diploid.txt \
      Gaur.Chr16.$sample.decode.annotation.diploid.txt \
      Gaur.Chr17.$sample.decode.annotation.diploid.txt \
      Gaur.Chr18.$sample.decode.annotation.diploid.txt \
      Gaur.Chr19.$sample.decode.annotation.diploid.txt \
      Gaur.Chr20.$sample.decode.annotation.diploid.txt \
      Gaur.Chr21.$sample.decode.annotation.diploid.txt \
      Gaur.Chr22.$sample.decode.annotation.diploid.txt \
      Gaur.Chr23.$sample.decode.annotation.diploid.txt \
      Gaur.Chr24.$sample.decode.annotation.diploid.txt \
      Gaur.Chr25.$sample.decode.annotation.diploid.txt \
      Gaur.Chr26.$sample.decode.annotation.diploid.txt \
      Gaur.Chr27.$sample.decode.annotation.diploid.txt \
      Gaur.Chr28.$sample.decode.annotation.diploid.txt \
      Gaur.Chr29.$sample.decode.annotation.diploid.txt > Gaur.AllChr.$sample.decode.annotation.diploid.txt
  rm Gaur.Chr*.$sample.decode.annotation.diploid.txt
  mv Gaur.AllChr.$sample.decode.annotation.diploid.txt /home/wlk579/2.0Indonesia_Bos_project/Indonesia_Bos_project/5.Analyses_King015/hmmix_new/$pop.annotation
done < ../$pop.pure.inds

