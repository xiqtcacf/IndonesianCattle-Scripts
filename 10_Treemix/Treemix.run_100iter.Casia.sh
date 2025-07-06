num_procs=150
num_jobs="\j"

FILE=treemix

for i in {0..10}
do
                for j in {1..100}
                do
                                 while (( ${num_jobs@P} >= num_procs )); do
                                        #echo "${num_jobs@P} "
                                        wait -n
                                 done
                nice -n10 treemix -i $FILE.frq.gz -m $i -o m_${i}/$FILE_m$i_${j} -root water_buffalo -bootstrap -k 1000  > m_${i}/treemix_${i}_log &
                done
done
