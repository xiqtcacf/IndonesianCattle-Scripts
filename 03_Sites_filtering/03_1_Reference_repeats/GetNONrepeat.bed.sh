python3 XiaodongLiuWritten.extSoftRepeatFA.py bosTau9.fasta.masked > bosTau9.fasta.masked.repeat.bed
/home/users/xi/software/bedtools2/bin/bedtools subtract -a bosTau9.fasta.fai.bed -b bosTau9.fasta.masked.repeat.bed > bosTau9.fasta.masked.NOrepeat.bed

python3 XiaodongLiuWritten.extSoftRepeatFA.py water_buffalo.reference.fna.masked > water_buffalo.reference.fna.masked.repeat.bed
/home/users/xi/software/bedtools2/bin/bedtools subtract -a water_buffalo.reference.fna.fai.bed -b water_buffalo.reference.fna.masked.repeat.bed > water_buffalo.reference.fna.masked.NOrepeat.bed
