while read breed; do
echo $breed

Rscript myhet.R ${breed}/*AC* > ${breed}/het.txt

paste <(grep homZebu ${breed}/het.txt ) <(grep homBanteng ${breed}/het.txt) <(grep Mix ${breed}/het.txt) <(grep -v 'Mix\|hom' ${breed}/het.txt) -d" " > ${breed}/het_4cols.txt

nice Rscript get_plot.R ${breed}  &

done < 5breeds.txt










