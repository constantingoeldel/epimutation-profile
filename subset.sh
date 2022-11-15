# !/bin/bash
for i in /mnt/nas/zhilin/others/MA3_new_total_original_methylome/methylome_*.txt;
do echo Filtering methylome for CGs in $i \
awk 'BEGIN{OFS="\t"}{if($0~/^1|^2|^3|^4|^5/){if($4=="CG") {print $1,$2,$2+1,$4,0,$3,$5,$6,$7,$8,$9}}}' $i 
> ./$i.change.bed; done

for i in methylome_*.txt.change.bed; do bedtools intersect -a $i -b /mnt/nas/zhilin/others/gbM_gene_anotation_extract_Arabidopsis.bed  -wa -wb -f 1
|bedtools groupby -i - -g 1,2,6,4,7,8,9,10,11 -c 12 -o collapse 
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' 
| sed '1i\seqnames\tstart\tstrand\tcontext\tcounts.methylated\tcounts.total\tposteriorMax\tstatus\trc.meth.lvl' - > 
./$i.within_CS-clock.txt; done