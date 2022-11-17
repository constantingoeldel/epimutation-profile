# !/bin/bash
cd /mnt/extStorage/constantin/epimutation-profile/methylome/MA3_new_total_original_methylome
for i in methylome_*.txt; \
do awk 'BEGIN{OFS="\t"}{if($0~/^1|^2|^3|^4|^5/){if($4=="CG") {print $1,$2,$2+1,$4,0,$3,$5,$6,$7,$8,$9}}}' $i \
> /mnt/extStorage/constantin/epimutation-profile/methylome/cg_bed/$i.change.bed; done
echo Done filtering CGs
cd /mnt/extStorage/constantin/epimutation-profile/methylome/cg_bed
for i in methylome_*.txt.change.bed; do bedtools intersect -a $i -b /mnt/extStorage/constantin/epimutation-profile/gbM_gene_anotation_extract_Arabidopsis.bed  -wa -wb -f 1 \
|bedtools groupby -i - -g 1,2,6,4,7,8,9,10,11 -c 12 -o collapse \
| awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
| sed '1i\seqnames\tstart\tstrand\tcontext\tcounts.methylated\tcounts.total\tposteriorMax\tstatus\trc.meth.lvl' - > \
/mnt/extStorage/constantin/epimutation-profile/methylome/within_gbM_genes/$i.within_CS-clock.txt; done