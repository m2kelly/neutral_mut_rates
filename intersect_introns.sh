tail -n +2 /home/maria/synon_mut_rates/auxillary/introns.bed | \
bedtools intersect -wa -wb -a stdin -b /home/maria/filter_transcripts/output/intron_cords.bed | \
uniq |\
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $10, $12}' \
> /home/maria/synon_mut_rates/auxillary/introns_genes.bed