# K4 multicov
samtools view -c -F 0x904 H3K4me3_Rep1.bam > H3K4me3_Rep1_total.txt
samtools view -c -F 0x904 H3K4me3_Rep2.bam > H3K4me3_Rep2_total.txt
bedtools multicov -bed RenLab-H3K4me3-heart-peaks.bed -bams H3K4me3_Rep1.bam H3K4me3_Rep2.bam > H3K4me3_peak_count.txt


#bivalent:H3K4me3 low & H3K27me3
grep -v -E '^(#|track|browser)' H3K4me3_low.bed | awk 'BEGIN{OFS="\t"} NF>=3{print $1,$2,$3}' | sort -k1,1 -k2,2n > K4_low.3col.sorted.bed
bedtools intersect -a K4_low.3col.sorted.bed -b K27.3col.sorted.bed -wo > K4_low_K27_overlap_wo.txt
awk 'BEGIN{OFS="\t"}{s=($2>$5?$2:$5);e=($3<$6?$3:$6);if(e-s>=200)print $1,s,e}' K4_low_K27_overlap_wo.txt | sort -k1,1 -k2,2n | bedtools merge -i - > overlap_abs200.merged_K4_low_K27.bed



