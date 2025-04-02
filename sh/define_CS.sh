#!/usr/bin/bash

f=$1

basef=$(basename $f "_R1_001.fastq.gz")
outfile=${basef}_trimmed.fastq.gz
report=${basef}_trimReport.txt

if [[ 1 == 1 ]]; then
zcat $f | cutadapt -g "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" -o $outfile -q 20 -m 25 -O 19 -j 20 - > $report # 52 Ts

bowtie2 --fast -N 1 -p 20 -x $HOME/hg19/hg19 -U $outfile -S ${basef}.sam
	
	# filter out unmapped -F 4 and multimapped -q 10 reads, -b output bam, -T reference sequence file, -h include header
samtools view -h -F 4 -q 10 -b -T $HOME/hg19/hg19.fa ${basef}.sam | samtools sort -o ${basef}_filtered.bam -
samtools index ${basef}_filtered.bam
	# convert to bed to use bedtools 
bedtools bamtobed -i ${basef}_filtered.bam > ${basef}_filtered.bed

	# get upstream 20 flanking regions
bedtools flank -i ${basef}_filtered.bed -g $HOME/hg19/human.hg19.genome -l 20 -r 0 -s > ${basef}_filtered_flankL20.bed

	# get flanking sequences and translate lower case to uppper case
bedtools getfasta -bed ${basef}_filtered_flankL20.bed -fi $HOME/hg19/hg19.fa -name -s | tr [a-z] [A-Z] > ${basef}_filtered_flankL20.fa

	# get flanking sequences that pass the polyT filter, 5 Ts out of 10 bases (or 10 consecutive Ts with 5 mismatch)
cat ${basef}_filtered_flankL20.fa | cutadapt -g "TTTTTTTTTT" -o ${basef}_filtered_flankL20_noT.fa -e 3 -O 10 -j 20 --discard-trimmed - > ${basef}_noT_report.txt

	# get the reads id of flanking sequence that pass the polyT filter
cat ${basef}_filtered_flankL20_noT.fa | grep ">" | sed 's/::/@/g' | cut -f1 -d@ | sed 's/>//g' > ${basef}_CS_read.txt
fi
	# get cleavage sites: only keep read without polyT in upstream, then switch strand so the bed has the same strand as gene, finally get the location of 3'end
join -1 1 -2 4 -a 1 <(sort -k1,1 ${basef}_CS_read.txt) <(sort -k4,4 ${basef}_filtered.bed) | awk 'BEGIN{OFS="\t"} {if($6=="+") $6="-"; else if($6=="-") $6="+"; print $2, $3, $4, $1, $5, $6}' | bedtools sort -i - > ${basef}_filteredCS.bed 
cat ${basef}_filteredCS.bed | awk 'BEGIN{OFS="\t"} {if($6=="-") $3=$2+1; else if($6=="+") $2=$3-1; print}' > ${basef}_CS.bed

	# get internal priming (IP) sites
removeRow.pl -q 3 -f 3 ${basef}_filtered.bed ${basef}_filteredCS.bed | awk 'BEGIN{OFS="\t"} {if($6=="+") $6="-"; else if($6=="-") $6="+"; print}' > ${basef}_filteredIP.bed
cat ${basef}_filteredIP.bed | awk 'BEGIN{OFS="\t"} {if($6=="-") $3=$2+1; else if($6=="+") $2=$3-1; print}' > ${basef}_IP.bed


