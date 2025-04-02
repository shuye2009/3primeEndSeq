#!/usr/bin/bash

#wd=/home/greenblattlab/shuyepu/Nujhat/3endseq/Dec21_2022/GREENBLATT
wd=/home/greenblattlab/shuyepu/Nujhat/3endseq/June04_2024/GREENBLATT

cd $wd

for f in *_R1_001.fastq.gz; do
	echo $f
#	submitjob -w 30 -c 1 -m 20 /home/greenblattlab/shuyepu/Nujhat/3endseq/scripts/define_CS.sh $f
done

#for f in *_filteredIP.bed; do basef=$(basename $f ".bed"); submitjob -w 2 -m 5 -c 1 bedtools bedtobam -g $HOME/hg19/human.hg19.genome -i $f \> ${basef}.bam; done

submitjob -w 20 -c 1 -m 20 /home/greenblattlab/shuyepu/Nujhat/3endseq/scripts/define_cluster.sh




