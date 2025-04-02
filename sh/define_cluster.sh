#!/usr/bin/bash
# tag2cluster.pl and tag2peak.pl are part of ctk-1.1.4

#tag2cluster.pl combined_IPs.bed IP_cluster.bed -big -s -maxgap 10 -collapse 0 -v
#awk 'BEGIN{OFS="\t"} {if($5>2) print}' IP_cluster.bed > IP_cluster_3.bed


#tag2cluster.pl combined_CSs.bed SC_cluster.bed -big -s -maxgap 10 -collapse 0 -v
#awk 'BEGIN{OFS="\t"} {if($5>2) print}' SC_cluster.bed > SC_cluster_3.bed

#tag2peak.pl combined_CSs.bed SC_peak_valleySeeking.bed -big -ss --valley-seeking --valley-depth 0.9 --out-boundary SC_peak_boundary.bed --dbkey hg19 -p 0.1 --multi-test -minPH 3 -maxPH -1 -gap -1 --prefix CSpeak -v

cat *_CS.bed | bedtools sort -i - | bedtools merge -s -d -1 -c 4,5,6 -o count,count,distinct -i - | awk 'BEGIN{OFS="\t"} {if($5>4 && $1 != "chrM") print}' | sort -k5,5 -nr > merged_CSs.bed

cat *_IP.bed | bedtools sort -i - | bedtools merge -s -d -1 -c 4,5,6 -o count,count,distinct -i - | awk 'BEGIN{OFS="\t"} {if($5>4 && $1 != "chrM") print}' | sort -k5,5 -nr > merged_IPs.bed

