#!/usr/bin/bash

#######################################################
#		Data download, quality control, mapping		  #
#######################################################

###### Build index for bovine genome
extract_splice_sites.py Bos_taurus.UMD3.1.92.gtf >bosgenome.ss
extract_exons.py Bos_taurus.UMD3.1.92.gtf >bosgenome.exon
hisat2-build –p 10 --ss genome.ss --exon genome.exon bosgenome.fa bosgenome

###### Get all id
x=$(cat SRR_Acc_List.txt)
for i in $x
do 
###### Download SRA data
prefetch "$i" -O ./ -X 100000000
fastq-dump --split-3 --gzip "$i".sra

###### Quality control
java -jar /Your/path/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 "$i"_1.fastq.gz "$i"_2.fastq.gz "$i"_1.clean.fq.gz "$i"_1_unpaired.fastq.gz "$i"_2.clean.fq.gz "$i"_2_unpaired.fastq.gz \
ILLUMINACLIP:/Your/path/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

###### mapping to genome 
hisat2  -p 10 -x /Bos/bosgenome -1 "$i"_1.clean.fq.gz -2 "$i"_2.clean.fq.gz -S "$i".sam
samtools view -bS "$i".sam -o "$i"_unsorted.bam
samtools sort "$i"_unsorted.bam  -o "$i"_sorted.bam

#######################################################
#				Gene quantification					  #
#######################################################

stringtie -p 20 -e -B -G ../Bos_taurus.UMD3.1.92.gtf -o ./"$i"/"$i".gtf  -A ./Expression/"$i".tsv "$i"_sorted.bam
awk '{print $1"\t"$8}' "$i".tsv|sed '1d'|sort -k1,1>"$i".fpkm
echo -e "ID\t"$i""|cat - "$i".fpkm >"$i".fpkm.txt
rm "$i".fpkm
done
########merge all expression together##############
paste *.fpkm.txt|awk '{printf("%s\t",$1);for(i=2;i<=NF;i+=2){printf("%s\t",$i)};print ""}'>All.gene.expression.fpkm.txt



#######################################################
#				 Detecting Splicing					  #
#######################################################
######transcript assembly
x=$(cat SRR_Acc_List.txt)
for i in $x
do 
stringtie -p 4  -G Bos_taurus.UMD3.1.92.gtf -o ./"$i"/"$i".gtf  -l "$i" "$i"_sorted.bam
done
find ./ -name "SRR*.gtf" > list.txt
cuffcompare -i list.txt -r ../Bos_taurus.UMD3.1.92.gtf -o All
cat ./*/*.gtf>All.gtf
cat ./*/*.tmap>All.tmap

##########Modify gtf style
y=$(sed '/^#/d' All.gtf|awk '{print $1}'|sort|uniq)
for i in $y
do
awk '{if($1=="'"$i"'")print}' All.gtf>"$i".gtf
awk '{if($3=="transcript")print}' "$i".gtf|awk 'BEGIN{FS = ";"}{for (f=1; f <= NF; f+=1) {if ($f ~ /transcript_id/) {printf $f"\t"}}}{for (g=1; g <= NF; g+=1){if($g~/FPKM/){print $g}}}'|sed 's/transcript_id \"//g' |sed 's/"//'>"$i".expression
perl extract_expression.pl "$i".gtf "$i".expression>"$i".convert.gtf
awk -F ";" '{if($4~/FPKM/) print $1";"$2";"$4";"$3";"$5";"$6;else print $1";"$2";"$3";"$NF";"$4";"$5";"$6}' "$i".convert.gtf|sed '$d'>"$i".FPKM.gtf
done
cat *.FPKM.gtf>All.FPKM.gtf

########## extract AS events
./ASprofile.b-1.0.4/extract-as All.gtf Bos_taurus.hdrs -r All.tmap Bos_taurus.UMD3.1.92.gtf>all.tmap.as
########## summarize events, per gene; also, create a catalog of non-redundant events 
perl ./ASprofile.b-1.0.4/summarize_as.pl All.gtf all.tmap.as -p bovine

########## calculate fpkm of events from transcripts in each sample
for i in $x
do
grep $i All.FPKM.gtf > $i.variants.gtf
extract-as-fpkm $i.variants.gtf Bos_taurus.hdrs bovine.as.nr -W 9 > $i.W9.fpkm
done

