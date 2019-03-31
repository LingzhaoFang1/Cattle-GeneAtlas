#############################################################################################################
#TF enrichment analysis for promoter regions ( -1500bp ~ +500 bp to TSS) of tissue-specific expressed genes (t-value: top 500) using meme software
for i in `ls *.fasta`
do
name=`echo $i | awk -F "." '{print $1}'`
meme $i -oc ${name} -neg all.gene.promoter.-1500+500.fasta -dna -mod zoops -objfun de -nmotifs 8 
cd ${name}
tomtom -oc TF -no-ssc -min-overlap 5 -dist pearson -evalue -thresh 10.0 -png -eps meme.html ../meme-5.0.1/db/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme ../meme-5.0.1/db/motif_databases/EUKARYOTE/jolma2013.meme ../meme-5.0.1/db/motif_databases/MOUSE/uniprobe_mouse.meme 
cd ../
done

#############################################################################################################
#Bar plot for methylation levels in tissues for tissue-specific expressed genes (t-value: top 500) 

# Filter the top genes occured in more than two tissues. 
for n in `ls *.t-value.txt.s`; do name=`echo $n |awk -F "." '{print $1}'`; sed -n '1,500p' $n > ${name}.top.10.txt; done 
cat *.top.10.txt | cut -f 1 |awk '{count[$1]++}END{for (i in count) if (count[i]<=2) print i}' > gene.less3.txt 
for p in `ls *.top.10.txt`; do grep -f gene.less3.txt $p> $p.uniq; done

#plot the barplot
for j in `ls *.top.10.txt.uniq`
do
 name=`echo $j | awk -F "." '{print $1}' ` 
 sort -k1,1 ${name}.top.10.txt.uniq > /top10methylation/methylation/top.10.txt.s
 for i in `ls -d *`
  do
  cd /tissues/${i}
  rm join.*
    for k in `ls *.promoter.-1500+500.bed.meth`
     do
     a=$(awk '{sum+=$5}END{print sum/NR}' ${k}) 
     awk -F "[:\t]" 'OFS="\t"{print $4,$9}' ${k} |sort -k 1 > join.${k}
     join -j 1 -o 1.1,1.2,2.2 join.${k} /top10methylation/methylation/top.10.txt.s | sort -g -k 3 |uniq > join.${k}.t-value
    done
   awk 'OFS="\t"{print $1"\t"$2}' join.*.promoter.-1500+500.bed.meth.t-value | sed "1s/^/gene\t${i}\n/" > /top10methylation/methylation/${i}.promoter.-1500+500.t.meth
  done
 sed -i "1s/^/gene\tvalue\n/" /top10methylation/methylation/top.10.txt.s 
 cd /top10methylation/
 Rscript /methylation_for_top_bar.R $name
done

###############################################################################################################
#Tissue specific methylation analysis using SMART software.
module load bedtools
cd /CGmethylation
bedtools unionbedg -i uterus_3842.bed Spleen_3886.bed Spleen_3842.bed Rumen_3886.bed rumen_3842.bed ovary_3842.bed MAM3886_1.bed MAM3842_1.bed Lung_3886.bed Lung_3842.bed liver_3886.bed Liver_3842.bed Ldorsi_3886.bed Ldorsi_3842.bed Kidney_3842.bed Ileum_3842.bed Heart_3842.bed CORTEX3886_1.bed CORTEX3842_1.bed blood_4254.bed Blood_4091.bed Adipose_3886.bed Adipose_3842.bed -header -names G1_1 G2_1 G2_2 G3_1 G3_2 G4_1 G5_1 G5_2 G6_1 G6_2 G7_1 G7_2 G8_1 G8_2 G9_1 G10_1 G11_1 G12_1 G12_2 G13_1 G13_2 G14_1 G14_2 -filler - > methylMatrix_3.txt
source /virtualenv-master/myVE/bin/activate
SMART methylMatrix_3.txt -t DeNovoDMR -n tissue_specific_no_placentasperm -o ./tissue_specific_no_placentasperm -MR 0.5 -AG 1.0 -MS 0.5 -CN 5 -SL 20 -PD 0.05 -AD 0.3 



