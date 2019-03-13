#!/bin/bash
# Requirement:
# AWK, Bedtools, bbmap, tophat(gtf_to_fasta)
# CPC2, Biopython, HMMER(Pfam-A dataset)

file=$1
BED="`dirname $0`/bed"
genome="/home/quanyi/genome/hg19/GRCh37.p13.genome.fa"
CPC2="/home/quanyi/app/CPC2-beta/bin/CPC2.py"
translation="/home/quanyi/app/bbmap/translate6frames.sh"
Pfam="/home/quanyi/app/PLAR/Pfam/Pfam-A.hmm"

# gunzip BED files
files=$(ls -1 $BED)
for i in files:
do
	if [ "${i##*.}"x = "gz"x ];then
		gunzip ${BED}/$i
	fi
done

# Grep all trnscripts
awk '$3=="transcript"' $file > transcript.gtf
num_trans=$(wc -l transcript.gtf|awk '{print $1}')
echo "Total $num_trans transcripts in the input GTF"

# get known lncRNA transcripts ID
#num_ref=$(wc -l $BED/NONCODE.list|awk '{print $1}')
#num_ref=$(wc -l $BED/gencode.v29lift37.long_noncoding_RNAs.list|awk '{print $1}')
num_ref=$(wc -l $BED/gencode.v19.long_noncoding_RNAs.list|awk '{print $1}')
echo "Load $num_ref known lncRNA transcripts from GENCODE"

####
# Grep known lncRNA transcripts from input GTF
####

awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($6 in a) print $0}' $BED/gencode.v19.long_noncoding_RNAs.list transcript.gtf known_lncRNA.gtf
#awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($6 in a) print $0}' $BED/gencode.v29lift37.long_noncoding_RNAs.list transcript.gtf known_lncRNA.gtf
#awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($6 in a) print $0}' $BED/NONCODE.list transcript.gtf > known_lncRNA.gtf
num_known=$(wc -l known_lncRNA.gtf|awk '{print $1}')

# filter out known lncRNA transcripts with FPKM<0.1
awk -F"\"" '$14>0.1' known_lncRNA.gtf > known_lncRNA_fpkm.1.gtf
num_fpkm_1=$(wc -l known_lncRNA_fpkm.1.gtf|awk '{print $1}')
low=$((num_known-num_fpkm_1))
echo "Found $num_known known lncRNA transcripts in the input GTF"
echo "$low out of $num_known known lncRNA transcripts have low (FPKM<0.1) expression"
# Get filtered known lncRNA list  
cut -f 4 -d '"' known_lncRNA_fpkm.1.gtf > known_lncRNA_fpkm.1.list
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' known_lncRNA_fpkm.1.list $file > known_lncRNA_f1.gtf

###
# De novo lncRNA filtering
###

# Remove all reference transcripts
grep -v ENST transcript.gtf > denovo.gtf
# get de novo transcripts ID
cut -f 4 -d '"' denovo.gtf > denovo.list
num_novo=$(wc -l denovo.list |awk '{print $1}')
echo "total $num_novo unannotated transcripts"
# get exon of all trnascripts
awk '{if($3=="exon"){print $0}}' $file |grep -v ENST > exon.gtf

# Sort transcripts by exon number
cut -f4 -d '"' exon.gtf |sort |uniq -c > exon.list
awk '{if($1>=2)print $2}' exon.list > multi_exon.list
awk '{if($1<2)print $2}' exon.list > single_exon.list
num_single=$(wc -l single_exon.list |awk '{print $1}')
num_multi=$(wc -l multi_exon.list |awk '{print $1}')
echo "Of these, $num_single are single-exon and $num_multi are multi-exon transcripts"

# Get single/multi-exon transcripts
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' single_exon.list denovo.gtf > single_exon.gtf
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon.list denovo.gtf > multi_exon.gtf

# filter 1: FPKM and length

# Remove FPKM<1 and length <200 bp single exon
awk -F"\"" '$8>1' single_exon.gtf |awk '{if(($5-$4)>=200){print $0}}'> single_exon_f1.gtf
num_single_f1=$(wc -l single_exon_f1.gtf |awk '{print $1}')
cut -f 4 -d '"' single_exon_f1.gtf >single_exon_f1.list
echo "$num_single_f1 out of $num_single single-exon transcripts have FPKM>1 and length>200nt"

# Remove FPKM<0.1 and length <200 bp multi exon
awk -F"\"" '$8>0.1' multi_exon.gtf |awk '{if(($5-$4)>=200){print $0}}'> multi_exon_f0.gtf
cut -f 4 -d '"' multi_exon_f0.gtf > multi_exon_f0.list
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f0.list exon.gtf > multi_exon_f0_exon.gtf
awk '{a[$12]+=($5-$4)}END{for(i in a){if(a[i]>=200){print i}}}' multi_exon_f0_exon.gtf |cut -f2 -d'"' > multi_exon_f1.list
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}'  multi_exon_f1.list multi_exon_f0.gtf > multi_exon_f1.gtf
num_multi_f1=$(wc -l multi_exon_f1.gtf |awk '{print $1}')
echo "$num_multi_f1 out of $num_multi multi-exon transcripts have FPKM>0.1 and length>200nt"

# filter 2: repeat/gap overlaps

# Remove single-exon transcripts overlapping with REPEAT sequences
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' single_exon_f1.list exon.gtf > single_exon_f1_exon.gtf
cut -f 1,4,5,7,9 single_exon_f1_exon.gtf |awk -v OFS="\t" '{print $1,$2,$3,$7,$8,$4}'|grep ^chr > single_exon_f1_exon.bed
intersectBed -a single_exon_f1_exon.bed -b $BED/hg19.repeat_gap.bed -v -s -f 0.5 > single_exon_f2_exon.bed
num_single_f2=$(wc -l single_exon_f2_exon.bed |awk '{print $1}')
single_ol_repeat=$((num_single_f1-num_single_f2))
echo "$single_ol_repeat out of $num_single_f1 single-exon transcripts have overlaps with repeat sequences"

# Remove multi-exon transcripts overlapping with REPEAT sequences
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f1.list exon.gtf > multi_exon_f1_exon.gtf
cut -f 1,4,5,7,9 multi_exon_f1_exon.gtf |awk -v OFS="\t" '{print $1,$2,$3,$7,$8,$4}'|grep ^chr > multi_exon_f1_exon.bed
intersectBed -a multi_exon_f1_exon.bed -b $BED/hg19.repeat_gap.bed -u -s -f 0.5 |cut -f 2 -d '"' |sort -u > multi_exon_f1_ol_repeat
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a){}else{print $0}}' multi_exon_f1_ol_repeat multi_exon_f1.gtf > multi_exon_f2.gtf
num_multi_f2=$(wc -l multi_exon_f2.gtf|awk '{print $1}')
multi_ol_repeat=$(wc -l multi_exon_f1_ol_repeat|awk '{print $1}')
echo "$multi_ol_repeat out of $num_multi_f1 mutli-exon transcripts have overlaps with repeat sequences"

# filter 3: CDS overlaps

# Remove single-exon transcripts overlapping with CDS
intersectBed -a single_exon_f2_exon.bed -b $BED/CDS.bed -v -s > single_exon_f3_exon.bed
num_single_f3=$(wc -l single_exon_f3_exon.bed |awk '{print $1}')
single_ol_cds=$((num_single_f2-num_single_f3))
echo "$single_ol_cds out of $num_single_f2 single-exon transcripts have overlaps with CDS"

# Remove multi-exon transcripts overlapping with CDS
cut -f 4 -d '"' multi_exon_f2.gtf > multi_exon_f2.list
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f2.list exon.gtf > multi_exon_f2_exon.gtf
cut -f 1,4,5,7,9 multi_exon_f2_exon.gtf |awk -v OFS="\t" '{print $1,$2,$3,$7,$8,$4}'|grep ^chr > multi_exon_f2_exon.bed
intersectBed -a multi_exon_f2_exon.bed -b $BED/CDS.bed -u -s |cut -f 2 -d '"' |sort -u > multi_exon_f2_ol_cds
#intersectBed -a multi_exon_f2_exon.bed -b $BED/CDS.bed -u -s -f 0.5|cut -f2 -d '"'|sort |uniq -c|awk -v OFS="\t" '{print $2,$1}' >  multi_exon_f2_ol_cds.list
#cut -f 2 -d '"' multi_exon_f2_exon.bed|sort|uniq -c|awk -v OFS="\t" '{print $2,$1}' >multi_exon_f2_exon.list
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a){}else{print $0}}' multi_exon_f2_ol_cds multi_exon_f2.gtf > multi_exon_f3.gtf
multi_ol_cds=$(wc -l multi_exon_f2_ol_cds|awk '{print $1}')
num_multi_f3=$(wc -l multi_exon_f3.gtf|awk '{print $1}')
echo "$multi_ol_cds out of $num_multi_f2 mutli-exon transcripts have overlaps with CDS"

# filter 4: close to genes

# Remove single-exon transcripts near (<2000bp) protein-coding genes
bedtools sort -i single_exon_f3_exon.bed > single_exon_f3_exon.bed2
closestBed -s -d -io -a $BED/hg19_pc_txSE.bed -b single_exon_f3_exon.bed2 > single_exon_f3_exon.dist
awk '$13>2000' single_exon_f3_exon.dist |cut -f 7-12 |sort -u > single_exon_f4_exon.bed
num_single_f4=$(wc -l single_exon_f4_exon.bed|awk '{print $1}')
single_ol_ud=$((num_single_f3-num_single_f4))
echo "$single_ol_ud out of $num_single_f3 single-exon transcripts are too close (2000bp) to protein-coding genes"

# Remove multi-exon transcripts near (<500bp) protein-coding genes
awk -v OFS="\t" '{print $1,$4,$5,$11,$12,$7}' multi_exon_f3.gtf |grep ^chr > multi_exon_f3.bed
bedtools sort -i multi_exon_f3.bed > multi_exon_f3.bed2
closestBed -s -d -io -a $BED/hg19_pc_txSE.bed -b multi_exon_f3.bed2 > multi_exon_f3.dist
awk '$13>500' multi_exon_f3.dist |cut -f 7-12 |sort -u > multi_exon_f4.bed
num_multi_f4=$(wc -l multi_exon_f4.bed|awk '{print $1}')
multi_ol_ud=$((num_multi_f3-num_multi_f4))
echo "$multi_ol_ud out of $num_multi_f3 multi-exon transcripts are too close (500bp) to protein-coding genes"

# filter 5: CPC2 and HMMER protein-coding potential tests

# Remove CPC2 coding poteintal transcripts for multi-exon
cut -f 2 -d '"' multi_exon_f4.bed > multi_exon_f4.list
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f4.list $file > multi_exon_f4.gtf
gtf_to_fasta multi_exon_f4.gtf $genome multi_exon_f4.fa
grep ^\> multi_exon_f4.fa|sed 's/^>*//g'|cut -f1,2 -d ' ' > multi_exon_f4.fa.list
echo "Running CPC2 on multi-exon transcripts:"
$CPC2 -i multi_exon_f4.fa -o multi_exon_f4.out
awk '{if ($8=="noncoding"){print $1}}' multi_exon_f4.out > multi_exon_f4.nc
awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $2}' multi_exon_f4.nc multi_exon_f4.fa.list > multi_exon_f4.nc.list
multi_cpc_nc=$(wc -l multi_exon_f4.nc.list|awk '{print $1}')
echo "$multi_cpc_nc are non-coding by CPC2"
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f4.nc.list multi_exon_f4.gtf > multi_exon_f4_cpc.gtf
gtf_to_fasta multi_exon_f4_cpc.gtf $genome multi_exon_f4_cpc.fa

# Remove HMMER coding poteintal transcripts for multi-exon
$translation in=multi_exon_f4_cpc.fa out=multi_exon_f4_cpc.aa frames=3
echo "Running HMMER on multi-exon transcripts:"
hmmscan --cpu 16 -o hmmer.log --noali -E 0.001 --tblout multi_exon_f4_hmm.txt $Pfam multi_exon_f4_cpc.aa
grep ^\> multi_exon_f4_cpc.fa|sed 's/^>*//g'|cut -f1,2 -d ' ' > multi_exon_f4_cpc.fa.list
awk '{print $3}' multi_exon_f4_hmm.txt | grep [0-9] |sort -n -u > multi_exon_f4_hmm.list
awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $2}' multi_exon_f4_hmm.list multi_exon_f4_cpc.fa.list > multi_exon_f4.nc.hmm.list
multi_hmm_nc=$(wc -l multi_exon_f4.nc.hmm.list|awk '{print $1}')
echo "$multi_hmm_nc are coding by HMMER"
awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a){}else{print $0}}' multi_exon_f4.nc.hmm.list multi_exon_f4.nc.list > multi_exon_f5.list
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f5.list multi_exon_f4.gtf > multi_exon_f5.gtf

# Remove CPC2 coding poteintal transcripts for single-exon
cut -f 2 -d '"' single_exon_f4_exon.bed > single_exon_f4.list
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' single_exon_f4.list $file > single_exon_f4.gtf
gtf_to_fasta single_exon_f4.gtf $genome single_exon_f4.fa
grep ^\> single_exon_f4.fa|sed 's/^>*//g'|cut -f1,2 -d ' ' > single_exon_f4.fa.list
echo "Running CPC2 on single-exon transcripts:"
$CPC2 -i single_exon_f4.fa -o single_exon_f4.out
awk '{if ($8=="noncoding"){print $1}}' single_exon_f4.out > single_exon_f4.nc
awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $2}' single_exon_f4.nc single_exon_f4.fa.list > single_exon_f4.nc.list
single_cpc_nc=$(wc -l single_exon_f4.nc.list|awk '{print $1}')
echo "$single_cpc_nc are non-coding by CPC2"
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' single_exon_f4.nc.list single_exon_f4.gtf > single_exon_f4_cpc.gtf
gtf_to_fasta single_exon_f4_cpc.gtf $genome single_exon_f4_cpc.fa

# Remove HMMER coding poteintal transcripts for single-exon
$translation in=single_exon_f4_cpc.fa out=single_exon_f4_cpc.aa frames=3
echo "Running HMMER on single-exon transcripts:"
hmmscan --cpu 16 -o hmmer.log --noali -E 0.001 --tblout single_exon_f4_hmm.txt $Pfam single_exon_f4_cpc.aa
grep ^\> single_exon_f4_cpc.fa|sed 's/^>*//g'|cut -f1,2 -d ' ' > single_exon_f4_cpc.fa.list
awk '{print $3}' single_exon_f4_hmm.txt | grep [0-9] |sort -n -u > single_exon_f4_hmm.list
awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $2}' single_exon_f4_hmm.list single_exon_f4_cpc.fa.list > single_exon_f4.nc.hmm.list
single_hmm_nc=$(wc -l single_exon_f4.nc.hmm.list|awk '{print $1}')
echo "$single_hmm_nc are coding by HMMER"
awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a){}else{print $0}}' single_exon_f4.nc.hmm.list single_exon_f4.nc.list > single_exon_f5.list
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' single_exon_f5.list single_exon_f4.gtf > single_exon_f5.gtf

# summary
num_single_f5=$(awk '$3=="transcript"' single_exon_f5.gtf |wc -l|awk '{print $1}')
num_multi_f5=$(awk '$3=="transcript"' multi_exon_f5.gtf |wc -l|awk '{print $1}')
total=$((num_single_f5+num_fpkm_1+num_multi_f5))
echo "Total $total final lncRNA, out of these:"
echo "$num_single_f5 are final single-exon lncRNA"
echo "$num_multi_f5 are final multi-exon lncRNA" 
echo "$num_fpkm_1 are known lncRNA"
cat known_lncRNA_f1.gtf single_exon_f5.gtf multi_exon_f5.gtf > final.gtf

# clean
rm known_lncRNA_fpkm* known_lncRNA.gtf hmmer.log single_exon.list single_exon.gtf transcript.gtf multi_exon.gtf exon.gtf denovo.gtf multi_exon.list denovo.list exon.list
rm multi_exon_f0* multi_exon_f1* multi_exon_f2* multi_exon_f3* multi_exon_f4* single_exon_f1* single_exon_f2* single_exon_f3* single_exon_f4*

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################