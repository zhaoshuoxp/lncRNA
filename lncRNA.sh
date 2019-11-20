#!/bin/bash
# Requirement:
# AWK, Bedtools, bbmap, tophat(gtf_to_fasta)
# CPC2, Biopython, HMMER(Pfam-A dataset)

BED=$(cd `dirname $0`; pwd)"/bed"
genome="/home/quanyi/genome/hg19/GRCh37.p13.genome.fa"
CPC2="/home/quanyi/app/CPC2-beta/bin/CPC2.py"
translation="/home/quanyi/app/bbmap/translate6frames.sh"
Pfam="/home/quanyi/app/PLAR/Pfam/Pfam-A.hmm"
s_fpkm=1
m_fpkm=0.1
r_fpkm=0.1

# gunzip BED files
files=$(ls -1 $BED)
for i in files:
do
    if [ "${i##*.}"x = "gz"x ];then
        gunzip ${BED}/$i
    fi
done

# help message
help(){
	cat <<-EOF
    Usage: lncRNA.sh <options> -g|-n <transcriptd.gtf> 

    ### INPUT: GTF file (stringtie output) ###
    ### python3/biopython/CPC2/bbmap/HMMER/bedtools/gtf_to_fasta required ###

    Options:
        -g indicate GENCODE reference was used for transcript assembly 
        -n indicate NONCODE reference was used for transcript assembly 
        -s single-exon FPKM cutoff (defualt:1)
        -m multi-exon FPKM cutoff (defualt:0.1)
        -r referenced transcript FPKM cutoff (defualt:0.1)
        -h Print this help message
EOF
	exit 0
}

# no ARGs error
if [ $# -lt 1 ];then
    help
    exit 1
fi

while getopts "hgns:m:r:" arg
do
    case $arg in
        s) s_fpkm=$OPTARG;;
        m) m_fpkm=$OPTARG;;
        r) r_fpkm=$OPTARG;;
        g) mod='gencode';;
        n) mod='noncode';;
        h) help;;
        ?) help
            exit 1;;
    esac
done

# shift ARGs to GTF
shift $(($OPTIND - 1))

prefix=$(basename $1 .gtf)


known_f1(){ # $1=mod $2=r_fpkm $3=input.gtf $4=prefix
    # filter out known lncRNA transcripts with FPKM<0.1

    if [ $1 = 'gencode' ];then
        # V19
        #awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($6 in a) print $0}' $BED/gencode.v19.long_noncoding_RNAs.list $2 known_lncRNA.gtf
        #num_ref=$(wc -l $BED/gencode.v19.long_noncoding_RNAs.list|awk '{print $1}')
        
        # V29 hg38lift19
        awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($6 in a) print $0}' $BED/gencode.v29lift37.long_noncoding_RNAs.list transcript.gtf > known_lncRNA.gtf
        awk -F"\"" '$14>"'"$2"'"' known_lncRNA.gtf > known_lncRNA_fpkm.1.gtf
        num_ref=$(wc -l $BED/gencode.v29lift37.long_noncoding_RNAs.list|awk '{print $1}')
        echo "Load $num_ref known lncRNA transcripts from GENCODE"
        
    elif [ $1 = 'noncode' ];then
        awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($6 in a) print $0}' $BED/NONCODE.list transcript.gtf > known_lncRNA.gtf
        awk -F"\"" '$12>"'"$2"'"' known_lncRNA.gtf > known_lncRNA_fpkm.1.gtf
        num_ref=$(wc -l $BED/NONCODE.list|awk '{print $1}')
        echo "Load $num_ref known lncRNA transcripts from NONCODE"
    fi     
    
    num_known=$(wc -l known_lncRNA.gtf|awk '{print $1}')
    echo "Found $num_known known lncRNA transcripts in the input GTF"
    num_fpkm_1=$(wc -l known_lncRNA_fpkm.1.gtf|awk '{print $1}')
    low=$((num_known-num_fpkm_1))
    echo "$low out of $num_known known lncRNA transcripts have low (FPKM<$2) expression"
    
    # Get filtered known lncRNA list  
    cut -f 4 -d '"' known_lncRNA_fpkm.1.gtf > known_lncRNA_fpkm.1.list
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' known_lncRNA_fpkm.1.list $3 > ${4}_known_lncRNA_f1.gtf

    rm known_lncRNA_fpkm.1.list known_lncRNA_fpkm.1.gtf known_lncRNA.gtf
}

split_sm(){ # $1=input.gtf $2=mod $3=for NONCODE known add filters
    if [ $3 ];then
        cut -f 4 -d '"' denovo.gtf > denovo.list
        num_novo=$(wc -l denovo.list |awk '{print $1}')
        echo "total $num_novo known NONCODE transcripts"
        # get exon of all trnascripts
        awk '{if($3=="exon"){print $0}}' $1 > exon.gtf
        
    else
        if [ $2 = 'gencode' ];then
            keyw='ENST'
        elif [ $2 = 'noncode' ];then
            keyw='NONHSAT'
        fi
    
        # Remove all reference transcripts
        grep -v $keyw transcript.gtf > denovo.gtf
        # get de novo transcripts ID
        cut -f 4 -d '"' denovo.gtf > denovo.list
        num_novo=$(wc -l denovo.list |awk '{print $1}')
        echo "total $num_novo unannotated transcripts"
    
        # get exon of all trnascripts
        awk '{if($3=="exon"){print $0}}' $1 |grep -v $keyw > exon.gtf
    fi
    
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
   
    rm single_exon.list multi_exon.list denovo.gtf exon.list denovo.list
}

denovo_f1(){ # $1=s_fpkm $2=m_fpkm $3=single_length $4=multi_length $5=fpkm_pos
    # Remove FPKM<1 and length <200 bp single exon
    awk -F"\"" '$"'"$5"'">"'"$1"'"' single_exon.gtf |awk '{if(($5-$4)>="'"$3"'"){print $0}}'> single_exon_f1.gtf
    num_single_f1=$(wc -l single_exon_f1.gtf |awk '{print $1}')
    cut -f 4 -d '"' single_exon_f1.gtf >single_exon_f1.list
    echo "$num_single_f1 out of $num_single single-exon transcripts have FPKM>$1 and length>$3nt"

    # Remove FPKM<0.1 and length <200 bp multi exon
    awk -F"\"" '$"'"$5"'">"'"$2"'"' multi_exon.gtf |awk '{if(($5-$4)>="'"$4"'"){print $0}}'> multi_exon_f0.gtf
    cut -f 4 -d '"' multi_exon_f0.gtf > multi_exon_f0.list
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f0.list exon.gtf > multi_exon_f0_exon.gtf
    awk '{a[$12]+=($5-$4)}END{for(i in a){if(a[i]>=200){print i}}}' multi_exon_f0_exon.gtf |cut -f2 -d'"' > multi_exon_f1.list
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}'  multi_exon_f1.list multi_exon_f0.gtf > multi_exon_f1.gtf
    num_multi_f1=$(wc -l multi_exon_f1.gtf |awk '{print $1}')
    echo "$num_multi_f1 out of $num_multi multi-exon transcripts have FPKM>$2 and length>$4nt"
    
    rm single_exon.gtf multi_exon.gtf multi_exon_f0.gtf multi_exon_f0.list multi_exon_f0_exon.gtf 
}

denovo_f2(){ # $1=overlap fraction
    # Remove single-exon transcripts overlapping with REPEAT sequences
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' single_exon_f1.list exon.gtf > single_exon_f1_exon.gtf
    cut -f 1,4,5,7,9 single_exon_f1_exon.gtf |awk -v OFS="\t" '{print $1,$2,$3,$7,$8,$4}'|grep ^chr > single_exon_f1_exon.bed
    intersectBed -a single_exon_f1_exon.bed -b $BED/hg19.repeat_gap.bed -v -s -f $1 > single_exon_f2_exon.bed
    num_single_f2=$(wc -l single_exon_f2_exon.bed |awk '{print $1}')
    single_ol_repeat=$((num_single_f1-num_single_f2))
    echo "$single_ol_repeat out of $num_single_f1 single-exon transcripts have overlaps with repeat sequences"

    # Remove multi-exon transcripts overlapping with REPEAT sequences
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f1.list exon.gtf > multi_exon_f1_exon.gtf
    cut -f 1,4,5,7,9 multi_exon_f1_exon.gtf |awk -v OFS="\t" '{print $1,$2,$3,$7,$8,$4}'|grep ^chr > multi_exon_f1_exon.bed
    echo '#' >multi_exon_f1_ol_repeat
    intersectBed -a multi_exon_f1_exon.bed -b $BED/hg19.repeat_gap.bed -u -s -f $1 |cut -f 2 -d '"' |sort -u >> multi_exon_f1_ol_repeat
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a){}else{print $0}}' multi_exon_f1_ol_repeat multi_exon_f1.gtf > multi_exon_f2.gtf
    num_multi_f2=$(wc -l multi_exon_f2.gtf|awk '{print $1}')
    multi_ol_repeat=$(wc -l multi_exon_f1_ol_repeat|awk '{print $1}')
    echo "$multi_ol_repeat out of $num_multi_f1 mutli-exon transcripts have overlaps with repeat sequences"

    rm single_exon_f1.list single_exon_f1_exon.gtf single_exon_f1_exon.bed multi_exon_f1.list multi_exon_f1_exon.gtf multi_exon_f1_exon.bed multi_exon_f1_ol_repeat single_exon_f1.gtf multi_exon_f1.gtf
}

denovo_f3(){ # $1=overlapping fraction
    # Remove single-exon transcripts overlapping with CDS
    intersectBed -a single_exon_f2_exon.bed -b $BED/CDS.bed -v -s -f $1> single_exon_f3_exon.bed
    num_single_f3=$(wc -l single_exon_f3_exon.bed |awk '{print $1}')
    single_ol_cds=$((num_single_f2-num_single_f3))
    echo "$single_ol_cds out of $num_single_f2 single-exon transcripts have overlaps with CDS"

    # Remove multi-exon transcripts overlapping with CDS
    cut -f 4 -d '"' multi_exon_f2.gtf > multi_exon_f2.list
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f2.list exon.gtf > multi_exon_f2_exon.gtf
    cut -f 1,4,5,7,9 multi_exon_f2_exon.gtf |awk -v OFS="\t" '{print $1,$2,$3,$7,$8,$4}'|grep ^chr > multi_exon_f2_exon.bed
    echo '#' >multi_exon_f2_ol_cds
    intersectBed -a multi_exon_f2_exon.bed -b $BED/CDS.bed -u -s -f $1|cut -f 2 -d '"' |sort -u >> multi_exon_f2_ol_cds
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a){}else{print $0}}' multi_exon_f2_ol_cds multi_exon_f2.gtf > multi_exon_f3.gtf
    multi_ol_cds=$(wc -l multi_exon_f2_ol_cds|awk '{print $1}')
    num_multi_f3=$(wc -l multi_exon_f3.gtf|awk '{print $1}')
    echo "$multi_ol_cds out of $num_multi_f2 mutli-exon transcripts have overlaps with CDS"
    
    rm single_exon_f2_exon.bed multi_exon_f2.gtf multi_exon_f2.list multi_exon_f2_exon.gtf multi_exon_f2_exon.bed multi_exon_f2_ol_cds
}

denovo_f4(){ # $1=single_dist $2=multi_dist
    # Remove single-exon transcripts near (<2000bp) protein-coding genes
    bedtools sort -i single_exon_f3_exon.bed > single_exon_f3_exon.bed2
    closestBed -s -d -a $BED/hg19_pc_txSE.bed -b single_exon_f3_exon.bed2 > single_exon_f3_exon.dist
    awk '$13<"'"$1"'"' single_exon_f3_exon.dist |awk -F"\"" '{print $2}'| sort -u > singe_dist_re.list
    awk -F"\"" 'NR==FNR{a[$0]}NR>FNR{if ($2 in a){}else{print $0}}' singe_dist_re.list single_exon_f3_exon.bed > single_exon_f4_exon.bed
    num_single_f4=$(wc -l single_exon_f4_exon.bed|awk '{print $1}')
    single_ol_ud=$((num_single_f3-num_single_f4))
    echo "$single_ol_ud out of $num_single_f3 single-exon transcripts are too close ($1bp) to protein-coding genes"

    # Remove multi-exon transcripts near (<500bp) protein-coding genes
    awk -v OFS="\t" '{print $1,$4,$5,$11,$12,$7}' multi_exon_f3.gtf |grep ^chr > multi_exon_f3.bed
    bedtools sort -i multi_exon_f3.bed > multi_exon_f3.bed2
    closestBed -s -d -a $BED/hg19_pc_txSE.bed -b multi_exon_f3.bed2 > multi_exon_f3.dist
    awk '$13<"'"$2"'"' multi_exon_f3.dist |awk -F"\"" '{print $2}'| sort -u > multi_dist_re.list
    awk -F"\"" 'NR==FNR{a[$0]}NR>FNR{if ($2 in a){}else{print $0}}' multi_dist_re.list multi_exon_f3.bed > multi_exon_f4.bed
    num_multi_f4=$(wc -l multi_exon_f4.bed|awk '{print $1}')
    multi_ol_ud=$((num_multi_f3-num_multi_f4))
    echo "$multi_ol_ud out of $num_multi_f3 multi-exon transcripts are too close ($2bp) to protein-coding genes"

    rm single_exon_f3_exon.bed single_exon_f3_exon.bed2 single_exon_f3_exon.dist multi_exon_f3.gtf multi_exon_f3.bed multi_exon_f3.bed2 multi_exon_f3.dist multi_dist_re.list singe_dist_re.list
}

cpc2_f(){ # $1=input.gtf
    # Remove CPC2 coding poteintal transcripts for single-exon
    cut -f 2 -d '"' single_exon_f4_exon.bed > single_exon_f4.list
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' single_exon_f4.list $1 > single_exon_f4.gtf
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

    # Remove CPC2 coding poteintal transcripts for multi-exon
    cut -f 2 -d '"' multi_exon_f4.bed > multi_exon_f4.list
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f4.list $1 > multi_exon_f4.gtf
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

    rm single_exon_f4_exon.bed single_exon_f4.list single_exon_f4.fa single_exon_f4.fa.list single_exon_f4.out single_exon_f4.nc single_exon_f4_cpc.gtf
    rm multi_exon_f4.bed multi_exon_f4.list multi_exon_f4.fa multi_exon_f4.fa.list multi_exon_f4.out multi_exon_f4.nc multi_exon_f4_cpc.gtf
}

hmmer_f(){ # $1=E-value $2=$prefix
    # Remove HMMER coding poteintal transcripts for single-exon
    $translation in=single_exon_f4_cpc.fa out=single_exon_f4_cpc.aa frames=3
    echo "Running HMMER on single-exon transcripts:"
    hmmscan --cpu 16 -o hmmer.log --noali -E $1 --tblout single_exon_f4_hmm.txt $Pfam single_exon_f4_cpc.aa
    grep ^\> single_exon_f4_cpc.fa|sed 's/^>*//g'|cut -f1,2 -d ' ' > single_exon_f4_cpc.fa.list
    awk '{print $3}' single_exon_f4_hmm.txt | grep [0-9] |sort -n -u > single_exon_f4_hmm.list
    echo '#' >single_exon_f4.nc.hmm.list
    awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $2}' single_exon_f4_hmm.list single_exon_f4_cpc.fa.list >> single_exon_f4.nc.hmm.list
    single_hmm_nc=$(wc -l single_exon_f4.nc.hmm.list|awk '{print $1-1}')
    echo "$single_hmm_nc are coding by HMMER"
    awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a){}else{print $0}}' single_exon_f4.nc.hmm.list single_exon_f4.nc.list > single_exon_f5.list
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' single_exon_f5.list single_exon_f4.gtf > ${2}_single_exon_f5.gtf

    # Remove HMMER coding poteintal transcripts for multi-exon
    $translation in=multi_exon_f4_cpc.fa out=multi_exon_f4_cpc.aa frames=3
    echo "Running HMMER on multi-exon transcripts:"
    hmmscan --cpu 16 -o hmmer.log --noali -E $1 --tblout multi_exon_f4_hmm.txt $Pfam multi_exon_f4_cpc.aa
    grep ^\> multi_exon_f4_cpc.fa|sed 's/^>*//g'|cut -f1,2 -d ' ' > multi_exon_f4_cpc.fa.list
    awk '{print $3}' multi_exon_f4_hmm.txt | grep [0-9] |sort -n -u > multi_exon_f4_hmm.list
    echo '#' >multi_exon_f4.nc.hmm.list
    awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $2}' multi_exon_f4_hmm.list multi_exon_f4_cpc.fa.list >> multi_exon_f4.nc.hmm.list
    multi_hmm_nc=$(wc -l multi_exon_f4.nc.hmm.list|awk '{print $1-1}')
    echo "$multi_hmm_nc are coding by HMMER"
    awk 'NR==FNR{a[$0]}NR>FNR{if ($1 in a){}else{print $0}}' multi_exon_f4.nc.hmm.list multi_exon_f4.nc.list > multi_exon_f5.list
    awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' multi_exon_f5.list multi_exon_f4.gtf > ${2}_multi_exon_f5.gtf

    rm single_exon_f4_cpc.fa single_exon_f4_cpc.aa hmmer.log single_exon_f4_hmm.txt single_exon_f4_cpc.fa.list single_exon_f4_hmm.list single_exon_f4.nc.hmm.list single_exon_f5.list single_exon_f4.gtf
    rm multi_exon_f4_cpc.fa multi_exon_f4_cpc.aa multi_exon_f4_hmm.txt multi_exon_f4_cpc.fa.list multi_exon_f4_hmm.list multi_exon_f4.nc.hmm.list multi_exon_f5.list multi_exon_f4.gtf
    rm single_exon_f4.nc.list single_exon_f4.fa.tlst single_exon_f4_cpc.fa.tlst multi_exon_f4.nc.list multi_exon_f4_cpc.fa.tlst multi_exon_f4.fa.tlst
}


main(){ # $input.gtf
    # Grep all trnscripts
    awk '$3=="transcript"' $1 > transcript.gtf
    num_trans=$(wc -l transcript.gtf|awk '{print $1}')
    echo "Total $num_trans transcripts in the input GTF"
    
    # # Filter out known lncRNA transcripts with FPKM<0.1
    if [ $mod ];then
        known_f1 $mod $r_fpkm $1 $prefix
        # Get de novo transcript
        split_sm $1 $mod
    else
        echo "Reference has to be either GENCODE (-g) or NONCODE (-n)!"
        exit 1
    fi
    
    # De novo lncRNA filtering
    # filter 1: FPKM and length
    denovo_f1 $s_fpkm $m_fpkm 200 200 8
    # filter 2: repeat/gap overlaps
    denovo_f2 0.5 # overlapping fraction
    # filter 3: CDS overlaps
    denovo_f3 1E-9 # overlapping fraction
    # filter 4: close to genes
    denovo_f4 2000 500 # max distance to genes
    # filter CPC2: coding potential of NT sequences
    cpc2_f $1
    # filter hmmscan: peptides alignment to Pfam-A
    hmmer_f 0.001 $prefix # E-value cutoff
    
    # summary
    num_single_f5=$(awk '$3=="transcript"' ${prefix}_single_exon_f5.gtf |wc -l|awk '{print $1}')
    num_multi_f5=$(awk '$3=="transcript"' ${prefix}_multi_exon_f5.gtf |wc -l|awk '{print $1}')
    
    if [ $mod = 'noncode' ];then
        echo ""
        echo "Addtional filters for known transcripts in NOCODE database need to be applied:"
        if [ ! -d non_add ];then
            mkdir non_add
        fi
        cd non_add
        awk '$3=="transcript"' ../${prefix}_known_lncRNA_f1.gtf > denovo.gtf
        
        split_sm ../${prefix}_known_lncRNA_f1.gtf $mod 'non_2nd'
        denovo_f1 $s_fpkm $m_fpkm 200 200 12
        denovo_f2 0.5 
        denovo_f3 1E-9 
        denovo_f4 2000 500 
        cpc2_f ../${prefix}_known_lncRNA_f1.gtf
        hmmer_f 0.001 $prefix 
        
        cat ${prefix}_single_exon_f5.gtf ${prefix}_multi_exon_f5.gtf > ${prefix}_final_known.gtf
        num_know=$(awk '$3=="transcript"' ${prefix}_final_known.gtf |wc -l |awk '{print $1}')
        echo "$num_know known NONCODE transcripts were kept after applying addtional filters"
        echo ""
        rm exon.gtf
        cd ..
        
        total=$((num_single_f5+num_know+num_multi_f5))
        echo "Total $total final lncRNA, out of these:"
        echo "$num_single_f5 are final single-exon lncRNA"
        echo "$num_multi_f5 are final multi-exon lncRNA" 
        echo "$num_know are known lncRNA"
        
        cat ./non_add/${prefix}_final_known.gtf ${prefix}_single_exon_f5.gtf ${prefix}_multi_exon_f5.gtf > ${prefix}_final.gtf

    else
        total=$((num_single_f5+num_fpkm_1+num_multi_f5))
        echo "Total $total final lncRNA, out of these:"
        echo "$num_single_f5 are final single-exon lncRNA"
        echo "$num_multi_f5 are final multi-exon lncRNA" 
        echo "$num_fpkm_1 are known lncRNA"
        cat ${prefix}_known_lncRNA_f1.gtf ${prefix}_single_exon_f5.gtf ${prefix}_multi_exon_f5.gtf > ${prefix}_final.gtf
    fi
    
    rm exon.gtf  transcript.gtf
}

main $1

# check running status
if [ $? -ne 0 ]; then
    help
    exit 1
else
    echo ""
    echo "Run succeed"
fi

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################