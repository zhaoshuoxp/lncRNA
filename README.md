# *de novo* lncRNA discovery pipeline
-----
This is a shell/awk re-writted [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) script for lncRNAs discovery. 

> Requirements:
> awk, Python, Biopython, bedtools, [bbmap](https://sourceforge.net/projects/bbmap/), gtf_to_fasta(tophat module), [CPC2.0beta](http://cpc2.cbi.pku.edu.cn/download.php), [HMMER](http://hmmer.org)([with Pfam-A dataset](https://pfam.xfam.org))



----

### Input

This script uses original output GTF file of [stringtie](https://ccb.jhu.edu/software/stringtie/). There should be a unique transcript ID in field 2, FPKM in field 4 for *de novo* or in field 7 for reference transcript of column 9 when you run stringtie with [GENCODE ref GTF](https://www.gencodegenes.org/human/release_19.html) guide.

> The position of ID or FPKM could change if you use other references, e.g. NONCODE. Change line43, line57 if necessary.

### How it works
The filtering principle and threshold are based on [PLAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4576741/) with several steps simplification and improvement.
#### **For known lncRNAs:**
1. Seperate reference lncRNA from input by matching the unique ID of  known lncRNAs in the reference to the input.
> [GRCh37.p13 long non-coding RNA gene annotation](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.long_noncoding_RNAs.gtf.gz) was used by default, can be changed on line37. !!!Use the SAME version of ref GTF annotation ([GRCh37.p13 comprehensive](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz)) for transcript assembly.!!!

2. Filter the known lncRNA by FPKM>1.
> The value can be changde on line43.

#### **For *de novo* lncRNAs:**
1. Get all *de novo* transcripts by removing all transcripts which have been assigned with a reference ID.
2. Sort the rest of transcripts into single- or multi- exon group.
3. Filter 1 - FPKM and length: transcipt length has to be >200bp and FPKM has to be >1 for single-exon or >0.1 for multi-exon to be kept.
> Thresholds are on line80 and 86.

4. Filter 2 - repeat/gap regions: Any transcrpt which has at least 1 exon overlaps with genomic repeat sequences or gap regions for >=50% fraction will be removed.
> line99 ans 107.

5. Filter 3 - CDS: Any transcrpt which overlaps with protein coding sequences on the same strand will be removed.
6. Filter 4 - distance to genes: For the single-exon transcripts which have distance < 2000bp to protein coding gene and <500 for multi-exon transcripts will be removed.
> line138 and 147.

7. Filter 5 - CPC2 and HMMER: Non-coding potential and protein coding potential of the transcripts will be calculated by CPC2.0beta and hmmscan with Pfam-A dataset. Transcripts passed both will be kept and could be consider as *de novo* lncRNAs. 

Then the script combines known and *de novo* lncRNAs together as final.gtf.

### Example

**1.** Reads mapping and transcripts assembly. [trans_assemble.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers/blob/master/trans_assemble.sh) can be used as an example:

```shell
    wget https://github.com/zhaoshuoxp/Pipelines-Wrappers/blob/master/trans_assemble.sh
    chmod 755 trans_assemble.sh
    ./trans_assemble.sh test_R1.fastq.gz test_R2.fastq.gz rf
```

The output test.gtf can be used for lncRNA discovery.
**2.** lncRNA filtering:

```shell
    git clone https://github.com/zhaoshuoxp/lncRNA
    cd lncRNA
    chmod 755 lncRNA.sh
    ./lncRNA.sh test.gtf
```

###  Output
All results will be store in current (./) directory. Log will be printed when running.

* final.gtf: combined final transcripts in GTF format.
* known_lncRNA_f1.gtf: FPKM>1 filtered reference lncRNA transcripts.
* multi_exon_f5.gtf: all filter passed multi-exon *de novo* lncRNA transcripts.
* single_exon_f5.gtf: all filter passed single-exon *de novo* lncRNA transcripts.

Further transcript deduplication could be performed if you merge multiple GTFs before or after running this pipeline:

```shell
    wget https://github.com/zhaoshuoxp/Converters/blob/master/GTF_rmdup.sh
    chmod 755 GTF_rmdup.sh
    ./GTF_rmdup.sh final.gtf final_uniq.gtf
```

> NOTE:[UCSC Genome Browser utility](http://hgdownload.soe.ucsc.edu/admin/exe/) gtfToGenePred and genePredToBed are required.



----

Author [@zhaoshuoxp](https://github.com/zhaoshuoxp)  
Mar 27 2019  