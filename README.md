# *de novo* lncRNA discovery pipeline
-----
This is a shell/awk re-written [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) script for lncRNA discovery. 

> Requirements:
> awk, Python, Biopython, bedtools, [bbmap](https://sourceforge.net/projects/bbmap/), gtf_to_fasta(tophat module), [CPC2.0beta](http://cpc2.cbi.pku.edu.cn/download.php), [HMMER](http://hmmer.org)([with Pfam-A dataset](https://pfam.xfam.org))

[![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu) [![LICENSE](https://img.shields.io/badge/license-Anti%20996-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)

----

### Input

This script uses original output GTF file of [stringtie](https://ccb.jhu.edu/software/stringtie/). There should be unique transcript IDs on field 2, FPKMs on field 4 for *de novo* or on field 7 for reference transcript of column if stringtie was ran with [GENCODE ref GTF](https://www.gencodegenes.org/human/release_19.html) guide.

### How it works

The filtering principle and threshold are based on [PLAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4576741/) with several steps simplification and improvement.

#### **For known lncRNAs:**

1. Separate reference lncRNA from input by matching the unique ID of  known lncRNAs of the reference.
> [GENCODE  release29 (GRCh37)](https://www.gencodegenes.org/human/release_29lift37.html) was used by default, can be changed on line 73-79.
>
> !NOTE: Use the SAME version of ref GTF annotation for transcript assembly.

2. Filter the known lncRNA by FPKM (default>1).

#### **For *de novo* lncRNAs:**

1. Get all *de novo* transcripts by removing the transcripts which have been assigned with reference ID.
2. Sort the rest of transcripts into single- or multi- exon groups.
3. Filter 1 - FPKM and length: transcript length has to be >200bp and FPKM has to be >1 for single-exon or >0.1 for multi-exon to be kept.
> Thresholds are on line 311.

4. Filter 2 - repeat/gap regions: Any transcripts which have at least 1 exon overlaps with genomic repeat sequences or gap regions for >=50% fraction will be removed.
> line 313.

5. Filter 3 - CDS: Any transcripts which overlap with protein coding sequences on the same strand will be removed.

>line 315.

6. Filter 4 - distance to genes: For the single-exon transcripts which have distance < 2000bp to protein coding gene and <500 for multi-exon transcripts will be removed.
> line 317.

7. Filter 5 - CPC2 and HMMER: Non-coding potential and protein coding potential of the transcripts will be calculated by *CPC2.0beta* and *hmmscan* with Pfam-A dataset. Transcripts passed both will be kept and could be consider as potential *de novo* lncRNAs. 

> line 321.

Then the script combines known and *de novo* lncRNAs together to a final.gtf.

### Options

help message can be shown by `lncRNA.sh -h`

```shell
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
```

### Example

**1.** Reads mapping and transcripts assembly. [trans_assemble.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#trans_assemblesh) can be used. i.e.

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/trans_assemble.sh
chmod 755 trans_assemble.sh
./trans_assemble.sh test_R1.fastq.gz test_R2.fastq.gz rf
```

The output test.gtf can be used for lncRNA discovery.

**2.** lncRNA filtering:

```shell
git clone https://github.com/zhaoshuoxp/lncRNA
cd lncRNA
chmod 755 lncRNA.sh
./lncRNA.sh -g test.gtf
```

###  Output

All results will be store in current (./) directory. Log will be printed during running.

* test_final.gtf: combined final transcripts in GTF format.
* test_known_lncRNA_f1.gtf: FPKM>1 filtered reference lncRNA transcripts.
* test_multi_exon_f5.gtf: all filter passed multi-exon *de novo* lncRNA transcripts.
* test_single_exon_f5.gtf: all filter passed single-exon *de novo* lncRNA transcripts.
* non_add: for NONCODE reference transcripts, the pipeline will put them into *de novo* filters and output the filtered transcripts in this subdirectory.

Further transcript deduplication could be performed if you merged multiple GTFs before or after running this pipeline:

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/GTF_rmdup.sh
chmod 755 GTF_rmdup.sh
./GTF_rmdup.sh final.gtf final_uniq.gtf
```

> NOTE:[UCSC Genome Browser utility](http://hgdownload.soe.ucsc.edu/admin/exe/) gtfToGenePred and genePredToBed are required.



----

Author [@zhaoshuoxp](https://github.com/zhaoshuoxp)  

Nov 8 2019  