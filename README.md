
# Overview

This repo contains the python souce code to run the method described in "Probabilistic method corrects previously uncharacterized Hi-C artifact" (Shen and Kingsford), which characterizes a previously uncharacterized Hi-C artifact and then design a probabilistic method to correct it. The method is combined with a standard Hi-C pipeline, HiC-Pro. The inputs are BAM files from mapping step of HiC-Pro, and the correction module (called Long Range Contact Correction), the outputs are new BAM files which can be used for Hi-C filtering step (proc_hic) of HiC-Pro. 


# Dependencies


To use these python source codes, you must have installed the following packages:

HiC-Pro: A standard pipeline for processing reads from Hi-C data. Please visit the webpage (http://nservant.github.io/HiC-Pro/) for the guide of installation and usage. 

Bowtie2: A widely used read alignment software. Please visit the webpage (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for the guide of installation and usage. 

Samtools: A suite of programs for dealing with high-throughput sequecing data. Please visit the webpage (http://www.htslib.org) for the guide of installation and usage. 

pysam: A python package for manipulating genomic data. Please visit the webpage (https://pysam.readthedocs.io/en/latest/index.html) for the guide of installation and usage. 

PyVCF: A python package for manipulating VCF files. Please visit the webpage (https://pyvcf.readthedocs.io/en/latest/index.html) for the guide of installation and usage.


# Usage

The whole Long Range Contact Correction module is divided into two parts. The first part is for estimating genomic distribution distance, base calling errrors and re-alignment ambiguous read pairs (read pairs with secondary alignments). The script used is LRCC_part_1.py. 

Here is an example: 

```
python LRCC_part_1.py -i ./example/input/ -r1 NHEK_test_R1_hg19.bwt2merged.bam -r2 NHEK_test_R2_hg19.bwt2merged.bam -o ./example/ -rs 10000 -k 2
```

The folder "./example/input" contains test BAM files. Each file has 5000 reads.

Option Tag | Description
----------------------- | -----------------------------
-i | The input path of BAM files
-r1 | BAM file from mapping step containing the first half read ends
-r2 | BAM file from mapping step containing the second half read ends
-o | The output path
-rs | The size of bins for calculating genomic distance distribution
-k | The number of alignments reported from Bowtie2 for each read

The second part is running probabilistic method to generate new BAM files. The script used for this part is LRCC_part_2.py. Note that in this script we use multiprocessing module to do parallelization. 

Here is an example:

```
python LRCC_part_2.py -i ./example/input/ -r1 *_R1_hg19.bwt2merged.bam -r2 *_R2_hg19.bwt2merged.bam -o1 ./example/ -o2 ./example/output/ -rs 10000 -m 1
```

Note that all the files in the input path that match the form '*_R1_hg19.bwt2merged.bam'/'*_R2_hg19.bwt2merged.bam' will be processed parallelly.

Option Tag | Description
----------------------- | -----------------------------
-i | The input path of BAM files, should be the same as part 1 
-r1 | BAM file from mapping step containing the first half read ends, regular expression can be used 
-r2 | BAM file from mapping step containing the second half read ends, regular expression can be used 
-o1 | The output path of part 1
-o2 | The output path of part 2
-rs | The size of bins for calculating genomic distance distribution, should be the same as part 1
-m | 0 or 1. 1 means choosing the location pair with highest posterior probability, 0 means sampling a location pair according to posterior probabilities
-gf | VCF file containing genetic variation information 

The option "-gf" can be used if genetic variation information are provided. Please visit the webpage (https://gatk.broadinstitute.org/hc/en-us) for the guide of genotype calling. 

The output folder './example/output/' contains two kinds of files. 
- NHEK_test_R1_hg19_new_multialign.bam/NHEK_test_R2_hg19_new_multialign.bam  These files are new alignment results of ambiguous pairs according to the probabilistic model (step 6 discussed in the paper). 
- NHEK_test_R1_hg19.bwt2merged.bam/NHEK_test_R2_hg19.bwt2merged.bam  These files are the combination of results in '*_multialign.bam' and the alignment of reliable read pairs (pairs with no secondary alignments). They can be used back to the filtering step of HiC-Pro. 
