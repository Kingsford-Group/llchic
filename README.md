
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

The whole Long Range Contact Correction module are divided into two parts. The first part is for estimating genomic distribution distance, base calling errrors and re-alignment ambiguous read pairs (read pairs with secondary alignments). The script used is LRCC_part_1.py. 

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

Note that all the files in the input path that match the form "*_R1_hg19.bwt2merged.bam"/"*_R2_hg19.bwt2merged.bam" will be processed parallelly.

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



This will generate distance matrix of HiC data and output persistence pairs and all the simplices appeared during the persistent homology. 



Here is an example:
```
python HiC_TDA.py -i example/input/RUES2_CM_combined_100000_iced_chr22.matrix -p example/output/ -o RUES2_CM_combined_100000_iced_chr22 -r 100000
```
Run this command line, and we will get three output files:

1: RUES2_CM_combined_100000_iced_chr22_distmat.txt
This file saves the distance matrix generated from original HiC contact matrix. 

2: RUES2_CM_combined_100000_iced_chr22_persisdiagram.txt
This file contains all the persistent diagram generated from persistent homology. 
- the first column: the dimension of a homology class
- the second column: birth time
- the third column: death time
- the fourth column: persitence pairs

3: RUES2_CM_combined_100000_iced_chr22_skeleton.txt
This file contains all the sinmplices generated from persistent homology. 
- the first column: a set of nodes representing one simplex
- the second column: the birth time of the simplex