
# Overview

This repo contains the python souce code to run the method described in "Probabilistic method corrects previously uncharacterized Hi-C artifact" (Shen and Kingsford), which characterizes a previously uncharacterized Hi-C artifact and then design a probabilistic method to correct it. The method is combined with a standard Hi-C pipeline, HiC-Pro. The inputs are BAM files from mapping step of HiC-Pro, and the outputs are new BAM files which can be used for Hi-C filtering step (proc_hic) of HiC-Pro. 


# Dependencies


To use these python source codes, you must have installed the following packages:

HiC-Pro: A standard pipeline for processing reads from Hi-C data. Please visit the webpage (http://nservant.github.io/HiC-Pro/) for the guide of installation and usage. 

Bowtie2: A widely used read alignment software. Please visit the webpage (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for the guide of installation and usage. 

Samtools: A suite of programs for dealing with high-throughput sequecing data. Please visit the webpage (http://www.htslib.org) for the guide of installation and usage. 

pysam: A python package for manipulating genomic data. Please visit the webpage (https://pysam.readthedocs.io/en/latest/index.html) for the guide of installation and usage. 

PyVCF: A python package for manipulating VCF files. Please visit the webpage (https://pyvcf.readthedocs.io/en/latest/index.html) for the guide of installation and usage.


# Usage
```
python HiC_TDA.py -i input_file_name -o output_file_name -p output_path -r resolution
```
This will generate distance matrix of HiC data and output persistence pairs and all the simplices appeared during the persistent homology. 

Option Tag | Description
----------------------- | -----------------------------
-i \<inputfile>| the input file, a normalized HiC contact matrix
-p \<path> | the directory containing output files 
-o | the name of output files
-r | The resolution of HiC data

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