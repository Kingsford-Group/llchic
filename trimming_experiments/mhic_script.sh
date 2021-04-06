#!/bin/bash

## mHi-C
## Author: Ye Zheng
## Contact: yezheng@stat.wisc.edu

#######################################################
## Script to call each step of mHi-C
## Take ring stage of Plasmodium falciparum for example
## Update May 2018
#######################################################

projectPath="/mnt/disk69/user/yihangs/HiC_correction_revision/PLOSComp_revise_1/test_mHiC/mHiC-master" ## path to save all the raw data, intermediate files and outputs.

## ************************************************
## step 0 - Download raw data - Example shown here.
## ************************************************

#echo "Start step 0 - downloading!"

## Download Plasmodium falciparum - TROPHOZOITES fastq files for both ends
#id="SRR1658689"
#fastqPath="$projectPath/fastqFiles/NHEK/trimmed_11"

#if [ ! -d "$fastqPath" ]; then
#    mkdir -p $fastqPath
#fi
#wget -r "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1215nnn/GSM1215593/suppl/GSM1215593_trimmedAndFiltered-TROPHOZOITES-XL-AGGG-L2_1.fastq.gz" -O $fastqPath/$id\_1.fastq.gz
#wget -r "ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1215nnn/GSM1215593/suppl/GSM1215593_trimmedAndFiltered-TROPHOZOITES-XL-AGGG-L2_2.fastq.gz" -O $fastqPath/$id\_2.fastq.gz

#gunzip $fastqPath/$id\_1.fastq.gz
#gunzip $fastqPath/$id\_2.fastq.gz

## remove end suffix to have matching id between two ends
#mv $fastqPath/$id\_1.fastq $fastqPath/$id\_1.tmp.fastq
#awk '{print $1}' $fastqPath/$id\_1.tmp.fastq > $fastqPath/$id\_1.fastq
#mv $fastqPath/$id\_2.fastq $fastqPath/$id\_2.tmp.fastq
#awk '{print $1}' $fastqPath/$id\_2.tmp.fastq > $fastqPath/$id\_2.fastq
#rm -rf $fastqPath/$id\_*tmp.fastq



## Download IMR90
# id="SRR1658591"
# sraDir="/projects/sratoolkit.2.8.2-1-centos_linux64"
# path="/projects/fastqFiles"
# mkdir -p $path

# ## tar -zxvf sratoolkit.tar.gz
# wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX764/SRX764954/SRR1658591/SRR1658591.sra -O $path/$id.sra
# $sraDir/bin/fastq-dump -F --split-files $path/$id.sra -O $path

## step 1-3: Can be run in parallel.

## ******************
## step 1: Alignment
## ******************
name="NHEK_bwa_uni_trimmed_26"
ref="$projectPath/bin/hg19.fa"
bwaDir="/opt/local/bin" ##"$projectPath/Softwares/bwa-0.5.9" ## need downloading bwa software
samtoolsDir="/opt/local/bin" ## "$projectPath/Softwares/samtools-1.3" ## need dowloading samtools software
fastqDir="$projectPath/fastqFiles/NHEK/trimmed_26"
resultsDir="$projectPath/result/NHEK/trimmed_26"
bin="$projectPath/bin"
nCores=20
seqLength=25
resolution=10000
saveFiles=0
cutsite="GATCGATC" ## for MboI ##"AAGCTAGCTT" for HindIII or (GATCGATC ******** ********) for multiple cutters.
## compile cutsite to trim chimeric reads
g++ -std=c++0x -o $bin/cutsite_trimming_mHiC $bin/cutsite_trimming_mHiC.cpp

## generate reference genome bwa index
$bwaDir/bwa index $ref

## alignment
echo "Start step 1 - alignment!"
bash s1_bwaAlignment.sh "$name" "$ref" "$bwaDir" "$samtoolsDir" "$fastqDir" "$resultsDir/s1" "$bin" "$nCores"  "$resultsDir/mHiC.summary_w${resolution}_s1" "$saveFiles" "$seqLength" "${cutsite[@]}"


## **************************
## step 2: Read ends pairing
## **************************
name="NHEK_bwa_uni_trimmed_26"
resultsDir="$projectPath/result/NHEK/trimmed_26"
resolution=10000

echo "Start step 2 - joining read ends!"
python3 s2_joinEnd.py -r1 ${resultsDir}/s1/${name}_1.sam -r2 ${resultsDir}/s1/${name}_2.sam -o ${resultsDir}/s2/${name}.sam -sf $resultsDir/mHiC.summary_w${resolution}_s2


## *********************************
## step 3: Valid fragment filtering
## *********************************
name="NHEK_bwa_uni_trimmed_26"
resultsDir="$projectPath/result/NHEK/trimmed_26"
bin="$projectPath/bin"
refrag="MboI_resfrag_hg19.bed" #restriction fragment file
resolution=10000
lowerBound=$((resolution * 2))
refragL=50 #$((seqLength * 2))
refragU=500

echo "Start step 3 - categorize read pairs!"
python3 s3_categorizePairs.py -f ${bin}/${refrag} -r ${resultsDir}/s2/${name}.sam -o ${resultsDir}/s3 -l $refragL -u $refragU -d $lowerBound -m "window" -b $resolution -sf $resultsDir/mHiC.summary_w${resolution}_s3

# ## In case, chrM is not needed in downstream analysis
# awk -v OFS="\t" '$2 != "chrM" && $7!="chrM" {print $0}' $validP >$validP.noChrM
# rm -rf $validP
# mv $validP.noChrM $validP


## ***************************************
## step 4 - Remove duplicates and binning.
## ***************************************
name="NHEK_bwa_uni_trimmed_26"
resultsDir="$projectPath/result/NHEK/trimmed_26"
resolution=10000
bin="$projectPath/bin"
validP="${resultsDir}/s3/${name}.validPairs"
validI="${resultsDir}/s4/${name}.validPairs"
minCount=1 #min contact counts allowed

normMethod="KR" #1. "ICE" 2. "KR" 3."None"
ICEmappFile="${bin}/pfal3D7.MboI.w${resolution}" ## mappability file for ICE method
ICEminMap=0.5 ## min mappability threshold for ICE method
ICEmaxIter=150 ## maximum number of iteration for ICE method
KRchromSizeFile=${bin}/chrom_hg19.sizes ## chromosome size file for KR method
KRsparsePerc=10 ## remove *% of sparse regions for KR method
splitByChrom=1
saveSplitContact=0
chrList=($(seq 1 22) X Y) ## If the chromosomes start with "chr", you can just put in the array of chromosome number, such as (1 2 3) or ($(seq 1 22) X Y). If not, the prefix and chromosome number should be both given, such as (chromosome1 chromosome2 chromosome3) or it can be generated by chrList=(chromosome{1..22}). Such chromosome names should be consistent with those in reference genome and the chromosome size file for KR normalization.

echo "Start step 4 - duplicates removal and binning!"
bash s4_bin.sh "$validP" "$validI" "$bin" "$resolution" "$minCount" "$normMethod" "$ICEmappFile" "$ICEminMap" "$ICEmaxIter" "whole" "$KRchromSizeFile" "$KRsparsePerc" "$resultsDir/mHiC.summary_w${resolution}_s4" "$splitByChrom" "$saveSplitContact" "${chrList[@]}"

## **********************
## step 5 - Build prior.
## **********************
name="NHEK_bwa_uni_trimmed_26"
bin="$projectPath/bin"
resultsDir="$projectPath/result/NHEK/trimmed_26"
validI="${resultsDir}/s4/${name}.validPairs"
resolution=10000
splineBin=150
lower=10000
priorName="uniPrior"
normMethod="KR" #"ICE" ##"KR"
chromSizeFile=${bin}/chrom_hg19.sizes
contactFile=$validI.binPairCount.uni

echo "Starts step 5 - prior construction based on uni-reads only!"
python3 $bin/createFitHiCFragments-fixedsize.py --chrLens "$chromSizeFile" --resolution "$resolution" --outFile "$resultsDir/$name.uni.fragments.mHiC"

python3 s5_prior.py -f ${resultsDir}/$name.uni.fragments.mHiC -i $contactFile -o ${resultsDir} -t $validI.binPairCount.uni.KRnorm.bias -b $splineBin -L $lower -r $resolution -p 2


## ************************************************************************************
## step 6 - Generative model to assign probability to multi-reads potential alignments.
## ************************************************************************************
name="NHEK_bwa_uni_trimmed_26"
resultsDir="$projectPath/result/NHEK/trimmed_26"
resolution=10000
prior="${resultsDir}/s5_prior.mhic"
multi="${resultsDir}/s4/${name}.validPairs.MULTI.binPair.multi"
multiKeys="$resultsDir/s4/${name}.validPairs.MULTI.binPair.multiKeys" 
uni="$resultsDir/s4/${name}.validPairs.binPairCount.uni"
filename="${name}.validPairs.binPair.multi"
threshold=0.5

echo "Starts step 6 - assign probability to multi-reads potential alignment positions !"

if [ ! -d ${resultsDir}/s4/sorttmp ]; then
    mkdir -p ${resultsDir}/s4/sorttmp
fi

awk -v OFS="_" '{print $2, $3, $4, $5}' $multi | sort -u -T ${resultsDir}/s4/sorttmp >$multiKeys

rm -rf ${resultsDir}/s4/sorttmp

python3 s6_em.py -p $prior -u $uni -m $multi -mk $multiKeys -t $threshold -o "${resultsDir}/s6" -f $filename


## Merge the multi-reads contact counts with the uni-reads contact counts.

filterT=0.5
multiOut=${resultsDir}/s6/${name}.validPairs.binPair.multi.mHiC

#if [ ! -d ${resultsDir}/s6/sorttmp ]; then
#    mkdir -p ${resultsDir}/s6/sorttmp
#fi

#awk -v OFS="\t" -v fT=$filterT '$6>fT {print $2, $3, $4, $5}' $multiOut | sort -T ${resultsDir}/s6/sorttmp | uniq -c | awk -v OFS="\t" '{print $2, $3, $4, $5, $1}' >$multiOut.binPairCount.multi ## get binPair Count for multi-reads

#rm -rf ${resultsDir}/s6/sorttmp
#mkdir ${resultsDir}/s6/sorttmp

#cat $uni $multiOut.binPairCount.multi | sort -k1,1V -k2,2n -k3,3V -k4,4n | awk -v OFS="\t" '{a[$1" "$2" "$3" "$4]+=$5}END{for (i in a) print i,a[i]}' | sort -k1,1V -k2,2n -k3,3V -k4,4n >${resultsDir}/s6/${name}.validPairs.binPairCount.uniMulti ## merged with uni-reads binpair count

#rm -rf ${resultsDir}/s6/sorttmp