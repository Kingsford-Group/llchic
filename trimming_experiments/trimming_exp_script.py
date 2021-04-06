### This script contain codes used for the trimming experiment
### Date: March, 2021

import csv,pysam,sys,pdb,copy,random
from Bio import SeqIO
import numpy as np

def compare_srrname(name_1,name_2):
    #return the index of the smaller one
    if(name_1==name_2):
        return 1
    srr_id_1,srr_num_1 = name_1.split(".")
    srr_id_2,srr_num_2 = name_2.split(".")
    if(srr_id_1<srr_id_2 or (srr_id_1==srr_id_2 and int(srr_num_1)<int(srr_num_2))):
        return 0
    return 2

def Next_iter(iter_list):
    try:
        return next(iter_list)
    except StopIteration:
        return 0

def find_sd_pair_2(pos_list_1,pos_list_2):
    sd = -1
    index_1 = 0
    index_2 = 0
    for i in range(len(pos_list_1)):
        for j in range(len(pos_list_2)):
            if(pos_list_1[i][0]==pos_list_2[j][0]):
                dist = abs(pos_list_1[i][1]-pos_list_2[j][1])
                if(sd==-1 or sd>dist):
                    sd = dist
                    index_1 = i
                    index_2 = j
    return (index_1,index_2)

def get_read_name(read):
    name = read.qname
    return name.split("/",1)[0]

### Generate uni-mapped read pairs
sam_file = "/mnt/disk51/user/yihangs/PLOSCOMP_revise/NHEK_mHiC/mHiC-master/result/NHEK/s2/NHEK_all.sam"
sam_pt = pysam.AlignmentFile(sam_file,"r")
r1 = Next_iter(sam_pt.fetch(until_eof=True))
r2 = Next_iter(sam_pt.fetch(until_eof=True))
output_path = "/mnt/disk69/user/yihangs/HiC_correction_revision/PLOSComp_revise_1/trim_HiC_data/"
output_1_file = output_path + "NHEK_bwa_uni_1.sam"
output_2_file = output_path + "NHEK_bwa_uni_2.sam"
output_1 = pysam.AlignmentFile(output_1_file,"w",template=sam_pt)
output_2 = pysam.AlignmentFile(output_2_file,"w",template=sam_pt)
with open("/mnt/disk51/user/yihangs/PLOSCOMP_revise/NHEK_mHiC/mHiC-master/result/NHEK/s3/NHEK_all.uniMulti","r")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        if(line[1]!="Uni"):
            continue
        #pdb.set_trace()
        name_comp = compare_srrname(line[0],r1.query_name)
        if(name_comp==0):
            print("what? How can that be!")
            print(line)
            break
        while(compare_srrname(line[0],r1.query_name)==2):
            r1 = Next_iter(sam_pt.fetch(until_eof=True))
            r2 = Next_iter(sam_pt.fetch(until_eof=True))
            assert r1.query_name==r2.query_name
        #check same
        assert compare_srrname(line[0],r1.query_name)==1
        output_1.write(r1)
        output_2.write(r2)

sam_pt.close()
output_1.close()
output_2.close()


### Generate the ground truth for the trimming experiment, i.e, read pairs that both BWA and bowtie2 uniquely map to the same location. 
file_1 = "NHEK_bwa_uni_1.bam"
file_2 = "NHEK_bwa_uni_2.bam"
bowtie_file_path = "/mnt/disk73/user/yihangs/HiC_artifact/origin_hicpro/NHEK/hic_results/bowtie_results/bwt2/sample/"
SRR_name = "SRR1658689"
bowtie_pt = pysam.AlignmentFile(bowtie_file_path+"SRR1658689_hg19.bwt2pairs.bam","rb")
br1 = Next_iter(bowtie_pt.fetch(until_eof=True))
br2 = Next_iter(bowtie_pt.fetch(until_eof=True))
output_1_file = "NHEK_uni_bwa_bowtie_diff_1.bam"
output_2_file = "NHEK_uni_bwa_bowtie_diff_2.bam"
same_count = 0
uni_same_count = 0
count = 0
reso = 10000
output_file_1 = open("NHEK_bwa_bowtie.binPair","w")
fwriter_1 = csv.writer(output_file_1,delimiter='\t')
output_file_2 = open("NHEK_bwa_bowtie_alluni.binPair","w")
fwriter_2 = csv.writer(output_file_2,delimiter='\t')
with pysam.AlignmentFile(file_1,"rb")as hr1, pysam.AlignmentFile(file_2,"rb")as hr2:
    #output_1 = pysam.AlignmentFile(output_1_file,"w",template=hr1)
    #output_2 = pysam.AlignmentFile(output_2_file,"w",template=hr2)
    for r1,r2 in zip(hr1.fetch(until_eof=True), hr2.fetch(until_eof=True)):
        count+=1
        if(count%100000==0):
            print("##"+str(count))
        if(r1.query_name.split(".")[0]>SRR_name):
            bowtie_pt.close()
            SRR_name = r1.query_name.split(".")[0]
            bowtie_pt = pysam.AlignmentFile(bowtie_file_path+SRR_name+"_hg19.bwt2pairs.bam","rb")
            br1 = Next_iter(bowtie_pt.fetch(until_eof=True))
            br2 = Next_iter(bowtie_pt.fetch(until_eof=True))
        while(compare_srrname(r1.query_name,br1.query_name)==2):
            br1 = Next_iter(bowtie_pt.fetch(until_eof=True))
            br2 = Next_iter(bowtie_pt.fetch(until_eof=True))
            if(type(br1)==type(0)):
                if(SRR_name=="SRR1658691"):
                    break
                else:
                    bowtie_pt.close()
                    SRR_name = "SRR" + str(int(SRR_name[3:])+1)
                    bowtie_pt = pysam.AlignmentFile(bowtie_file_path+SRR_name+"_hg19.bwt2pairs.bam","rb")
                    br1 = Next_iter(bowtie_pt.fetch(until_eof=True))
                    br2 = Next_iter(bowtie_pt.fetch(until_eof=True))
        if(type(br1)==type(0) and SRR_name=="SRR1658691"):
            break
        name_comp = compare_srrname(r1.query_name,br1.query_name)
        if(name_comp==0):
            continue
        if(r1.reference_name==br1.reference_name and abs(r1.reference_start-br1.reference_start)<100 and r2.reference_name==br2.reference_name and abs(r2.reference_start-br2.reference_start)<100):
            same_count+=1
            fwriter_1.writerow([r1.query_name,r1.reference_name,str(int((r1.reference_start+1)/reso)*reso+int(reso/2)),r2.reference_name,str(int((r2.reference_start+1)/reso)*reso+int(reso/2))])
            if(br1.has_tag('XS')==0 and br2.has_tag('XS')==0):
                uni_same_count+=1
                fwriter_2.writerow([r1.query_name,r1.reference_name,str(int((r1.reference_start+1)/reso)*reso+int(reso/2)),r2.reference_name,str(int((r2.reference_start+1)/reso)*reso+int(reso/2))])
            continue
        if(r1.reference_name==br2.reference_name and abs(r1.reference_start-br2.reference_start)<100 and r2.reference_name==br1.reference_name and abs(r2.reference_start-br1.reference_start)<100):
            same_count+=1
            fwriter_1.writerow([r1.query_name,r1.reference_name,str(int((r1.reference_start+1)/reso)*reso+int(reso/2)),r2.reference_name,str(int((r2.reference_start+1)/reso)*reso+int(reso/2))])
            if(br1.has_tag('XS')==0 and br2.has_tag('XS')==0):
                fwriter_2.writerow([r1.query_name,r1.reference_name,str(int((r1.reference_start+1)/reso)*reso+int(reso/2)),r2.reference_name,str(int((r2.reference_start+1)/reso)*reso+int(reso/2))])
                uni_same_count+=1
            continue
        #output_1.write(br1)
        #output_2.write(br2)
        
bowtie_pt.close()
output_file_1.close()
output_file_2.close()


### Generate trimmed reads
fastq_file_1 = "NHEK_bwa_uni_1.fastq"
fastq_file_2 = "NHEK_bwa_uni_2.fastq"
trim_unit = 55
t_len_1 = 46
t_len_2 = 36
output_file_1 = "NHEK_bwa_uni_trimmed_" + str(trim_unit) + "_1.fastq"
output_file_2 = "NHEK_bwa_uni_trimmed_" + str(trim_unit) + "_2.fastq"
output_handle_1 = open(output_file_1,"w")
output_handle_2 = open(output_file_2,"w")
#output_sequences_1 = []
#output_sequences_2 = []
counter = 0
for fastq_1,fastq_2 in zip(SeqIO.parse(fastq_file_1, "fastq"),SeqIO.parse(fastq_file_2, "fastq")):
    if(fastq_1.name.split('/')[0]!=fastq_2.name.split('/')[0]):
        print("fastq name error!")
        sys.exit(1)
    counter+=1
    if(counter%10000==0):
        print(counter)
    phred_info_1 = fastq_1.letter_annotations['phred_quality']
    fastq_1.letter_annotations = {}
    phred_info_2 = fastq_2.letter_annotations['phred_quality']
    fastq_2.letter_annotations = {}
    fastq_1.seq = fastq_1.seq[:t_len_1]
    fastq_1.letter_annotations['phred_quality'] = phred_info_1[:t_len_1]
    fastq_2.seq = fastq_2.seq[:t_len_2]
    fastq_2.letter_annotations['phred_quality'] = phred_info_2[:t_len_2]
    SeqIO.write(fastq_1,output_handle_1,"fastq")
    SeqIO.write(fastq_2,output_handle_2,"fastq")
    #output_sequences_1.append(fastq_1)
    #output_sequences_2.append(fastq_2)

output_handle_1.close()
output_handle_2.close()


### Obtain multi-mapped trimmed read pairs, these are the pairs used for calculating allocation accuracy
NHEK_all_bwa_bowtie = {}
count = 0
with open("NHEK_bwa_bowtie_alluni.binPair","r")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        count+=1
        NHEK_all_bwa_bowtie[line[0]] = [line[1],int(line[2]),line[3],int(line[4])]
        if(count%1000000==0):
            print("##"+str(count))


trim_num = str(55)
output_multitouni = open("NHEK_bwa_uni_trimmed_"+trim_num+".compare.MultitoUniList","w")
fwriter_1 = csv.writer(output_multitouni,delimiter='\t')
output_multi = open("NHEK_bwa_uni_trimmed_"+trim_num+".compare.MultiList","w")
fwriter_2 = csv.writer(output_multi,delimiter='\t')

with open("/mnt/disk78/user/yihangs/PLOSComp_revise_1/trimming_exp/mHiC-master/result/NHEK/trimmed_"+trim_num+"/s3/NHEK_bwa_uni_trimmed_"+trim_num+".uniMulti","r")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        if(line[1]=='MultiToUni'):
            if(line[0] in NHEK_all_bwa_bowtie.keys()):
                fwriter_1.writerow(line)
        elif(line[1]=='Multi'):
            if(line[0] in NHEK_all_bwa_bowtie.keys()):
                fwriter_2.writerow(line)

output_multitouni.close()
output_multi.close()



### Calcaulte the allocation accuracy of each algorithm
trim_num = str(55)
NHEK_all_bwa_bowtie = {}
count = 0
###Read the ground truth
with open("NHEK_bwa_bowtie_alluni.binPair","r")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        count+=1
        NHEK_all_bwa_bowtie[line[0]] = [line[1],int(line[2]),line[3],int(line[4])]
        if(count%1000000==0):
            print("##"+str(count))

###Read the information of multi-mapped trimmed read pairs
multitouni_set = set()
with open("NHEK_bwa_uni_trimmed_"+trim_num+".compare.MultitoUniList","r")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        multitouni_set.add(line[0])

multi_set = set()
with open("NHEK_bwa_uni_trimmed_"+trim_num+".compare.MultiList","r")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        multi_set.add(line[0])

multi_set = multi_set | multitouni_set

###mHi-C allocation accuracy
right_count_intra = 0
right_count_inter = 0
wrong_count_intra = 0
wrong_count_inter = 0
#wrong_list = {}
with open("/mnt/disk78/user/yihangs/PLOSComp_revise_1/trimming_exp/mHiC-master/result/NHEK/trimmed_"+trim_num+"/s4/NHEK_bwa_uni_trimmed_"+trim_num+".validPairs.binPair.uni","r")as f:
#with open("/mnt/disk69/user/yihangs/HiC_correction_revision/PLOSComp_revise_1/test_mHiC/mHiC-master/result/NHEK/trimmed_"+trim_num+"/s4/NHEK_bwa_uni_trimmed_"+trim_num+".validPairs.binPair.uni","r")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        if(line[0] in multi_set):
            ref_1,pos_1,ref_2,pos_2 = NHEK_all_bwa_bowtie[line[0]]
            if(line[1]==ref_1 and int(line[2])==pos_1 and line[3]==ref_2 and int(line[4])==pos_2):
                if(ref_1==ref_2):
                    right_count_intra+=1
                else:
                    right_count_inter+=1
            elif(line[1]==ref_2 and int(line[2])==pos_2 and line[3]==ref_1 and int(line[4])==pos_1):
                if(ref_1==ref_2):
                    right_count_intra+=1
                else:
                    right_count_inter+=1
            else:
                if(ref_1==ref_2):
                    wrong_count_intra+=1
                else:
                    wrong_count_inter+=1
                #wrong_list[line[0]] = line
                #print(NHEK_all_bwa_bowtie[line[0]])
                #print(line)
                #if(wrong_count>=100):
                    #break
            if((right_count_intra+right_count_inter)%100000==0):
                print("##"+str(right_count_intra+right_count_inter))
                print(line[0])

tmp_container = ["SRR",[],[]]
with open("/mnt/disk78/user/yihangs/PLOSComp_revise_1/trimming_exp/mHiC-master/result/NHEK/trimmed_"+trim_num+"/s6/NHEK_bwa_uni_trimmed_"+trim_num+".validPairs.binPair.multi.mHiC","r")as f:
#with open("/mnt/disk69/user/yihangs/HiC_correction_revision/PLOSComp_revise_1/test_mHiC/mHiC-master/result/NHEK/trimmed_"+trim_num+"/s6/NHEK_bwa_uni_trimmed_"+trim_num+".validPairs.binPair.multi.mHiC","r")as f:
    freader = csv.reader(f,delimiter='\t')
    for line in freader:
        if(line[0]==tmp_container[0]):
            tmp_container[1].append([line[1],int(line[2]),line[3],int(line[4])])
            tmp_container[2].append(float(line[5]))
        else:
            if(len(tmp_container[2])>0):
                #pdb.set_trace()
                tmp_container[2] = np.array(tmp_container[2])
                max_idx = np.random.choice(np.where(tmp_container[2]==np.max(tmp_container[2]))[0])
                chr_1,pos_pred_1,chr_2,pos_pred_2 = tmp_container[1][max_idx]
                ref_1,pos_1,ref_2,pos_2 = NHEK_all_bwa_bowtie[tmp_container[0]]
                if(chr_1==ref_1 and pos_pred_1==pos_1 and chr_2==ref_2 and pos_pred_2==pos_2):
                    if(ref_1==ref_2):
                        right_count_intra+=1
                    else:
                        right_count_inter+=1
                elif(chr_1==ref_2 and pos_pred_1==pos_2 and chr_2==ref_1 and pos_pred_2==pos_1):
                    if(ref_1==ref_2):
                        right_count_intra+=1
                    else:
                        right_count_inter+=1
                else:
                    if(ref_1==ref_2):
                        wrong_count_intra+=1
                    else:
                        wrong_count_inter+=1
                if((right_count_intra+right_count_inter)%100000==0):
                    print("##"+str(right_count_intra+right_count_inter))
                    print(line[0])
            if(line[0] in multi_set):
                tmp_container[0] = line[0]
                tmp_container[1] = [[line[1],int(line[2]),line[3],int(line[4])]]
                tmp_container[2] = [float(line[5])]
            else:
                tmp_container[0] = 'SRR'
                tmp_container[1] = []
                tmp_container[2] = []


if(len(tmp_container[2])>0):
    tmp_container[2] = np.array(tmp_container[2])
    max_idx = np.random.choice(np.where(tmp_container[2]==np.max(tmp_container[2]))[0])
    chr_1,pos_pred_1,chr_2,pos_pred_2 = tmp_container[1][max_idx]
    ref_1,pos_1,ref_2,pos_2 = NHEK_all_bwa_bowtie[tmp_container[0]]
    if(chr_1==ref_1 and pos_pred_1==pos_1 and chr_2==ref_2 and pos_pred_2==pos_2):
        if(ref_1==ref_2):
            right_count_intra+=1
        else:
            right_count_inter+=1
    elif(chr_1==ref_2 and pos_pred_1==pos_2 and chr_2==ref_1 and pos_pred_2==pos_1):
        if(ref_1==ref_2):
            right_count_intra+=1
        else:
            right_count_inter+=1
    else:
        if(ref_1==ref_2):
            wrong_count_intra+=1
        else:
            wrong_count_inter+=1


###LLC-HiC allocation accuracy
right_count_llc_intra = 0
right_count_llc_inter = 0
wrong_count_llc_intra = 0
wrong_count_llc_inter = 0
reso=10000
file_1 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/realign_result/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R1_hg19.bwt2merged.bam"
file_2 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/realign_result/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R2_hg19.bwt2merged.bam"
with pysam.AlignmentFile(file_1,"rb")as hr1, pysam.AlignmentFile(file_2,"rb")as hr2:
    for r1,r2 in zip(hr1.fetch(until_eof=True), hr2.fetch(until_eof=True)):
        if(r1.query_name.split("/")[0]!=r2.query_name.split("/")[0]):
            print("query name error!")
        if(r1.query_name.split("/")[0] in multi_set):
            if((right_count_llc_intra+right_count_llc_inter)%100000==0):
                print("##"+str(right_count_llc_intra+right_count_llc_inter))
                print(r1.query_name.split("/")[0])
            ref_1,pos_1,ref_2,pos_2 = NHEK_all_bwa_bowtie[r1.query_name.split("/")[0]]
            if(r1.reference_name==ref_1 and r2.reference_name==ref_2 and abs(pos_1-r1.reference_start-1)<reso/2 and abs(pos_2-r2.reference_start-1)<reso/2):
                if(ref_1==ref_2):
                    right_count_llc_intra+=1
                else:
                    right_count_llc_inter+=1
                continue
            if(r1.reference_name==ref_2 and r2.reference_name==ref_1 and abs(pos_2-r1.reference_start-1)<reso/2 and abs(pos_1-r2.reference_start-1)<reso/2):
                if(ref_1==ref_2):
                    right_count_llc_intra+=1
                else:
                    right_count_llc_inter+=1
                continue
            if(ref_1==ref_2):
                wrong_count_llc_intra+=1
            else:
                wrong_count_llc_inter+=1


###AlignerChoose allocation accuracy
right_count_bowtie_intra = 0
right_count_bowtie_inter = 0
wrong_count_bowtie_intra = 0
wrong_count_bowtie_inter = 0
reso=10000
file_1 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/trimmed_"+trim_num+"/hic_results/bowtie_results/bwt2/sample/NHEK_bwa_uni_trimmed_"+trim_num+"_R1_hg19.bwt2merged.bam"
file_2 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/trimmed_"+trim_num+"/hic_results/bowtie_results/bwt2/sample/NHEK_bwa_uni_trimmed_"+trim_num+"_R2_hg19.bwt2merged.bam"
with pysam.AlignmentFile(file_1,"rb")as hr1, pysam.AlignmentFile(file_2,"rb")as hr2:
    for r1,r2 in zip(hr1.fetch(until_eof=True), hr2.fetch(until_eof=True)):
        if(r1.query_name.split("/")[0]!=r2.query_name.split("/")[0]):
            print("query name error!")
        if(r1.query_name.split("/")[0] in multi_set):
            if((right_count_bowtie_intra+right_count_bowtie_inter)%100000==0):
                print("##"+str(right_count_bowtie_intra+right_count_bowtie_inter))
                print(r1.query_name.split("/")[0])
            ref_1,pos_1,ref_2,pos_2 = NHEK_all_bwa_bowtie[r1.query_name.split("/")[0]]
            if(r1.reference_name==ref_1 and r2.reference_name==ref_2 and abs(pos_1-r1.reference_start-1)<reso/2 and abs(pos_2-r2.reference_start-1)<reso/2):
                if(ref_1==ref_2):
                    right_count_bowtie_intra+=1
                else:
                    right_count_bowtie_inter+=1
                #right_count_bowtie+=1
                continue
            if(r1.reference_name==ref_2 and r2.reference_name==ref_1 and abs(pos_2-r1.reference_start-1)<reso/2 and abs(pos_1-r2.reference_start-1)<reso/2):
                if(ref_1==ref_2):
                    right_count_bowtie_intra+=1
                else:
                    right_count_bowtie_inter+=1
                #right_count_bowtie+=1
                continue
            if(ref_1==ref_2):
                wrong_count_bowtie_intra+=1
            else:
                wrong_count_bowtie_inter+=1


### SDChoose allocation accuracy
multi_align_file_1 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/LLCHiC/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R1_hg19_multialign.bam"
multi_align_file_2 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/LLCHiC/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R2_hg19_multialign.bam"

remap_file_1 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/LLCHiC/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R1_hg19_multireads.bwt2glob.bam"
remap_file_2 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/LLCHiC/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R2_hg19_multireads.bwt2glob.bam"

new_multi_align_file_1 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/shortest_dist/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R1_hg19_new_multialign.bam"
new_multi_align_file_2 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/shortest_dist/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R2_hg19_new_multialign.bam"

rp1 = pysam.AlignmentFile(remap_file_1,"rb")
rp2 = pysam.AlignmentFile(remap_file_2,"rb")
next_read_1 = Next_iter(rp1.fetch(until_eof=True))
next_read_2 = Next_iter(rp2.fetch(until_eof=True))

counter = 0

with pysam.AlignmentFile(multi_align_file_1,"rb")as hr1, pysam.AlignmentFile(multi_align_file_2,"rb")as hr2:
    output_1 = pysam.AlignmentFile(new_multi_align_file_1,"wb",template=hr1)
    output_2 = pysam.AlignmentFile(new_multi_align_file_2,"wb",template=hr2)
    for r1,r2 in zip(hr1.fetch(until_eof=True), hr2.fetch(until_eof=True)):
        counter+=1
        reads_list_1 = [copy.deepcopy(r1)]
        pos_info_1 = [[r1.reference_name,r1.reference_start+1]]
        reads_list_2 = [copy.deepcopy(r2)]
        pos_info_2 = [[r2.reference_name,r2.reference_start+1]]
        if(r1.has_tag('XS')):
            while True:
                if(type(next_read_1)!=type(r1)):
                    break
                if(next_read_1.query_name!=r1.query_name):
                    break
                reads_list_1.append(copy.deepcopy(next_read_1))
                pos_info_1.append([next_read_1.reference_name,next_read_1.reference_start+1])
                next_read_1 = Next_iter(rp1.fetch(until_eof=True))
        if(r2.has_tag("XS")):
            while True:
                if(type(next_read_2)!=type(r2)):
                    break
                if(next_read_2.query_name!=r2.query_name):
                    break
                reads_list_2.append(copy.deepcopy(next_read_2))
                pos_info_2.append([next_read_2.reference_name,next_read_2.reference_start+1])
                next_read_2 = Next_iter(rp2.fetch(until_eof=True))
        ind_1,ind_2 = find_sd_pair_2(pos_info_1,pos_info_2)
        output_1.write(reads_list_1[ind_1])
        output_2.write(reads_list_2[ind_2])
        if(counter % 100000 == 0):
            print("##",counter)


hr1.close()
hr2.close()
rp1.close()
rp2.close()
output_1.close()
output_2.close()


output_1_file = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/shortest_dist/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R1.bam"
output_2_file = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/shortest_dist/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R2.bam"
R1file = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/trimmed_"+trim_num+"/hic_results/bowtie_results/bwt2/sample/NHEK_bwa_uni_trimmed_"+trim_num+"_R1_hg19.bwt2merged.bam"
R2file = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/trimmed_"+trim_num+"/hic_results/bowtie_results/bwt2/sample/NHEK_bwa_uni_trimmed_"+trim_num+"_R2_hg19.bwt2merged.bam"

nr1 = pysam.AlignmentFile(new_multi_align_file_1,"rb")
nr2 = pysam.AlignmentFile(new_multi_align_file_2,"rb")
new_read_1 = Next_iter(nr1.fetch(until_eof=True))
new_read_2 = Next_iter(nr2.fetch(until_eof=True))

counter = 0
with pysam.AlignmentFile(R1file,"rb")as hr1, pysam.AlignmentFile(R2file,"rb")as hr2:
    output_1 = pysam.AlignmentFile(output_1_file,"wb",template=hr1)
    output_2 = pysam.AlignmentFile(output_2_file,"wb",template=hr2)
    for r1,r2 in zip(hr1.fetch(until_eof=True), hr2.fetch(until_eof=True)):
        counter+=1
        if(type(new_read_1)==type(r1) and new_read_1.query_name==r1.query_name):
            output_1.write(new_read_1)
            new_read_1 = Next_iter(nr1.fetch(until_eof=True))
        else:
            output_1.write(r1)
        if(type(new_read_2)==type(r2) and new_read_2.query_name==r2.query_name):
            output_2.write(new_read_2)
            new_read_2 = Next_iter(nr2.fetch(until_eof=True))
        else:
            output_2.write(r2)
        if(counter % 1000000 == 0):
            print("##",counter)

hr1.close()
hr2.close()
nr1.close()
nr2.close()
output_1.close()
output_2.close()

right_count_sd_bowtie_intra = 0
right_count_sd_bowtie_inter = 0
wrong_count_sd_bowtie_intra = 0
wrong_count_sd_bowtie_inter = 0
reso=10000
file_1 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/shortest_dist/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R1.bam"
file_2 = "/mnt/disk76/user/yihangs/PLOSCOMP_revise/trimming_experiment/LLCHiC/NHEK/shortest_dist/trimmed_"+trim_num+"/NHEK_bwa_uni_trimmed_"+trim_num+"_R2.bam"
with pysam.AlignmentFile(file_1,"rb")as hr1, pysam.AlignmentFile(file_2,"rb")as hr2:
    for r1,r2 in zip(hr1.fetch(until_eof=True), hr2.fetch(until_eof=True)):
        if(r1.query_name.split("/")[0]!=r2.query_name.split("/")[0]):
            print("query name error!")
        if(r1.query_name.split("/")[0] in multi_set):
            if((right_count_sd_bowtie_intra+right_count_sd_bowtie_inter)%100000==0):
                print("##"+str(right_count_sd_bowtie_intra+right_count_sd_bowtie_inter))
                print(r1.query_name.split("/")[0])
            ref_1,pos_1,ref_2,pos_2 = NHEK_all_bwa_bowtie[r1.query_name.split("/")[0]]
            if(r1.reference_name==ref_1 and r2.reference_name==ref_2 and abs(pos_1-r1.reference_start-1)<reso/2 and abs(pos_2-r2.reference_start-1)<reso/2):
                if(ref_1==ref_2):
                    right_count_sd_bowtie_intra+=1
                else:
                    right_count_sd_bowtie_inter+=1
                continue
            if(r1.reference_name==ref_2 and r2.reference_name==ref_1 and abs(pos_2-r1.reference_start-1)<reso/2 and abs(pos_1-r2.reference_start-1)<reso/2):
                if(ref_1==ref_2):
                    right_count_sd_bowtie_intra+=1
                else:
                    right_count_sd_bowtie_inter+=1
                continue
            if(ref_1==ref_2):
                wrong_count_sd_bowtie_intra+=1
            else:
                wrong_count_sd_bowtie_inter+=1


### Random allocation accuracy
right_count_rd_bwa_intra = 0
right_count_rd_bwa_inter = 0
wrong_count_rd_bwa_intra = 0
wrong_count_rd_bwa_inter = 0
reso = 10000
#file_input = "/mnt/disk69/user/yihangs/HiC_correction_revision/PLOSComp_revise_1/test_mHiC/mHiC-master/result/NHEK/trimmed_"+trim_num+"/s2/NHEK_bwa_uni_trimmed_"+trim_num+".sam"
file_input = "/mnt/disk78/user/yihangs/PLOSComp_revise_1/trimming_exp/mHiC-master/result/NHEK/trimmed_"+trim_num+"/s2/NHEK_bwa_uni_trimmed_"+trim_num+".sam"
with pysam.AlignmentFile(file_input,"r")as hr1:
    while True:
        r1 = Next_iter(hr1.fetch(until_eof=True))
        r2 = Next_iter(hr1.fetch(until_eof=True))
        if(type(r1)==type(0)):
            break
        if(r1.query_name in multi_set):
            if((right_count_rd_bwa_intra+right_count_rd_bwa_inter)%100000==0):
                print("##"+str(right_count_rd_bwa_intra+right_count_rd_bwa_inter))
                print(r1.query_name)
            ref_1,pos_1,ref_2,pos_2 = NHEK_all_bwa_bowtie[r1.query_name]
            chr_1 = r1.reference_name
            pos_pred_1 = r1.reference_start + 1
            chr_2 = r2.reference_name
            pos_pred_2 = r2.reference_start + 1
            if r1.has_tag("XA"):
                alter_pos_1 = r1.tags[-1][1].split(";")
                if(len(alter_pos_1[-1])==0):
                    alter_pos_1 = alter_pos_1[:-1]
                rand_num = random.randint(0,len(alter_pos_1))
                if(rand_num>0):
                    pos_tmp = alter_pos_1[rand_num-1].split(",")[:2]
                    chr_1 = pos_tmp[0]
                    pos_pred_1 = abs(int(pos_tmp[1]))
            if r2.has_tag("XA"):
                alter_pos_2 = r2.tags[-1][1].split(";")
                if(len(alter_pos_2[-1])==0):
                    alter_pos_2 = alter_pos_2[:-1]
                rand_num = random.randint(0,len(alter_pos_2))
                if(rand_num>0):
                    pos_tmp = alter_pos_2[rand_num-1].split(",")[:2]
                    chr_2 = pos_tmp[0]
                    pos_pred_2 = abs(int(pos_tmp[1]))
            if(chr_1==ref_1 and abs(pos_pred_1-pos_1)<reso/2 and chr_2==ref_2 and abs(pos_pred_2-pos_2)<reso/2):
                if(ref_1==ref_2):
                    right_count_rd_bwa_intra+=1
                else:
                    right_count_rd_bwa_inter+=1
                continue
            if(chr_1==ref_2 and abs(pos_pred_1-pos_2)<reso/2 and chr_2==ref_1 and abs(pos_pred_2-pos_1)<reso/2):
                if(ref_1==ref_2):
                    right_count_rd_bwa_intra+=1
                else:
                    right_count_rd_bwa_inter+=1
                continue
            if(ref_1==ref_2):
                wrong_count_rd_bwa_intra+=1
            else:
                wrong_count_rd_bwa_inter+=1