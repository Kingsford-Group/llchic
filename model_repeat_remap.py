import csv
import pdb
import time
import json
import numpy as np
import random
import os
from os import listdir
from os.path import isfile,join
import sys
import pysam
import copy
import vcf
import pickle

cigar_op = {"M","I","D","N","S","H","P","=","X"}
cigar_op_q_only = {"I","S"}
cigar_op_r = {"D","N","M","=","X"}
cigar_op_q = {"M","I","S","=","X"}
cigar_op_r_only = {"D","N"}
cigar_op_no_r_q = {"H","P"}
cigar_op_q_r = {"M","=","X"}
base_set = {"A","T","C","G"}
base_pair = {"A":"T","C":"G","T":"A","G":"C","N":"N"}
num_set = {"0","1","2","3","4","5","6","7","8","9"}


class read_info(object):
    def __init__(self,read_list,inter_read,pysam_form=0):
        if(pysam_form==0):
            self.name = read_list[0]
            self.flag = int(read_list[1])
            self.chr = read_list[2]
            self.pos = int(read_list[3])
            self.mapq = int(read_list[4])
            self.cigar = read_list[5]
            self.inter_read = inter_read
            self.seq = read_list[9]
            self.seq_len = len(self.seq)
            self.qual = read_list[10]
            self.tags = {}
            for i in range(11,len(read_list)):
                self.tags[read_list[i][:2]] = read_list[i][5:]
        else:
            self.name = read_list.query_name
            self.flag = read_list.flag
            self.chr = read_list.reference_name
            self.pos = read_list.reference_start+1
            self.mapq = read_list.mapq
            self.cigar = read_list.cigarstring
            self.inter_read = inter_read
            self.seq = read_list.query_sequence
            self.seq_len = read_list.query_length
            self.qual = read_list.qual
            self.tags = {}
            for tag_name,val in read_list.tags:
                self.tags[tag_name] = val
    def decompose_cigar(self):
        self.cigar_list = []
        num_begin = 0
        for i in range(len(self.cigar)):
            if(self.cigar[i] in cigar_op):
                self.cigar_list.append([self.cigar[i],int(self.cigar[num_begin:i])])
                num_begin = i+1
    def decompose_MD_tag(self):
        self.MD_list = []
        num_begin = 0
        for i in range(len(self.tags["MD"])):
            if(self.tags["MD"][i] not in num_set):
                if(self.tags["MD"][i-1] in num_set):
                    self.MD_list.append(["M",int(self.tags["MD"][num_begin:i])])
                    if(self.tags["MD"][i] in base_set):
                        self.MD_list.append([self.tags["MD"][i],1])
                num_begin = i+1
        if(self.tags["MD"][-1] not in base_set):
            self.MD_list.append(["M",int(self.tags["MD"][num_begin:])])
    #insertion base include insertion and soft-clipped, this function must be used after decompose cigar
    def get_insert_base(self):
        self.insert_base = []
        seq_point = 0
        for info in self.cigar_list:
            if(info[0] in {"M","=","X"}):
                seq_point+=info[1]
            elif(info[0] in {"I","S"}):
                self.insert_base+=[[i,self.seq[i],self.qual[i]] for i in range(seq_point,seq_point+info[1])]
                seq_point+=info[1]
    #get mismatch base from MD tag, this function must be used after decompose MD tag and get insert base
    def get_mismatch_base(self):
        self.mismatch_base = []
        insert_base_index = 0
        seq_point = 0
        for info in self.MD_list:
            #print(seq_point)
            if(info[0]=="M"):
                seq_point+=info[1]
            for i in range(insert_base_index,len(self.insert_base)):
                if(self.insert_base[i][0]<=seq_point):
                    seq_point+=1
                    if(i==len(self.insert_base)-1):
                        insert_base_index = i+1
                        break
                else:
                    insert_base_index = i
                    break
            if(info[0] in base_set):
                self.mismatch_base.append([seq_point,self.seq[seq_point],self.qual[seq_point],info[0]])
                seq_point+=1
            #print(seq_point)
    #check if the sum of bases in M,I,S,=,X is equal to the sequence length
    def check_cigar_decomp_correct(self):
        len_from_cigar=sum(info[1] for info in self.cigar_list if info[0] in {"M","I","S","=","X"})
        if(len_from_cigar!=self.seq_len):
            print("Error:sequence length form cigar is not equal to the true length!")
            exit()
    #make sure that the qual is ASCII of base quality plus 33
    def convert_qual_to_num(self):
        self.qual_val = [ord(one_char) for one_char in self.qual]
    def get_delete_base_num(self):
        delete_base_num = 0
        for info in self.cigar_list:
            if(info[0]=="D"):
                delete_base_num+=info[1]
        return delete_base_num
    def search_gt(self,prob_val,gt_file={}):
        if_find = 0
        for r_pos in range(self.pos,self.r_pos_largest):
            if(r_pos in gt_file.data[self.chr].keys()):
                if_find = 1
                break
        if(if_find):
            while(True):
                one_record = gt_file.data[self.chr][r_pos].record
                prob_val = prob_val*sum([one_record.samples[i]["AD"][0] for i in range(len(one_record.samples))])/one_record.INFO["DP"]
                if(gt_file.data[self.chr][r_pos].next_chr==self.chr and gt_file.data[self.chr][r_pos].next_loc<self.r_pos_largest):
                    r_pos = gt_file.data[self.chr][r_pos].next_loc
                else:
                    break
        return prob_val   
    def get_correspond_q(self,r_pos_val):
        for tag,length,tmp_pos_q,tmp_pos_r in self.q_r_pos_list:
            if(tag=='I'):
                continue
            if(tag=='M' and r_pos_val<tmp_pos_r+length):
                return ('M',r_pos_val-tmp_pos_r+tmp_pos_q)
            if(tag=='D' and r_pos_val<tmp_pos_r+length):
                return ('D',tmp_pos_q)  
    #This function must be used after getting mismatch base
    def map_prob(self,basequal_prob={},if_gt=0,gt_file={},gt_record={}):
        delete_base_qual = int(np.mean(self.qual_val))
        qual_prob = np.array([basequal_prob[str(i)] for i in self.qual_val])
        if(if_gt==0 or self.chr not in gt_file.data.keys() or self.r_pos_largest<=gt_record[self.chr][self.pos-1]):
            prob_val = 1
            for info in self.insert_base:
                prob_val = prob_val*qual_prob[info[0]]/4.0
                qual_prob[[info[0]]] = 0.0
                #prob_val = prob_val*basequal_prob[str(ord(info[2]))]
            for info in self.mismatch_base:
                prob_val = prob_val*qual_prob[info[0]]/3.0
                qual_prob[[info[0]]] = 0.0
                #prob_val = prob_val*basequal_prob[str(ord(info[2]))]
            prob_val = prob_val*np.prod(1-qual_prob)
            prob_val = prob_val*pow(basequal_prob[str(delete_base_qual)],self.get_delete_base_num())
            return prob_val
        
        prob_val = 1
        for v_type,q_pos,r_pos,v_alt,v_ref in self.variants:
            alt_find = 0
            if(r_pos in gt_file.data[self.chr].keys()):
                one_record = gt_file.data[self.chr][r_pos].record
                for alt in range(len(one_record.ALT)):
                    if(str(one_record.ALT[alt])==v_alt and v_ref==one_record.REF):
                        alt_find = 1
                        prob_val = prob_val*sum([one_record.samples[i]["AD"][alt+1] for i in range(len(one_record.samples))])/one_record.INFO["DP"]
                        for i in range(len(v_alt)):
                            qual_prob[q_pos-self.pos+i] = 0.0
                        break
            if(alt_find==0):
                if(v_type=="W"):
                    base_qual = self.qual[q_pos-self.pos]
                    prob_val = prob_val*basequal_prob[str(ord(base_qual))]/3.0
                    qual_prob[q_pos-self.pos] = 0.0
                elif(v_type=="I"):
                    for i in range(1,len(v_alt)):
                        base_qual = self.qual[q_pos-self.pos+i]
                        prob_val = prob_val*basequal_prob[str(ord(base_qual))]/4.0
                        qual_prob[q_pos-self.pos+i] = 0.0
                elif(v_type=="D"):
                    prob_val = prob_val*pow(basequal_prob[str(delete_base_qual)],len(v_ref)-1)
        r_pos = self.pos-1
        while(True):
            r_pos = int(gt_record[self.chr][r_pos])
            if(r_pos>=self.r_pos_largest):
                break
            if(r_pos in self.variants_keys):
                continue
            one_record = gt_file.data[self.chr][r_pos].record
            prob_val = prob_val*sum([one_record.samples[i]["AD"][0] for i in range(len(one_record.samples))])/one_record.INFO["DP"]
            _,q_pos = self.get_correspond_q(r_pos)
            qual_prob[q_pos-self.pos] = 0.0
        prob_val = prob_val*np.prod(1-qual_prob)
        return prob_val          

    def get_q_r_pos(self):
        tmp_pos_q = self.pos
        tmp_pos_r = self.pos
        self.q_r_pos_list = []
        for tag,length in self.cigar_list:
            if(tag in cigar_op_r and tag in cigar_op_q):
                self.q_r_pos_list.append([tag,length,tmp_pos_q,tmp_pos_r])
                tmp_pos_q+=length
                tmp_pos_r+=length
            elif(tag in cigar_op_q):
                self.q_r_pos_list.append([tag,length,tmp_pos_q,tmp_pos_r-1])
                tmp_pos_q+=length
            elif(tag in cigar_op_r):
                self.q_r_pos_list.append([tag,length,tmp_pos_q-1,tmp_pos_r])
                tmp_pos_r+=length
            else:
                self.q_r_pos_list.append([tag,length,tmp_pos_q-1,tmp_pos_r-1])
        self.r_pos_largest = tmp_pos_r
        self.mismatch_ref_pos = []
        for read_pos,_,_,_ in self.mismatch_base:
            tmp_q_pos = self.pos + read_pos
            for tag,length,q_pos,r_pos in self.q_r_pos_list:
                if(tag in cigar_op_q_r and q_pos<=tmp_q_pos and q_pos + length>tmp_q_pos):
                    self.mismatch_ref_pos.append(tmp_q_pos-q_pos+r_pos)
    #we only consider insertion, deletion and snps
    def get_variants(self):
        self.variants = []
        self.variants_keys = []
        #mismatch
        for i in range(len(self.mismatch_ref_pos)):
            self.variants_keys.append(self.mismatch_ref_pos[i])
            self.variants.append(["W",self.pos+self.mismatch_base[i][0],self.mismatch_ref_pos[i],self.mismatch_base[i][1],self.mismatch_base[i][3]])
        #insertion and deletion
        tmp_md_point = 0
        for tag,length,pos_q,pos_r in self.q_r_pos_list:
            if(tag=='I'):
                if(pos_q==self.pos):
                    continue
                self.variants_keys.append(pos_r)
                if(pos_r not in self.mismatch_ref_pos):
                    self.variants.append(["I",pos_q-1,pos_r,self.seq[pos_q-self.pos-1:pos_q-self.pos+length],self.seq[pos_q-self.pos-1]])
                else:
                    self.variants.append(["I",pos_q-1,pos_r,self.seq[pos_q-self.pos-1:pos_q-self.pos+length],self.mismatch_base[self.mismatch_ref_pos.index(pos_r)][3]])
            elif(tag=='D'):
                if(pos_q==self.pos):
                    continue
                while(True):
                    if(self.tags["MD"][tmp_md_point]=="^"):
                        break
                    tmp_md_point+=1
                tmp_md_point+=1
                tmp_begin_point = tmp_md_point
                while(True):
                    if(tmp_md_point==len(self.tags["MD"])):
                        break
                    if(self.tags["MD"][tmp_md_point] in num_set):
                        break
                    tmp_md_point+=1
                self.variants_keys.append(pos_r-1)
                if(pos_r-1 not in self.mismatch_ref_pos):
                    self.variants.append(["D",pos_q,pos_r-1,self.seq[pos_q-self.pos],self.seq[pos_q-self.pos]+self.tags["MD"][tmp_begin_point:tmp_md_point]])
                else:
                    self.variants.append(["D",pos_q,pos_r-1,self.seq[pos_q-self.pos],self.mismatch_base[self.mismatch_ref_pos.index(pos_r-1)][3]+self.tags["MD"][tmp_begin_point:tmp_md_point]])
        self.variants_keys = set(self.variants_keys)
    def basic_op(self):
        self.decompose_cigar()
        self.check_cigar_decomp_correct()
        self.decompose_MD_tag()
        self.get_insert_base()
        self.get_mismatch_base()
        self.convert_qual_to_num()
        self.get_q_r_pos()
        self.get_variants()

class remap_reads(object):
    def __init__(self,reads_list):
        self.reads_list = reads_list
        self.reads_num = len(reads_list)
    def reads_list_to_info(self,pysam_form=0):
        self.reads_info = [read_info(read,None,pysam_form) for read in self.reads_list]
    def reads_map_prob(self,basequal_prob={},if_gt=0,gt_file={},gt_record={}):
        if(len(self.reads_info)==1):
            self.maps_prob = np.array([1],dtype=np.float64)
            return
        self.maps_prob = []
        for read in self.reads_info:
            read.basic_op()
            self.maps_prob.append(read.map_prob(basequal_prob=basequal_prob,if_gt=if_gt,gt_file=gt_file,gt_record=gt_record))
        if(sum(self.maps_prob)==0):
            self.maps_prob = np.array(self.maps_prob,dtype=np.float64)
            self.maps_prob[:] = 1/self.maps_prob.size
            return
            #print("attention: " + str(len(self.maps_prob)))
        self.maps_prob = np.array(self.maps_prob,dtype=np.float64)/sum(self.maps_prob)
        return 

class ref_genome(object):
    def __init__(self,file):
        self.ref_genome_set = {}
        for record in SeqIO.parse(file,"fasta"):
            self.ref_genome_set[record.id] = record
    def get_sub_string(self,chr_name,start,end):
        return str(self.ref_genome_set[chr_name].seq[start:end]).upper()

class basequal_to_prob_independent(object):
    def __init__(self,ori_qual_set={}):
        if(ori_qual_set=={}):
            self.qual_set = {i:[0,0] for i in range(33,127)}
        else:
            self.qual_set = ori_qual_set
        self.checkpoint = 0
    def update_qual_set(self,read):
        mismatch_dict = {mismatch_base[0]:mismatch_base[1:] for mismatch_base in read.mismatch_base}
        insert_dict = {insert_base[0]:insert_base[1:] for insert_base in read.insert_base}
        self.checkpoint+=1
        for i in range(len(read.qual)):
            if(i in mismatch_dict.keys()):
                if(mismatch_dict[i][1]!=read.qual[i]):
                    print("Error: quality does not match!")
                    exit()
                self.qual_set[ord(read.qual[i])][0]+=1
                self.qual_set[ord(read.qual[i])][1]+=1
            elif(i in insert_dict.keys()):
                if(insert_dict[i][1]!=read.qual[i]):
                    print("Error: quality does not match!")
                    exit()
                self.qual_set[ord(read.qual[i])][0]+=1
                self.qual_set[ord(read.qual[i])][1]+=1
            else:
                self.qual_set[ord(read.qual[i])][1]+=1
        if(self.checkpoint%100000==0):
            print("reads finish: %d" % self.checkpoint)
    #this function must be used after all the update quality set are done
    def get_qual_prob(self):
        self.prob_set = {}
        for i in self.qual_set.keys():
            #print(i,type(i))
            if(self.qual_set[i][1]==0):
                self.prob_set[i] = pow(10,-(int(i)-33)/10)
            else:
                self.prob_set[i] = self.qual_set[i][0]/self.qual_set[i][1]
    def write_to_file(self,filename):
        with open(filename,"w")as f:
            json.dump(self.qual_set,f)

class record_info(object):
    def __init__(self,record,back_chr='',back_loc=0,next_chr='',next_loc=0):
        self.back_chr = back_chr
        self.back_loc = back_loc
        self.next_chr = next_chr
        self.next_loc = next_loc
        self.record = record
    def update_back(self,back_chr,back_loc):
        self.back_chr = back_chr
        self.back_loc = back_loc
    def update_next(self,next_chr,next_loc):
        self.next_chr = next_chr
        self.next_loc = next_loc

class genotype_info(object):
    def __init__(self):
        self.size = 0
        self.data = {}
    def update_record(self,record):
        if(record.CHROM not in self.data.keys()):
            self.data[record.CHROM] = {record.POS:record_info(record)}
        else:
            self.data[record.CHROM][record.POS] = record_info(record)
        self.size+=1
    def update_back(self,chr_num,loc,back_chr,back_loc):
        self.data[chr_num][loc].update_back(back_chr,back_loc)
    def update_next(self,chr_num,loc,next_chr,next_loc):
        self.data[chr_num][loc].update_next(next_chr,next_loc)

def get_posterior_inter(reads_1,reads_2,dist_prob,region_size,if_max=0):
    inter_prob = np.zeros((reads_1.reads_num,reads_2.reads_num))
    for i in range(reads_1.reads_num):
        for j in range(reads_2.reads_num):
            dist_bet_two = abs(reads_1.reads_info[i].pos-reads_2.reads_info[j].pos)
            #print(dist_bet_two)
            if(reads_1.reads_info[i].chr==reads_2.reads_info[j].chr):
                inter_prob[i,j] = reads_1.maps_prob[i]*reads_2.maps_prob[j]*get_dist_prob(dist_bet_two,dist_prob,region_size)
    if(np.sum(inter_prob)==0):
        return (reads_1.reads_list[0],reads_2.reads_list[0])
    inter_prob = inter_prob/np.sum(inter_prob)
    #print(inter_prob)
    index_1,index_2 = choose_inter(inter_prob,if_max)
    #print(index_1,index_2)
    return (reads_1.reads_list[index_1],reads_2.reads_list[index_2])

def get_dist_prob(dist_val,dist_prob,region_size):
    #print(dist_prob[int(dist_val/region_size)])
    if(int(dist_val/region_size)>len(dist_prob)-1):
        return dist_prob[-1]
    else:
        return dist_prob[int(dist_val/region_size)]

def choose_inter(inter_prob,if_max=0):
    pos = -1
    prob_sum = 0
    inter_prob_flat = inter_prob.flatten()
    if(if_max==1):
        pos = inter_prob_flat.argmax()
        return (int(pos/inter_prob.shape[1]),pos%inter_prob.shape[1])
    rnd_num = random.uniform(0,1)
    #print(rnd_num)
    for i in range(len(inter_prob_flat)):
        if(prob_sum+inter_prob_flat[i]>=rnd_num):
            pos = i
            break
        else:
            prob_sum+=inter_prob_flat[i]
    if(pos==-1):
        pos = len(inter_prob_flat)-1
    return (int(pos/inter_prob.shape[1]),pos%inter_prob.shape[1])

def formalize_read_pair(pair_list,name):
    read_1,read_2 = pair_list
    new_read_1 = read_1.copy()
    new_read_2 = read_2.copy()
    new_read_1[0] = name
    new_read_2[0] = name
    new_read_1[7] = new_read_2[3]
    new_read_2[7] = new_read_1[3]
    if(new_read_1[2]==new_read_2[2]):
        new_read_1[6] = '='
        new_read_2[6] = '='
    else:
        new_read_1[6] = new_read_2[2]
        new_read_2[6] = new_read_1[2]
    return (new_read_1,new_read_2)

def Next_iter(iter_list):
    try:
        return next(iter_list)
    except StopIteration:
        return 0

def get_pair_seq(seq):
    return ('').join([base_pair[base] for base in seq])

def reverse_seq(seq):
    return get_pair_seq(seq)[::-1]

'''
def get_region_prob(dist_list,region_size=50000):
    t0 = time.time()
    reg_num = int(np.max(dist_list)/region_size)
    freq_list = [len(np.intersect1d(np.where(dist_list>=i*region_size)[0],np.where(dist_list<(i+1)*region_size)[0])) for i in range(reg_num+1)]
    print(time.time()-t0)
    return (np.array(freq_list),np.array(freq_list)/sum(freq_list))
'''

def generate_region_prob(dist_list,region_size=50000):
    #t0 = time.time()
    dist_list.sort()
    freq_list = []
    region_pos = 1
    region_start = 0
    dist_list_len = len(dist_list)
    i=0
    while(i<dist_list_len):
        if(dist_list[i]>=region_pos*region_size):
            freq_list.append(i-region_start)
            region_pos+=1
            region_start=i
        else:
            i+=1
    freq_list.append(dist_list_len-region_start)
    #print(time.time()-t0)
    return (np.array(freq_list,dtype=np.float64),np.array(freq_list,dtype=np.float64)/sum(freq_list))

def get_read_name(read):
    name = read.qname
    return name.split("/",1)[0]

def read_json_file(file_name):
    with open(file_name)as f:
        data = json.load(f)
    return data

def write_json_file(data,file_name):
    with open(file_name,"w")as f:
        json.dump(data,f)
    return

def new_map_read_processing(read_pair):
    new_pair = []
    for read in read_pair:
        if((read.flag & 0x100)==0x100):
            read.flag = read.flag-256
        if(read.mapping_quality==255):
            #print(read.to_string())
            read.mapping_quality = 0
        if(read.has_tag('XS')):
            if(read.get_tag("AS")<=read.get_tag("XS")):
                read.set_tags([tag for tag in read.tags if tag[0]!="XS"])
        new_pair.append(read)
    return new_pair

def mapping_back(input_path,input_file,output_path,align_num=2):
    file_name = input_file.split(".")[0]
    cmd = 'samtools bam2fq ' + input_path + input_file + ' > ' + output_path + file_name + '.fastq'
    os.system(cmd)
    cmd = 'bowtie2 --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder --un '+output_path+file_name+\
    '.bwt2glob.unmap.fastq --rg-id BMG --rg SM:merge --phred33-quals -p 2 -k ' + str(align_num) + ' -x /mnt/disk53/user/nsauerwa/HiCPro_stuff/bt2idx_hg19/hg19 -U '+\
    output_path+file_name+'.fastq 2> '+output_path+'bowtie_'+file_name+'_global.log | samtools view -F 4 -bS - > '+output_path+file_name+'.bwt2glob.bam'
    os.system(cmd)
    return

def rm_file(file_name):
    if(os.path.exists(file_name)):
        cmd = 'rm ' + file_name
        os.system(cmd)
    return

def find_all_files(path):
    file_list = [f for f in listdir(path) if isfile(join(path,f))]
    return file_list

def create_new_dir(path_name):
    path_list = path_name.split("/")
    tmp_path = ""
    for tmp_dir in path_list:
        tmp_path = tmp_path + tmp_dir + "/"
        if not os.path.exists(tmp_path):
            os.makedirs(tmp_path)
    return 

def read_chr_len(file_name):
    chr_len_dict = {}
    with open(file_name,"rt")as f:
        freader = csv.reader(f,delimiter='\t')
        for line in freader:
            chr_len_dict[line[0]] = int(line[1])
    return chr_len_dict


def get_vcf_info(file_name,chr_len_dict):
    vcf_reader = vcf.Reader(filename=file_name)
    back_chr = 'chr1'
    back_loc = 0
    count = 0
    vcf_record_dict = {}
    for chr_name in chr_len_dict.keys():
        vcf_record_dict[chr_name] = np.zeros(chr_len_dict[chr_name]+1)
        vcf_record_dict[chr_name][:] = chr_len_dict[chr_name]+1
    vcf_info = genotype_info()
    for record in vcf_reader:
        '''
        if(record.INFO["DP"]<10 or record.FILTER!=[]):
            continue
        '''
        if(record.FILTER!=[]):
            continue
        if(record.CHROM!=back_chr):
            #vcf_record_dict[back_chr][back_loc:] = chr_len_dict[back_chr]+1
            vcf_record_dict[record.CHROM][:record.POS] = record.POS
        else:
            vcf_record_dict[record.CHROM][back_loc:record.POS] = record.POS
        vcf_info.update_record(record)
        vcf_info.update_back(record.CHROM,record.POS,back_chr,back_loc)
        if(count>0):
            vcf_info.update_next(back_chr,back_loc,record.CHROM,record.POS)
        back_chr = record.CHROM
        back_loc = record.POS
        count+=1
        #print(count)
        if(count % 10000==0):
            print("## " + str(count))
    return (vcf_info,vcf_record_dict)

def if_unique(read):
    if not read.is_unmapped and read.has_tag('AS'):
        if(read.has_tag('XS')):
            if(read.get_tag('AS')>read.get_tag('XS')):
                return True
        else:
            return True
    return False