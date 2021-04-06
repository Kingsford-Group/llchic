from model_repeat_remap import *
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',type=str,help='the path of input files')
    parser.add_argument('-r1',type=str,help='the first input file')
    parser.add_argument('-r2',type=str,help='the second input file')
    parser.add_argument('-o',type=str,help='the output path')
    parser.add_argument('-rs',type=int,help='the size of each region in distance distribution')
    parser.add_argument('-k',type=int,help='report up to <int> aligns per read')
    parser.add_argument('-m',type=int,help='decide if we pick the pair with largest posterior probability')
    parser.add_argument('-gf',type=str,default='',help='the genotype info file')
    args=parser.parse_args()
    return args

args = get_args()
#print(args.i,args.r1,args.r2,args.o)

SRA_id = args.r1.split(".")[0].split("_")[0]

R1file = args.i + args.r1
R2file = args.i + args.r2

multi_align_file_1 = args.o + args.r1.split(".")[0] + "_multialign.bam"
multi_align_file_2 = args.o + args.r2.split(".")[0] + "_multialign.bam"

multi_read_file_1 = args.o + args.r1.split(".")[0] + "_multireads.bam"
multi_read_file_2 = args.o + args.r2.split(".")[0] + "_multireads.bam"

output_basequal = args.o + SRA_id + "_basequal_set.txt"

if_comb_gt = (args.gf!='')

time.sleep(30)
create_new_dir(args.o)

ori_qual_file = None

if(ori_qual_file!=None):
    basequal_prob = basequal_to_prob_independent(read_json_file(ori_qual_file))
else:
    basequal_prob = basequal_to_prob_independent()

reads_counter = 0
print("extract multi aligned reads...")
t0 = time.time()
dist_list = []
with pysam.AlignmentFile(R1file,"rb")as hr1, pysam.AlignmentFile(R2file,"rb")as hr2:
    output_1 = pysam.AlignmentFile(multi_align_file_1,"wb",template=hr1)
    output_2 = pysam.AlignmentFile(multi_align_file_2,"wb",template=hr2)
    output_multi_read_1 = pysam.AlignmentFile(multi_read_file_1,"wb",template=hr1)
    output_multi_read_2 = pysam.AlignmentFile(multi_read_file_2,"wb",template=hr2)
    for r1,r2 in zip(hr1.fetch(until_eof=True), hr2.fetch(until_eof=True)):
        reads_counter+=1
        if(reads_counter % 1000000 == 0):
            print("##",reads_counter)
        if(get_read_name(r1)!=get_read_name(r2)):
            print("Forward and reverse reads not paired. Check that BAM files have the same read names and are sorted.")
            sys.exit(1)
        if(r1.is_unmapped == True or r2.is_unmapped == True or if_not_unmap(r1) == False or if_not_unmap(r2) == False):
            continue
        #if(r1.is_unmapped == True or r2.is_unmapped == True or r1.reference_name!=r2.reference_name or if_not_unmap(r1) == False or if_not_unmap(r2) == False):
            #continue
        if(r1.has_tag('XS') or r2.has_tag('XS')):
            output_1.write(r1)
            output_2.write(r2)
            if(r1.has_tag('XS')):
                output_multi_read_1.write(r1)
            if(r2.has_tag('XS')):
                output_multi_read_2.write(r2)
        else:
            #unqiue mapping update base quality set, update distance list
            read_1 = read_info(r1,None,pysam_form=1)
            read_2 = read_info(r2,None,pysam_form=1)
            read_1.basic_op()
            read_2.basic_op()
            basequal_prob.update_qual_set(read_1)
            basequal_prob.update_qual_set(read_2)
            dist_list.append(abs(r1.reference_start-r2.reference_start))

hr1.close()
hr2.close()
output_1.close()
output_2.close()
output_multi_read_1.close()
output_multi_read_2.close()

basequal_prob.write_to_file(output_basequal)
np.save(args.o + SRA_id + "_dist_list.npy",np.array(dist_list))
dist_freq,dist_prob = generate_region_prob(np.array(dist_list),region_size=args.rs)
np.save(args.o + SRA_id + "_dist_freq.npy",dist_freq)
np.save(args.o + SRA_id + "_dist_prob.npy",dist_prob)
print("extract multi aligned reads finish, time duration: " + str(time.time()-t0))



print("mapping back...")
t0 = time.time()
mapping_back(args.o,args.r1.split(".")[0] + "_multireads.bam",args.o,args.k)
mapping_back(args.o,args.r2.split(".")[0] + "_multireads.bam",args.o,args.k)
print("mapping back finish, tinme duration: " + str(time.time()-t0))

'''

print("mapping back maize...")
t0 = time.time()
mapping_back_maize(args.o,args.r1.split(".")[0] + "_multireads.bam",args.o,args.k)
mapping_back_maize(args.o,args.r2.split(".")[0] + "_multireads.bam",args.o,args.k)
print("mapping back finish, tinme duration: " + str(time.time()-t0))



remap_file_1 = multi_read_file_1.split(".")[0] + ".bwt2glob.bam"
remap_file_2 = multi_read_file_2.split(".")[0] + ".bwt2glob.bam"

new_multi_align_file_1 = args.o + args.r1.split(".")[0] + "_new_multialign.bam"
new_multi_align_file_2 = args.o + args.r2.split(".")[0] + "_new_multialign.bam"


dist_prob = np.load(args.o + SRA_id + "_dist_prob.npy")
basequal_prob = basequal_to_prob_independent(read_json_file(output_basequal))
basequal_prob.get_qual_prob()

gt_info = {}
gt_record_info = {}
if(if_comb_gt):
    gt_info,gt_record_info = get_vcf_info(args.gf,read_chr_len("/mnt/disk75/user/yihangs/data/commen_data/NCBI_hg19_GRCh37_chr_len.txt"))

print("calculate posterior probability and re-choose read pairs...")
t0 = time.time()

counter = 0

rp1 = pysam.AlignmentFile(remap_file_1,"rb")
rp2 = pysam.AlignmentFile(remap_file_2,"rb")
next_read_1 = Next_iter(rp1.fetch(until_eof=True))
next_read_2 = Next_iter(rp2.fetch(until_eof=True))
with pysam.AlignmentFile(multi_align_file_1,"rb")as hr1, pysam.AlignmentFile(multi_align_file_2,"rb")as hr2:
    output_1 = pysam.AlignmentFile(new_multi_align_file_1,"wb",template=hr1)
    output_2 = pysam.AlignmentFile(new_multi_align_file_2,"wb",template=hr2)
    for r1,r2 in zip(hr1.fetch(until_eof=True), hr2.fetch(until_eof=True)):
        #pdb.set_trace()
        counter+=1
        #print(counter)
        reads_list_1 = [copy.deepcopy(r1)]
        reads_list_2 = [copy.deepcopy(r2)]
        if(get_read_name(r1)!=get_read_name(r2)):
            print("Forward and reverse reads not paired. Check that BAM files have the same read names and are sorted.")
            sys.exit(1)
        if(r1.has_tag('XS')):
            while(1):
                if(type(next_read_1)!=type(r1)):
                    break
                if(get_read_name(next_read_1)!=get_read_name(r1)):
                    break
                if(next_read_1.reference_start==r1.reference_start):
                    next_read_1 = Next_iter(rp1.fetch(until_eof=True))
                    continue
                reads_list_1.append(copy.deepcopy(next_read_1))
                next_read_1 = Next_iter(rp1.fetch(until_eof=True))
        #pdb.set_trace()
        if(r2.has_tag("XS")):
            while(1):
                if(type(next_read_2)!=type(r2)):
                    break
                if(get_read_name(next_read_2)!=get_read_name(r2)):
                    break
                if(next_read_2.reference_start==r2.reference_start):
                    next_read_2 = Next_iter(rp2.fetch(until_eof=True))
                    continue
                reads_list_2.append(copy.deepcopy(next_read_2))
                next_read_2 = Next_iter(rp2.fetch(until_eof=True))
        #pdb.set_trace()
        if(len(reads_list_1)==1 and len(reads_list_2)==1):
            output_1.write(r1)
            output_2.write(r2)
            continue
        reads_1 = remap_reads(reads_list_1)
        reads_1.reads_list_to_info(pysam_form=1)
        reads_1.reads_map_prob(basequal_prob=basequal_prob.prob_set,if_gt=if_comb_gt,gt_file=gt_info,gt_record=gt_record_info)
        reads_2 = remap_reads(reads_list_2)
        reads_2.reads_list_to_info(pysam_form=1)
        reads_2.reads_map_prob(basequal_prob=basequal_prob.prob_set,if_gt=if_comb_gt,gt_file=gt_info,gt_record=gt_record_info)
        inter_pair = get_posterior_inter(reads_1,reads_2,dist_prob,args.rs,if_max=args.m)
        inter_pair = new_map_read_processing(inter_pair)
        #pdb.set_trace()
        if(inter_pair[0].reference_id==r1.reference_id and inter_pair[0].reference_start==r1.reference_start and inter_pair[1].reference_id==r2.reference_id and inter_pair[1].reference_start==r2.reference_start):
            output_1.write(r1)
            output_2.write(r2)
        else:
            output_1.write(inter_pair[0])
            output_2.write(inter_pair[1])
        if(counter % 100000 == 0):
            print("##",counter)

hr1.close()
hr2.close()
rp1.close()
rp2.close()
output_1.close()
output_2.close()

print("calculate posterior probability and re-choose read pairs finish, time duration: " + str(time.time()-t0))

output_1_file = args.o + args.r1
output_2_file = args.o + args.r2

print("create new merged bam file...")
t0 = time.time()

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
        if(type(new_read_1)==type(r1) and get_read_name(new_read_1)==get_read_name(r1)):
            output_1.write(new_read_1)
            new_read_1 = Next_iter(nr1.fetch(until_eof=True))
        else:
            output_1.write(r1)
        if(type(new_read_2)==type(r2) and get_read_name(new_read_2)==get_read_name(r2)):
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

print("create new merged bam file finish, time duration: " + str(time.time()-t0))

'''