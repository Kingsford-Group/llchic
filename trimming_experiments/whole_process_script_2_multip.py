from model_repeat_remap import *
import argparse
import multiprocessing
import glob

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',type=str,help='the path of input files')
    parser.add_argument('-r1',type=str,help='the first input file')
    parser.add_argument('-r2',type=str,help='the second input file')
    parser.add_argument('-o1',type=str,help='the output path for intermediate results')
    parser.add_argument('-o2',type=str,help='the output path for final results')
    parser.add_argument('-rs',type=int,help='the size of each region in distance distribution')
    #parser.add_argument('-k',type=int,help='report up to <int> aligns per read')
    parser.add_argument('-m',type=int,help='decide if we pick the pair with largest posterior probability')
    parser.add_argument('-gf',type=str,default='',help='the genotype info file')
    args=parser.parse_args()
    return args


args = get_args()

R1file_list = [file_name.split('/')[-1] for file_name in glob.glob(args.i + args.r1)]
R1file_list.sort()
R2file_list = [file_name.split('/')[-1] for file_name in glob.glob(args.i + args.r2)]
R2file_list.sort()
print(R1file_list)
print(R2file_list)

if_comb_gt = (args.gf!='')

gt_info = {}
gt_record_info = {}
if(if_comb_gt):
    gt_info,gt_record_info = get_vcf_info(args.gf,read_chr_len("/mnt/disk75/user/yihangs/data/commen_data/NCBI_hg19_GRCh37_chr_len.txt"))

#pdb.set_trace()

def sub_process(process_number):
    global args
    global gt_info
    global if_comb_gt
    global gt_record_info
    global R1file_list
    global R2file_list
    r1_file = R1file_list[process_number]
    R1file = args.i + r1_file
    print(R1file)
    r2_file = R2file_list[process_number]
    R2file = args.i + r2_file
    print(R2file)
    SRA_id = r1_file.split(".")[0].split("_")[0]
    if(r2_file.split(".")[0].split("_")[0]!=SRA_id):
        print("Error: Two input files are not compatible!")
        sys.exit(1)

    multi_align_file_1 = args.o1 + r1_file.split(".")[0] + "_multialign.bam"
    multi_align_file_2 = args.o1 + r2_file.split(".")[0] + "_multialign.bam"

    multi_read_file_1 = args.o1 + r1_file.split(".")[0] + "_multireads.bam"
    multi_read_file_2 = args.o1 + r2_file.split(".")[0] + "_multireads.bam"

    output_basequal = args.o1 + SRA_id + "_basequal_set.txt"   

    remap_file_1 = multi_read_file_1.split(".")[0] + ".bwt2glob.bam"
    remap_file_2 = multi_read_file_2.split(".")[0] + ".bwt2glob.bam"

    new_multi_align_file_1 = args.o2 + r1_file.split(".")[0] + "_new_multialign.bam"
    new_multi_align_file_2 = args.o2 + r2_file.split(".")[0] + "_new_multialign.bam"

    dist_prob = np.load(args.o1 + SRA_id + "_dist_prob.npy")
    basequal_prob = basequal_to_prob_independent(read_json_file(output_basequal))
    basequal_prob.get_qual_prob()

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
            #reads_list_1 = []
            #reads_list_2 = []
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
            #if(len(reads_list_1)==0 and len(reads_list_2)==0):
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

    output_1_file = args.o2 + r1_file
    output_2_file = args.o2 + r2_file

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

def main_process():
    global R1file_list
    
    process_list = []
    for i in range(len(R1file_list)):
        tmp_process = multiprocessing.Process(target=sub_process,args=(i,))
        process_list.append(tmp_process)
    
    for process in process_list:
        process.start()
    for process in process_list:
        process.join()

if __name__ == "__main__":
    #print(gt_info.size)
    main_process()