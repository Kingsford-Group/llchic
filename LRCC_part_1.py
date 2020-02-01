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
        if(r1.is_unmapped == True or r2.is_unmapped == True or r1.reference_name!=r2.reference_name or if_unique(r1) == False or if_unique(r2) == False):
            continue
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
dist_freq,dist_prob = generate_region_prob(np.array(dist_list),region_size=args.rs)
np.save(args.o + SRA_id + "_dist_list.npy",np.array(dist_list)) 
np.save(args.o + SRA_id + "_dist_freq.npy",dist_freq)
np.save(args.o + SRA_id + "_dist_prob.npy",dist_prob)
print("extract multi aligned reads finish, time duration: " + str(time.time()-t0))



print("mapping back...")
t0 = time.time()
mapping_back(args.o,args.r1.split(".")[0] + "_multireads.bam",args.o,args.k)
mapping_back(args.o,args.r2.split(".")[0] + "_multireads.bam",args.o,args.k)
print("mapping back finish, tinme duration: " + str(time.time()-t0))

