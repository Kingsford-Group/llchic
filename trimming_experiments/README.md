


This repo contains the python souce code to run the method described in "Probabilistic method corrects previously uncharacterized Hi-C artifact" (Shen and Kingsford), which characterizes a previously uncharacterized Hi-C artifact and then design a probabilistic method to correct it. The method is combined with a standard Hi-C pipeline, HiC-Pro. The inputs are BAM files from mapping step of HiC-Pro, and the correction module (called Long Range Contact Correction), the outputs are new BAM files which can be used for Hi-C filtering step (proc_hic) of HiC-Pro. 


