#%%
import os
import subprocess
import sys
import multiprocessing as mp

working_path = sys.argv[-2]
softmasked_genome = sys.argv[-1]


def Exonerate_Run(fasta_file):

    file_name = fasta_file.split("fasta_chunk_")
    prefix_name, number_name = file_name[0], file_name[1]
    output = prefix_name + "_Exo_" + number_name + ".gff"

    subprocess.run("echo 'Job for file {}'".format(number_name), shell=True, stdout=subprocess.DEVNULL)
    subprocess.run("exonerate --model protein2genome --bestn 3 --showtargetgff T --showvulgar F --maxintron 500000 --fsmmemory 16000 --softmasktarget T {} {} > {}".format(fasta_file, softmasked_genome,output), shell=True, stdout=subprocess.DEVNULL)


proteome_files = [file for file in os.listdir(working_path) if ".fasta_chunk_" in file]


print("NUMBER OF CORES AVAILABLE: {}".format(mp.cpu_count()))


pool = mp.Pool(mp.cpu_count())

# map pool to function
results = pool.map(Exonerate_Run, [proteome_file for proteome_file in proteome_files][:32])

pool.close()