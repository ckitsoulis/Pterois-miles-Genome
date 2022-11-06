#%%
import pandas as pd
from Bio import SeqIO
from datetime import datetime
import sys
import os

'''
This scripts produces the Paralog-sequences HOG files (FASTA).

inputs:
1) PATH to file created by "single_orthoX.py" --> HOGs_ManyCopyOrth.tsv
2) PATH to proteome files (FASTA)

outputs:
1) Directory with Paralog-sequences HOGs (FASTA)
'''

file_PATH = "/home/christos/Documents/master_bioinfo/Thesis_PT/Ortho_Pmiles/HOGs_ManyCopyOrth.tsv"
path2fasta = ""

#os.mkdir("ParalogOrthogroups")

#=============================================================================================================
HPGs = pd.read_csv(file_PATH, sep = "\t", header = 0, index_col = None, low_memory = False)
HPGs = HPGs.fillna("")

print("[{}] Import paralog HOGs ===> DONE".format(datetime.now()))
#=============================================================================================================
#%%
HPGs
#%%
#=============================================================================================================
print("Start selecting paralogs and writing to HOG files (FASTA):")
# Parse a fasta file and return a dictionary with headers as keys and sequences as values
def parse_fasta(file_path):
    HEADERS = []
    SEQUENCES = []
    for record in SeqIO.parse(file_path, 'fasta'):
        HEADERS.append(record.description)
        SEQUENCES.append(str(record.seq))
    
    return dict(zip(HEADERS, SEQUENCES))


def SequencesPerSpecies(a_dict, species):

    species_name = species.replace("_","")

    species_file = parse_fasta(path2fasta+"/"+species+".fasta")

    for HOG, Sheader in list(a_dict.items()):

        if len(Sheader) == 0:
            continue
        
        else:
            
            with open("ParalogOrthogroups"+"/"+HOG+".fasta", "a+") as HOG_file:
                
                if "," in Sheader:
                    for seq in Sheader.split(", "):
                        gene_name = seq.replace("_","")
                
                        HOG_file.write(">"+species_name+"_"+gene_name+"\n")
                        HOG_file.write(species_file[seq]+"\n")
                
                else:
                    gene_name = Sheader.replace("_","")

                    HOG_file.write(">"+species_name+"_"+gene_name+"\n")
                    HOG_file.write(species_file[Sheader]+"\n")

Species = HPGs.columns[1:].to_list()

for sp in Species:
    species_dict = dict(zip(HPGs["HOG"], HPGs[sp]))

    print("[{}] Selecting sequences for: {}".format(datetime.now(),sp))
    SequencesPerSpecies(species_dict, species = sp)
#=============================================================================================================

#%%
import pandas as pd

c = pd.read_csv("/home/christos/Documents/master_bioinfo/Thesis_PT/Ortho_Pmiles/HOG_counts.tsv", sep="\t", header=0, index_col=0)
c
#%%
HOGs = pd.read_csv("/home/christos/Documents/master_bioinfo/Thesis_PT/Ortho_Pmiles/N0.tsv", header = 0, sep = "\t", index_col = None, low_memory=False)
HOGs = HOGs.drop(['OG', 'Gene Tree Parent Clade'], axis = 1)
HOGs["HOG"] = HOGs["HOG"].str.replace("N0.","")
HOGs

#%%
HOGs[HOGs["HOG"] == "HOG0028402"]