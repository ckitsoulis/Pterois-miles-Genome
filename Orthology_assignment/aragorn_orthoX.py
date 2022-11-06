#%%
import pandas as pd
import numpy as np
from Bio import SeqIO
import sys
import os

'''
This script counts the number of genes in every species per Hierarchical Orthogroup (HOG), as produced by OrthoFinder in N0.tsv file, and
creates the HOG files (FASTA) for Single-Copy-Orthologs with a user-defined representation. The second group of files can be used for multiple
alignment to create later a superalignment and phylogenomic tree inferrence.

inputs:
1) full PATH to N0.tsv file
2) full PATH to proteome files (FASTA)

outputs:
1) sub-table that contains only the paralog-sequences HOGs
2) count matrix for use by CAFE
3) Directory with Single-Copy-Sequences HOGs (FASTA)
'''

#working_PATH = sys.argv[-1]
#specimen = sys.argv[-2]
#path2fasta = sys.argv[-3]

working_PATH ="/home/christos/Documents/master_bioinfo/Thesis_PT/Ortho_Pmiles/"
path2fasta = ""

#os.mkdir("SingleCopyOrthogroups")

#=============================================================================================================
# Open HOGs file, drop 2nd & 3rd column
HOGs = pd.read_csv(working_PATH+"N0_clean.tsv", header = 0, sep = "\t", index_col = None, low_memory=False)
HOGs = HOGs.drop(['OG', 'Gene Tree Parent Clade'], axis = 1)
HOGs["HOG"] = HOGs["HOG"].str.replace("N0.","")
HOGs

print("Import N0.tsv ===> DONE", flush=True)
#=============================================================================================================

#=============================================================================================================
# Count the number of genes in each HOG per species
def GeneCounts (x, Df):

    tmp = np.array(Df[x].str.count(","))
    masked_array = (tmp >= 0)

    tmp[masked_array] += 1
    
    tmp2 = np.nan_to_num(tmp)

    Df[x] = tmp2

# Drop rows (standing for HOGs) where a selected species is not represented 
def SelectSpecies(df, specimen):
    return df.dropna(axis=0, subset=[specimen]).reset_index(drop = True)


HOGs_SpeciesSpecific = SelectSpecies(HOGs, "Pterois_miles")

HOGs_SpeciesSpecific
print("Make species-specific dataframe ===> DONE")
#=============================================================================================================

#=============================================================================================================
# List the name of species
Species = HOGs_SpeciesSpecific.columns[1:].to_list()

# Copy dataframes with genes to produce count matrices later on
HOGs_count = HOGs.copy()
HOGs_count_SpeciesSpecific = HOGs_SpeciesSpecific.copy()

for sp in Species:

    GeneCounts(sp, HOGs_count_SpeciesSpecific)
    GeneCounts(sp, HOGs_count)

HOGs_count.to_csv("HOG_counts.tsv", sep="\t", header=True, index=False)
print("Save count-matrix as 'HOG_counts.tsv' ===> DONE")
#=============================================================================================================

#=============================================================================================================
#Keep the rows that have at most one representation per species in each HOG
id = np.where((np.array(HOGs_count_SpeciesSpecific[Species]) <= 1).all(axis=1))[0]

s = HOGs_count_SpeciesSpecific[Species].loc[id].apply(sum, axis = 1)

threshold_HOGs = np.where(s >= 43)[0]
len(threshold_HOGs) # number of Single-Copy-Ortholog HOGs with a representation for at least 90% of total species

masked_ids = id[threshold_HOGs]

# Dataframe with Single-Copy-Orthologs per HOG
SingleCopyOrth = HOGs_SpeciesSpecific.loc[masked_ids].fillna("")
SingleCopyOrth

ManyCopyOrth = HOGs.loc[HOGs["HOG"].isin(list(SingleCopyOrth["HOG"]))==False].fillna("")
ManyCopyOrth.to_csv("HOGs_ManyCopyOrth.tsv", sep = "\t", header = True, index = False)
print("Save dataframe with multiple genes per species for each HOG ===> DONE")
#=============================================================================================================

#=============================================================================================================
print("Start selecting single-copy-orthologs and writing to HOG files (FASTA):")
# Parse a fasta file and return a dictionary with headers as keys and sequences as values
def parse_fasta(file_path):
    HEADERS = []
    SEQUENCES = []
    for record in SeqIO.parse(file_path, 'fasta'):
        HEADERS.append(record.description)
        SEQUENCES.append(str(record.seq))
    
    return dict(zip(HEADERS, SEQUENCES))


# Find the requested sequences in a fasta file and write them in the appropriate HOG file
def SequencesPerSpecies(a_dict, species):

    species_file = parse_fasta(path2fasta+"/"+species+".fasta")

    for HOG, Sheader in list(a_dict.items()):

        if len(Sheader) > 0:

            species_name = species.replace("_","")
            gene_name = Sheader.replace("_","")

            with open("SingleCopyOrthogroups"+"/"+HOG+".fasta", "a+") as HOG_file:
                HOG_file.write(">"+species_name+"_"+gene_name+"\n")
                HOG_file.write(species_file[Sheader]+"\n")
        
        else:
            pass

for sp in Species:
    species_dict = dict(zip(SingleCopyOrth["HOG"], SingleCopyOrth[sp]))

    print("Selecting sequences for: {}".format(sp))
    SequencesPerSpecies(species_dict, species = sp)
#=============================================================================================================

#%%
import pandas as pd
import numpy as np


working_PATH ="/home/christos/Documents/master_bioinfo/Thesis_PT/Ortho_Pmiles/"

#HOGs = pd.read_csv(working_PATH+"HOG_counts.tsv", header = 0, sep = "\t", index_col = 0, low_memory=False)
HOGs = pd.read_csv(working_PATH+"Vpap/N0_hogCounts_BEST-10_filtTaxon-19_clean_READY.tsv", header = 0, sep = "\t", index_col = 0, low_memory=False)
HOGs
#%%

sp = list(HOGs.columns)[1:]

HOGs[sp].max(axis=1) - HOGs[sp].min(axis=1)
#%%

tmp = HOGs.loc[(HOGs == 0).astype(int).sum(axis=1) < 11]

#%%
tmp.to_csv("HOG_counts_filt10.tsv", sep = "\t", header = True, index = True)
