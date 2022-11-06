import imp
from Bio import SeqIO
from datetime import datetime
import os

alignments_PATH = ""
path2fasta = ""

os.mkdir(alignments_PATH+"Corrected_MultipleAlignments/")

Species_global = [fasta.replace(".fasta", "") for fasta in os.listdir(path2fasta) if fasta.endswith(".fasta")]


def parse_fasta(file_path):
    HEADERS = []
    SEQUENCES = []
    for record in SeqIO.parse(file_path, 'fasta'):
        HEADERS.append((record.description).split("_")[0])
        SEQUENCES.append(str(record.seq))
    
    return dict(zip(HEADERS, SEQUENCES))


for aln_file in os.listdir(alignments_PATH):

    if aln_file.endswith(".fasta"):
        print("[{}] Checking file {}".format(datetime.now(), aln_file))
    
        Species_local = dict(zip(map(lambda x: x.replace("_",""),Species_global), [""]*len(Species_global)))

        aln_dict = parse_fasta(alignments_PATH+aln_file)

        if len(aln_dict.keys()) == len(Species_global):
            print("[{}] File {} ===> OK".format(datetime.now(), aln_file))
            
            with open(alignments_PATH+"Corrected_MultipleAlignments/{}".format(aln_file), "a+") as new_file:
                for header, seq in aln_dict.items():
                    new_file.write(">"+header+"\n")
                    new_file.write(seq+"\n")
            
            continue

        else:
            
            print("[{}] Correcting file {}".format(datetime.now(), aln_file))
            average_length = int(sum( map(len, aln_dict.values()) ) / len(aln_dict.values()))
            
            for sp in Species_local.keys():
                
                if sp in aln_dict.keys():
                    Species_local[sp] = aln_dict[sp]
                
                else:
                    Species_local[sp] = "-"*average_length

            with open(alignments_PATH+"Corrected_MultipleAlignments/{}".format(aln_file), "a+") as new_file:
                
                for index, seq in enumerate(Species_local.values()):

                    new_file.write(">"+Species_global[index]+"\n")
                    new_file.write(seq+"\n")
        
    else:
        continue


