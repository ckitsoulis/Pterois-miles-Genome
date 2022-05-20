from Bio import SeqIO
import os

unknown_TE = "TE_libraries/Unknown_EDTA.fasta"
unordered_TE = "TE_libraries/UnOrdered_EDTA.fasta"
classified_TE = "TE_libraries/Classified_EDTA.fasta"

def parse_fasta(file_path):
    HEADERS = []
    SEQUENCES = []
    for record in SeqIO.parse(file_path, 'fasta'):
        HEADERS.append(">"+ record.description)
        SEQUENCES.append(str(record.seq))
    
    return HEADERS, SEQUENCES

Headers, Sequences = parse_fasta("TE_libraries/lionfish_EDTA_TEfamilies.lib")

#Split EDTA library to Classified TEs, Unordered TEs and Unclassified TEs

with open(unknown_TE, "w") as file1, open(classified_TE, "w") as file2, open(unordered_TE, "w") as file3:

    for number, header in enumerate(Headers):

        if "#Unknown" in header:
            file1.write(header + "\n")
            file1.write(Sequences[number] + "\n")
        
        elif "/unknown" in header:
            file3.write(header + "\n")
            file3.write(Sequences[number] + "\n")

        else:
            file2.write(header + "\n")
            file2.write(Sequences[number] + "\n")

if os.path.getsize("TE_libraries/Unknown_EDTA.fasta") > 0:
    print(".... Unknown_EDTA.fasta is written ....")
else:
    print("Something went wrong with File 1!")

if os.path.getsize("TE_libraries/Classified_EDTA.fasta") > 0:
    print(".... Classified_EDTA.fasta is written ....")
else:
    print("Something went wrong with File 2!")

if os.path.getsize("TE_libraries/UnOrdered_EDTA.fasta") > 0:
    print(".... UnOrdered_EDTA.fasta is written ....")
else:
    print("Something went wrong with File 3!")