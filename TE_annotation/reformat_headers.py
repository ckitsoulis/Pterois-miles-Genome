#%%
from library_split import parse_fasta
import os

reformation_dict = {"DTT": "TcMar-Mariner", "DTA": "hAT-hobo", 
"DTM": "Mutator", "DTE": "Merlin", "DTR": "Transib", "DTP": "P", 
"DTB": "PiggyBac", "DTH": "PIF-Harbinger", "DTC": "CACTA", "DYC": "Crypton", 
"DHH": "Helitron", "DMM": "Maverick", "polinton": "DNA/Maverick", 
"TIR/Tc1_Mariner": "DNA/TcMar-Tc1", "DIRS": "LTR/DIRS1"}

classified_h, classified_s = parse_fasta("TE_libraries/Classified_EDTA.fasta")

re_formatted = classified_h.copy()

for i, h in enumerate(re_formatted):

    for k,v in reformation_dict.items():
        if k in h:
            re_formatted[i] = h.replace(k,v)

classified_reformat = list(map(lambda st: str.replace(st, "MITE", "DNA"), re_formatted))

list(zip(classified_h, classified_reformat))[:10]

#%%

with open("TE_libraries/ReFormatted_Classified_EDTA.fasta", "w") as file1:
    for number, header in enumerate(classified_reformat):
        file1.write(header + "\n")
        file1.write(classified_s[number] + "\n")

if os.path.getsize("TE_libraries/ReFormatted_Classified_EDTA.fasta") > 0:
    print(".... ReFormatted_Classified_EDTA.fasta is written ....")
else:
    print("Something went wrong with File!")