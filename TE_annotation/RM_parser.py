from sys import argv
import pandas as pd
import numpy as np
import itertools
import os

"""
This scripts parses the BED file that came from .out file, produced by RepeatMasker, and returns an up-to-date categorisation of TEs.
"""

#input_file = sys.argv[-2]
#assembly_size = sys.argv[-1]

file_out = "/christos/Documents/Thesis/Genome_Annotation/TE_Annotation/BEDfiles/RM_output.merged.bed"
RM = pd.read_csv(file_out, sep="\t", header=None)
assembly_pt = 902353306

def RM_parser(RM_out, assembly):

    print("===> STARTING . . .")
    # name the columns of bed file
    names = ["Contig", "Start","End", "RepeatCl", "ID"]
    RM_out.columns = names

    # replace nan values in the column "ID" with the word "empty" to use afterwards for elements count
    RM_out["ID"] = RM_out["ID"].replace(np.nan, "empty")

    # change naming given by previous classification tools to new ones corresponding to up-to-date classification naming of TEs
    mapping_dict_regex = {
        "DNA/hAT-hobo":"TIR/hAT",
        "DNA/Mutator": "TIR/Mutator",
        "DNA/CACTA": "TIR/CACTA", 
        "DNA/PIF-Harbinger": "TIR/PIF-Harbinger",
        "DNA/Helitron": "Helitron/Helitron",
        "DNA/nMITE": "DNA", 
        "DNA/TcMar-Mariner": "TIR/Tc1-Mariner",
        "DNA/TcMar-Tc1": "TIR/Tc1-Mariner", 
        "DNA/MITE": "DNA", 
        "DNA/TcMar": "TIR/Tc1-Mariner", 
        "PLE": "PLE/Penelope", 
        "DNA/TcMar/nMITE": "TIR/Tc1-Mariner", 
        "LTR/BEL": "LTR/BEL-Pao", 
        "DNA/P": "TIR/P", 
        "LINE/Penelope": "PLE/Penelope", 
        "LTR/DIRS1": "DIRS/DIRS1", 
        "DNA/PiggyBac": "TIR/PiggyBac", 
        "DNA/Maverick": "Maverick/Maverick"
        }

    # fix some older irregularities lost in the previous steps
    mapping_dict2_regex = {
        "Unknown": "Unclassified",
        "DNA/":"DNA",
        "/nMITE":""
        }

    RM_out["RepeatCl"].replace(mapping_dict_regex, inplace=True, regex = True)
    RM_out["RepeatCl"].replace(mapping_dict2_regex, inplace=True, regex = True)

    # add extra column with the length of each element for percentage calculation
    RM_out["length"] = (RM_out["End"] - RM_out["Start"]) + 1

    masked_genome = sum(RM_out["length"])
    percentage = round((masked_genome/assembly) * 100, ndigits = 2)

    print("Bases masked: {} bp ({}%)".format(masked_genome, percentage))

    ## ====================================== NESTED FUNCTION ====================================== ##
    # Split the reformatted .bed file in two subfiles: 1) Only clean (individual) defined TEs, 2) Overlapping-complex TEs
    def split_dataframe(a_Df):

        clean_TEs = a_Df[~a_Df["RepeatCl"].str.contains(",")]
        merged_TEs = a_Df[a_Df["RepeatCl"].str.contains(",")]

        return clean_TEs, merged_TEs
    
    ## ============================================================================================= ##

    ## ====================================== NESTED FUNCTION ====================================== ##
    # Count the number of elements taking into account specific assumptions
    # IF pre-defined distinct elements merged by bedtools, their number is calculated based on given ID
    # IF there was no ID in an element and is replaced with "empty" is count as indvidual element
    # IF there is only one ID/element, is counted that way

    def count_elements(a_Df):
        IDs1, IDs2, IDs3 = [], 0, []
        
        for id in list(a_Df["ID"]):
            
            if "," in id:
                IDs3.append([int(i) for i in id.split(",")])

            else:
                if id == "empty":
                    IDs2 += 1
                else:
                    IDs1.append(int(id))

        return len(set(IDs1)) + IDs2 + len(set(list(itertools.chain(*IDs3))))

    ## ============================================================================================= ##

    ## ====================================== NESTED FUNCTION ====================================== ##
    # Calculate the total length occupied by a specific class of TEs
    def length_calc(a_Df):

        return sum(a_Df["length"])
    ## ============================================================================================= ##
        
    CLEAN, MERGED = split_dataframe(RM_out)


    RetroElements = [
        "LTR", "LTR/BEL-Pao", "LTR/Copia", "LTR/Gypsy", "LTR/ERV", 
        "DIRS", "DIRS/DIRS1", "DIRS/Ngaro", 
        "PLE/Penelope", 
        "LINE", "LINE/R2", "LINE/RTE", "LINE/Jockey", "LINE/L1",
        "SINE", "SINE/tRNA", "SINE/7L", "SINE/5S"
        ]

    DNAtransposons = [
        "TIR","TIR/Tc1-Mariner", "TIR/hAT", "TIR/Mutator", "TIR/Merlin", "TIR/Transib", "TIR/P", "TIR/PiggyBac", "TIR/PIF-Harbinger", "TIR/CACTA",
        "Helitron/Helitron", 
        "Maverick/Maverick"
        ]

    Others = ["Unclassified", "Simple_repeat", "Low_complexity"]

    ## ============================================================================================= ##
    # DNA transposons sub-table
    DNA_TRANSPOSONS = CLEAN[CLEAN["RepeatCl"].str.contains("DNA|TIR|Helitron|Maverick", regex=True)]

    ## ============================================================================================= ##
    # Reteroelements sub-table
    RETROELEMENTS = CLEAN[CLEAN["RepeatCl"].str.contains("LTR|SINE|LINE|DIRS|PLE", regex=True)]

    ## ============================================================================================= ##
    # Match overlapping elements to 3 categories (1.overlapping Retroelements, 2.overlapping DNA transposons, 3.Complex overlapping elements)
    conv_tmp = []
    for mixed in list(MERGED["RepeatCl"].unique()):

        if ("DNA" in mixed) or ("TIR" in mixed) or ("Helitron" in mixed) or ("Maverick" in mixed) :

            if ("LTR" in mixed) or ("SINE" in mixed) or ("LINE" in mixed) or ("PLE" in mixed) or ("DIRS" in mixed):
                conv_tmp.append("Complex")
                continue

            if ("Unclassified" in mixed) or ("Simple_repeat" in mixed) or ("Low_complexity" in mixed):
                conv_tmp.append("DNA transposons")
                continue

            else:
                if ("DNA" not in mixed) and ("Helitron" not in mixed) and ("Maverick" not in mixed):
                    conv_tmp.append("TIR elements")
                    continue
                else:
                    conv_tmp.append("DNA transposons")
                    continue
        
        else:

            if ("Unclassified" in mixed) or ("Simple_repeat" in mixed) or ("Low_complexity" in mixed):

                if ("LTR" in mixed) or ("SINE" in mixed) or ("LINE" in mixed) or ("PLE" in mixed) or ("DIRS" in mixed):
                    conv_tmp.append("Retroelements")
                    continue
                else:
                    conv_tmp.append("Complex")
                    continue
            else:
                if ("SINE" not in mixed) and ("LINE" not in mixed) and ("PLE" not in mixed) and ("DIRS" not in mixed):
                    conv_tmp.append("LTR elements")
                    continue
                else:
                    conv_tmp.append("Retroelements")
                    continue

    tmp = dict(zip(list(MERGED["RepeatCl"].unique()), conv_tmp))

    new = MERGED.replace(tmp)


    MERGED_DNA = new[new["RepeatCl"] == "DNA transposons"]
    MERGED_TIR = new[new["RepeatCl"] == "TIR elements"]
    MERGED_RETRO = new[new["RepeatCl"] == "Retroelements"]
    MERGED_LTR = new[new["RepeatCl"] == "LTR elements"]
    COMPLEX = new[new["RepeatCl"] == "Complex"]

    ## ============================================================================================= ##
    # export matching table as .TSV ("complex_elements_categorisation.tsv")
    tmp_df = pd.DataFrame(list(zip(list(MERGED["RepeatCl"].unique()), conv_tmp)), 
    columns = ["Overlapping elements", "Categorisation"])

    tmp_df.to_csv("complex_elements_categorisation.tsv", sep = "\t", index = False)

    if os.path.exists(os.getcwd() + "/complex_elements_categorisation.tsv"):
        print("====> 'complex_elements_categorisation.tsv' is written")
    else:
        print("Something went wrong with supplementary file!")

    ## ============================================================================================= ##
    # Calculate number of elements, their occupied length and percentage per TE category
    numbers_of_elements, lenth_occupied = [], []

    numbers_of_elements.append(count_elements(RETROELEMENTS) + len(MERGED_RETRO) + len(MERGED_LTR))
    lenth_occupied.append(length_calc(RETROELEMENTS) + length_calc(MERGED_RETRO) + length_calc(MERGED_LTR))

    for query1 in RetroElements:
        subRetro = RETROELEMENTS[RETROELEMENTS["RepeatCl"].str.contains(query1, regex=True)]
        if query1 == "LTR":
            numbers_of_elements.append(count_elements(subRetro) + len(MERGED_LTR))
            lenth_occupied.append(length_calc(subRetro) + length_calc(MERGED_LTR))
        else:
            numbers_of_elements.append(count_elements(subRetro))
            lenth_occupied.append(length_calc(subRetro))

    ## Only Retro-overlapping
    numbers_of_elements.append(len(MERGED_RETRO))
    lenth_occupied.append(length_calc(MERGED_RETRO))

    numbers_of_elements.append(count_elements(DNA_TRANSPOSONS) + len(MERGED_DNA) + len(MERGED_TIR))
    lenth_occupied.append(length_calc(DNA_TRANSPOSONS) + length_calc(MERGED_DNA) + length_calc(MERGED_TIR))
    for query2 in DNAtransposons:
        subDNA = DNA_TRANSPOSONS[DNA_TRANSPOSONS["RepeatCl"].str.contains(query2, regex=True)]
        if query2 == "TIR":
            numbers_of_elements.append(count_elements(subDNA) + len(MERGED_TIR))
            lenth_occupied.append(length_calc(subDNA) + length_calc(MERGED_TIR))
        else:
            numbers_of_elements.append(count_elements(subDNA))
            lenth_occupied.append(length_calc(subDNA))

    ## Only DNA-overlapping
    numbers_of_elements.append(len(MERGED_DNA))
    lenth_occupied.append(length_calc(MERGED_DNA))

    for query3 in Others:
        subOthers = CLEAN[CLEAN["RepeatCl"].str.contains(query3)]
        numbers_of_elements.append(len(subOthers))
        lenth_occupied.append(length_calc(subOthers))

    # DNA-Retro-Others-
    numbers_of_elements.append(len(COMPLEX))
    lenth_occupied.append(length_calc(COMPLEX))

    ## ============================================================================================= ##
    # Export statistics table as .TSV (TE_annotation_results.tsv)

    table_names = [
        "Retroelements (Class I)", "LTR", 
        "BEL-Pao", "Copia", "Gypsy", "ERV", 
        "DIRS", "DIRS1", "Ngaro",
        "Penelope",
        "LINE", "R2", "RTE", "Jockey", "L1",
        "SINE", "tRNA", "7L", "5S",
        "Overlapping Retroelements",
        "DNA transposons (Class II)", "TIR",
        "Tc1-Mariner", "hAT", "Mutator", "Merlin", "Transib", "P", "PiggyBac", "PIF-Harbinger", "CACTA",
        "Helitron", "Maverick", 
        "Overlapping DNA transposons",
        "Unclassified", "Simple repeats", "Low complexity", 
        "Complex overlapping elements"
        ]

    column_names = ["Transposable elements", "number of elements", "length occupied (bp)", "percentage of sequence (%)"]

    percentages = list(np.around((np.array(lenth_occupied)/assembly)*100, 2))

    results = pd.DataFrame(list(zip(table_names,numbers_of_elements,lenth_occupied, percentages)), 
    columns = column_names)

    results.to_csv("TE_annotation_results.tsv", sep = "\t", index = False)

    if os.path.exists(os.getcwd() + "/TE_annotation_results.tsv"):
        print("====> 'TE_annotation_results.tsv' is written")
        print("===> END !")
    else:
        print("Something went wrong with annotation file!")

    ## ============================================================================================= ##

#%%
RM_parser(RM, assembly_pt)
