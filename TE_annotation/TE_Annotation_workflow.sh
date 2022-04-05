# ================================================================================================================================================================================== #
# ============================================================================= TE ANNOTATION WORKFLOW ============================================================================= #
# ================================================================================================================================================================================== #
# Christos Kitsoulis, 2021-2022
# Genomics & Bioinformatics Group at IMBBC, HCMR

# ======================================================================================================================================
## STEP 1: Install EDTA pipeline & Run on genome assembly for de novo TE families prediction and classification -- EDTA --
# INPUT: Assembly
# DEPENDENCIES: Check DeepTE documentation
# parameters: species = Others, sensitive = 1 (use repeatmodeler extension), anno = 1 (for whole genome TE annotation), threads

ASSEMBLY=/home1/christos_k/pterois_miles/Assembly_lionfish/GENOME_ASSEMBLY/lionfish_assembly.fasta
EDTA_PATH=/home1/christos_k/software/EDTA

perl $EDTA_PATH/EDTA.pl --genome $ASSEMBLY --species Others --sensitive 1 --anno 1 --threads 32
# ======================================================================================================================================


# ======================================================================================================================================
## STEP 2: Separate TE.lib to 3 sub-libraries: 1) Classified families, 2) UnOrdered families, 3) UnClassified elements -- library_split.py --
# INPUT: lionfish_EDTA_TEfamilies.lib
# DEPENDENCIES: Bio (biopython), os
# OUTPUT: Unknown_EDTA.fasta, UnOrdered_EDTA.fasta, Classified_EDTA.fasta

python library_split.py

mkdir libraries/

mv Unknown_EDTA.fasta UnOrdered_EDTA.fasta Classified_EDTA.fasta libraries/
# ======================================================================================================================================


# ======================================================================================================================================
## STEP 3: Install DeepTE & Reclassify sub-libraries (2) & (3) with DeepTE -- DeepTE --
# INPUT : Unknown_EDTA.fasta, UnOrdered_EDTA.fasta, results_dir, tmp_dir
# DEPENDENCIES: Check DeepTE documentation
# parameters: m (model) = M (for metazoa), sp (species) = M (for metazoans), prop_thr = 0.8 (probability classification threshold)
# OUTPUTS:

DEEPTE_PATH=/home1/christos_k/software/DeepTE

mkdir DeepTE_UnOrdered_temp2 DeepTE_UnOrdered2 DeepTE_Unknown_temp2 DeepTE_Unknown2

# rename .fasta file to ReClassified_Unknown_DeepTE08.fasta and mv to libraries dir

echo "Starting classification for unclassified"

$DEEPTE_PATH/DeepTE.py -d DeepTE_Unknown_temp2 -o DeepTE_Unknown2 -i libraries/Unknown_EDTA.fasta -m M -sp M -prop_thr 0.8

echo "Starting classification for unordered TEs"

$DEEPTE_PATH/DeepTE.py -d DeepTE_UnOrdered_temp2 -o DeepTE_UnOrdered2 -i libraries/UnOrdered_EDTA.fasta -m M -sp M -prop_thr 0.8

# rename .fasta file to ReOrdered_DeepTE08.fasta and mv to libraries dir
# ======================================================================================================================================


# ======================================================================================================================================
## STEP 4: ReFormat headers -- reformat_headers.py & bash --
# INPUT: Classified_EDTA.fasta, ReClassified_Unknown_DeepTE08.fasta, ReOrdered_DeepTE08.fasta

python reformat_headers.py #for classified lib

# custom reformat with bash commands in FormatIDs.sh
# ======================================================================================================================================


# ======================================================================================================================================
## STEP 5: Concatenate libraries -- bash --
# INPUT: sub-libs after header simplificatiion step

cat ReFormatted_Classified.fasta ReClassified_DeepTE08_Curated.fasta ReOrdered_DeepTe_Curated.fasta > TE_lib_EDTA_DeepTE_customized.fasta
# ======================================================================================================================================


# ======================================================================================================================================
## STEP 6: Align TE families in the genome and soft-mask the corresponding regions -- RepeatMasker --
# INPUT: TE_lib_EDTA_DeepTE_customized.fasta, $ASSEMBLY
# parameters: -s (sensitive search), -xsmall (soft-masking, lowercase), -pa (number of processors to run in parallel), -gff (to return GFF file)
# OUTPUT: .tbl (stats), .gff, .out

mkdir RMasker_EDTA_DeepTE

RepeatMasker -s -xsmall -lib libraries/TE_lib_EDTA_DeepTE_customized.lib -dir RMasker_EDTA_DeepTE -pa 10 -gff $ASSEMBLY
# ======================================================================================================================================


# ======================================================================================================================================
## STEP 7: Parse .out file to calculate TE families representation based on up-to-date TE classification -- bash & BEDTOOLS & RM_parser.py --
# INPUT: lionfish_assembly.fasta.out
# OUTPUT: .bed file, .csv file (up-to-date TE classification stats)

# produce BED file from .out with bash command:
awk 'BEGIN{OFS="\t"}{if(NR>3) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6,$7,$10,".",strand,$11,$15}}' lionfish_assembly.fasta.out > RM_output.bed

# sort produced BED file
sort -k1,1 -k2,2n RM_output.bed > RM_output.sorted.bed

# merge overlapping intervals in sorted file
bedtools merge -i RM_output.sorted.bed -c 7,8 -o distinct,distinct -delim "," > RM_output.merged.bed

# parse sorted .bed file
python RM_parser.py
# ======================================================================================================================================
