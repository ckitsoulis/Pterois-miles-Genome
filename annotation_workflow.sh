#!/bin/bash

#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
## #SBATCH --mem-per-cpu=4000
#SBATCH --job-name="S"
#SBATCH --output=s.output
#SBATCH --mail-user=chriskits@hotmail.com
#SBATCH --mail-type=ALL

# conda activate CupCake

MASKED_GENOME=/home1/christos_k/pterois_miles/Genome_annotation/gene_prediction/Workflow_Run2/lionfish_assembly_MASKED.fasta
UNIPROT_FISH=/home1/christos_k/pterois_miles/Genome_annotation/gene_prediction/UniProt_24FishProteomes.fasta

## ===== STEP 1 =====
### Map collapsed transcripts to softmasked genome assembly -- GMAP -- ###

# Create DB from genome assembly
# input: $MASKED_GENOME, $TRANSCRIPTS
# output: .GFF3

#gmap_build -d DevilFirefish_s lionfish_assembly_MASKED.fasta -t 18

#  Map long reads to assembly

#gmap -D /home1/christos_k/software/miniconda3/envs/CupCake/share/DevilFirefish_s -d DevilFirefish_s -f 3 -n 0 -x 50 -t 20 -B 4 --cross-species --gff3-add-separators=0 $TRANSCRIPTS > PacBioVsGenome_gmap.gff3 2> $

# conda activate MIKADO

## ===== STEP 2 =====
### Generate a BLAST DB for protein homology -- DIAMOND -- ###

#diamond makedb --in $UNIPROT_FISH --db UniProt_24FishProteomes

## ===== STEP 3 =====
### Generate a conﬁguration ﬁle -- MIKADO -- ###

#mikado configure -t 20 --list Mikado_input.txt --reference $MASKED_GENOME --mode permissive --scoring mammalian.yaml  --copy-scoring mammalian.yaml -bt $UNIPROT_FISH configuration.yaml

## ===== STEP 4 =====
### Mikado prepare step ###

#mikado prepare --json-conf configuration.yaml

## ===== STEP 5 =====
### Generate homology-based evidence against dmnd DB which produced on STEP 2 ###

#diamond blastx --query Mikado_results/mikado_prepared.fasta --max-target-seqs 24 --sensitive --index-chunks 1 --threads 20 --db UniProt_24FishProteomes.dmnd --evalue 1e-6 --outfmt 5 --out Mikado_results/mikado.diamond.xml

## ===== STEP 6 =====
### Generate ORF predictions of mikado_prepared.fasta transcripts -- TRANSCODER -- ###

#TransDecoder.LongOrfs -t Mikado_results/mikado_prepared.fasta

#TransDecoder.Predict -t Mikado_results/mikado_prepared.fasta

## ===== STEP 7 =====
### Merge all information and generate final output consensus gene set

#mikado serialise --procs 1 --json-conf configuration.yaml --xml ./Mikado_results/mikado.diamond.xml --orfs ./Transdecoder_results/mikado_prepared.fasta.transdecoder.gff3 --blast_targets $UNIPROT_FISH --transcripts ./Mikado_results/mikado_prepared.fasta

#mikado pick --procs 20 --json-conf configuration.yaml

# conda actiivate PyEnv

## ===== STEP 8 =====

#python ../Annotation_workflow/select_mik_train.py -f 0.5 -e 2 mikado.loci.metrics.tsv mikado.loci.gff3

# conda activate BRAKER2

## ===== STEP 9 =====

### First round of training Augustus --AUGUSTUS -- ###

#export AUGUSTUS_CONFIG_PATH=/home1/christos_k/software/miniconda3/envs/BRAKER2/config

#mkdir aug_training/

#autoAugTrain.pl --trainingset=training.gff3 --genome=$MASKED_GENOME --species=lionfish --workingdir=aug_training --optrounds=2 --cpus=20 --verbose

# conda activate Agat/PyEnv

## ===== STEP 10 =====
### Convert gff3 to gtf [(-T flag) and keep only coding transfrag(-C flag), highlighting any error (-E flag)] -- AGAT -- ###

#gffread --gtf -C -o mikado.loci.gtf mikado.loci.gff3
#agat_convert_sp_gff2gtf.pl --gff mikado.loci.gff3 -o mikado.loci.gtf

#python ../Annotation_workflow/gtfToHintsMik.py mikado.loci.gtf

# conda activate EXONERATE/PyEnv

## ===== STEP 11 =====
### Generate spliced protein alignments from a well annotated species -- EXONERATE -- ###
# Download proteome from Ensembl and split to smaller fasta files
# Multithread exonerate (x2 batch)

#mkdir exonerate_results
#fastasplit -f Oryzias_latipes_UniProt_180422.fasta -o Exonerate_results/ -c 90

#python exonerate_parallel1.py ~/pterois_miles/Genome_annotation/gene_prediction/Annotation_workflow/exonerate_results $MASKED_GENOME

#python exonerate_parallel2.py ~/pterois_miles/Genome_annotation/gene_prediction/Annotation_workflow/exonerate_results $MASKED_GENOME

#python exonerate_parallel2.py ~/pterois_miles/Genome_annotation/gene_prediction/Annotation_workflow/exonerate_results $MASKED_GENOME

#========================================================================================================================================
#fastasplit -f Argyrosomous_regius_vpapadog.fasta -o Exonerate_results2/ -c 132

#cd Exonerate_results2/

#input_PATTERN=fasta_chunk_0
#output_PATTERN=fasta_chunk_

#for f in Argyrosomous_regius_vpapadog.fasta_chunk_0*;
#do
#        mv "$f" "$(echo "$f" | sed s/$input_PATTERN/$output_PATTERN/)";
#done


#python exonerate_parallel.py ~/pterois_miles/Genome_annotation/gene_prediction/Workflow_Run2/Exonerate_results2 $MASKED_GENOME

#python exonerate_parallel2.py ~/pterois_miles/Genome_annotation/gene_prediction/Workflow_Run2/Exonerate_results2 $MASKED_GENOME

#python exonerate_parallel3.py ~/pterois_miles/Genome_annotation/gene_prediction/Workflow_Run2/Exonerate_results2 $MASKED_GENOME

#python exonerate_parallel4.py ~/pterois_miles/Genome_annotation/gene_prediction/Workflow_Run2/Exonerate_results2 $MASKED_GENOME
#=======================================================================================================================================

## ===== STEP 12 =====
### Concatenate all gff files to one and filter out the alignment blocks/extract a hint file using Ferdi's script -- Bash/Python -- ###
# in exonerate_results/

#cat Exonerate_results/*.gff > Oryzias_latipes.UniProt._Exo.gff
#cat Exonerate_results2/*.gff > Argyrosomous_regius.vpapadog._Exo.gff

#grep -E "^scaffold|^contig" Oryzias_latipes.UniProt._Exo.gff  > Oryzias_latipes.UniProt._Exo_clean.gff
#grep -E "^scaffold|^contig" Argyrosomous_regius.vpapadog._Exo.gff > Argyrosomous_regius.vpapadog._Exo_clean.gff

#python ../Annotation_workflow/exoToHints.py Oryzias_latipes.UniProt._Exo_clean.gff
#python ~/pterois_miles/Genome_annotation/gene_prediction/Annotation_workflow/exoToHints.py Argyrosomous_regius.vpapadog._Exo_clean.gff

#sort -k1,1 -k2,2n Areguis_Olatipes_Gaceleatus._Exo_clean.exh.gff > Areguis_Olatipes_Gaceleatus._Exo_clean_sorted.exh.gff
#cat Argyrosomous_regius.vpapadog._Exo_clean.exh.gff Oryzias_latipes.UniProt._Exo_clean.exh.gff Gasterosteus_aculeatus.ensembl._Exo_clean.exh.gff > Areguis_Olatipes_Gaceleatus._Exo_clean.exh.gff

# merge all hint-evidence
#cat mikado.loci.exh.gff Oryzias_latipes.UniProt._Exo_clean.exh.gff > Pmiles2_merged_hints.gff
#cat mikado.loci.exh.gff Areguis_Olatipes_Gaceleatus._Exo_clean_sorted.exh.gff > Pmiles_merged_hints.all.gff

## ===== STEP 13 =====
### Run Augustus for gene prediction using hint evidence and ab-initio predictions -- AUGUSTUS -- ###

#export AUGUSTUS_CONFIG_PATH=/home1/christos_k/software/miniconda3/envs/BRAKER2/config

#augustus --uniqueGeneId=true --progress=true --gff3=on --species=lionfish --hintsfile=Pmiles_merged_hints.all.gff --extrinsicCfgFile=extrinsic.Chris2.E.W.P.cfg --allow_hinted_splicesites=atac --alternatives-from-evidence=false $MASKED_GENOME > Pmiles.aug.final.out

# extract protein sequences from augustus
#agat_convert_sp_gxf2gxf.pl -g mikado.loci.gff3 -o mikado.loci.AGAT.gff3

#gffread mikado.loci.AGAT.gff3 -V -w mikado.AGAT.genes.fasta -g $MASKED_GENOME

# evaluate completeness of predicted genes with BUSCO
#busco -i mikado.AGAT.genes.fasta -m genome -o Pmiles_Mikado -l actinopterygii_odb10 --cpu 20 --augustus

## ===== STEP 14 =====
### Merge outputs from MIKADO & AUGUSTUS to one single gene set -- PASA -- ##

#/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/bin/seqclean mikado.AGAT.genes.fasta

#/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c alignAssembly2.config -C -R -g $MASKED_GENOME -t mikado.AGAT.genes.fasta.clean -T -u mikado.AGAT.genes.fasta --ALIGNERS blat --CPU 20

# remove non-canonical terms from Augustus GFF3 file -- AGAT --
#agat_convert_sp_gxf2gxf.pl -g Pmiles.aug.final.out.gff3 -o Pmiles.aug.final.AGAT.out.gff3

# extract gene sequences
#gffread Pmiles.aug.final.AGAT.out.gff3 -V -w Pmiles.aug.final.AGAT.out.genes.fasta -g $MASKED_GENOME

# evaluate completeness of genes with BUSCO
#busco -i Pmiles.aug.final.AGAT.out.genes.fasta -m genome -o Pmiles_Augustus -l actinopterygii_odb10 --cpu 20 --augustus

#/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly2.config -g $MASKED_GENOME -P Pmiles.aug.final.AGAT.out.gff3

## ===== STEP 15 =====
### Update PASA DB with Mikado predictions -- PASA --

#/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c autoCompare2.config -A -g $MASKED_GENOME -t mikado.AGAT.genes.fasta.clean --CPU 20

#agat_sp_statistics.pl --gff $FILE -g $MASKED_GENOME --output _Statistics_.txt

## ===== STEP 16 - ROUND 2 =====
# Update PASA round 2 -- PASA --

#/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly2.config -g $MASKED_GENOME -P PmilesPASAdb.gene_structures_post_PASA_updates.24554.gff3

#/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c autoCompare2.config -A -g $MASKED_GENOME -t mikado.AGAT.genes.fasta.clean --CPU 20

# extract protein sequences from GFF file
#gffread -E PmilesPASAdb.gene_structures_post_PASA_updates.8106.gff3 -g $MASKED_GENOME -y Pmiles_PASAMikAug_rnd2_prot.fasta

# evaluate protein set with BUSCO 
#busco -m proteins --cpu 20 -i Pmiles_PASAMikAug_rnd2_prot.fasta -o Pmiles_rnd2 -l actinopterygii_odb10

## ===== STEP 17 =====
### Filter predicted gene models

# clean GFF file from possible identical isoforms
#agat_convert_sp_gxf2gxf.pl -g PmilesPASAdb.gene_structures_post_PASA_updates.8106.gff3 -o Pmiles_PASAMikAug_rnd2.AGAT.gff3

# identify potential gene models with in-frame STOP codons
#gffread -E Pmiles_PASAMikAug_rnd2.AGAT.gff3 -g $MASKED_GENOME -V -H -o Pmiles_PASAMikAug_rnd2.simp.gff3

#========================================================================================================================
# filter out gene models with identified in-frame STOP codons                                                           #
# grep -o "'.*'" gene_filt2.output | sed "s/'//g" > gene_models.inframeSTOP.txt                                         #
#                                                                                                                       #
# grep -v -f gene_models.inframeSTOP.txt Pmiles_PASAMikAug_rnd2.AGAT.gff3 > Pmiles_PASAMikAug_rnd2.AGAT.noSTOP.gff3     #
#========================================================================================================================

#==========================================================================================================================
# find gene models overlapping TE (min overlap = 0.50) and exclude them                                                   #
#bedtools intersect -a Pmiles_PASAMikAug_rnd2.AGAT.noSTOP.gff3 -b lionfish_assembly.fasta.out.gff -wa -f 0.50             #
#                                                                                                                         #
#grep "gene" bedtools_intersect.output | awk -F"Name=" '{print $2}' > gene_models.overlappingTE.txt                       #
#                                                                                                                         #
#grep -v -f gene_models.overlappingTE.txt Pmiles_PASAMikAug_rnd2.AGAT.noSTOP.gff3 > Pmiles_PASAMikAug_rnd2.AGAT.filt.gff3 #
#==========================================================================================================================

# extract protein sequences from GFF file
#gffread -E Pmiles_PASAMikAug_rnd2.AGAT.filt.gff3 -g $MASKED_GENOME -y Pmiles_final.proteins.fasta

# evaluate protein set with BUSCO
#busco -m proteins --cpu 20 -i Pmiles_final.proteins.fasta -o Pmiles_filt -l actinopterygii_odb10

# rename GFF file based on Ensembl-like nomenclature
#agat_sp_manage_IDs.pl -f Pmiles_PASAMikAug_rnd2.AGAT.filt.gff3 --ensembl --prefix PMIL --type_dependent --tair -o Pmiles_annotation_v050522.gff3

# calculate statistics
#agat_sp_statistics.pl --gff Pmiles_annotation_v050522.gff3 -g $MASKED_GENOME --output Statistics_final.txt

# extract protein sequences from GFF file
#gffread -E Pmiles_annotation_v050522.gff3 -g $MASKED_GENOME -y Pmiles_prot_v050522.fasta

#====================================================================== EXTRA =======================================================================
#grep -v -f inFrame2.txt Pmiles_annotation_v050522.gff3 > Pmiles_annotation_v090522.gff3

#agat_sp_manage_IDs.pl -f Pmiles_annotation_v090522.gff3 --ensembl --prefix PMIL --type_dependent --tair -o Pmiles_v090522.gff3

#gffread -E Pmiles_v090522.gff3 -g $MASKED_GENOME -y Pmiles_v090522.pep.fasta

agat_sp_statistics.pl --gff Pmiles_v090522.gff3 -g $MASKED_GENOME --output Statistics_100522.txt
