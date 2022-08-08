# ================================================================================================================================================================================== #
# ============================================================================ GENE PREDICTION WORKFLOW ============================================================================ #
# ================================================================================================================================================================================== #
# Christos Kitsoulis, 2021-2022
# Genomics & Bioinformatics Group at IMBBC, HCMR

MASKED_GENOME=/home/christos/pterois_miles/Genome_annotation/gene_prediction/pm_assembly_MASKED.fasta
UNIPROT_FISH=/home/christos/pterois_miles/Genome_annotation/gene_prediction/UniProt_24FishProteomes.fasta
HQ_TRANSCRIPTS=/home/christos/pterois_miles/Genome_annotation/gene_prediction/hq_transcripts

## ===== STEP 1 =====
### Map high-quality consensus transcripts (Iso-seq) to soft-masked genome assembly -- GMAP -- 
# Create DB from genome assembly
# INPUT: $MASKED_GENOME, $HQ_TRANSCRIPTS
# parameters: -t (threads), -d (database), -f samse (SAM format), -n 0 (single alignment), --max-intronlength-ends (Max length for one internal intron),
# --cross-species (more sensitive search for canonical splicing), -z (cDNA direction), -D (database directory)
# output: .GFF3

gmap_build -d DevilFirefish $MASKED_GENOME -t 20

# Map long reads to assembly
gmap -D /home/christos/software/miniconda3/envs/CupCake/share/DevilFirefish -d DevilFirefish_s -f samse -n 0 -t 20 --cross-species --max-intronlength-ends 200000 -z sense_force $TRANSCRIPTS > hq_transcripts.sam 2> log_file.log

# Convert SAM to BAM + SORT
samtools view -bS hq_transcripts.sam > hq_transcripts.bam
samtools sort hq_transcripts.bam > hq_transcripts.sorted.bam
samtools index hq_transcripts.sorted.bam

### Collapse hq transcripts to produce a non-redundant hq set -- cDNA_Cupcake --
# parameters: --dun-merge-5-shorter (Don't collapse shorter 5' transcripts), --cpus
python /home/christos/software/cDNA_Cupcake/cupcake/tofu/collapse_isoforms_by_sam.py --input hq_transcripts.fasta  --bam hq_transcripts.sorted.bam --dun-merge-5-shorter -o hq_isoforms

## ===== STEP 2 =====
### Generate a BLAST DB for protein homology -- DIAMOND -- 
# INPUT: database (.fasta) of fish proteomes

diamond makedb --in $UNIPROT_FISH --db UniProt_24FishProteomes

## ===== STEP 3 =====
### Generate a conﬁguration ﬁle -- MIKADO -- 
# INPUT: Mikado_input.txt (list of trascripts evidence), mammalian.yaml (scoring parameters file, default)

mikado configure -t 20 --list Mikado_input.txt --reference $MASKED_GENOME --mode permissive --scoring mammalian.yaml  --copy-scoring mammalian.yaml -bt $UNIPROT_FISH configuration.yaml

## ===== STEP 4 =====
### Mikado prepare step -- MIKADO --
# INPUT: configuration file (.yaml) produced in STEP 3

mikado prepare --json-conf configuration.yaml

## ===== STEP 5 =====
### Generate homology-based evidence against diamond DB which produced on STEP 2 -- DIAMOND --
# INPUT: mikado_prepared.fasta (mikado proposed evidence set)
# parameters: --max-target-seqs (maximum number of sequences to consider), --sensitive (sensitive mode), -db (diamond db), --evalue (cut-off), --outfmtt 5
# output: mikado.diamond.xml

diamond blastx --query mikado_prepared.fasta --max-target-seqs 24 --sensitive --index-chunks 1 --threads 20 --db UniProt_24FishProteomes.dmnd --evalue 1e-6 --outfmt 5 --out mikado.diamond.xml

## ===== STEP 6 =====
### Generate ORF predictions of mikado_prepared.fasta transcripts -- TRANSCODER --

TransDecoder.LongOrfs -t mikado_prepared.fasta

TransDecoder.Predict -t mikado_prepared.fasta

## ===== STEP 7 =====
### Merge all information-evidence to generate final consensus gene set -- MIKADO --
# INPUT: configuration file (.yaml), mikado-diamond (.xml), ORF predictions (.gff3), fish proteomes db (.fasta), mikado transcripts (_prepared.fasta)
# parameters: --procs (processors)

mikado serialise --procs 1 --json-conf configuration.yaml --xml ./Mikado_results/mikado.diamond.xml --orfs ./Transdecoder_results/mikado_prepared.fasta.transdecoder.gff3 --blast_targets $UNIPROT_FISH --transcripts ./Mikado_results/mikado_prepared.fasta

mikado pick --procs 20 --json-conf configuration.yaml

## ===== STEP 8 =====
### Select a well annotated subset of genes to train AUGUSTUS later on
# INPUT: mikado annotation & metrics files
# parameters: -bc 0.5 (blast coverage threshold), -e (number of exons)

python select_training.py -bc 0.5 -e 2 mikado.loci.metrics.tsv mikado.loci.gff3

## ===== STEP 9 =====
### First round of training Augustus --AUGUSTUS --
# INPUT: training set (.gff3), masked genome
# parameters: --species (predefined or custom), --optrounds 2 (optimization rounds), --cpus (threads)

export AUGUSTUS_CONFIG_PATH=/home/christos/software/miniconda3/envs/BRAKER/config

mkdir aug_training/

autoAugTrain.pl --trainingset=training.gff3 --genome=$MASKED_GENOME --species=lionfish --workingdir=aug_training --optrounds=2 --cpus=20 --verbose

## ===== STEP 10 =====
### Convert gff3 to gtf keeping only coding regions and generate species-specific exon hints -- AGAT, Python --

agat_convert_sp_gff2gtf.pl --gff mikado.loci.gff3 -o mikado.loci.gtf

## ===== STEP 11 =====
### Generate spliced protein alignments to masked genome, from well annotated species -- EXONERATE --
# Download proteome from Ensembl and split to smaller fasta files
Species1_proteome=/home/christos/pterois_miles/Genome_annotation/gene_prediction/Oryzias_latipes.fasta
Species2_proteome=/home/christos/pterois_miles/Genome_annotation/gene_prediction/Argyrosomous_regius.fasta
Species3_proteome=/home/christos/pterois_miles/Genome_annotation/gene_prediction/Gasterosteus_aculeatus.fasta

mkdir exonerate_results/

fastasplit -f $Species1_proteome -o Exonerate_results/ -c 90 #same for Species 2 & 3

#========================================================================================================================================
# rename chunk files
input_PATTERN=fasta_chunk_0
output_PATTERN=fasta_chunk_

for f in species.fasta_chunk_0*;
do
        mv "$f" "$(echo "$f" | sed s/$input_PATTERN/$output_PATTERN/)";
done
#========================================================================================================================================

python exonerate_parallel1.py ~/pterois_miles/Genome_annotation/gene_prediction/exonerate_results $MASKED_GENOME

python exonerate_parallel2.py ~/pterois_miles/Genome_annotation/gene_prediction/exonerate_results $MASKED_GENOME

python exonerate_parallel3.py ~/pterois_miles/Genome_annotation/gene_prediction/exonerate_results $MASKED_GENOME

#========================================================================================================================================

## ===== STEP 12 =====
### Concatenate all gff files to one and filter out the alignment blocks/extract hints -- Bash/Python --
# in exonerate_results/

cat exonerate_results/*.gff > Oryzias_latipes._Exo.gff
cat exonerate_results2/*.gff > Argyrosomous_regius._Exo.gff
cat exonerate_results3/*.gff > Gasterosteus_aculeatus._Exon.gff

grep -E "^scaffold|^contig" Oryzias_latipes._Exo.gff  > Oryzias_latipes._Exo_clean.gff
grep -E "^scaffold|^contig" Argyrosomous_regius._Exo.gff > Argyrosomous_regius._Exo_clean.gff
grep -E "^scaffold|^contig" Gasterosteus_aculeatus._Exon.gff > Gasterosteus_aculeatus._Exo_clean.gff

# sort .exh.gff files
sort -k1,1 -k2,2n Species._Exo_clean.exh.gff > Species._Exo_clean_sorted.exh.gff
cat *._Exo_clean.exh.gff > Areguis_Olatipes_Gaceleatus._Exo_clean.exh.gff

# merge all hint-evidence
cat mikado.loci.exh.gff Areguis_Olatipes_Gaceleatus._Exo_clean_sorted.exh.gff > Pmiles_merged_hints.all.gff

## ===== STEP 13 =====
### Run Augustus for gene prediction using hint evidence and ab-initio predictions -- AUGUSTUS --
# INPUT: hints file (Pmiles_merged_hints.all.gff), configuration file (.cfg, custom or default), masked genome
# parameters: --gff3=on (produce GFF3 file), --species-lionfish (custom species for initialize training params), --progress=true (show progress),
# --alternatives-from-evidence=false (not report alternative transcripts when they are suggested by hints), --allow_hinted_splicesites=atac (allows Augustus to predict the (rare) introns)

export AUGUSTUS_CONFIG_PATH=/home/christos/software/miniconda3/envs/BRAKER/config

augustus --uniqueGeneId=true --progress=true --gff3=on --species=lionfish --hintsfile=Pmiles_merged_hints.all.gff --extrinsicCfgFile=extrinsic.Chris2.E.W.P.cfg --allow_hinted_splicesites=atac --alternatives-from-evidence=false $MASKED_GENOME > Pmiles.aug.final.out

# extract genes' sequences from mikado
agat_convert_sp_gxf2gxf.pl -g mikado.loci.gff3 -o mikado.loci.AGAT.gff3

gffread mikado.loci.AGAT.gff3 -V -w mikado.AGAT.genes.fasta -g $MASKED_GENOME

# evaluate completeness of predicted genes with BUSCO
busco -i mikado.AGAT.genes.fasta -m genome -o Pmiles_Mikado -l actinopterygii_odb10 --cpu 20 --augustus

## ===== STEP 14 =====
### Merge outputs from both MIKADO & AUGUSTUS to one single (consensus )gene set -- PASA --

# clean any spurious gene sequences
/home/christos/software/miniconda3/envs/PASA/opt/pasa-2.4.1/bin/seqclean mikado.AGAT.genes.fasta

# align mikado evidence genes to masked genome with PASA
/home/christos/software/miniconda3/envs/PASA/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c alignAssembly2.config -C -R -g $MASKED_GENOME -t mikado.AGAT.genes.fasta.clean -T -u mikado.AGAT.genes.fasta --ALIGNERS blat --CPU 20

# remove non-canonical terms from Augustus GFF3 file -- AGAT --
agat_convert_sp_gxf2gxf.pl -g Pmiles.aug.final.out.gff3 -o Pmiles.aug.final.AGAT.out.gff3

# extract gene sequences from augustus predicted gene set
gffread Pmiles.aug.final.AGAT.out.gff3 -V -w Pmiles.aug.final.AGAT.out.genes.fasta -g $MASKED_GENOME

# evaluate completeness of genes with BUSCO
busco -i Pmiles.aug.final.AGAT.out.genes.fasta -m genome -o Pmiles_Augustus -l actinopterygii_odb10 --cpu 20 --augustus

# load augustus predictions to initial gene set
/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly2.config -g $MASKED_GENOME -P Pmiles.aug.final.AGAT.out.gff3

## ===== STEP 15 =====
### Update PASA DB with Mikado predictions -- PASA --

/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c autoCompare2.config -A -g $MASKED_GENOME -t mikado.AGAT.genes.fasta.clean --CPU 20

# calculate statistics of produced gene set
agat_sp_statistics.pl --gff $FILE -g $MASKED_GENOME --output _Statistics_.txt

## ===== STEP 16 - ROUND 2 =====
# Update PASA round 2 -- PASA --

# load gene set from round 1 to initial PASA db
/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly2.config -g $MASKED_GENOME -P PmilesPASAdb.gene_structures_post_PASA_updates.24554.gff3

# update PASA db with mikado evidence 
/home1/christos_k/software/miniconda3/envs/PASA/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c autoCompare2.config -A -g $MASKED_GENOME -t mikado.AGAT.genes.fasta.clean --CPU 20

# extract protein sequences from GFF file
gffread -E PmilesPASAdb.gene_structures_post_PASA_updates.8106.gff3 -g $MASKED_GENOME -y Pmiles_PASAMikAug_rnd2_prot.fasta

# evaluate protein set with BUSCO 
busco -m proteins --cpu 20 -i Pmiles_PASAMikAug_rnd2_prot.fasta -o Pmiles_rnd2 -l actinopterygii_odb10

## ===== STEP 17 =====
### Filter predicted gene models -- gffread/Agat --

# clean GFF file from possible identical isoforms
agat_convert_sp_gxf2gxf.pl -g PmilesPASAdb.gene_structures_post_PASA_updates.8106.gff3 -o Pmiles_PASAMikAug_rnd2.AGAT.gff3

# identify potential gene models with in-frame STOP codons
gffread -E Pmiles_PASAMikAug_rnd2.AGAT.gff3 -g $MASKED_GENOME -V -H -o Pmiles_PASAMikAug_rnd2.simp.gff3

#========================================================================================================================
# filter out gene models with identified in-frame STOP codons                                                           
grep -o "'.*'" gene_filt2.output | sed "s/'//g" > gene_models.inframeSTOP.txt                                         
                                                                                                                       
grep -v -f gene_models.inframeSTOP.txt Pmiles_PASAMikAug_rnd2.AGAT.gff3 > Pmiles_PASAMikAug_rnd2.AGAT.noSTOP.gff3     
#========================================================================================================================

#========================================================================================================================
# find gene models overlapping TE (min overlap = 0.50) and exclude them                                                   
bedtools intersect -a Pmiles_PASAMikAug_rnd2.AGAT.noSTOP.gff3 -b lionfish_assembly.fasta.out.gff -wa -f 0.50             
                                                                                                                         
grep "gene" bedtools_intersect.output | awk -F"Name=" '{print $2}' > gene_models.overlappingTE.txt                       
                                                                                                                         
grep -v -f gene_models.overlappingTE.txt Pmiles_PASAMikAug_rnd2.AGAT.noSTOP.gff3 > Pmiles_PASAMikAug_rnd2.AGAT.filt.gff3 
#========================================================================================================================

# rename GFF file based on Ensembl-like nomenclature
agat_sp_manage_IDs.pl -f Pmiles_PASAMikAug_rnd2.AGAT.filt.gff3 --ensembl --prefix PMIL --type_dependent --tair -o Pmiles_annotation_v090522.gff3

# calculate statistics
agat_sp_statistics.pl --gff Pmiles_annotation_v090522.gff3 -g $MASKED_GENOME --output Statistics_final.txt

# extract protein sequences from GFF file
gffread -E Pmiles_annotation_v090522.gff3 -g $MASKED_GENOME -y Pmiles_v090522.fasta

# evaluate protein set with BUSCO
busco -m proteins --cpu 20 -i Pmiles_v090522.fasta -o Pmiles_filtered -l actinopterygii_odb10
