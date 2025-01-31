hello-world.nf

#!/usr/bin/env nextflow

/*
 * Use echo to print 'Hello World!' to standard out
 */

####################################### STAR INDEXING ##########################################################

#!/usr/bin/env bash

##For Galaxy use (SGE)
##qsub -cwd -pe threads 8 -j yes -o /slipstream/home/mmariani/projects/hhv6_time_course/logs ./index_star.bash
##qsub -cwd -pe threads 8 -j yes -o /slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs ./index_star.bash
##11/09/2019 (human) and 02/21/2021 (virus)
##qsub -cwd -pe threads 20 -j yes -o /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs /slipstream/home/mmariani/scripts/vzv_interactions_scripts/vzv_rna_seq_index_star_mm.bash

##For VACC use (PBS-Torque):
##PBS -q poolmemq
##PBS -l nodes=1:ppn=16,mem=64gb,vmem=64gb
##PBS -l walltime=03:00:00
##PBS -N star_index_job
##PBS -j oe   

########################################### NOTES ##############################################################
################################################################################################################

##This script will do indexing for STAR
##note that:
##"I believe you don't need --sjdbGTFtagExonParentTranscript Parent 
##to index GTF files 
##https://www.biostars.org/p/300066/

##To combine GTF of hg38 with vzv:
##gffread hhv3.gff3 -T -o hhv3.gtf
##cat ucsc_hg38.gtf hhv3.gtf > ucsc_hg38_with_vzv.gtf 

#https://www.biostars.org/p/93883/
#This is quiet close to what you have asked:
#https://groups.google.com/forum/#!topic/rna-star/J6qH9JCysZw
#Basically, sjdbOverhang should be set as readlength -1. 
#So if you have 75 bp read then it should be set to 74. 
#Whereas, alignSJDBoverhangMin option ignores the alignment 
#with a small spilce overhangs. I use the default settings for this parameter.

##For larger genome, e.g. hg38, or combined hg38+virus genome:

##If using gff3 for example with virus ref:
##From Seth: "I’ve used 7 for this value with VZV in the past. 
##From the manual for small genomes, ‘this needs to be scaled down, 
##with a typical value of min(14, log2(GenomeLength)/2 - 1). 
##For example, for 1 megaBase genome, this is equal to 9, 
##for 100 kiloBase genome, this is equal to 7.’" 

##The below parameters should work for indexing small gff3 files (or gtf files, e.g. after 
##converting .gff3 to .gtf with gffread tool) such as our herpesvirus genomes, 
##but need to set CDS as exon if there are no exons present,
##or could create lines with "exon" corresponding to the "CDS": in column 3,
##or if there are a handful of exon lines already present then could
##change the few "exons" in column 3 to "CDS", or could change all "CDS"
##in column 3 to "exon". Make sure to note what you decide to do when altering
##gff/gtf/gff3 files and indexing with STAR.

##For newest HHV6a genome for example, gffread produces gtf with 
##only "CDS" so I changed --sjdbGTFfeatureExon exon to
##--sjdbGTFfeatureExon CDS when indexing with STAR

##This works and will do counts by gene for a gtf file
##With no "exon" in 3rd column - instead has "CDS":
##--sjdbGTFfeatureExon CDS \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id

##This works and will do counts by gene for a gff3 file
##With no "exon" in 3rd column - instead has "CDS":
##--sjdbGTFfeatureExon CDS \
##--sjdbGTFtagExonParentTranscript Parent \
##--sjdbGTFtagExonParentGene gene \

##Can also consider setting sjdbOverhang (based on read size)
##--sjdbOverhang 49
##star_path="/gpfs1/home/s/f/sfrietze/programs/STAR/bin/Linux_x86_64/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/"
##star_gtf_path="/gpfs1/home/s/f/sfrietze/references/hhv1/hhv1_gffread.gtf"
##star_fasta_path="/gpfs1/home/s/f/sfrietze/references/hhv1/hhv1.fasta"

########################################### CODE ##################################################################################
###################################################################################################################################

star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"

######################################### Index Small virus genomes (e.g. VZV) ####################################################
###################################################################################################################################

##HHV6 Newest
##gffread hhv6_newest.gff3 -T -o hhv6_newest.gtf
##star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/star"
##star_fasta_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest_mm_cds_to_exon.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.gtf"
##star_gff3_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.gff3"

##HHV6_GFP
##star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/star"
##star_fasta_path="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/hhv6_gfp.fasta"
##star_gff3_path="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/hhv6_gfp.gff3"

##VZV
##gffread vzv.gff3 -T -o vzv.gtf
##September? 2019 attempt to index small viral genomes with STAR:
##star_ref_dir="/slipstream/home/mmariani/references/hhv3/star/"
##star_fasta_path="/slipstream/home/mmariani/references/hhv3/hhv3.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.exon.removed.gtf"
##star_gff3_path="/slipstream/home/mmariani/references/hhv3/hhv3.gff3"
##star_gff3_path="/slipstream/home/mmariani/references/hhv3/hhv3.exon.removed.gff3"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.single.exon.changed.to.cds.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.all.changed.to.exon.gtf"
##Here (above) I changed the single "exon" (3rd column) line in the gtf file created 
##from the genbank gff3 to "CDS", thus all 3rd column entries are "CDS" in the gtf 
##file then I indexed with STAR below

##02/21/2021 need to reindec for newer version of star
##dont know which one I used above, so I will go with 
##vzv.gtf which I think was created from 
##the gffread program, first.

##star_ref_dir="/slipstream/home/mmariani/references/vzv/star/"
##star_fasta_path="/slipstream/home/mmariani/references/vzv/vzv.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv.gtf"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 20 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--genomeSAindexNbases 7 \
##--sjdbGTFfile $star_gtf_path \
##--sjdbGTFfeatureExon exon \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id

##--readFilesCommand zcat \

######################################### Index Hg38 and VZV all exons gtf ##############################################################################
###################################################################################################################################

##star_ref_dir="/slipstream/home/mmariani/references/ucsc_hg38_canonical_with_vzv/star"
##star_fasta_path="/slipstream/home/mmariani/references/ucsc_hg38_canonical_with_vzv/ucsc_hg38_canonical_with_vzv.fa"
##star_gtf_path="/slipstream/home/mmariani/references/ucsc_hg38_canonical_with_vzv/ucsc_hg38_canonical_with_vzv_all_exons.gtf"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 20 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--genomeSAindexNbases 14 \
##--sjdbGTFfile $star_gtf_path \
##--sjdbGTFfeatureExon exon \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id

##--readFilesCommand zcat \

######################################### Index Small virus genomes (e.g. VZV) ####################################################
###################################################################################################################################

##HHV6 Newest
##gffread hhv6_newest.gff3 -T -o hhv6_newest.gtf
##star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/star"
##star_fasta_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest_mm_cds_to_exon.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.gtf"
##star_gff3_path="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest.gff3"

##HHV6_GFP
##star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##star_ref_dir="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/star"
##star_fasta_path="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/hhv6_gfp.fasta"
##star_gff3_path="/slipstream/home/mmariani/references/hhv6/hhv6_gfp/hhv6_gfp.gff3"

##VZV
##gffread vzv.gff3 -T -o vzv.gtf
##September? 2019 attempt to index small viral genomes with STAR:
##star_ref_dir="/slipstream/home/mmariani/references/hhv3/star/"
##star_fasta_path="/slipstream/home/mmariani/references/hhv3/hhv3.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.exon.removed.gtf"
##star_gff3_path="/slipstream/home/mmariani/references/hhv3/hhv3.gff3"
##star_gff3_path="/slipstream/home/mmariani/references/hhv3/hhv3.exon.removed.gff3"
##star_gtf_path="/slipstream/home/mmariani/references/hhv3/hhv3.single.exon.changed.to.cds.gtf"
##star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv.all.changed.to.exon.gtf"
##Here (above) I changed the single "exon" (3rd column) line in the gtf file created 
##from the genbank gff3 to "CDS", thus all 3rd column entries are "CDS" in the gtf 
##file then I indexed with STAR below

##02/24/2021 Now try with all "exon" or "CDS" as the third column

##star_ref_dir="/slipstream/home/mmariani/references/vzv/star_all_exon/"
##star_fasta_path="/slipstream/home/mmariani/references/vzv/vzv.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv_all_exon.gtf"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 20 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--genomeSAindexNbases 7 \
##--sjdbGTFfile $star_gtf_path \
##--sjdbGTFfeatureExon exon \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id

##--readFilesCommand zcat \

##02/24/2021 Now try with all "CDS" as the third column

##star_ref_dir="/slipstream/home/mmariani/references/vzv/star_all_cds/"
##star_fasta_path="/slipstream/home/mmariani/references/vzv/vzv.fasta"
##star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv.single.exon.changed.to.cds.gtf"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 20 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--genomeSAindexNbases 7 \
##--sjdbGTFfile $star_gtf_path \
##--sjdbGTFfeatureExon CDS \
##--sjdbGTFtagExonParentTranscript transcript_id \
##--sjdbGTFtagExonParentGene gene_id
##
####--readFilesCommand zcat \

##02/24/2021
##Now try with the new depledge gtf that we were working 
##with for the scRNA data

star_ref_dir="/slipstream/home/mmariani/references/vzv/vzv_depledge/star"
star_fasta_path="/slipstream/home/mmariani/references/vzv/vzv_depledge/vzv_depledge_cellranger.fa"
star_gtf_path="/slipstream/home/mmariani/references/vzv/vzv_depledge/vzv_depledge_adjusted_cellranger.gtf"

$star_path \
--runMode genomeGenerate \
--runThreadN 20 \
--genomeDir $star_ref_dir \
--genomeFastaFiles $star_fasta_path \
--genomeSAindexNbases 7 \
--sjdbGTFfile $star_gtf_path \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript transcript_id \
--sjdbGTFtagExonParentGene gene_id

##--readFilesCommand zcat \

######################################### Index Hg38 ##############################################################################
###################################################################################################################################

##11/09/2019, Hg38 indexing:
##Use the prebuilt ucsc reference files for hg38 to build 
##a new STAR genome, for newest version of STAR:

##star_fasta_path="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa"
##star_ref_dir="/slipstream/home/mmariani/references/ucsc_hg38_star"
##star_gtf_path="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
##
##cd "/slipstream/home/mmariani/references/ucsc_hg38_star"
##
##$star_path \
##--runMode genomeGenerate \
##--runThreadN 16 \
##--genomeDir $star_ref_dir \
##--genomeFastaFiles $star_fasta_path \
##--sjdbGTFfile $star_gtf_path \
##> "star.index.hg38.11092019.logs.txt"

####################################### BOWTIE2 INDEX and SORTING ##############################

####################################### STAR ALIGNMENT #########################################
#!/usr/bin/env bash

##For Galaxy use (SGE):
##qsub -cwd -pe threads 8 -j yes -o /slipstream/home/mmariani/projects/hhv6_time_course/logs ./star_script_mm.bash
##qsub -cwd -pe threads 8 -j yes -o /slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs ./star_script_mm.bash
##11/09/2019 (human) and 02/21/2021 (virus)
##qsub -cwd -pe threads 16 -j yes -o /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs /slipstream/home/mmariani/scripts/vzv_interactions_scripts/vzv_rna_seq_star_script_mm.bash

##For VACC use (PBS-Torque):
##PBS -q poolmemq
##PBS -l nodes=1:ppn=16,mem=128gb,vmem=128gb
##PBS -l walltime=30:00:00
##PBS -N hhv6_rna_seq
##PBS -j oe
##PBS -o /gpfs1/home/s/f/sfrietze/mike_m/hhv6_rna_project/logs

##PBS -M michael.mariani@uvm.edu
##PBS -m bea

##spack load STAR

##This script I have smoothed out in August 2019 and pairs well the the 
##index_star.bash script, it can be used for single-end or paired-end star data
##there is a featureCounts block added to both SE and PE routines, and 
##bamCoverage routine and/or samtools stats routines (along with multiqc) 
##to compile can also be added, or run separately

star_path="/slipstream/home/mmariani/programs/miniconda3/bin/STAR"
##mode="single_end_h"
mode="single_end_v"
num_threads=16
t_value="exon"
g_value="gene_id"

##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_with_vzv"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_then_vzv"
##output_dir="/slipstream/home/mmariani/projects/hhv6_time_course/output_hhv6"
##output_dir="/slipstream/home/mmariani/projects/hhv6_time_course/output_hg38"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_alignment"

##ref_dir="/slipstream/home/mmariani/references/ucsc_hg38_with_vzv/star"
##ref_dir_h="/slipstream/home/mmariani/references/ucsc_hg38_star"
##ref_dir_v="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/star"
##ref_dir_v="/slipstream/home/mmariani/references/hhv3/star"

##gtf_file="/slipstream/home/mmariani/references/ucsc_hg38_with_vzv/ucsc_hg38_with_vzv.gtf"
##gtf_file_h="/slipstream/home/mmariani/references/ucsc_hg38_with_vzv/ucsc_hg38.gtf"
##gtf_file_v="/slipstream/home/mmariani/references/hhv6/hhv6a_newest/hhv6a_newest_mm_cds_to_exon.gtf"
##gtf_file_v="/slipstream/home/mmariani/references/hhv3/hhv3.single.exon.changed.to.cds.gtf"
##gtf_file_v="/slipstream/home/mmariani/references/hhv3/hhv3.all.changed.to.exon.gtf"

############################### Single end alignment: human #################################################################
#######################################################################################################################

if [ $mode == "single_end_h" ]
then

output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_alignment"
ref_dir_h="/slipstream/home/mmariani/references/ucsc_hg38_star"
gtf_file_h="/slipstream/home/mmariani/references/ucsc_hg38_with_vzv/ucsc_hg38.gtf"

while IFS=$'\t' read -r -a my_array
do

cd $output_dir

base_name=$(basename "${my_array[0]}" ".fastq.gz")

$star_path \
--runThreadN $num_threads \
--genomeDir $ref_dir_h \
--readFilesIn "${my_array[0]}" \
--readFilesCommand zcat \
--outFileNamePrefix $output_dir"/"$base_name".hg38." \
--outSAMtype BAM SortedByCoordinate \
> $base_name".log.txt"

file_name=$output_dir"/"$base_name".hg38.Aligned.sortedByCoord.out.bam"
samtools index $file_name
samtools idxstats $file_name > $output_dir"/"$(basename $file_name ".bam")".idxstats"
samtools stats $file_name > $output_dir"/"$(basename $file_name ".bam")".stats"
samtools flagstat $file_name > $output_dir"/"$(basename $file_name ".bam")".flagstat"

##--outReadsUnmapped Fastx
##feature_counts_output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts"

##~/programs/miniconda3/bin/featureCounts \
##-T $num_threads \
##-t $t_value \
##-g $g_value \
##-a $gtf_file_h \
##-o $feature_counts_output_dir"/"$base_name".hg38.counts.txt" \
##$output_dir"/"$base_name".hg38.Aligned.sortedByCoord.out.bam"

##Can do a second round of alignment here (e.g. subtractive alignment)
##$star_path \
##--runThreadN 8 \
##--genomeDir $ref_dir_v \
##--readFilesIn $output_dir"/"$base_name".hUnmapped.out.mate1" \
##--outFileNamePrefix $output_dir"/"$base_name."v" \
##--outSAMtype BAM SortedByCoordinate 
##
##~/programs/miniconda3/bin/featureCounts \
##-T 8 \
##-t exon \
##-g gene_id \
##-a $gtf_file \
##-o $output_dir"/"$base_name".v.counts.txt" \
##$output_dir"/"$base_name".v.Aligned.sortedByCoord.out.bam"

done < /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/files/vzv_rna_files.txt

fi

############################### Single end alignment: virus #################################################################
#######################################################################################################################

if [ $mode == "single_end_v" ]
then

##output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv_star_and_counts_02212021"
##output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv_all_exon_star_and_counts_02242021"
##output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_and_vzv_star_and_counts_02242021"
output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_and_vzv_star_cellranger_and_counts_02242021"
##ref_dir_v="/slipstream/home/mmariani/references/vzv/star"
##ref_dir_v="/slipstream/home/mmariani/references/vzv/star_all_exon"
ref_dir_v="/slipstream/home/mmariani/references/vzv/vzv_depledge/star"
##gtf_file_v="/slipstream/home/mmariani/references/vzv/vzv.gtf"
##gtf_file_v="/slipstream/home/mmariani/references/vzv/vzv_all_exon.gtf"
gtf_file_v="/slipstream/home/mmariani/references/vzv/vzv_depledge/vzv_depledge_adjusted_cellranger.gtf"

while IFS=$'\t' read -r -a my_array
do

cd $output_dir

base_name=$(basename "${my_array[0]}" ".fastq.gz")

$star_path \
--runThreadN $num_threads \
--genomeDir $ref_dir_v \
--genomeSAindexNbases 7 \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript transcript_id \
--sjdbGTFtagExonParentGene gene_id \
--readFilesIn "${my_array[0]}" \
--readFilesCommand zcat \
--outFileNamePrefix $output_dir"/"$base_name".vzv." \
> $base_name".log.txt"

##--outSAMtype BAM SortedByCoordinate \

file_name=$output_dir"/"$base_name".vzv.Aligned.sortedByCoord.out.bam"
samtools index $file_name
samtools idxstats $file_name > $output_dir"/"$(basename $file_name ".bam")".idxstats"
samtools stats $file_name > $output_dir"/"$(basename $file_name ".bam")".stats"
samtools flagstat $file_name > $output_dir"/"$(basename $file_name ".bam")".flagstat"

##--outReadsUnmapped Fastx
##feature_counts_output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts"

##~/programs/miniconda3/bin/featureCounts \
##-T $num_threads \
##-t $t_value \
##-g $g_value \
##-a $gtf_file_h \
##-o $feature_counts_output_dir"/"$base_name".hg38.counts.txt" \
##$output_dir"/"$base_name".hg38.Aligned.sortedByCoord.out.bam"

##Can do a second round of alignment here (e.g. subtractive alignment)
##$star_path \
##--runThreadN 8 \
##--genomeDir $ref_dir_v \
##--readFilesIn $output_dir"/"$base_name".hUnmapped.out.mate1" \
##--outFileNamePrefix $output_dir"/"$base_name."v" \
##--outSAMtype BAM SortedByCoordinate 
##
##~/programs/miniconda3/bin/featureCounts \
##-T 8 \
##-t exon \
##-g gene_id \
##-a $gtf_file \
##-o $output_dir"/"$base_name".v.counts.txt" \
##$output_dir"/"$base_name".v.Aligned.sortedByCoord.out.bam"

done < /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/files/vzv_rna_files.txt

fi

######################################### Paired end alignment ########################################################
#######################################################################################################################

####Testing:
##while read -ra array
##do
##echo "${array[0]} ${array[1]}"
####echo $file1 $file2
##done < /slipstream/home/mmariani/projects/hhv6_time_course/files/hhv6_rna_files.txt

##if [ $mode == "paired_end" ]
##then
##
##while read -ra array
##do
##
##cd $output_dir
##
##base_name=$(basename "${array[0]}" ".fastq.gz")
##
##$star_path \
##--runThreadN $num_threads \
##--genomeDir $ref_dir_v \
##--readFilesIn "${array[0]}" "${array[1]}" \
##--readFilesCommand zcat \
##--outFileNamePrefix $output_dir"/"$base_name".pe.v." \
##--outSAMtype BAM SortedByCoordinate \
##> $base_name".log.txt"
##
####--outReadsUnmapped Fastx
##feature_counts_output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts"
##~/programs/miniconda3/bin/featureCounts \
##-T $num_threads \
##-t $t_value \
##-g $g_value \
##-a $gtf_file_v \
##-o $feature_counts_output_dir"/"$base_name".pe.v.counts.txt" \
##$output_dir"/"$base_name".pe.v.Aligned.sortedByCoord.out.bam"
##
####Can do a second round of alignment here (e.g. subtractive alignment)
####$star_path \
####--runThreadN 8 \
####--genomeDir $ref_dir_v \
####--readFilesIn $output_dir"/"$base_name".hUnmapped.out.mate1" \
####--outFileNamePrefix $output_dir"/"$base_name."v" \
####--outSAMtype BAM SortedByCoordinate 
####
####~/programs/miniconda3/bin/featureCounts \
####-T 8 \
####-t exon \
####-g gene_id \
####-a $gtf_file \
####-o $output_dir"/"$base_name".v.counts.txt" \
####$output_dir"/"$base_name".v.Aligned.sortedByCoord.out.bam"
##
##done < /slipstream/home/mmariani/projects/hhv6_time_course/files/hhv6_rna_files.txt
##
##fi


##########################################################################################################

#!/usr/bin/env bash

##Mike Mariani UVM 2019-2021

##11/09/2019 - 02/24/2021
##qsub -cwd -pe threads 20 -j yes -o /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/logs /slipstream/home/mmariani/scripts/vzv_interactions_scripts/vzv_rna_seq_feature_counts_mm.bash

##Let's run feature counts on the 
##Hg38 aligned STAR files
##featureCounts is part of the Subread package
##I installed it via bioconda

#http://bioinf.wehi.edu.au/featureCounts/

##input_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_alignment"
##output_dir="/slipstream/home/mmariani/projects/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_star_and_counts_11092019/hg38_star_featureCounts"
##gtf_path="/slipstream/home/mmariani/references/hg38_ucsc_prebuilt/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
##
##cd $input_dir
##
##/slipstream/home/mmariani/programs/miniconda3/bin/featureCounts \
##-T 20 \
##-t exon \
##-g gene_id \
##-a $gtf_path \
##-o $output_dir"/"vzv_mock_rna_re-seq_star_feature_counts_combined_11092019.txt \
##H4_R4_S4_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##H5_R7_S5_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##H6_R8_S6_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV1_R1_S1_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV2_R2_S2_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV3_R3_S3_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam

##input_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv_star_and_counts_02212021"
##output_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_vzv_star_and_counts_02212021"
##gtf_path="/slipstream/home/mmariani/references/vzv/vzv.gtf"
##
##cd $input_dir
##
##/slipstream/home/mmariani/programs/miniconda3/bin/featureCounts \
##-T 20 \
##-t CDS \
##-g gene_id \
##-a $gtf_path \
##-o $output_dir"/"vzv_mock_rna_re-seq_star_feature_counts_combined_02242021.txt \
##H4_R4_S4_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##H5_R7_S5_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##H6_R8_S6_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV1_R1_S1_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV2_R2_S2_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam \
##HV3_R3_S3_L002_R1_001.hg38.Aligned.sortedByCoord.out.bam
##
########################################## VIRUS #######################################################
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

input_dir="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_and_vzv_star_cellranger_and_counts_02242021"
vzv_path_1=$input_dir"/HV1_R1_S1_L002_R1_001.vzv.Aligned.out.sorted.bam"
vzv_path_2=$input_dir"/HV2_R2_S2_L002_R1_001.vzv.Aligned.out.sorted.bam"
vzv_path_3=$input_dir"/HV3_R3_S3_L002_R1_001.vzv.Aligned.out.sorted.bam"
mock_path_1=$input_dir"/H4_R4_S4_L002_R1_001.vzv.Aligned.out.sorted.bam"
mock_path_2=$input_dir"/H5_R7_S5_L002_R1_001.vzv.Aligned.out.sorted.bam"
mock_path_3=$input_dir"/H6_R8_S6_L002_R1_001.vzv.Aligned.out.sorted.bam"

genes_file="/slipstream/home/mmariani/references/vzv_10x_depledge/FW__current_VZV_annotation/gtf_dumas_adjusted_cellranger.gtf"

cd /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_rna_seq_hfl_resequence_08072019/output_hg38_and_vzv_star_cellranger_and_counts_02242021

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $vzv_path_1 ".bam")".counts" \
$vzv_path_1

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $vzv_path_2 ".bam")".counts" \
$vzv_path_2

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $vzv_path_3 ".bam")".counts" \
$vzv_path_3

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $mock_path_1 ".bam")".counts" \
$mock_path_1

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $mock_path_2 ".bam")".counts" \
$mock_path_2

/slipstream/home/mmariani/programs/subread-1.6.4-Linux-x86_64/bin/featureCounts \
-s 2 \
-T 16 \
-a $genes_file \
-o $(basename $mock_path_3 ".bam")".counts" \
$mock_path_3

##htseq-count \
##-s reverse \
##$vzv_path_1 \
##$genes_file \
##> $(basename $vzv_path_1 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$vzv_path_2 \
##$genes_file \
##> $(basename $vzv_path_2 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$vzv_path_3 \
##$genes_file \
##> $(basename $vzv_path_3 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$mock_path_1 \
##$genes_file \
##> $(basename $mock_path_1 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$mock_path_2 \
##$genes_file \
##> $(basename $mock_path_2 ".bam")".counts"
##
##htseq-count \
##-s reverse \
##$mock_path_3 \
##$genes_file \
##> $(basename $mock_path_3 ".bam")".counts"




process sayHello {

    output:
        stdout

    script:
    """
    echo 'Hello World!'
    """
}

workflow {

    // emit a greeting
    sayHello()
}
