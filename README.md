# smallRNAseq

#screen
# conda envi

srun --pty --partition=scu-cpu --mem=200G --cpus-per-task=20 bash -i 

# NXF_OPTS='-Xms1g -Xmx4g'

nextflow run nf-core/smrnaseq -profile singularity --input samplesheet.csv --genome 'GRCh37' --mirtrace_species 'hsa' --protocol 'custom' --outdir <OUTDIR> --mirgenedb TRUE 


########mirDB
nextflow run nf-core/smrnaseq -r 2.2.0 -profile singularity --input sample_sheet_k562_small_RNA.csv --mirtrace_species 'hsa' --protocol 'custom' --outdir k562_z8_ko --mirgenedb TRUE -params-file params.yml --mirgenedb_mature /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_MirGeneDB_mature.fa --mirgenedb_gff /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_miRNA_MirGeneDB.gff --mirgenedb_hairpin /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_MirGeneDB_precusor.fa -resume --igenomes_ignore True --fasta /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/GRCh38.primary_assembly.genome.fa --bt_indices /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/



#### mirBASE

nextflow run nf-core/smrnaseq -r 2.2.0 -profile singularity --input sample_sheet_k562_small_RNA.csv --mirtrace_species 'hsa' --protocol 'custom' --outdir k562_z8_ko -params-file params.yml --mature /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/mature_mirBASE_all.fa --mirna_gtf /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_mirBase.gff3 --hairpin /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hairpin_mirBASE_all.fa -resume --igenomes_ignore True --fasta /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/GRCh38.primary_assembly.genome.fa


#########    dev version is used 

(env_nf) [sor4003@scu-node038 run6_diff_expression_iN_ESC_H9_K562_zswim8_ko]$ nextflow run nf-core/smrnaseq -r dev -profile singularity --input sample_sheet_k562_small_RNA.csv --mirtrace_species 'hsa' --protocol 'custom' --outdir k562_z8_ko_MIRGeneDB -params-file params.yml --mature /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_MirGeneDB_mature.fa --mirna_gtf /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_miRNA_MirGeneDB.gff --hairpin /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_MirGeneDB_precusor.fa -resume --igenomes_ignore True --fasta /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/GRCh38.primary_assembly.genome.fa


########## extracting fasta file 
grep -A 1 "hsa" mature_mirBASE_all.fa | grep -v "^--$" > hsa_mature_mirBase.fa


######
(env_nf) [sor4003@scu-node023 Kleaveland-SR-14763_2023_07_28]$ nextflow run nf-core/smrnaseq -r dev -profile singularity --input trimmed_sample_sheet_k562_small_RNA.csv --mirtrace_species 'hsa' --protocol 'custom' --outdir k562_z8_ko_mirbase_hsa_subsetted_cutadat_timmed -params-file params.yml --mature /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_mature_mirBase.fa --mirna_gtf /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_mirBase.gff3 --hairpin /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_hairpin_mirBase.fa -resume --igenomes_ignore True --fasta /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/GRCh38.primary_assembly.genome.fa --trim_fastq false --clip_r1 0 --three_prime_clip_r1 0



####### cassava naming conviton 

file.fastq.gz ##### file.trim.fastq.gz >>> is wrong 


###### contamination check using human genome as --cdna 
nextflow run nf-core/smrnaseq -r dev -profile singularity --input trimmed_sample_sheet_k562_small_RNA.csv --mirtrace_species 'hsa' --protocol 'custom' --outdir k562_z8_ko_mirbase_hsa_subsetted_cutadat_trimmed_genome_as_contmination -params-file params.yml --mature /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_mature_mirBase.fa --mirna_gtf /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_mirBase.gff3 --hairpin /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_hairpin_mirBase.fa -resume --igenomes_ignore True --fasta /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/GRCh38.primary_assembly.genome.fa --trim_fastq false --clip_r1 0 --three_prime_clip_r1 0 --cdna /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/smallRNA_contamination_fasta/GRCh38.primary_assembly.genome.fa


##### mapping into human genome with bowtie index 
####https://scilifelab.github.io/courses/rnaseq/labs/smallRNA-lab
location-/home/sor4003/store_sor4003/2_star_genome_index_nexflow/small_RNA_genomes/smallRNA_contamination_fasta

bowtie2-build GRCh38.primary_assembly.genome.fa GRCh38.primary_assembly.genome_index

