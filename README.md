# smallRNAseq

nextflow run nf-core/smrnaseq -profile singularity --input samplesheet.csv --genome 'GRCh37' --mirtrace_species 'hsa' --protocol 'custom' --outdir <OUTDIR> --mirgenedb TRUE 


########mirDB
nextflow run nf-core/smrnaseq -r 2.2.0 -profile singularity --input sample_sheet_k562_small_RNA.csv --mirtrace_species 'hsa' --protocol 'custom' --outdir k562_z8_ko --mirgenedb TRUE -params-file params.yml --mirgenedb_mature /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_MirGeneDB_mature.fa --mirgenedb_gff /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_miRNA_MirGeneDB.gff --mirgenedb_hairpin /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/hsa_MirGeneDB_precusor.fa -resume --igenomes_ignore True --fasta /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/GRCh38.primary_assembly.genome.fa --bt_indices /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/MirGeneDB_hsa/



#### mirBASE
