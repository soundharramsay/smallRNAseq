# smallRNAseq conda activate bowtie2

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

bash script 
#!/bin/bash
#wtie2 index name and output directory
index_base_name="/home/sor4003/store_sor4003/2_star_genome_index_nexflow/small_RNA_genomes/smallRNA_contamination_fasta/GRCh38.primary_assembly.genome_index"
output_dir="5_K562_bowtie2_map_hsa_genome"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# List of input FASTQ files (replace these with your file names)
fastq_files=("trim_E10_1_Z8KO_k562_sample4_S4_L002_R1_001.fastq.gz"
"trim_E13_1_z8ko_sample5_S5_L002_R1_001.fastq.gz"
"trim_E13_2_z8ko_sample6_S6_L002_R1_001.fastq.gz"
"trim_E13_3_z8ko_sample7_S7_L002_R1_001.fastq.gz"
"trim_nt1_k562_sample1_S1_L002_R1_001.fastq.gz"
"trim_nt2_k562_sample2_S2_L002_R1_001.fastq.gz"
"trim_nt3_k562_sample3_S3_L002_R1_001.fastq.gz")

# Loop over each FASTQ file
for fastq_file in "${fastq_files[@]}"; do
    filename=$(basename "$fastq_file")
    sam_file="${output_dir}/${filename%.fastq.gz}.sam"
    bam_file="${output_dir}/${filename%.fastq.gz}.bam"
    sorted_bam_file="${output_dir}/${filename%.fastq.gz}.sorted.bam"


    ################################################################################################################ things to do in analysis 
    #### are zswim8 microRNA are isomirs ----- mature -- +1-mature +2-mature , -1-mature and -2-mature ----- only for zswim8 sensitive microRNA 
    #### to ur microRNA reference add ADAR microRNA 
    #### what is the mismatch your mapping allowing on mismatch reads 

    # Run Bowtie2 alignment
    bowtie2 -q --very-sensitive-local -x "$index_base_name" -U "$fastq_file" -S "$sam_file"

    # Convert SAM to BAM and sort
    samtools view -S -b "$sam_file" | samtools sort -o "$sorted_bam_file"

    # Optional: Index the sorted BAM file
    samtools index "$sorted_bam_file"
done

######### feature counts ------- includes reads with mismatch 

#!/bin/bash
#input SAM files directory
input_dir="/home/sor4003/store_sor4003/RNAseq_results_fastq/run7_k562_small_RNA_Seq/Kleaveland-SR-14763_2023_07_28/5_K562_bowtie2_map_hsa_genome"

# Define the desired output directory name
output_dir_name="my_miRNA_counts"

# Define the path to the miRBase GTF file
mirbase_gtf="/home/sor4003/store_sor4003/2_star_genome_index_nexflow/small_RNA_genomes/mirBASE/hsa_mirBase.gff3"
# Create the output directory
output_dir="/home/sor4003/store_sor4003/RNAseq_results_fastq/run7_k562_small_RNA_Seq/Kleaveland-SR-14763_2023_07_28/5_K562_bowtie2_map_hsa_genome/$output_dir_name"
mkdir -p "$output_dir"

# Loop through all SAM files in the input directory
for samfile in "$input_dir"/*.sam; do
    # Get the base name of the SAM file without extension
    base_name="$(basename "$samfile" .sam)"

    # Define the output file name
    outfile="$output_dir/$base_name.featureCounts.txt"

    # Run featureCounts
    featureCounts -t miRNA -g Name -O -s 1 -M -a "$mirbase_gtf" -o "$outfile" "$samfile"

    # Print a message indicating completion for this file
    echo "Counting completed for $samfile"
done


###################################################### Feature_count_no_mismatch#################################################

#!/bin/bash
#input SAM files directory
input_dir="/home/sor4003/store_sor4003/RNAseq_results_fastq/run7_k562_small_RNA_Seq/Kleaveland-SR-14763_2023_07_28/5_K562_bowtie2_map_hsa_genome"

# Define the desired output directory name
output_dir_name="my_miRNA_counts_no_mismtach_featurecount"

# Define the path to the miRBase GTF file
mirbase_gtf="/home/sor4003/store_sor4003/2_star_genome_index_nexflow/small_RNA_genomes/mirBASE/hsa_mirBase.gff3"
# Create the output directory
output_dir="/home/sor4003/store_sor4003/RNAseq_results_fastq/run7_k562_small_RNA_Seq/Kleaveland-SR-14763_2023_07_28/5_K562_bowtie2_map_hsa_genome/my_miRNA_counts_no_mismtach_featurecount/$output_dir_name"
mkdir -p "$output_dir"

# Loop through all SAM files in the input directory
for samfile in "$input_dir"/*.sam; do
    # Get the base name of the SAM file without extension
    base_name="$(basename "$samfile" .sam)"

    # Define the output file name
    outfile="$output_dir/$base_name.featureCounts.txt"

    # Run featureCounts
    featureCounts -t miRNA -g Name -O -s 1 -M -a "$mirbase_gtf" -o "$outfile" "$samfile"

    # Print a message indicating completion for this file
    echo "Counting completed for $samfile"
done








#### remove first line from text file 
for file in E10_1_Z8KO_k562_sample4_S4_L002_R1_001.featureCounts.txt E13_3_z8ko_sample7_S7_L002_R1_001.featureCounts.txt nt3_k562_sample3_S3_L002_R1_001.featureCounts.txt E13_1_z8ko_sample5_S5_L002_R1_001.featureCounts.txt nt1_k562_sample1_S1_L002_R1_001.featureCounts.txt E13_2_z8ko_sample6_S6_L002_R1_001.featureCounts.txt nt2_k562_sample2_S2_L002_R1_001.featureCounts.txt; do
    sed -i '1d' "$file"
done

####### adding starting pattern of the file name as seven column name 
for file in E10_1_Z8KO_k562_sample4_S4_L002_R1_001.featureCounts.txt E13_3_z8ko_sample7_S7_L002_R1_001.featureCounts.txt nt3_k562_sample3_S3_L002_R1_001.featureCounts.txt E13_1_z8ko_sample5_S5_L002_R1_001.featureCounts.txt nt1_k562_sample1_S1_L002_R1_001.featureCounts.txt E13_2_z8ko_sample6_S6_L002_R1_001.featureCounts.txt nt2_k562_sample2_S2_L002_R1_001.featureCounts.txt; do     new_col_name="$(echo "$file" | cut -d'_' -f1,2)";     awk -v new_col_name="$new_col_name" 'BEGIN {FS=OFS="\t"} NR==1 {$7=new_col_name} {print}' "$file" > temp && mv temp "$file"; done




##### done on DEC18 

####### counting mir7-5p in original and trimmed.fastq files 

for file in trim_*.fastq.gz; do
    output_file="mir7_5p_18nt_subset_$(basename "$file" .fastq.gz).txt"
    zgrep -i TGGAAGACTAGTGATTTT "$file" > "$output_file"
done
#### extracting len and number of lines 

# Create a single output file
output_file="merged.output.count.txt"

# Iterate over all *.txt files
for file in *.txt; do
    # Extract total number of lines and length of lines
    total_lines=$(wc -l < "$file")
    line_length=$(awk '{ if (length > max) max = length } END { print max }' "$file")
    
    # Append the information to the merged output file
    echo "Input File: $file" >> "$output_file"
    echo "Total Number of Lines: $total_lines" >> "$output_file"
    echo "Maximum Line Length: $line_length" >> "$output_file"
    echo "" >> "$output_file"  # Add a separator between file entries
done
############# line length frequency #### good to execute as .sh script

output_file="merged_line_length_frequencies.txt"

for file in mir7_5p*.txt; do
    echo "=== $file ===" >> "$output_file"
    awk '{ lengths[length]++ } END { for (i in lengths) print i, lengths[i] }' "$file" >> "$output_file"
    echo "" >> "$output_file"  # Add a separator between file sections
done





