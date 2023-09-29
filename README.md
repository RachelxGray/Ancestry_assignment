# Ancestry_assignment

## Mitochondrial Analysis 
Jensen & Miller 2018 - "Genotype calling was performed using SAMTOOLS mpileup (with the settings -gd -q 15 -t DP) with consensus genotypes determined using BCFTOOLS call (Li et al. 2009) and exported in vcf format. Filtering of the vcf files was done using VCFTOOLS (Danecek et al. 2011), requiring a minimum read depth of 6×."

bcftools mpileup -f lgeorge.mtgenome --threads 8 --annotate AD,DP --bam-list refs.bamlist | bcftools call --ploidy 1 -m -Ov --format-fields gq -o mt_genome_refs_allHist_Sep23_final.vcf

### bcftools mpileup
Generate VCF or BCF containing genotype likelihoods for one or multiple alignment (BAM or CRAM) files.
* -f: reference sequence. Supplying this option will turn on left-alignment and normalization
* --threads: Use multithreading with INT worker threads
* --annotate: Comma-separated list of FORMAT and INFO tags to output.
  * AD: Allelic depth (Number=R,Type=Integer)
  * DP: Number of high-quality bases (Number=1,Type=Integer)
* --bam-list: List of input alignment files, one file per line

### bcftools call
This command replaces the former bcftools view caller. Some of the original functionality has been temporarily lost in the process of transition under htslib, but will be added back on popular demand. 
* --ploidy: predefined ploidy
*  -m: alternative model for multiallelic and rare-variant calling designed to overcome known limitations in -c calling mode
*   -Ov: output type (v=vcf)
*   --format-fields: comma-separated list of FORMAT fields to output for each sample. Currently GQ and GP fields are supported. For convenience, the fields can be given as lower case letters.
*   -o: output vcf name 
