# Ancestry_assignment

## Mitochondrial Analysis 

### Submit dSQ script to remove hard and soft calls
```
module load dSQ
dSQ --job-file clipping_dSQ.txt --mem-per-cpu 10g -t 1-00:00:00 -J pileup --mail-type ALL
sbatch dsq*.sh
```
### Submit dSQ script to remove sites with > 10 mismatches 
```
module load dSQ
dSQ --job-file bamtools_dSQ.txt --mem-per-cpu 10g -t 1-00:00:00 -J pileup --mail-type ALL
sbatch dsq*.sh
```
### Call Genotypes 
```
bcftools mpileup -f /gpfs/gibbs/pi/caccone/ao566/genome/lgeorge.mtgenome.fasta --threads 8 --annotate AD,DP -Ou --bam-list mismatches.bamlist | bcftools call --ploidy 1 -Ov --format-fields gq -m -o Galaponly_bcftools_mt_genome_ploidy1_rmclipping_maxmismatch10.vcf
```
Generate VCF or BCF containing genotype likelihoods for one or multiple alignment (BAM or CRAM) files.
* -f: reference sequence. Supplying this option will turn on left-alignment and normalization
* --threads: Use multithreading with INT worker threads
* --annotate: Comma-separated list of FORMAT and INFO tags to output.
  * AD: Allelic depth (Number=R,Type=Integer)
  * DP: Number of high-quality bases (Number=1,Type=Integer)
* --bam-list: List of input alignment files, one file per line

This command replaces the former bcftools view caller. Some of the original functionality has been temporarily lost in the process of transition under htslib, but will be added back on popular demand. 
* --ploidy: predefined ploidy
*  -m: alternative model for multiallelic and rare-variant calling designed to overcome known limitations in -c calling mode
*   -Ov: output type (v=vcf)
*   --format-fields: comma-separated list of FORMAT fields to output for each sample. Currently GQ and GP fields are supported. For convenience, the fields can be given as lower case letters.
*   -o: output vcf name 
```
## Check for missing data and remove individuals with > 50% missing
vcftools --vcf Galaponly_bcftools_mt_genome_ploidy1_rmclipping_maxmismatch10.vcf --missing-indv --out missing_indvs

bcftools view -S keep.txt Galaponly_bcftools_mt_genome_ploidy1_rmclipping_maxmismatch10.vcf > missingremoved_Galaponly_bcftools_mt_genome_ploidy1_rmclipping_maxmismatch10_Jan25check.vcf
```

### Convert vcf to nexus 
```
python vcf2phylip.py --input missingremoved_Galaponly_bcftools_mt_genome_ploidy1_rmclipping_maxmismatch10.vcf --phylip-disable –nexus -m 1
```

## Nuclear Analysis
### Genomic Chunks 
The Galapagos giant tortoise genome consists of very long chromosomes. For computational and time efficieny the genome can be split into genomic chunks of a defined size, prior to analysis. 
```
#Load BEDTools:
module load BEDTools/2.30.0-GCCcore-10.2.0

#Retrieve the first two columns of the indexed genome file:7
cut -f 1,2 /home/rg974/palmer_scratch/CheloAbing2/CheloAbing2.fasta.fai > CheloAbing2.genome.txt

#Create a BED file divided by 1-Mb windows:
bedtools makewindows -g CheloAbing2.genome.txt -w 1000000 > CheloAbing2.1MBfragments.bed

#Convert the BED file into the format required by ANGSD (just adding a unit to the second column):
awk '{print $1,$2+1,$3}' CheloAbing2.1MBfragments.bed | sed 's/ /\t/g' - > CheloAbing2.1MBfragments.txt

# in excel have to manually convert CheloAbing2.1MBfragments.txt to the correct format for angsd -rf flag
# save as sites_flag.txt on cluster
```
### Genotype calling in bcftools  
Since the genome is split into chunks we will use a wrapper script (saved as vcf_step1.py) which will create a text file that runs the genotype calling on each desired chunk, dsq is then used to run these as seperate jobs in parallel on the cluster. 
```
#!/usr/bin/python

import sys 

# get scaffold list from input file #
scaffold_f_name = sys.argv[1]
scaffold_f = open (scaffold_f_name, "r")
scaffold_list = scaffold_f.readlines()
scaffold_list = [x.strip() for x in scaffold_list]

with open("vcf_ref_SNPS_contigs.txt", "w") as out_f:
	for scaffold in scaffold_list:
		out_f.write(f"module load BCFtools BEDTools; bcftools mpileup -f /home/rg974/palmer_scratch/Updated_Floreana_Ref_Analyses/nuclear_analyses/Genome_fragments/CheloAbing2.fasta -r {scaffold} --annotate AD,DP -Ou --bam-list above4x_coverage.bamlist | bcftools call --skip-variants indels -mv -Ov -f GQ -o ./output/{scaffold}.vcf; bedtools intersect -header -a ./output/{scaffold}.vcf -b /gpfs/gibbs/pi/caccone/ao566/genome/repeatmasker.sorted.merged.complement.bed > ./output/{scaffold}.noREPS.vcf; rm ./output/{scaffold}.vcf; bgzip ./output/{scaffold}.noREPS.vcf;\n")
out_f.close()
```
Then run the text file in a seperate bash script: 
```
python vcf_step1.py /home/rg974/palmer_scratch/Updated_Floreana_Ref_Analyses/nuclear_analyses/Genome_fragments/sites_1millionbp.txt
module load dSQ
dSQ --job-file vcf_ref_SNPS_contigs.txt --mem-per-cpu 8g -t 1-0:00:00 -J pileup --mail-type ALL
```
And submit the dsq file. 

### Genotype filtering in vcftools

### Principal Components Analysis (PCA) in Plink 

### Principal Components Analysis (PCA) in EMU

### Private SNP Analysis in bcftools 

### Genome-wide heterozygosity in PIXY

### AdmixFrog

