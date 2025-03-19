# Ancestry_assignment

## Mitochondrial Analysis 

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

## Nuclear Analysis - Genotype Calling
### Run genotype calling script: 
```
module load dSQ
dSQ --job-file bcftools_calling_100kb_dSQ.txt --mem-per-cpu 8g -t 1-0:00:00 -J pileup --mail-type ALL
```
### Make and Run vcftools script filtering for depth of 4 and minimum genotype quality of 18. 
```
for x in *.vcf.gz; do 
base_name=$(basename "$x")
name_without_extension="${base_name%.vcf.gz}"
echo "vcftools --gzvcf $x --minDP 4 --minGQ 18 --recode  --stdout | bgzip > ./depth4_GQ18/filtered_minDP4GQ18_${name_without_extension}.recode.vcf.gz" >>  vcftools_filterdepth_dSQtxt; 
done
module load dSQ
dSQ --job-file vcftools_filterdepth_dSQ.txt --mem-per-cpu 8g -t 1-0:00:00 -J pileup --mail-type ALL
```
### Concat the genome chunks
```
ls *.gz > bcftools_concat.txt
bcftools concat --file-list bcftools_concat.txt --threads 8 -O v -o - | bgzip -c > refshybs_minDP4GQ18_v2.1.vcf.gz
```

### Diversity Statistics (Nucleotide Diversity and Genome-Wide Heterozygosity)
## Filter concatted vcf 
```
vcftools --gzvcf refshybs_minDP4GQ18_v2.1.vcf.gz --max-missing 0.9 --recode --stdout  | bgzip -c > ./missing10/missing10_refshybs_minDP4GQ18_v2.1.vcf.gz

gunzip missing10_refshybs_minDP4GQ18_v2.1.vcf.gz
vcftools --vcf missing10_refshybs_minDP4GQ18_v2.1.vcf.gz --site-mean-depth --out missing10
bgzip missing10_refshybs_minDP4GQ18_v2.1.vcf.gz

Rscript mean.R

vcftools --gzvcf missing10_refshybs_minDP6GQ18_v2.1.vcf.gz --max-meanDP 14.02 --recode --stdout  | bgzip -c > ./max_DP/maxDP_missing10_refshybs_minDP6GQ18_v2.1.vcf.gz
```

## Nucleotide Diversity in PIXY 
```
pixy --stats pi \
--vcf maxDP_missing10_refshybs_minDP6GQ18_v2.1.vcf.gz \
--populations popfile.txt \
--window_size 10000 \
--n_cores 8
```
### Genome-wide heterozygosity
```
vcftools --gzvcf maxDP_missing10_refshybs_minDP6GQ18_v2.1.vcf.gz --het
```

## Principal Components Analysis (PCA) with smartSNP
### Filter concatted vcf 
```
# filter to only ESP, FLO, VW, Putative Hybs
bcftools view -S refshybs.txt refshybs_minDP4GQ18_v2.1.vcf.gz > ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
bgzip ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
# retain only biallelic sites 
vcftools --gzvcf ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf --max-alleles 2 --min-alleles 2 --recode --stdout | bgzip > ./biallelic/biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
# allow 10% missing data
vcftools --gzvcf biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz --max-missing 0.9 --recode --stdout  | bgzip -c > ./missing10/missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz
# max DP, one SD above mean
gunzip missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz
vcftools --vcf missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf --site-mean-depth --out missing10
bgzip missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz
Rscript mean.R
vcftools --gzvcf missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz.gz --max-meanDP 21.06 --recode --stdout  | bgzip -c > ./max_DP/maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz
# mac of 1
vcftools --gzvcf maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz --mac 1 --recode --stdout | bgzip -c > ./mac1/mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz
# find LD sites
gunzip mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz
awk 'BEGIN{OFS="\t"} !/#/ {sub(/\./, $1"_"$2, $3)}1' mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf  > annotated_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
bgzip mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
bgzip annotated_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
plink --vcf annotated_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz --indep-pairwise 50 5 0.5 --double-id --allow-extra-chr --out Keep_loci
sed 's/_/\t/' Keep_loci.prune.in > Keep_loci.txt
# filter to keep sites not in LD
vcftools --gzvcf mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz --positions Keep_loci.txt --recode --out LDpruned_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf 
bgzip LDpruned_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
```
### Prepare vcf for smartSNP format
```
gunzip LDpruned_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz
awk 'BEGIN{OFS="\t"} !/#/ {sub(/\./, $1"_"$2, $3)}1' LDpruned_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf  > annotated_LDpruned_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
bgzip LDpruned_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
bgzip annotated_LDpruned_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf
plink --vcf annotated_LDpruned_mac1_maxDP21_missing10_biallelic_ReducedRefsHybs_refshybs_minDP4GQ18_v2.1.vcf.gz --make-bed --out 10RefsHybs --allow-extra-chr 0 --double-id
plink --bfile 10RefsHybs --recode A-transpose --out 10RefsHybs --allow-extra-chr 0 --double-id
```
### Run smartSNP PCA with Projection 
```
library(smartsnp)
my_groups <- c(1:140)
my_ancient <- c(1:96)
numSamples = nrow(read.table("10RefsHybs.fam"))
head(numSamples)
pca <- smart_pca(snp_data = "RefsHybs_genotypeMatrix.traw", sample_group = my_groups, sample_project = my_ancient,missing_value = NA, scaling="none",pc_project = c(1, 2))
# extract eigens and other data
eigen <- pca$pca.eigenvalues # extract PCA eigenvalues
load <- pca$pca$pca.snp_loadings # extract principal coefficients (SNP loadings)
pca <- pca$pca.sample_coordinates
# save data for plotting
write.table(eigen,"eigenvals.txt")
write.table(coord,"eigenvec.txt")
write.table(load,"loads.txt")
```
### Visualise PCA


## Private SNP Analysis in bcftools 
### Hybrids 
```
filename="hybs.txt"

for x in $(cat "$filename")
do
bcftools view -s $x /home/rg974/palmer_scratch/Updated_Floreana_Ref_Analyses/Chapter3_analyses/bcftools_calling/version2.1/unequal_indep/depth4_GQ18/biallelic/missing10/max_DP/mac1/RefsHybs_alone/new_10missing/mac1_maxDP21_missing10_ReducedRefsHybs_biallelic_refshybs_minDP4GQ18_v2.1.vcf.gz  > $x.vcf
bgzip $x.vcf
tabix $x.vcf.gz
bcftools merge $x.vcf.gz /home/rg974/palmer_scratch/Updated_Floreana_Ref_Analyses/Chapter3_analyses/bcftools_calling/version2.1/unequal_indep/depth4_GQ18/biallelic/missing10/max_DP/mac1/RefsHybs_alone/new_10missing/Refs_mac1_maxDP21_missing10_ReducedRefsHybs_biallelic_refshybs_minDP4GQ18_v2.1.vcf.gz > "$x"_Refs_.vcf
bgzip "$x"_Refs_.vcf
tabix "$x"_Refs_.vcf
bcftools view -s $x -x "$x"_Refs_.vcf.gz > "$x"_private_snps.vcf
bcftools stats "$x"_private_snps.vcf | awk '/^SN/' > "$x"_psnps.txt
rm $x*.vcf
rm $x*.gz
rm $x*.tbi
printf "Completed, %s\n" $x
done
```
### References 
```
filename="refs.txt"

for x in $(cat "$filename")
do
bcftools view -s $x -x /home/rg974/palmer_scratch/Updated_Floreana_Ref_Analyses/Chapter3_analyses/bcftools_calling/version2.1/unequal_indep/depth4_GQ18/biallelic/missing10/max_DP/mac1/RefsHybs_alone/new_10missing/Refs_mac1_maxDP21_missing10_ReducedRefsHybs_biallelic_refshybs_minDP4GQ18_v2.1.vcf.gz > "$x"_private_snps.vcf
bcftools stats "$x"_private_snps.vcf | awk '/^SN/' > "$x"_psnps.txt
rm $x*.vcf
printf "Completed, %s\n" $x
done
```

## AdmixFrog
### Create Input File 
```
module load Python/3.8.6-GCCcore-10.2.0
module load XZ/5.2.5-GCCcore-10.2.0

admixfrog-ref --out AIMS_ESP_FLO_PBR.xz --vcf-ref /home/rg974/palmer_scratch/Updated_Floreana_Ref_Analyses/Chapter3_analyses/bcftools_calling/version2.1/unequal_indep/depth4_GQ18/biallelic/max_DP/AIMs_FLO_ESP_VW_maxDP_RefsHybs_biallelic_minDP4GQ18_v2.1.vcf.recode.vcf.gz \
    --state-file pops2.yaml \
    --states FLO ESP VW \
    --chroms chr01,chr02,chr03,chr04,chr05,chr06,chr07,chr08,chr09,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,scaffold0001,scaffold0027 

```

### Run admixfrog
```
module load dSQ
dSQ --job-file admix_frog.txt --mem-per-cpu 8g -t 1-0:00:00 -J pileup --mail-type ALL
```
### Extract bin results in R for genome-wide proportions 
```
dat <- read.csv("***.bin")
head(dat)

VW <- dat$VW
ESP<- dat$ESP
FLO <- dat$FLO
VWFLO <- dat$VWFLO
VWESP <- dat$VWESP
FLOESP <- dat$FLOESP


## VW
VW_results <- numeric(4481)
for(i in 1:length(VW)) {
  VW_results[i] <- VW[i]*500000
}
VW_SUM <- sum(VW_results)
VW_PROP <- (VW_SUM/(4481*500000))*100

## ESP
ESP_results <- numeric(4481)
for(i in 1:length(ESP)) {
  ESP_results[i] <- ESP[i]*500000
}
ESP_SUM <- sum(ESP_results)
ESP_PROP <- (ESP_SUM/(4481*500000))*100

## FLO
FLO_results <- numeric(4481)
for(i in 1:length(FLO)) {
  FLO_results[i] <- FLO[i]*500000
}
FLO_SUM <- sum(FLO_results)
FLO_PROP <- (FLO_SUM/(4481*500000))*100

## VWESP
VWESP_results <- numeric(4481)
for(i in 1:length(VWESP)) {
  VWESP_results[i] <- VWESP[i]*500000
}
VWESP_SUM <- sum(VWESP_results)
VWESP_PROP <- (VWESP_SUM/(4481*500000))*100

## VWFLO
VWFLO_results <- numeric(4481)
for(i in 1:length(VWFLO)) {
  VWFLO_results[i] <- VWFLO[i]*500000
}
VWFLO_SUM <- sum(VWFLO_results)
VWFLO_PROP <- (VWFLO_SUM/(4481*500000))*100

## ESPFLO
FLOESP_results <- numeric(4481)
for(i in 1:length(FLOESP)) {
  FLOESP_results[i] <- FLOESP[i]*500000
}
FLOESP_SUM <- sum(FLOESP_results)
FLOESP_PROP <- (FLOESP_SUM/(4481*500000))*100

check <- VW_PROP + ESP_PROP + FLO_PROP + VWESP_PROP + VWFLO_PROP + FLOESP_PROP

VW_PROP
ESP_PROP
FLO_PROP
VWESP_PROP
VWFLO_PROP
FLOESP_PROP
```
