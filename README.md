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
?smart_pca
pca_PC12 <- smart_pca(snp_data = "10RefsHybs.traw", sample_group = my_groups,missing_value = NA, scaling="none",program_svd = "bootSVD",sample_project = my_ancient,pc_project = c(1,2))
pca_PC34 <- smart_pca(snp_data = "10RefsHybs.traw", sample_group = my_groups,missing_value = NA, scaling="none",program_svd = "bootSVD",sample_project = my_ancient,pc_project = c(3,4))
# extract eigens and other data
eigen <- pca_PC12$pca.eigenvalues # extract PCA eigenvalues
load <- pca_PC12$pca$pca.snp_loadings # extract principal coefficients (SNP loadings)
pca_PC12 <- pca_PC12$pca.sample_coordinates
pca_PC34 <- pca_PC34$pca.sample_coordinates

eigenvalues <- c(130428.6877, 41236.67048,17697.919817,6742.161622,5714.661940, 5286.143612,5056.339411,4577.852233,4342.846216,4033.359383,3826.748154,3747.234795, 3711.567199,3638.707716,3536.139638,3471.33905,3412.336317,3378.152617,3350.169572, 3326.825417,3202.100253,3144.111055,3108.952277,3097.608666,3021.843927,2926.888843,2900.836395,2764.814801,2645.211651,2616.156115,2553.674668,2509.080031,2409.254483,2310.344832,2275.209167,2222.696024,2215.368914,2070.505338,1834.3073272,1734.286243,1632.9293045,1562.1737595,1403.0073478,1220.2198550)

total_variance <- sum(eigenvalues)
variance_explained <- sum(eigenvalues) / total_variance * 100
PC1 = (130428.6877 / total_variance) * 100 ## 41.02854
PC2 = (41236.67048 / total_variance) * 100 ## 12.97169
PC3 = (17697.919817 / total_variance) * 100 ## 5.567179
PC4 = (6742.161622 / total_variance) * 100 ## 2.120861

# save data for plotting
write.table(eigen,"Project_44_eigenvals.txt")
write.table(pca_PC12,"Project_PC12_eigenvec.txt")
write.table(pca_PC34,"Project_PC34_eigenvec.txt")
```
### Run smartSNP PCA NO Projection 
```
library(smartsnp)
my_groups <- c(1:140)
my_ancient <- c(1:96)
numSamples = nrow(read.table("10RefsHybs.fam"))
head(numSamples)
?smart_pca
pca <- smart_pca(snp_data = "10RefsHybs.traw", sample_group = my_groups,missing_value = NA, scaling="none",program_svd = "bootSVD")
# extract eigens and other data
eigen <- pca$pca.eigenvalues # extract PCA eigenvalues
load <- pca$pca$pca.snp_loadings # extract principal coefficients (SNP loadings)
pca <- pca$pca.sample_coordinates

eigenvalues <- c(120116.62124,20700.09474,8522.932741,7705.185363,5758.841694,3579.48966,2930.646962,2625.221834,2194.865296,2165.193960,2034.628699,1943.568941,1897.813399,1723.3798446,1696.0344557,1674.9075883,1586.191324,1570.2171842,1536.1059070,1514.3372783,1469.5757164,1447.0736982,1419.886123,1416.1533840,1389.5326206,1371.2621222,1349.3218639,1335.4450745,1322.458613,1299.4168746,1290.0103049,1273.0990080,1269.3338072,1248.6370930,1232.2238244,1222.6163262,1209.9564860,1194.0109487,1184.1151192,1174.8272562,1166.5631165,1160.5085892,1151.8567133,1143.6126727,1129.9704257,1117.7535049,1111.4437921,1098.9798330,1094.3025717,1079.4664752,1068.3041109,    1055.1213702,1052.593518,1042.1387423,1040.1859420,1031.0811637,1020.6920604,1014.9710644,1013.1223778,1006.0560868,999.0382129,988.0367368,982.3132510,978.6355836,972.6734185,966.8485105,960.2781054,949.9614376,945.775898,941.5899851,934.9257116,923.2682102,913.2771400,909.3573570,904.8394389,900.7939172,897.6750669,891.4419062,886.1399177,879.3298606,875.8289703,871.5094448,862.4720342,853.5201752,850.6473558,844.4597809,833.4441935,832.287579,830.6390172,821.2977640,809.505973,803.4395254,798.500496,794.7866491,785.5697809,771.7860407,765.9281747,759.8616216,754.1716589,749.5781806,743.3349818,741.1502003,732.2300744,727.4785802,716.2792300,714.8788150,704.505299,702.8504698,694.7064856,691.5827127,689.7288409,667.0318508,664.0722357,655.6487994,642.3291367,632.7382972,617.191738,610.5830642,591.1298110,574.2207439,570.0503956,566.876740,543.7035718,535.2596511,522.260545,509.8278047,492.0476524,484.9361971,458.0872655,445.1846043,434.1950663,428.9811555,408.1157529,401.1552636,380.1668673,327.5952664,310.2162058,302.9949617,242.4340907,194.5157110)

total_variance <- sum(eigenvalues)
variance_explained <- sum(eigenvalues) / total_variance * 100
PC1 = (120116.62124 / total_variance) * 100 ## 40.27223
PC2 = (20700.09474/ total_variance) * 100 ## 6.940246
PC3 = (8522.932741/ total_variance) * 100 ## 2.857535
PC4 = (7705.185363 / total_variance) * 100 ## 2.583364

# save data for plotting
write.table(eigen,"NoProject_140_eigenvals.txt")
write.table(pca,"NoProject_eigenvec.txt")
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
### Visualise in R as a stacked Plot 
```
library(ggplot2)
library(tidyverse)
data <- read.csv("Refs_plot.csv")
head(data)

data_long <- data %>% 
  pivot_longer(cols = FLO:VWESP, names_to = "condition", values_to = "value")
head(data_long)

# Custom colours
custom_colours <- c("FLO" = "#DB9D85", "ESP" = "#B1AF64", "VW" = "#87AEDF", "FLOVW" = "indianred4", "FLOESP" = "antiquewhite2", "VWESP" = "deepskyblue4")
data_long$X <- factor(data_long$X, levels = data$X)

png("Refs_AdmixFrog_stacked_barplot.png", width = 1754, height = 1240, res = 300)
# Stacked + percent plot
ggplot(data_long, aes(fill = condition, y = value, x = X)) + 
  geom_bar(position = "fill", stat = "identity",width = 0.9) +
  scale_fill_manual(values = custom_colours) +
  labs(title = "", y = "Average Genome Wide Posteriod Probability", x = "Sample ID") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

dev.off()
```
