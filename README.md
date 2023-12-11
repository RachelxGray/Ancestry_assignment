# Ancestry_assignment

## Mitochondrial Analysis 
Key papers: 
* 


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

## ANGSD Filtering
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
# save as sites_flag.txt
file="/home/rg974/palmer_scratch/ANGSD_filtering/Contig_approach/Smaller_Contigs/files/tests/sites_flag.txt"
while IFS= read -r line; do echo "module load angsd; angsd -b All_refs_hybs_hist.bamlist -gl 2 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -doGlf 2 -doBcf 1 -doIBS 1 -doCounts 1 -doCov 1 -makeMatrix 1 -doGeno 4 -doPost 1 -postCutoff 0.95 -rmTrans 1 -skipTriallelic 1 -out /home/rg974/palmer_scratch/ANGSD_filtering/Contig_approach/Smaller_Contigs/files/tranversions_output/Refs_${line} -nThreads 10 -minInd 100 -minMaf 0.05 -r ${line}"  >>  angsd_filter.txt; done < "$file"
module load dSQ
dsq --job-file angsd_filter.txt --mem-per-cpu 100g -t 1-00:00:00 -J angsd_filter --mail-type ALL
```
# ANGSD general 
* -b: bam files to be read 
* -r: list of sites to be analysed
 ## data filtering
* -uniqueOnly 1: uniquely mapping reads
* -remove_bads 1: keep reads not labelled as bad 
* -only_proper_pairs 1: considering only proper pairs
* -trim 0: without trimming
* -C 50: reduces the effect of reads with excessive mismatches,
* -baq 1: computes base alignment quality used to rule out false SNPs close to INDELS.
* -minInd 100 : use only sites with data from at least half? individuals 
* -setMinDepth 5: minimum total depth
* -setMaxDepth 5: maximum total depth
## Genotype likelihoods 
* -gl 1 : Estimate genotype likelihoods (1 = SAMtools method, preffered when data uncertainty is high)
* -doGlf 2: associated with gl. output beagle likelihood file (beagle.gz)
## Call Genotypes 
* -doGeno 4: Call genotypes (4 = write the called genotype directly: eg AA,AC etc)
* -doPost 1: estimate the posterior genotype probability based on the allele frequency as a prior
* -doMajorMinor 1: Infer major and minor from GL
* -postCutoff 0.95: associated with doPost. Call only a genotype with a posterior above this threshold. For instance, we can set as missing genotypes when their (highest) genotype posterior probability is below 0.95 (what you set this to depends on teh depth of samples) 
## Estimate Allele Frequencies (estimate how many copies of different alleles (two in case of biallelic variants) we observe in our sample (across all sequenced individuals). However with low depth data direct counting of individually assigned genotypes can lead to biased allele frequencies.) 
* doMaf 1: Frequency (fixed major and minor)
* doMajorMinor 4: Infer major and minor from GL. When working with multiple populations -doMajorMinor 1 is not suitable as you need to ensure that each site is polarised in the same way among different populations. Therefore, in these conditions, it is better to use option 4 or 5. Option 4 is reference allele as major.
* -ref: reference sequence
## SNP Calling 
May be interested in looking at allele frequencies only for sites that are actually variable in our sample. Therefore we want to perform a SNP calling.
There are two main ways to call SNPs using ANGSD with these options:
*  -minMaf         0.05      (Remove sites with MAF below)
*  -SNP_pval       1e-6       (Remove sites with a pvalue larger)
* -rmTrans 1: remove transitions
* -skipTriallelic 1: Remove sites with a pvalue lower
Therefore we can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.

# prepare the beagle file for ngsLD
zcat Refs_Contig3709:1-9486.mafs.gz  | head -n 1  > header_.txt

for f in *.beagle.gz ; do
    zcat $f | tail -n -1 | gzip >> $f.tail
done
ls -lR ./*.tail | wc -l
ls -lR ./*.beagle.gz | wc -l

zcat *.tail |cat > concat_all.beagle
wc -l concat_all.beagle

awk '$1!="marker"' concat_all.beagle >  filtererd_final_concat_all.beagle
wc -l filtererd_final_concat_all.beagle

#manually copy and paste header.txt as hearer of the filtererd_final_concat_all.beagle file 
gzip  filtererd_final_concat_all.beagle

# get into correct format for ngsLD 
zcat filtererd_final_concat_all.beagle.gz | cut -f 4- |sed '1d'| gzip  > formatted_filtererd_final_concat_all.beagle.gz

# prepare the mafs file for ngsLD
zcat Refs_Contig3709:1-9486.mafs.gz  | head -n 1  > header_mafs.txt

for f in *.mafs.gz ; do
    zcat $f | tail -n -1 | gzip >> $f.postail
done
ls -lR ./*.mafs.gz | wc -l
ls -lR ./*.postail | wc -l

zcat *.postail |cat > concat_all.pos

head concat_all.pos

awk '$1!="chromo"' concat_all.pos >  filtererd_final_concat_all.pos

zcat filtererd_final_concat_all.pos.gz | cut -f 1,2 | sed 's/:/_/g'|sed '1d'| gzip > formatted_filtererd_final_concat_all.pos.gz
zcat formatted_filtererd_final_concat_all.pos.gz | sort -k 1 -V | gzip > sorted_formatted_filtererd_final_concat_all.pos.gz


## Linkage Disequilibrium 
# Install dependencies 
module load miniconda 
conda activate satmools
wget https://ftp.gnu.org/gnu/gsl/gsl-2.6.tar.gz
tar -zxvf gsl-2.6.tar.gz
mkdir gsl_install_26
cd gsl_install_26
./configure --prefix=/home/rg974/palmer_scratch/Software/gsl_install_26
make
make check
make install
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/rg974/palmer_scratch/Software/gsl_install_26/lib
# Install ngSLD
git clone https://github.com/fgvieira/ngsLD.git
cd ngsLD
make
ngsLD/ngsLD

# Running ngsLD

/home/rg974/palmer_scratch/Software/ngsLD/ngsLD \
--geno /home/rg974/palmer_scratch/ANGSD_filtering/Contig_approach/Smaller_Contigs/files/output_reduced/analysis/formatted_filtererd_final_concat_all.beagle.gz \
--pos /home/rg974/palmer_scratch/ANGSD_filtering/Contig_approach/Smaller_Contigs/files/output_reduced/analysis/sorted_formatted_filtererd_final_concat_all.pos.gz \
--probs \
--n_ind 69 \
--max_kb_dist 0 \
--n_sites 5888 \
--n_threads 1 \
--out ANGSD_PCA_subsampled.ld





