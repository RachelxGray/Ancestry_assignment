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





