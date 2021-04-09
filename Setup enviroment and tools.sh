### Install necessary programs and set up environment


### ------------------------- IndexGenome -----------------------
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=8:mem=16gb

module load samtools

PATHtoPicard=/rds/general/applications/picard/2.6.0
PATHtoGENOME=/rds/general/user/ll3515/home/SCRATCH/ax4/Reference_genome_hg38/Gencode_GRCh38 # Genome which the RNA-seq was mapped to

samtools faidx GRCh38.primary_assembly.genome.fa
java -jar $PATHtoPicard/picard.jar CreateSequenceDictionary R= GRCh38.primary_assembly.genome.fa O= GRCh38.primary_assembly.genome.dict



### ------------------------- Download GATK java script -----------------------
https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk
#GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

# extract file
tar xjf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

### ------------------------- Download VEP -----------------------

# load your conda environment
module load anaconda3/personal
source activate apollo

#Download and install VEP, its dependencies:
conda install -qy -c conda-forge -c bioconda -c defaults ensembl-vep==102.0 htslib==1.10.2 bcftools==1.10.2 ucsc-liftover==377

#Download VEP's offline cache for GRCh38, and the reference FASTA:
mkdir -p $HOME/.vep/homo_sapiens/102_GRCh38/
rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-102/variation/indexed_vep_cache/homo_sapiens_vep_102_GRCh38.tar.gz $HOME/.vep/
tar -zxf $HOME/.vep/homo_sapiens_vep_102_GRCh38.tar.gz -C $HOME/.vep/
rsync -avr --progress rsync://ftp.ensembl.org/ensembl/pub/release-102/fasta/homo_sapiens/dna_index/ $HOME/.vep/homo_sapiens/102_GRCh38/

#Test running VEP in offline mode on a GRCh38 VCF:
vep --species homo_sapiens --assembly GRCh38 \
    --offline --no_progress --no_stats --sift b \
    --ccds --uniprot --hgvs --symbol --numbers --domains \
    --gene_phenotype --canonical --protein --biotype --tsl \
    --pubmed --variant_class --shift_hgvs 1 --check_existing \
    --total_length --allele_number --no_escape --xref_refseq \
    --failed 1 --vcf --minimal --flag_pick_allele \
    --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir $HOME/.vep \
    --fasta $HOME/.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
    --input_file homo_sapiens_GRCh38.vcf \
    --output_file homo_sapiens_GRCh38.vep.vcf --polyphen b \
    --af --af_1kg --af_esp --regulatory



### ------------------------- Download vcf2maf -----------------------
#https://github.com/mskcc/vcf2maf

export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*
perl vcf2maf.pl --man
perl maf2maf.pl --man

# Set any default paths and constants in vcf2maf.pl
# line 16: conda environment name to VEP binary
# line 17: GRCh38
# line 18: GRCh38

16: my ( $vep_path, $vep_data, $vep_forks, $buffer_size, $any_allele, $inhibit_vep, $online ) = ( "$ENV{HOME}/anaconda3/envs/apollo/bin", "$ENV{HOME}/.vep", 4, 5000, 0, 0, 0 );
17: my ( $ref_fasta, $filter_vcf ) = ( "$ENV{HOME}/.vep/homo_sapiens/102_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz", "" );
18: my ( $species, $ncbi_build, $cache_version, $maf_center, $retain_info, $retain_fmt, $min_hom_vaf, $max_filter_ac ) = ( "homo_sapiens", "GRCh38", "", ".", "", "", 0.7, 10 );

# load your conda environment
module load anaconda3/personal
source activate apollo

# Test the script
perl vcf2maf.pl --input-vcf tests/test_grch38.vcf --output-maf tests/test_grch38.vep.maf
