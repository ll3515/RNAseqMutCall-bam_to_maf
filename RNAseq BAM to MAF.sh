#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=90:ompthreads=90:mem=1080gb

#@ Load modules
module load samtools/1.3.1
module load java
module load picard
module load gatk/4.0
module load perl
module load anaconda3/personal
source activate apollo # your conda envionment


#@ Set PATH parameters
PATHtoPicard=/rds/general/applications/picard/2.6.0
PATHtoGENOME=/rds/general/user/ll3515/home/SCRATCH/ax4/Reference_genome_hg38/Gencode_GRCh38 # Genome which the RNA-seq was mapped to
PATHtoGATK=/rds/general/user/ll3515/home/SCRATCH/ax4/tools/GenomeAnalysisTK
PATHtoTMP=/rds/general/user/ll3515/home/DATA/MM_Nel/ # temporary directory for SplitNCigarReads intermediate outputs
PATHtoSAMPLE=/rds/general/user/ll3515/home/DATA/MM_Nel # path to RNA .bam directory
PATHtoVEP=/rds/general/user/ll3515/home/SCRATCH/ax4/tools/vcf2maf

###@ Input file >> STAR mapped bam file (Aligned.out.bam)

# Define input sample prefix:
S=a14869_15

### ------------------------- SortBam -----------------------

cd $PATHtoSAMPLE

samtools sort ${S}Aligned.out.bam -o ${S}Aligned.out.sorted.bam       ### Note input file


### ------------------------- MarkDuplicates -----------------------

java -XX:ParallelGCThreads=60 -jar $PATHtoPicard/picard.jar MarkDuplicates \
      I=${S}Aligned.out.sorted.bam \
      O=${S}marked_duplicates.bam \
      M=${S}marked_dup_metrics.txt \
      MAX_RECORDS_IN_RAM=500


### ------------------------- AddReadGroups -----------------------

java -jar $PATHtoPicard/picard.jar AddOrReplaceReadGroups \
      I=${S}marked_duplicates.bam \
      O=${S}marked_duplicates_readGroups.bam \
      RGID=ID \
      RGLB=libUnstranded \
      RGPL=illumina \
      RGPU=HJW2GDSXY \
      RGSM=MM

#      RGID (String)	Read Group ID Default value: 1. This option can be set to 'null' to clear the default value.
#      RGLB (String)	Read Group library Required.
#      RGPL (String)	Read Group platform (e.g. illumina, solid) Required.
#      RGPU (String)	Read Group platform unit (eg. run barcode) Required.
#      RGSM (String)	Read Group sample name Required.


### ------------------------- Index bam file -----------------------

samtools index ${S}marked_duplicates_readGroups.bam ${S}marked_duplicates_readGroups.bai


### ------------------------- SplitNCigarReads java -----------------------

java -Djava.io.tmpdir=$PATHtoTMP -jar $PATHtoGATK/GenomeAnalysisTK3.8.jar -T SplitNCigarReads \
                  -R $PATHtoGENOME/GRCh38.primary_assembly.genome.fa -I ${S}marked_duplicates_readGroups.bam \
                  -o ${S}split.bam \
                  -rf ReassignOneMappingQuality \
                  -RMQF 255 -RMQT 60 \
                  -U ALLOW_N_CIGAR_READS \
                  --log_to_file /rds/general/user/ll3515/home/DATA/MM_Nel/log.out


### ------------------------- Uniquely aligned reads -----------------------
#The variant calling is done on the uniquely aligned reads only in order to reduce the number of
#false positive variants called:

(samtools view -H ${S}split.bam; samtools view ${S}split.bam| grep -w 'NH:i:1') \
  | samtools view -Sb -  > ${S}split.uniq.bam

samtools index ${S}split.uniq.bam


### ------------------------- Variant calling -----------------------

java -jar $PATHtoGATK/GenomeAnalysisTK3.8.jar -T HaplotypeCaller \
                  -R $PATHtoGENOME/GRCh38.primary_assembly.genome.fa -I ${S}split.uniq.bam \
                  -dontUseSoftClippedBases \
                  -stand_call_conf 20.0 \
                  -o ${S}.gatk.vcf.gz

# unqip gz
gzip -d ${S}.gatk.vcf.gz

### ------------------------- Variant annotation with VEP -----------------------

# Change chromosome annotation from chr to number
awk '{gsub(/^chr/,""); print}' ${S}.gatk.vcf > ${S}.gatk.no_chr.vcf

#perl $PATHtoVEP/vcf2maf.pl --input-vcf ${S}.gatk.no_chr.vcf --output-maf ${S}.gatk.no_chr.vep.maf
perl $PATHtoVEP/vcf2maf_custom.pl --input-vcf ${S}.gatk.no_chr.vcf --output-maf ${S}.gatk.no_chr.vep.maf
