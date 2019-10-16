## This script align the raw fastq reads to reference genome usin BWA MEM

#!/bin/bash

## Run bwa mem
/apps/bwa/0.7.9a/bwa mem \
-t 8 \
$path/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
$path/${sample}.R1.fq \
$path/${sample}.R2.fq \
> $path/${sample}.aln.sam

## Convert sams to bams
/apps/samtools/1.2/bin/samtools view -bS \
$path/${sample}.aln.sam \
> $path/${sample}.aln.bam

## Sort and indexing bam file
/apps/samtools/1.2/bin/samtools sort -m 16000000000 \
$path/${sample}.aln.bam \
$path/${sample}.aln.sorted

/apps/samtools/1.2/bin/samtools index \
$path/${sample}.aln.sorted.bam

## Fix read group information
/apps/java/jdk1.6.0_45/bin/java -Xmx64g \
-jar /apps/picard-tools/1.130/picard.jar AddOrReplaceReadGroups \
I=$path/${sample}.aln.sorted.bam \
O=$path/${sample}.aln.sorted.RGfixed.bam \
RGID=${sample} \
RGLB=${sample} \
RGPL=illumina \
RGPU=NNNNNNNN \
RGSM=${sample} \
SORT_ORDER=coordinate \
CREATE_INDEX=true \
TMP_DIR=./ \
VALIDATION_STRINGENCY=LENIENT

## Remove duplicates
/apps/java/jdk1.6.0_45/bin/java -Xmx64g \
-jar /apps/picard-tools/1.130/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=true \
INPUT=$path/${sample}.aln.sorted.RGfixed.bam \
OUTPUT=$path/${sample}.aln.sorted.RGfixed.rmdup.bam \
METRICS_FILE=$path/${sample}.metrics \
ASSUME_SORTED=true \
CREATE_INDEX=true \
TMP_DIR=./ \
VALIDATION_STRINGENCY=LENIENT

## Reorder bam file
java -Xmx64G -jar /apps/picard-tools/1.130/picard.jar \
ReorderSam \
I=$path/${sample}.aln.sorted.RGfixed.rmdup.bam \
O=$path/${sample}.aln.sorted.RGfixed.rmdup_reordered.bam \
R=$path/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa 

/apps/samtools/1.2/bin/samtools index \
$path/${sample}.aln.sorted.RGfixed.rmdup_reordered.bam

## Run realigner target creator
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $path/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path/${sample}.aln.sorted.RGfixed.rmdup_reordered.bam \
-known $path/1000G_indels_for_realignment.b37.vcf \
-o $path/${sample}.forIndelRealigner.intervals \
-L $path/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-filterRNC \
-filterMBQ  \
-filterNoBases

## Run Indel aligner
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R $path/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path/${sample}.aln.sorted.RGfixed.rmdup_reordered.bam \
-targetIntervals $path/${sample}.forIndelRealigner.intervals \
-known $path/1000G_indels_for_realignment.b37.vcf \
-L $path/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-o $path/${sample}.forIndelRealigner.bam

## Run picard Sortsam
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/picard-tools/1.130/picard.jar SortSam \
VALIDATION_STRINGENCY=LENIENT \
MAX_RECORDS_IN_RAM= 10000000 \
SORT_ORDER=coordinate \
I=$path/${sample}.forIndelRealigner.bam \
O=$path/${sample}.forIndelRealigner.sorted.bam \
CREATE_INDEX=true \
TMP_DIR= $path

## Run BaseRecalibrator
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $path/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path/${sample}.forIndelRealigner.sorted.bam \
-o $path/${sample}.forIndelRealigner.sorted.bam.table \
-L $path/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-filterRNC -filterMBQ  -filterNoBases \
-knownSites $path/dbsnp_132.b37.vcf \
-knownSites $path/1000G_indels_for_realignment.b37.vcf

## Run PrintReads
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T PrintReads \
-R $path/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path/${sample}.forIndelRealigner.sorted.bam \
-o $path/${sample}.recaliberated.bam \
-BQSR $path/${sample}.forIndelRealigner.sorted.bam.table \
-L $path/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-filterRNC -filterMBQ  -filterNoBases

## run DepthOfCoverage
/apps/java/jdk1.7.0_80/bin/java -Xmx64G -jar /apps/gatk3/3.4-0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R $path/Homo_sapiens.GRCh37.75.dna.primary_assembly_reordered.fa \
-I $path/${sample}.recaliberated.bam \
-o $path/${sample}.DepthOfCoverage \
-L $path/SureSelect_Exome_Capture_All_Tracks_v5_161011_mod.bed \
-ct 4 -ct 6 -ct 10 \
-pt readgroup

## Get WGS metrics
/apps/java/jdk1.6.0_45/bin/java -Xmx64g \
-jar /apps/picard-tools/1.130/picard.jar CollectWgsMetrics \
INPUT=$path/${sample}.recaliberated.bam \
OUTPUT=$path/${sample}.WGSmetrics \
REFERENCE_SEQUENCE=$path/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
TMP_DIR=./ \
MAX_RECORDS_IN_RAM=5000000 \
VALIDATION_STRINGENCY=LENIENT \
INCLUDE_BQ_HISTOGRAM=true
