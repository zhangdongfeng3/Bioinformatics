nohup bwa index -p chr1 Z



nohup bwa mem -t 20 /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B S1_1.clean.fq.gz S1_2.clean.fq.gz | samtools view -h -bF 12 -o S1.bam &
nohup bwa mem -t 20 /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B S2_1.clean.fq.gz S2_2.clean.fq.gz | samtools view -h -bF 12 -o S2.bam &

nohup STAR --runThreadN 10 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn S1_1.clean.fq.gz S1_2.clean.fq.gz --outFileNamePrefix S1 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 --readFilesCommand zcat &
nohup STAR --runThreadN 10 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn S2_1.clean.fq.gz S2_2.clean.fq.gz --outFileNamePrefix S2 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 --readFilesCommand zcat &

nohup samtools sort -@ 20 S1.bam > S1.sort.bam &
gatk ValidateSamFile -I S1.bam -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa
nohup gatk SortSam --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" -I S1.bam -O S1.sort.bam -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa -SO coordinate --CREATE_INDEX true &
nohup gatk SortSam --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" -I S2.bam -O S2.sort.bam -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa -SO coordinate --CREATE_INDEX true &


nohup samtools sort -@ 20 S2.bam > S2.sort.bam &

 

nohup gatk AddOrReplaceReadGroups --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" -I S1.sort.bam -O S1.sort.add.bam -ID 1 -LB group1 -PL Illumina -PU run0 -SM S1 &

nohup gatk AddOrReplaceReadGroups -I S2.sort.bam -O S2.sort.add.bam -ID 2 -LB group2 -PL Illumina -PU run1 -SM S2 &

nohup gatk MarkDuplicates --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" --REMOVE_DUPLICATES true -I S1.sort.add.bam -O S1.sort.add.rm.bam -M 1.metrics &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=12" --REMOVE_DUPLICATES true -I S2.sort.add.bam -O S2.sort.add.rm.bam -M 2.metrics &

*****
nohup gatk BaseRecalibrator
-R hg19.fa
-I ChrALL.100.sam.dedup.realn.07.bam
-knownSites dbsnp_137.hg19.vcf
-o Chr.table

gatk GatherBQSRReports \
    -I ${sep=' -I ' input_bqsr_reports} \
    -O ${output_report_filename}

gatk ApplyBQSR \
    -R ${ref_fasta} \
    -I ${input_bam} \
    -O ${output_bam_basename}.bam \
    -bqsr ${recalibration_report} \
    --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --create-output-bam-md5 \
    --use-original-qualities
*******

samtools index S1.sort.add.rm.bam

samtools index S2.sort.add.rm.bam

nohup gatk --java-options "-Xmx4G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=16" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I S1.sort.add.rm.bam \
 --native-pair-hmm-threads 16 \
 -O S1.gvcf -ERC GVCF -stand-call-conf 10 &
 
nohup gatk --java-options "-Xmx4G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=16" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I S2.sort.add.bam \
 --native-pair-hmm-threads 16 \
 -O S2.gvcf -ERC GVCF -stand-call-conf 10 &


RNASEQ

STAR --runThreadN 5 --runMode genomeGenerate --genomeDir B73V5S --genomeFastaFiles B73V5B.fa --sjdbGTFfile Zea_mays.Zm-B73-REFERENCE-NAM-5.0.51.chr.gtf --sjdbOverhang 149

nohup STAR --runThreadN 10 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn A_1.clean.fq.gz  A_2.clean.fq.gz --outFileNamePrefix A --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat &

nohup STAR --runThreadN 5 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn B_1.clean.fq.gz  B_2.clean.fq.gz --outFileNamePrefix B --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat &

nohup STAR --runThreadN 5 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn C_1.clean.fq.gz  C_2.clean.fq.gz --outFileNamePrefix C --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand gunzip &

nohup STAR --runThreadN 5 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn D_1.clean.fq.gz  D_2.clean.fq.gz --outFileNamePrefix D --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand gunzip &

nohup STAR --runThreadN 5 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn E_1.clean.fq.gz  E_2.clean.fq.gz --outFileNamePrefix E --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat &

bwa mem -t 20 /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B A_1.clean.fq.gz  A_2.clean.fq.gz | samtools view -h -bF 12 -o A.bam

STAR --twopassMode Basic 
--quantMode TranscriptomeSAM GeneCounts 
--runThreadN 6 
--genomeDir index_dir 
--alignIntronMin 20 
--alignIntronMax 50000 
--outSAMtype BAM SortedByCoordinate 
--sjdbOverhang 149 
--outSAMattrRGline ID:sample SM:sample PL:ILLUMINA 
--outFilterMismatchNmax 2 
--outSJfilterReads Unique 
--outSAMmultNmax 1 
--outFileNamePrefix out_prefix 
--outSAMmapqUnique 60 
--readFilesCommand gunzip -c 
--readFilesIn seq1.fq.gz seq2.fq.gz


nohup gatk AddOrReplaceReadGroups -I clean/AAligned.sortedByCoord.out.bam -O A.sort.add.bam -ID 3 -LB group3 -PL Illumina -PU run0 -SM A &

nohup gatk AddOrReplaceReadGroups -I clean/BAligned.sortedByCoord.out.bam -O B.sort.add.bam -ID 4 -LB group4 -PL Illumina -PU run0 -SM B &

nohup gatk AddOrReplaceReadGroups -I clean/CAligned.sortedByCoord.out.bam -O C.sort.add.bam -ID 5 -LB group5 -PL Illumina -PU run0 -SM C &

nohup gatk AddOrReplaceReadGroups -I clean/DAligned.sortedByCoord.out.bam -O D.sort.add.bam -ID 6 -LB group6 -PL Illumina -PU run0 -SM D &

nohup gatk AddOrReplaceReadGroups -I clean/EAligned.sortedByCoord.out.bam -O E.sort.add.bam -ID 7 -LB group7 -PL Illumina -PU run0 -SM E &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=2" --REMOVE_DUPLICATES true -I A.sort.add.bam -O A.sort.add.rm.bam -M 1.metrics &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=2" --REMOVE_DUPLICATES true -I B.sort.add.bam -O B.sort.add.rm.bam -M 2.metrics &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=2" --REMOVE_DUPLICATES true -I C.sort.add.bam -O C.sort.add.rm.bam -M 3.metrics &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=2" --REMOVE_DUPLICATES true -I D.sort.add.bam -O D.sort.add.rm.bam -M 4.metrics &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=2" --REMOVE_DUPLICATES true -I E.sort.add.bam -O E.sort.add.rm.bam -M 5.metrics &

nohup gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" SplitNCigarReads \
	-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
	-I B.sort.add.bam \
	-O B.sort.add.bam_split.bam &
	
nohup gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" SplitNCigarReads \
	-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
	-I C.sort.add.bam \
	-O C.sort.add.bam_split.bam &
	
nohup gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" SplitNCigarReads \
	-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
	-I D.sort.add.bam \
	-O D.sort.add.bam_split.bam &
nohup gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" SplitNCigarReads \
	-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
	-I E.sort.add.bam \
	-O E.sort.add.bam_split.bam &

nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I A.sort.add_split.bam \
 --native-pair-hmm-threads 2 \
 -O A.gvcf -ERC GVCF -stand-call-conf 10 &
 
nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I B.sort.add.bam_split.bam \
 --native-pair-hmm-threads 2 \
 -O B.gvcf -ERC GVCF -stand-call-conf 10 &

nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I C.sort.add.bam_split.bam \
 --native-pair-hmm-threads 2 \
 -O C.gvcf -ERC GVCF -stand-call-conf 10 &

nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I D.sort.add.bam_split.bam \
 --native-pair-hmm-threads 2 \
 -O D.gvcf -ERC GVCF -stand-call-conf 10 &

nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I E.sort.add.bam_split.bam \
 --native-pair-hmm-threads 2 \
 -O E.gvcf -ERC GVCF -stand-call-conf 10 &



smallRNASEQ

nohup STAR --runThreadN 10 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn clean_data/N_4_clean_total.fa.gz --outFileNamePrefix N4 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat &

nohup STAR --runThreadN 5 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn clean_data/N_6_clean_total.fa.gz --outFileNamePrefix N6 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat &

nohup STAR --runThreadN 5 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn clean_data/N_11_clean_total.fa.gz --outFileNamePrefix N11 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat &

nohup STAR --runThreadN 5 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn clean_data/N_12_clean_total.fa.gz --outFileNamePrefix N12 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat &

nohup STAR --runThreadN 5 --genomeDir /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5S --readFilesIn clean_data/N_14_clean_total.fa.gz --outFileNamePrefix N14 --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 5 --quantMode TranscriptomeSAM GeneCounts --readFilesCommand zcat &

nohup gatk AddOrReplaceReadGroups -I N4Aligned.sortedByCoord.out.bam -O N4.sort.add.bam -ID 31 -LB group31 -PL Illumina -PU run0 -SM N4 &

nohup gatk AddOrReplaceReadGroups -I N6Aligned.sortedByCoord.out.bam -O N6.sort.add.bam -ID 41 -LB group41 -PL Illumina -PU run0 -SM N6 &

nohup gatk AddOrReplaceReadGroups -I N11Aligned.sortedByCoord.out.bam -O N11.sort.add.bam -ID 51 -LB group51 -PL Illumina -PU run0 -SM N11 &

nohup gatk AddOrReplaceReadGroups -I N12Aligned.sortedByCoord.out.bam -O N12.sort.add.bam -ID 61 -LB group61 -PL Illumina -PU run0 -SM N12 &

nohup gatk AddOrReplaceReadGroups -I N14Aligned.sortedByCoord.out.bam -O N14.sort.add.bam -ID 71 -LB group71 -PL Illumina -PU run0 -SM N14 &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=12" --REMOVE_DUPLICATES true -I N4.sort.add.bam -O N4.sort.add.rm.bam -M N1.metrics &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=12" --REMOVE_DUPLICATES true -I N6.sort.add.bam -O N6.sort.add.rm.bam -M N2.metrics &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=12" --REMOVE_DUPLICATES true -I N11.sort.add.bam -O N11.sort.add.rm.bam -M N3.metrics &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=12" --REMOVE_DUPLICATES true -I N12.sort.add.bam -O N12.sort.add.rm.bam -M N4.metrics &

nohup gatk MarkDuplicates --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=12" --REMOVE_DUPLICATES true -I N14.sort.add.bam -O N14.sort.add.rm.bam -M N5.metrics &

nohup gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" SplitNCigarReads \
	-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
	-I B.sort.add.bam \
	-O B.sort.add.bam_split.bam &
	
nohup gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" SplitNCigarReads \
	-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
	-I B.sort.add.bam \
	-O B.sort.add.bam_split.bam &

nohup gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" SplitNCigarReads \
	-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
	-I B.sort.add.bam \
	-O B.sort.add.bam_split.bam &

nohup gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" SplitNCigarReads \
	-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
	-I B.sort.add.bam \
	-O B.sort.add.bam_split.bam &

nohup gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" SplitNCigarReads \
	-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
	-I B.sort.add.bam \
	-O B.sort.add.bam_split.bam &

nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I N4.sort.add.bam \
 --native-pair-hmm-threads 2 \
 -O N4.gvcf -ERC GVCF -stand-call-conf 10 &
 
nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I N6.sort.add.bam \
 --native-pair-hmm-threads 2 \
 -O N6.gvcf -ERC GVCF -stand-call-conf 10 &

nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I N11.sort.add.bam \
 --native-pair-hmm-threads 2 \
 -O N11.gvcf -ERC GVCF -stand-call-conf 10 &

nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I N12.sort.add.bam \
 --native-pair-hmm-threads 2 \
 -O N12.gvcf -ERC GVCF -stand-call-conf 10 &

nohup gatk --java-options "-Xmx1G \
 -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
 -R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
 -I N14.sort.add.bam \
 --native-pair-hmm-threads 2 \
 -O N14.gvcf -ERC GVCF -stand-call-conf 10 &



nohup gatk --java-options "-Xmx4G \
-XX:+UseParallelGC -XX:ParallelGCThreads=16" CombineGVCFs \
-R /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
--variant /data3/yanglab3/chengchen/dengyang/Feng/chen/wgs/cleandata/S1/S1.gvcf \
--variant /data3/yanglab3/chengchen/dengyang/Feng/chen/wgs/cleandata/S2/S2.gvcf \
--variant A.gvcf --variant B.gvcf \
--variant C.gvcf --variant D.gvcf --variant E.gvcf \
--variant /data3/yanglab3/chengchen/dengyang/Feng/chen/smallRNA/N4.gvcf \
--variant /data3/yanglab3/chengchen/dengyang/Feng/chen/smallRNA/N6.gvcf \
--variant /data3/yanglab3/chengchen/dengyang/Feng/chen/smallRNA/N11.gvcf \
--variant /data3/yanglab3/chengchen/dengyang/Feng/chen/smallRNA/N12.gvcf \
--variant /data3/yanglab3/chengchen/dengyang/Feng/chen/smallRNA/N14.gvcf -O gatksnp.gvcf &

nohup gatk --java-options "-Xmx4G \
-XX:+UseParallelGC -XX:ParallelGCThreads=16" GenotypeGVCFs \
-R R498_Final_Version2.soft.fasta --variant gatksnp.gvcf \
-O gatksnp.vcf &


####gvcf files

nohup bcftools mpileup -O v --threads 20 \
-f /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa \
N4.sort.add.bam \
N6.sort.add.bam \
N11.sort.add.bam \
N12.sort.add.bam \
N14.sort.add.bam -o smallrna.vcf &
nohup bcftools call --threads 20 -c -v smallrna.vcf -o smallrnav.vcf &
# Carry over all INFO and FORMAT annotations except FORMAT/GT
    bcftools annotate -a src.bcf -c INFO smallrnav.vcf


nohup freebayes-parallel <(fasta_generate_regions.py /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa.fai 100000) \
40 -f /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa N4.sort.add.bam N6.sort.add.bam N11.sort.add.bam N12.sort.add.bam N14.sort.add.bam > Nfree.vcf &

nohup freebayes -f /data3/yanglab3/chengchen/dengyang/Feng/chen/ref/B73V5B.fa N4.sort.add.bam N6.sort.add.bam N11.sort.add.bam N12.sort.add.bam N14.sort.add.bam > Nfree.vcf &






























