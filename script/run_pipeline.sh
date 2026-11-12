#!/bin/bash
set -e

mkdir -p data/raw data/processed script results

ref_file=data/raw/chr20.fa.gz
ref_seq20=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
ref_pro=data/processed/chr20.fa
BAM=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam
BAI=data/raw/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam.bai
region=chr20:10000000-12000000
BAM300=data/raw/NA12878_chr20_10-12Mb_300x.bam
BAM30=data/raw/NA12878_chr20_10-12Mb_30x.bam
collated=data/processed/NA12878_30x_collated.bam
r1=data/processed/NA12878_R1.fastq.gz
r2=data/processed/NA12878_R2.fastq.gz
singletons=data/processed/NA12878_singletons.fastq.gz
unsorted=data/processed/NA12878_aligned_unsorted.bam
QCreport=results/NA12878_R1_fastqc.html
fixmate=data/processed/NA12878_fixmate.bam
sorted=data/processed/NA12878_sorted.bam
markdup=data/processed/NA12878_markdup.bam
markdupindex=data/processed/NA12878_markdup.bam.bai
raw_vcf=results/NA12878_raw.vcf.gz
raw_filtered_vcf=results/NA12878_raw_filtered.vcf.gz
truth_url=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38
truth_vcf=data/raw/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
truth_bed=data/raw/HG001_GRCh38_1_22_v4.2.1_benchmark.bed
truth_region=data/processed/truth_region.vcf.gz
ours_region=data/processed/PASS_confident_region.vcf.gz
truth_norm=data/processed/truth_norm.vcf.gz
ours_norm=data/processed/ours_norm.vcf.gz
comparedir=results/compare

if [ ! -f $ref_file ]; then
 wget $ref_seq20 -P data/raw
fi

if [ !  -f $ref_pro ]; then
 zcat $ref_file > $ref_pro
fi

if [ ! -f ${ref_pro}.bwt ]; then
 bwa index $ref_pro
fi

if [ ! -f ${ref_pro}.fai ]; then
 samtools faidx $ref_pro
fi

if [ ! -f $BAI ]; then
  wget -c --tries=5 --waitretry=30 ${BAM}.bai -P data/raw
fi

if [ ! -f $BAM300 ]; then
 samtools view -b -X $BAM $BAI $region -o $BAM300
fi

if [ ! -f $BAM30 ]; then
 samtools view -b -s 42.1  $BAM300 -o $BAM30
fi

if [ ! -f $collated ]; then
 samtools collate -o $collated $BAM30
fi

if [ ! -f $r1 ]; then
 samtools fastq -1 $r1 -2 $r2 -s $singletons -0 /dev/null -n $collated
fi

if [ ! -f $QCreport ]; then
 fastqc $r1 $r2 --outdir results/ -t 2
fi

if [ ! -f $unsorted ]; then
 bwa mem -t 6 -R '@RG\tID:NA12878\tSM:NA12878\tPL:ILLUMINA' $ref_pro  $r1 $r2 2> results/bwa_mem.log | samtools view -b -o $unsorted
fi

if [ ! -f $fixmate ]; then
 samtools fixmate -m $unsorted  $fixmate
fi

if [ ! -f $sorted ]; then
 samtools sort -@ 4 -o $sorted $fixmate
fi

if [ ! -f $markdup ]; then
 samtools markdup  -s $sorted  $markdup
fi

if [ ! -f $markdupindex ]; then
 samtools index  $markdup
fi

if [ ! -f $raw_vcf ]; then
 bcftools mileup -f $ref_pro -r $region -q 20 -Q 20 -a AD,DP -Ou $markdup | bcftools view  -mv -Oz -o $raw_vcf
 bcftools index -f $raw_vcf
fi

if [ ! -f $raw_filtered_vcf ]; then
 bcftools filter -s lowQual -e 'QUAL<30 || INFO/DP<10' -Oz -o $raw_filtered_vcf $raw_vcf
 bcftools index -f $raw_filtered_vcf 
fi

for f in  HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi HG001_GRCh38_1_22_v4.2.1_benchmark.bed; do
 if [ ! -f data/raw/$f ]; then
  wget -c --tries=5 -P data/raw ${truth_url}/$f
 fi
done

if [ ! -f $conf_bed ]; then
 awk -v s0=$region_start -v e0=$region_end 'BEGIN{OFS="\t"} $1=="chr20" && $3>s0 && $2<e0 {s=($2<s0)?s0:$2; e=($3>e0)?e0:$3; print $1,s,e}' $truth_bed > $conf_bed
fi

if [ ! -f $truth_region ]; then
 bcftools view -R $conf_bed -Oz -o $truth_region $truth_vcf
 bcftools index -f $truth_region
fi

if [ ! -f $ours_region ]; then
 bcftools view -f PASS -R $conf_bed -Oz -o $ours_region $raw_filtered_vcf
 bcftools index -f $ours_region
fi

if [ ! -f $truth_norms ]; then
 bcftools norm -f $ref_pro -m -any -Oz -o $truth_norms $truth_region
 bcftools index -f $truth_norms
fi

if [ ! -f $ours_norms ]; then
 bcftools norm -f $ref_pro -m -any -Oz -o $ours_norms $ours_region
 bcftools index -f $ours_norms
fi

if [ ! -f $comparedir/0002.vcf ]; then
 bcftools isec -p $comparedir $truth_norms $ours_norms
fi

{

	echo -e "Type\tCorrect(TP)\tMissed(FN)\tMade_up(FP)"
for t in snps indels; do

 tp=$(bcftools view -H -v $t $comparedir/0002.vcf 2>/dev/null | wc -l)
 fn=$(bcftools view -H -v $t $comparedir/0000.vcf 2>/dev/null | wc -l)
 fp=$(bcftools view -H -v $t $comparedir/0001.vcf 2>/dev/null | wc -l)

 echo -e "$t\t$tp\t$fn\t$fp"
done
} > results/benchmark_summary.txt

cat results/benchmark_summary.txt
