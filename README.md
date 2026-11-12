# Germline Variant Calling Pipeline: NA12878 (chr20) Benchmarked Against GIAB

A Bash pipeline that calls small germline variants (SNVs and indels) from Illumina short reads and measures its own accuracy against the Genome in a Bottle (GIAB) truth set.

This was a learning project: the aim was to build and understand every step of a clinical-style variant calling workflow myself, from a reference genome and raw reads through to a benchmarked variant list.

---

## Overview

The pipeline analyses the GIAB reference sample **NA12878 (HG001)**. To fit the storage and compute of a laptop, it works on a **2 Mb region of chromosome 20 (chr20:10,000,000–12,000,000, GRCh38)** at **~30x** coverage.

Reads were taken from GIAB's 300x alignment, downsampled to 30x, converted back to FASTQ and **re-aligned from scratch** with BWA-MEM. After quality control, duplicate marking and variant calling with bcftools, the calls were compared with the GIAB v4.2.1 truth set inside GIAB's high-confidence regions.

**Result:** the pipeline found **3,059 of 3,086** true variants (**recall 99.1%**), and only **4** of its calls were false (**precision 99.9%**). SNVs reached **99.6% recall**; indels were harder at **96.5%**, mostly because many sit in repetitive DNA.

---

## Data

| Data | Source | Details |
|---|---|---|
| Reference genome | UCSC hg38 (`chr20.fa.gz`) | GRCh38 chromosome 20 only; 64,444,167 bp; uses `chr20` naming |
| Sequencing reads | GIAB, NA12878 / HG001, Illumina HiSeq 300x, aligned to GRCh38 (Novoalign) | Only the target region was fetched from the remote BAM using its index |
| Truth set | GIAB NISTv4.2.1, GRCh38 | Benchmark VCF + high-confidence regions (BED) for chromosomes 1–22 |

**NA12878** is one of the most thoroughly characterised human genomes. GIAB publishes a "truth set" of its variants, so any pipeline run on its reads can be scored, the standard way clinical and research labs validate variant calling pipelines.

**About the reference:** 0.8% of chr20 (499,586 bp) is unknown sequence (`N`), and 53.5% of the known sequence is soft-masked repeats (lowercase in the FASTA).

---

## Pipeline

```
Reference (UCSC chr20) ──► bwa index + samtools faidx ─────────────────┐
                                                                       │
GIAB 300x BAM (remote) ──► fetch chr20:10–12 Mb ──► keep 10% (~30x)    │
        ──► collate ──► FASTQ (R1, R2) ──► FastQC                      │
                              │                                        │
                              └──────────► bwa mem ◄───────────────────┘
                                              │
                         fixmate ──► sort ──► markdup ──► index
                                              │
                         bcftools mpileup + call ──► raw VCF
                                              │
                         bcftools filter (QUAL<30 or DP<10 → lowQual)
                                              │
     GIAB truth VCF + confident BED ──► restrict both to confident regions
                                              │
                         bcftools norm (left-align, split multiallelic)
                                              │
                         bcftools isec ──► TP / FN / FP ──► summary table
```

### Steps

| # | Step | Tool | Why |
|---|---|---|---|
| 1 | Download and index the reference | `wget`, `bwa index`, `samtools faidx` | BWA needs its index to place reads; bcftools needs `.fai` to look up reference bases |
| 2 | Fetch the target region from GIAB's BAM | `samtools view -X` with a local `.bai` | Downloads only ~350 MB instead of the full ~500 GB file |
| 3 | Downsample 300x → ~30x | `samtools view -s 42.1` | 30x is typical for clinical genomes; the seed (42) makes it reproducible |
| 4 | Convert back to paired FASTQ | `samtools collate`, `samtools fastq` | Discards GIAB's alignment so the reads can be re-aligned independently |
| 5 | Read quality control | `fastqc` | Checks read quality before alignment |
| 6 | Align to the reference | `bwa mem` (with read group `SM:NA12878`) | Places each read pair on chr20 |
| 7 | Prepare the alignment | `samtools fixmate -m`, `sort`, `markdup`, `index` | Adds mate information, sorts by position, marks PCR duplicates |
| 8 | Call variants | `bcftools mpileup -q 20 -Q 20 -a AD,DP` + `bcftools call -mv` | Ignores poorly mapped reads (MAPQ < 20) and low-quality bases |
| 9 | Filter | `bcftools filter -s lowQual -e 'QUAL<30 \|\| INFO/DP<10'` | Labels low-confidence calls instead of deleting them |
| 10 | Restrict to GIAB confident regions | `awk`, `bcftools view -R` | Scores only where GIAB is confident of the true answer |
| 11 | Normalise both call sets | `bcftools norm -m -any` | Left-aligns indels and splits multiallelic sites so equivalent variants match |
| 12 | Compare with the truth | `bcftools isec` | Sorts every variant into true positive, false negative or false positive |

---

## Results

### Read and alignment QC

| Metric | Value |
|---|---|
| Reads fetched (region, 300x) | 3,963,323 |
| Reads after downsampling (~30x) | 397,597 |
| Read pairs used for alignment | 198,257 (1,083 reads without a partner in the region were set aside) |
| Read length | 148 bp, paired-end |
| Reads mapped | **99.57%** |
| Properly paired | **98.92%** |
| Duplicate rate | **0.09%** (342 reads) |
| Split (supplementary) alignments | 489 |

The duplicate rate is low partly because random downsampling from 300x to 30x removes most copies of duplicated fragments.

### Variant calls (whole 2 Mb region)

| | Raw calls | After filtering (PASS) |
|---|---|---|
| SNVs | 2,851 | 2,812 |
| Indels | 590 | 563 |
| **Total** | **3,441** | **3,375** |

Of the 66 calls labelled `lowQual`, 62 were heterozygous. Weak or false calls usually affect only some reads at a position, so they tend to look heterozygous. Indels were filtered about three times more often than SNVs (4.6% vs 1.4%).

### Benchmark against GIAB v4.2.1

GIAB's high-confidence regions cover **1,952,263 bp (97.6%)** of the target region. Only PASS variants inside these regions were scored.

| Type | True positives | False negatives (missed) | False positives (false calls) | Precision | Recall | F1 |
|---|---|---|---|---|---|---|
| SNVs | 2,642 | 12 | 2 | 99.92% | 99.55% | 99.7% |
| Indels | 417 | 15 | 2 | 99.52% | 96.53% | 98.0% |
| **All** | **3,059** | **27** | **4** | **99.87%** | **99.13%** | **99.5%** |

- **Precision** = TP / (TP + FP): of the variants the pipeline reported, the share that are real.
- **Recall** = TP / (TP + FN): of the real variants, the share the pipeline found.

The summary table is written automatically to `results/benchmark_summary.txt`.

### Interpretation

- **SNVs** are called with very high accuracy on both measures.
- **Indels are the main weakness**: they make up only 14% of true variants but more than half of the missed ones. Many indels lie in tandem repeats (for example, a deletion of one `AAAT` unit in a long `AAAT` repeat), where reads are harder to align and repeat lengths harder to count.
- **Errors cluster in difficult DNA.** The non-confident 2.4% of the region contained 10% of the PASS calls (337 of 3,375). GIAB excludes these regions because the true answer is uncertain, and this is where short-read pipelines are least reliable.
- **Normalisation mattered.** `bcftools norm` rewrote 391 of the pipeline's records (mostly indels in repeats). Without it, many correct indels would have been scored as errors only because they were written differently from the truth set.

---

## Repository structure

```
variant-calling-pipeline/
├── script/
│   └── run_pipeline.sh          # the whole pipeline, from download to benchmark
├── results/
│   ├── benchmark_summary.txt    # TP / FN / FP table
│   ├── NA12878_R1_fastqc.html   # FastQC reports
│   ├── NA12878_R2_fastqc.html
│   └── bwa_mem.log              # alignment log
└── README.md
```

Raw and intermediate data (reference, BAM, FASTQ and VCF files, the GIAB truth set) are not stored in the repository because of their size. The script downloads or regenerates all of them.

---

## How to run

### Requirements

Tested with:

| Tool | Version |
|---|---|
| bwa | 0.7.19 |
| samtools | 1.24 |
| bcftools | 1.24 |
| FastQC | 0.13.0 |

```bash
conda create -n variant_calling -c conda-forge -c bioconda bwa samtools bcftools fastqc
conda activate variant_calling
```

### Run

From the project root folder:

```bash
bash script/run_pipeline.sh
```

Every step checks whether its output already exists and skips it if so, so the script can be safely re-run. A first run downloads about 520 MB (reference, region of the GIAB BAM and its index, truth set) and takes roughly 20–30 minutes, most of it fetching reads from the GIAB server.

---

## Limitations

- **Small region:** 2 Mb of one chromosome. Accuracy on the whole genome, especially in difficult regions, would be lower.
- **chr20-only reference:** reads that truly belong elsewhere can be forced onto chr20. Aligning to the full genome would be more accurate.
- **Reads were re-extracted from an existing alignment** rather than taken from raw sequencer output, so reads that GIAB's aligner could not place near this region were not included.
- **Simple hard filters** (QUAL and depth only); production pipelines use richer filters or models.
- **Comparison by `bcftools isec`** matches variants on position and alleles only, not genotype. Standard benchmarking tools (`hap.py`, `rtg vcfeval`) also check genotypes and handle complex indel representations; scores with those tools may be slightly lower.
- **A harmless header warning** (`MQ should be declared as Type=Float`) appears with this combination of bcftools/htslib versions and does not affect results.

---

## Possible next steps

- Benchmark with **hap.py** or **rtg vcfeval** for genotype-aware scores
- **Annotate** variants with VEP or SnpEff, with ClinVar and gnomAD
- Compare callers: **GATK HaplotypeCaller** and **DeepVariant** versus bcftools, especially for indels
- Use the kept 300x data to test how **coverage** (for example 10x, 30x, 60x) affects accuracy
- Scale up to **all of chromosome 20** or the whole genome
- Re-write the pipeline in **Nextflow** (or run nf-core/sarek) for parallel, resumable runs

---

## Acknowledgements

Data from the **Genome in a Bottle Consortium** (NIST) and the **UCSC Genome Browser**.
