---
layout: post
title: "Building an RNA-seq Pipeline for Rare Diseases on AWS HealthOmics"
date: 2025-08-28
categories: bioinformatics
tags: [bioinformatics, rnaseq, nextflow, aws, rare-disease, pipeline]
---

Over the past couple of years at the [UCLA Nelson Lab](https://nelsonlab.ucla.edu/), I built and maintained a comprehensive RNA-seq analysis pipeline designed specifically for rare disease diagnostics. The pipeline, [nl-rna-seq_wf](https://github.com/uclanelsonlab/nl-rna-seq_wf), is written in Nextflow (DSL2) and runs on [AWS HealthOmics](https://aws.amazon.com/healthomics/) — Amazon's managed service for genomics workflows. In this post I want to walk through the architecture, explain why each tool was chosen, and share some lessons learned from running production genomics in the cloud.

## Why AWS HealthOmics?

AWS HealthOmics is a purpose-built service for storing, querying, and analyzing genomic data at scale. We chose it for several reasons:

- **Managed infrastructure**: No need to provision or manage EC2 instances, auto-scaling groups, or job schedulers. HealthOmics handles compute orchestration for you.
- **Native Nextflow support**: HealthOmics supports Nextflow workflows natively as a "private workflow", so we could keep the same DSL2 codebase we developed locally.
- **ECR integration**: All Docker images are stored in AWS Elastic Container Registry (ECR), which HealthOmics pulls from directly, avoiding Docker Hub rate limits and ensuring reproducibility.
- **S3 storage**: Input FASTQ files, reference genomes, and STAR indices live in S3 buckets. Outputs are written to `/mnt/workflow/pubdir`, which HealthOmics automatically syncs back to S3.
- **HIPAA compliance**: Since we process clinical patient data from programs like the [Undiagnosed Diseases Network (UDN)](https://undiagnosed.hms.harvard.edu/), HealthOmics provides a compliant environment out of the box.

The main trade-off is that HealthOmics has its own constraints — fixed compute tiers, specific Nextflow version requirements, and some limitations on how processes can communicate. But for a production clinical pipeline, the reliability and compliance benefits far outweigh these.

## Pipeline Overview

The pipeline takes paired-end RNA-seq FASTQ files and produces a comprehensive set of outputs for rare disease analysis:

```
FASTQ R1/R2
    │
    ├── Fastp (QC + trimming)
    │     ├── BWA → rRNA contamination check
    │     └── BWA → Globin RNA contamination check
    │
    ├── STAR (splice-aware alignment)
    │     ├── Samtools Index
    │     └── Sambamba MarkDup
    │           ├── FeatureCounts (gene counts)
    │           ├── RNA-SeQC (QC metrics)
    │           ├── Qualimap (QC metrics)
    │           ├── IRFinder (intron retention)
    │           ├── Mosdepth (coverage analysis)
    │           │     └── Bedtools (CDS filtering)
    │           │           └── DeepVariant (variant calling)
    │           ├── BAM2SJ (splice junctions)
    │           └── Samtools CRAM (compression)
    │                 └── Mosdepth BED (XBP1 + MT coverage)
    │
    ├── Kallisto (transcript quantification)
    │
    └── MultiQC (aggregated report)
```

## The Tools and Why They're There

### Fastp — Quality Control and Trimming

[Fastp](https://github.com/OpenGene/fastp) is the first step: adapter trimming, quality filtering, and length filtering. We chose fastp over alternatives like Trimmomatic or Cutadapt because it's significantly faster (written in C++), produces clean HTML/JSON reports, and handles paired-end reads natively in a single pass. The JSON output feeds directly into MultiQC for aggregated reporting.

### BWA — Contamination Detection (rRNA and Globin RNA)

Before alignment, we check for two common sources of contamination in blood-derived RNA samples:

- **rRNA contamination**: Ribosomal RNA should be depleted during library prep, but incomplete depletion wastes sequencing capacity. We align trimmed reads against an rRNA reference using [BWA-MEM](https://github.com/lh3/bwa) and compute flagstats to get the percentage of reads mapping to rRNA.
- **Globin RNA contamination**: Blood samples are dominated by hemoglobin transcripts. If globin depletion wasn't performed (or failed), a large fraction of reads will be globin. Same approach: align against a globin reference and check flagstats.

High contamination rates in either category are a red flag that the sample may need to be re-prepped or that results should be interpreted with caution.

### STAR — Splice-Aware Alignment

[STAR](https://github.com/alexdobin/STAR) is the gold standard for RNA-seq alignment. It's splice-aware, meaning it can align reads that span exon-exon junctions — critical for RNA data. STAR also outputs gene-level read counts directly (`ReadsPerGene.out.tab`) and detects novel splice junctions.

One interesting aspect of our pipeline is that we maintain **six different STAR indices** for read lengths of 69, 75, 100, 120, 150, and 151 bp. This is because STAR's `sjdbOverhang` parameter should ideally be `read_length - 1`, and our lab processes samples from multiple sequencing platforms and protocols. The pipeline selects the appropriate index based on the detected read length, ensuring optimal alignment regardless of the input data source.

STAR is by far the most resource-intensive step, requiring 192 GB of RAM and 48 CPUs — driven by the size of the human genome index in memory.

### Sambamba — Duplicate Marking

[Sambamba](https://github.com/biod/sambamba) marks PCR duplicates in the aligned BAM. We chose Sambamba over Picard's MarkDuplicates because it's multi-threaded and significantly faster on large files, which matters when you're paying per minute on cloud compute. The marked BAM is the starting point for all downstream analyses.

### Subread FeatureCounts — Gene-Level Quantification

[FeatureCounts](https://subread.sourceforge.net/) from the Subread package assigns aligned reads to genomic features (genes) based on GENCODE annotations. It produces a count matrix that feeds into downstream differential expression tools like [DESeq2](https://bioconductor.org/packages/DESeq2/) or aberrant expression detection with [OUTRIDER](https://bioconductor.org/packages/OUTRIDER/). FeatureCounts is fast, handles multi-mapping reads well, and works at the exon level which gives us flexibility.

### Kallisto — Transcript-Level Quantification

[Kallisto](https://pachterlab.github.io/kallisto/) provides transcript-level abundance estimates using pseudoalignment — it doesn't perform traditional alignment at all, instead using k-mer matching against a transcriptome index. This makes it extremely fast and complementary to STAR+FeatureCounts. We use Kallisto's transcript-level quantification for analyses that need isoform resolution, while FeatureCounts handles gene-level counting.

### RNA-SeQC — Quality Assessment

[RNA-SeQC](https://github.com/getzlab/rnaseqc) from the Getz Lab at the Broad Institute computes a comprehensive set of quality metrics: mapping rates, rRNA rates, exonic/intronic/intergenic ratios, 3'/5' coverage bias, GC bias, and gene-level TPM values. These metrics are essential for identifying problematic samples before they contaminate downstream analyses, especially when building cohort-level models for rare disease detection.

### Qualimap — Additional QC

[Qualimap](http://qualimap.conesalab.org/) provides additional RNA-seq-specific quality metrics and visualizations, complementing RNA-SeQC. Having multiple QC tools gives us a more complete picture of sample quality.

### IRFinder — Intron Retention Detection

[IRFinder](https://github.com/williamritchie/IRFinder) detects intron retention events, where introns that should be spliced out are instead retained in the mature mRNA. Intron retention is increasingly recognized as a mechanism in rare diseases — it can lead to premature stop codons, nonsense-mediated decay, or altered protein function. This is particularly relevant for our work since aberrant splicing is one of the key things RNA-seq can reveal that exome/genome sequencing alone cannot.

### Mosdepth — Coverage Analysis

[Mosdepth](https://github.com/brentp/mosdepth) calculates sequencing depth and coverage quickly from BAM/CRAM files. We run it in two modes:

1. **Genome-wide per-base coverage**: This output feeds into Bedtools to identify coding regions with adequate coverage for variant calling.
2. **Targeted BED regions**: We specifically check coverage over the **XBP1** gene and **mitochondrial** regions. XBP1 is relevant because its unconventional splicing is a biomarker for certain conditions, and mitochondrial coverage is important for mitochondrial disease diagnostics — a significant subset of rare diseases.

### Bedtools — CDS Region Filtering

[Bedtools](https://bedtools.readthedocs.io/) takes the per-base coverage from Mosdepth and intersects it with GENCODE CDS (coding sequence) annotations, filtering for regions that meet a minimum coverage threshold. The resulting BED file defines the genomic regions where we have enough data to confidently call variants. This ensures DeepVariant only runs on well-covered coding regions, improving both accuracy and efficiency.

### DeepVariant — Variant Calling

[DeepVariant](https://github.com/google/deepvariant) is Google's deep learning-based variant caller. We use a **custom model trained on RNA-seq data** rather than the default WGS/WES models — this is important because RNA-seq has fundamentally different characteristics (splice junctions, allele-specific expression, variable coverage across exons). The pipeline accepts custom model files (`model_data`, `model_index`, `model_meta`, `model_info`) as parameters.

Variant calling from RNA-seq complements WGS/WES by providing expression-aware variant detection. A variant found in both DNA and RNA data, especially with biallelic expression, gives higher confidence for pathogenicity assessment.

### BAM2SJ — Splice Junction Analysis

BAM2SJ reconstructs splice junctions from the aligned BAM, producing a table of junction coordinates, read support, and strand information. This output is used downstream for aberrant splicing detection with tools like [FRASER](https://bioconductor.org/packages/FRASER/), which identifies samples with statistically unusual splicing patterns compared to a control cohort — a powerful approach for finding disease-causing splice variants.

### Samtools — CRAM Compression

After all analyses are complete, we convert the final BAM to [CRAM](https://www.ga4gh.org/cram/) format using samtools. CRAM files are typically 40-60% smaller than BAM files, which significantly reduces S3 storage costs when you're processing hundreds of samples.

### MultiQC — Aggregated Reporting

[MultiQC](https://multiqc.info/) collects outputs from Fastp, STAR, FeatureCounts, RNA-SeQC, Qualimap, Kallisto, and the contamination checks into a single interactive HTML report. This is invaluable for quickly assessing sample quality and identifying issues across a batch of samples.

## The Rare Disease Context

This pipeline was built to support the [Undiagnosed Diseases Network (UDN)](https://undiagnosed.hms.harvard.edu/) and similar rare disease programs. RNA-seq adds a critical layer of evidence beyond DNA sequencing:

- **Aberrant expression**: Genes with significantly reduced expression may harbor regulatory or deep intronic variants not visible on exome.
- **Aberrant splicing**: Novel or increased usage of cryptic splice sites can point to pathogenic variants affecting splicing machinery.
- **Intron retention**: Retained introns can indicate splicing defects caused by variants in splice regions.
- **Allele-specific expression**: Monoallelic expression of a heterozygous variant suggests the other allele may be silenced or degraded.

The combination of STAR alignment, IRFinder, BAM2SJ, and the downstream OUTRIDER/FRASER analyses provides a comprehensive view of these RNA-level effects.

## Lessons Learned

**Pre-build multiple STAR indices.** Different sequencing platforms produce different read lengths. Having indices ready for common lengths avoids pipeline failures and re-processing.

**Always check contamination first.** A 30% rRNA contamination rate means you effectively sequenced 30% less of your transcriptome. Catching this early saves time on interpretation.

**Coverage-guided variant calling matters.** Running DeepVariant on the entire genome wastes compute and produces low-confidence calls in poorly covered regions. The Mosdepth → Bedtools → DeepVariant chain focuses resources where the data supports confident calling.

**CRAM saves real money.** At hundreds of samples per year, the storage savings from CRAM over BAM are significant on S3.

**AWS HealthOmics simplifies compliance but adds constraints.** You trade flexibility for managed infrastructure and compliance. Worth it for clinical workflows, but be prepared to work within the platform's limitations.

---

The pipeline is open source and available at [github.com/uclanelsonlab/nl-rna-seq_wf](https://github.com/uclanelsonlab/nl-rna-seq_wf). Contributions and feedback are welcome.
