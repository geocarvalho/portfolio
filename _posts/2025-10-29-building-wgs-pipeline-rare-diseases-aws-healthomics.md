---
layout: post
title: "Building a WGS Pipeline for Rare Diseases on AWS HealthOmics"
date: 2025-10-29
categories: bioinformatics
tags: [bioinformatics, wgs, nextflow, aws, rare-disease, pipeline]
---

Alongside the [RNA-seq pipeline](/bioinformatics/2025/08/28/building-rnaseq-pipeline-rare-diseases-aws-healthomics.html) I built at the [UCLA Nelson Lab](https://uclanelsonlab.github.io/), I also developed a whole genome sequencing (WGS) pipeline for germline short-read data: [nl-wgs_wf](https://github.com/uclanelsonlab/nl-wgs_wf). Like its RNA-seq counterpart, it's written in Nextflow (DSL2), runs on [AWS HealthOmics](https://aws.amazon.com/healthomics/), and is designed for rare disease diagnostics. This post covers the architecture, tool choices, and design decisions behind it.

## Pipeline Architecture

One key design goal was flexibility in input formats. Clinical labs send data in different states — sometimes raw FASTQ, sometimes already-aligned BAM or CRAM files. Instead of forcing a single entry point, the pipeline uses **subworkflows** that converge into a shared analysis path:

```
Inputs
  ├── FASTQ R1/R2 → Fastp → BWA-MEM2 → Sort
  ├── BAM ──────────────────────────── → Sort
  └── CRAM ─────── → CRAM-to-BAM ──── → Sort
                                          │
                              ┌────────────┘
                              ▼
                     Picard MarkDuplicates
                     Samtools Index
                              │
              ┌───────────────┼───────────────────┐
              ▼               ▼                   ▼
         QC & Metrics    Variant Calling    SV & Repeat Analysis
              │               │                   │
              ▼               ▼                   ▼
        ┌─────────┐    ┌───────────┐    ┌─────────────────┐
        │ Picard   │    │DeepVariant│    │ Manta (SVs)     │
        │ Qualimap │    │    ↓      │    │ CNVpytor (CNVs) │
        │ Mosdepth │    │ AutoMap   │    │ ExpansionHunter │
        │ MultiQC  │    │ (ROH)    │    │ EH Denovo       │
        └─────────┘    │    ↓      │    └─────────────────┘
                        │ BCFtools  │
                        │    ↓      │
                        │ HapCUT2   │
                        │ (phasing) │
                        └───────────┘
                              │
                              ▼
                        BAM → CRAM
```

This subworkflow architecture means you can re-analyze existing alignments without re-running the expensive alignment step, and CRAM files received from external labs can enter the pipeline directly.

## Running on AWS HealthOmics

The pipeline is structured to run as a private workflow on AWS HealthOmics. Deployment is straightforward:

```bash
cd /path/to/nl-wgs_wf/
rm nl-wgs_wf.zip; zip -r nl-wgs_wf.zip *
```

You import the zip file into HealthOmics along with `parameters.json` for the parameter definitions, and use `run_parameters.json` to specify S3 paths when launching runs. All Docker images live in AWS ECR, and outputs go to `/mnt/workflow/pubdir` which syncs back to S3.

The same benefits from the RNA-seq pipeline apply here: managed compute, HIPAA compliance, no infrastructure to maintain, and native Nextflow support.

## The Tools and Why They're There

### Fastp — Quality Control

[Fastp](https://github.com/OpenGene/fastp) handles adapter trimming and quality filtering on the raw FASTQ input. Same reasons as in the RNA-seq pipeline: it's fast, single-pass, and produces JSON output that feeds directly into MultiQC.

### BWA-MEM2 — Alignment

[BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) is the SIMD-accelerated successor to BWA-MEM. For WGS data with billions of reads, the speed improvement over the original BWA-MEM is significant — roughly 2-3x faster on modern hardware. It produces identical alignments to BWA-MEM but takes better advantage of CPU vector instructions. Read group information is added during alignment so samples are properly identified throughout the pipeline.

### Picard — Duplicate Marking and Metrics

We use [Picard](https://broadinstitute.github.io/picard/) for multiple purposes:

- **MarkDuplicates**: Flags PCR duplicates. For WGS (unlike our RNA-seq pipeline where we used Sambamba), we chose Picard here because it integrates better with the GATK/DeepVariant ecosystem and we already need Picard for metrics collection.
- **CollectMultipleMetrics**: Alignment summary, insert size distribution, quality score distribution, and more.
- **CollectWgsMetrics**: WGS-specific metrics including coverage uniformity, mean coverage, and the percentage of bases reaching various depth thresholds.

These metrics are critical for assessing whether a genome has sufficient quality and coverage for clinical interpretation.

### Qualimap — BAM Quality Control

[Qualimap](http://qualimap.conesalab.org/) provides additional BAM-level QC including coverage across chromosomes, GC content distribution, and mapping quality statistics. It complements Picard metrics with different visualizations and is particularly useful for spotting coverage anomalies.

### DeepVariant — SNV/Indel Calling

[DeepVariant](https://github.com/google/deepvariant) is our primary variant caller for SNVs and small indels. Unlike traditional callers like GATK HaplotypeCaller, DeepVariant uses a deep learning model trained on truth sets to evaluate evidence for variants. It consistently performs well in [precisionFDA Truth Challenges](https://precision.fda.gov/challenges) and produces both VCF and gVCF outputs.

For WGS we use the standard WGS model (unlike the RNA-seq pipeline where we needed a custom model). DeepVariant is the most resource-intensive step at 192 GB RAM and 48 CPUs, but produces high-quality variant calls.

### AutoMap — Runs of Homozygosity

[AutoMap](https://github.com/mquinodo/AutoMap) detects runs of homozygosity (ROH) from the DeepVariant VCF output. ROH regions are important in rare disease because:

- They can indicate **consanguinity**, which increases the likelihood of autosomal recessive conditions.
- Homozygous pathogenic variants in ROH regions are strong candidates for disease causation.
- Large ROH blocks may point to **uniparental disomy** (UPD).

AutoMap outputs both TSV files with ROH coordinates and PDF visualizations across the genome.

### HapCUT2 — Haplotype Phasing

[HapCUT2](https://github.com/vibansal/HapCUT2) phases heterozygous variants into haplotype blocks. The phasing workflow is:

1. **BCFtools** filters the DeepVariant VCF to retain only diploid genotypes (0/0, 0/1, 1/1).
2. **HapCUT2 extractHAIRS** extracts haplotype-informative reads from the BAM.
3. **HapCUT2** assembles these into phased haplotype blocks.

Phasing is valuable for rare disease because it tells you which variants are on the same chromosome (cis) versus opposite chromosomes (trans). For compound heterozygous variants in a recessive gene, you need to confirm the two variants are on different alleles — phasing provides this evidence directly from the sequencing data.

### Manta — Structural Variant Detection

[Manta](https://github.com/Illumina/manta) detects structural variants (SVs) including large deletions, duplications, inversions, and insertions. Manta is fast, well-validated, and specifically designed for germline analysis. It uses split-read and paired-end evidence to call SVs with high sensitivity.

SVs are particularly important in rare diseases because they can disrupt genes in ways not captured by SNV callers — a large deletion removing an entire exon, an inversion disrupting a gene, or a translocation creating a fusion.

### CNVpytor — Copy Number Variant Analysis

[CNVpytor](https://github.com/abyzovlab/CNVpytor) detects copy number variants (CNVs) from read depth signals. It analyzes coverage across the genome at multiple bin sizes to identify regions with significantly more or fewer reads than expected. CNVpytor also produces Manhattan plots for visual inspection.

CNVpytor complements Manta: Manta excels at breakpoint-resolved SVs from paired-end/split-read evidence, while CNVpytor catches larger CNVs that may not have clear breakpoint signatures but show clear read depth changes.

### ExpansionHunter and ExpansionHunterDenovo — Repeat Expansions

Repeat expansion disorders are an important class of rare diseases (Huntington's disease, Fragile X, various ataxias, ALS). We use two complementary tools:

- [ExpansionHunter](https://github.com/Illumina/ExpansionHunter) genotypes known repeat expansion loci from a curated variant catalog. It can detect expansions beyond the read length by analyzing spanning reads, flanking reads, and in-repeat reads. This is the targeted approach — you look at known disease loci.
- [ExpansionHunterDenovo](https://github.com/Illumina/ExpansionHunterDenovo) takes the untargeted approach: it profiles the genome for any locus showing evidence of repeat expansion, even if it's not in a known catalog. This is powerful for discovering novel repeat expansions in undiagnosed patients.

### Mosdepth — Coverage Analysis

[Mosdepth](https://github.com/brentp/mosdepth) runs on the CRAM output to calculate coverage over specific BED regions, particularly **mitochondrial** regions. Mitochondrial coverage in WGS data reflects the mitochondrial copy number, which can be clinically relevant, and adequate coverage is needed for mitochondrial variant calling.

### MultiQC — Aggregated Reporting

[MultiQC](https://multiqc.info/) aggregates reports from Fastp, Picard (multiple metrics + WGS metrics), Qualimap, DeepVariant, and Mosdepth into a single HTML report for quick sample-level QC assessment.

## RNA-seq vs WGS: Complementary Pipelines

These two pipelines were designed to work together for rare disease diagnostics:

| Aspect | RNA-seq Pipeline | WGS Pipeline |
|--------|-----------------|--------------|
| **Primary input** | FASTQ only | FASTQ, BAM, or CRAM |
| **Aligner** | STAR (splice-aware) | BWA-MEM2 (linear) |
| **Variant calling** | DeepVariant (custom RNA model, CDS-restricted) | DeepVariant (standard WGS model, genome-wide) |
| **Unique analyses** | Intron retention, splice junctions, transcript quantification, rRNA/globin contamination | Structural variants, CNVs, repeat expansions, haplotype phasing, runs of homozygosity |
| **Key strength** | Detects expression and splicing anomalies | Detects all variant types genome-wide |

A variant found in WGS gains evidence when RNA-seq shows it's expressed and doesn't cause nonsense-mediated decay. Conversely, an aberrant splicing event found in RNA-seq points to where to look in the WGS data for the causal variant.

## Lessons Learned

**Support multiple input formats from day one.** Clinical collaborators send data in whatever format they have. The subworkflow architecture that accepts FASTQ, BAM, and CRAM saved us from constant format-conversion requests.

**Combine targeted and untargeted repeat analysis.** ExpansionHunter catches known loci reliably, but ExpansionHunterDenovo has found novel expansions that would have been missed otherwise.

**Phasing is underused in clinical genomics.** Adding HapCUT2 was a late addition but immediately proved its value for resolving compound heterozygotes without requiring parental samples.

**Complement SV callers.** No single SV caller catches everything. Manta + CNVpytor gives better coverage across the SV size spectrum than either alone.

**Runs of homozygosity are a quick win.** AutoMap runs fast and immediately flags consanguinity and candidate regions for recessive disease — information that's useful even before looking at individual variants.

---

The pipeline is open source at [github.com/uclanelsonlab/nl-wgs_wf](https://github.com/uclanelsonlab/nl-wgs_wf). Contributions and feedback are welcome.
