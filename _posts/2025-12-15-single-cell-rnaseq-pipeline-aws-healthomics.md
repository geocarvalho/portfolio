---
layout: post
title: "Building a Single-Cell RNA-seq Pipeline on AWS HealthOmics"
date: 2025-12-15
categories: bioinformatics
tags: [bioinformatics, single-cell, rnaseq, nextflow, aws, pipeline, 10x-genomics]
---

After building the [bulk RNA-seq](/bioinformatics/2025/08/28/building-rnaseq-pipeline-rare-diseases-aws-healthomics.html) and [WGS](/bioinformatics/2025/10/29/building-wgs-pipeline-rare-diseases-aws-healthomics.html) pipelines at the [UCLA Nelson Lab](https://uclanelsonlab.github.io/), the next step was single-cell. The [nl-snrna-seq_wf](https://github.com/uclanelsonlab/nl-snrna-seq_wf) pipeline processes 10x Genomics single-cell and single-nucleus RNA-seq data on [AWS HealthOmics](https://aws.amazon.com/healthomics/), following the same Nextflow + Docker + ECR pattern as the other pipelines.

Single-cell RNA-seq is a different beast from bulk RNA-seq. Instead of measuring the average expression across millions of cells, it captures the transcriptome of individual cells. This is powerful for rare disease research because it reveals cell-type-specific expression patterns, identifies rare cell populations, and can expose disease mechanisms that are diluted away in bulk data.

## Pipeline Overview

The pipeline is intentionally lean compared to the bulk RNA-seq and WGS pipelines. Cell Ranger, the 10x Genomics analysis suite, handles the heavy lifting — alignment, cell calling, UMI deduplication, and gene-barcode matrix generation. The pipeline wraps it with QC and storage optimization:

```
FASTQ R1/R2
    │
    ├── FastP (QC assessment)
    │
    ├── Cell Ranger Count
    │     ├── Alignment to transcriptome
    │     ├── Cell calling & UMI counting
    │     ├── Gene-barcode matrices
    │     ├── Loupe Browser file
    │     └── BAM (optional)
    │           └── Samtools BAM → CRAM
    │
    └── MultiQC (aggregated report)
```

## Why Keep It Simple?

The bulk RNA-seq pipeline has 15+ tools. This one has four. That's by design.

Cell Ranger is an opinionated, end-to-end tool from 10x Genomics. It handles alignment (using STAR internally), barcode error correction, UMI counting, cell/empty droplet classification, and secondary analysis (clustering, t-SNE, UMAP) all in one step. Trying to replicate or replace parts of this with individual tools (like STARsolo + custom barcode handling) adds complexity without clear benefits for our use case.

The downstream single-cell analysis — normalization, clustering, differential expression, cell type annotation — happens interactively in R (Seurat) or Python (Scanpy) by the researchers, not in the pipeline. The pipeline's job is to go from raw FASTQ to clean count matrices reliably and reproducibly.

## The Tools

### FastP — Quality Assessment

[FastP](https://github.com/OpenGene/fastp) runs in **QC-only mode** — no trimming is applied. This is important: Cell Ranger expects untrimmed reads with the full 10x barcode and UMI structure intact. Trimming could damage the barcode sequence in Read 1 or the UMI, breaking cell assignment. FastP is here purely to assess read quality and generate reports that feed into MultiQC.

### Cell Ranger Count — The Core

[Cell Ranger](https://support.10xgenomics.com/) `count` is the primary analysis step. Given paired FASTQ files and a transcriptome reference, it:

1. **Extracts barcodes and UMIs** from Read 1 (the structure depends on the chemistry version).
2. **Aligns** Read 2 to the transcriptome using an embedded STAR aligner.
3. **Corrects barcode errors** using a whitelist of valid barcodes for the specific chemistry.
4. **Deduplicates UMIs** to remove PCR duplicates at the molecular level.
5. **Calls cells** — distinguishes real cells from empty droplets using an algorithm that looks at the total UMI count distribution.
6. **Generates count matrices** — both raw (all barcodes) and filtered (cells only), in HDF5 format.
7. **Runs secondary analysis** — PCA, clustering, t-SNE/UMAP, and differential expression.
8. **Creates a Loupe file** for interactive visualization in the 10x Loupe Browser.

The pipeline supports several key parameters:

- **`chemistry`**: Auto-detected by default, but can be forced to a specific 10x chemistry version (SC3Pv3, SC5P-R2, etc.) when auto-detection fails or for non-standard libraries.
- **`expect_cells`**: The expected number of recovered cells (default 5000). This guides the cell calling algorithm — setting it too low might filter out real cells, too high might include empty droplets.
- **`create_bam`**: When `true`, Cell Ranger produces a BAM file with cell barcode and UMI tags, which the pipeline then converts to CRAM. When `false` (default), alignment files are skipped entirely, saving significant compute time and storage.

Cell Ranger is the most resource-intensive step, requiring 32 CPUs and 128 GB of RAM.

### Samtools — BAM to CRAM Conversion

When `create_bam` is enabled, [Samtools](http://www.htslib.org/) converts the Cell Ranger BAM to CRAM format using the reference genome extracted from the transcriptome directory. This provides 40-60% file size reduction. The original BAM is **not** published to S3 — only the CRAM makes it to the output, saving storage costs.

This is optional because many downstream single-cell analyses only need the count matrices, not the alignments. But when you need alignments — for example, to run [Velocyto](http://velocyto.org/) for RNA velocity analysis, or to inspect specific loci in IGV — having the CRAM available avoids re-running Cell Ranger.

### MultiQC — Aggregated Reporting

[MultiQC](https://multiqc.info/) collects FastP and Cell Ranger metrics into a single HTML report. The Cell Ranger web summary is already comprehensive, but MultiQC provides a standardized view that's consistent with our bulk RNA-seq and WGS pipeline reports, making it easier to assess quality across different assay types.

## AWS HealthOmics Deployment

Same pattern as the other pipelines — zip the repository, import as a private workflow, and launch runs with `run_parameters.json` pointing to S3 paths:

```bash
cd /path/to/nl-snrna-seq_wf/
rm nl-snrna-seq_wf.zip; zip -r nl-snrna-seq_wf.zip *
```

Docker images live in ECR, and outputs go to `/mnt/workflow/pubdir` → S3. The transcriptome reference (e.g., `refdata-gex-GRCh38-2024-A`) also lives in S3 and is mounted into the Cell Ranger process.

## Output Structure

```
output/
├── QC/
│   ├── <sample>_fastp.html
│   ├── <sample>_fastp.json
│   └── multiqc_report.html
└── CELLRANGER/
    ├── web_summary.html
    ├── metrics_summary.csv
    ├── filtered_feature_bc_matrix.h5
    ├── raw_feature_bc_matrix.h5
    ├── molecule_info.h5
    ├── cloupe.cloupe
    ├── <sample>.cram (if create_bam: true)
    ├── <sample>.cram.crai (if create_bam: true)
    └── analysis/
```

The key outputs for downstream analysis:

- **`filtered_feature_bc_matrix.h5`**: The filtered count matrix — this is what goes into Seurat or Scanpy.
- **`raw_feature_bc_matrix.h5`**: The raw matrix with all barcodes, useful for custom cell calling with tools like [CellBender](https://github.com/broadinstitute/CellBender) or [EmptyDrops](https://bioconductor.org/packages/DropletUtils/).
- **`molecule_info.h5`**: Per-molecule information needed for sample aggregation (`cellranger aggr`) across multiple samples.
- **`cloupe.cloupe`**: For quick interactive exploration in the 10x Loupe Browser before diving into code-based analysis.

## How It Fits With the Other Pipelines

The three pipelines at the Nelson Lab form a complementary toolkit for rare disease genomics:

| Pipeline | Assay | What It Reveals |
|----------|-------|-----------------|
| [nl-wgs_wf](https://github.com/uclanelsonlab/nl-wgs_wf) | Whole Genome Sequencing | All variant types genome-wide: SNVs, SVs, CNVs, repeat expansions |
| [nl-rna-seq_wf](https://github.com/uclanelsonlab/nl-rna-seq_wf) | Bulk RNA-seq | Aberrant expression, aberrant splicing, intron retention across tissue |
| [nl-snrna-seq_wf](https://github.com/uclanelsonlab/nl-snrna-seq_wf) | Single-cell/nucleus RNA-seq | Cell-type-specific expression, rare cell populations, cellular heterogeneity |

A pathogenic variant found in WGS might show reduced expression in bulk RNA-seq, and single-cell data can reveal whether that effect is cell-type-specific — for example, only affecting motor neurons but not glial cells. This multi-omics approach is increasingly important for solving the hardest undiagnosed cases.

## Lessons Learned

**Don't fight Cell Ranger.** For 10x data, Cell Ranger is the path of least resistance. The alternative (STARsolo + custom scripts) can work but requires careful handling of barcode whitelists, chemistry detection, and cell calling — all of which Cell Ranger handles robustly.

**Skip the BAM by default.** Most single-cell workflows only need count matrices. Making BAM generation optional (`create_bam: false`) cuts runtime and storage significantly. Only enable it when you actually need alignments.

**Don't trim 10x reads.** This is a common mistake. The barcode + UMI structure in Read 1 must be intact for Cell Ranger. FastP runs in assessment-only mode for this reason.

**Keep the pipeline focused.** The temptation is to add Seurat/Scanpy analysis into the pipeline. Resist it — single-cell analysis is inherently interactive and exploratory. The pipeline should produce clean, reproducible count matrices and let researchers drive the analysis.

---

The pipeline is open source at [github.com/uclanelsonlab/nl-snrna-seq_wf](https://github.com/uclanelsonlab/nl-snrna-seq_wf). Contributions and feedback are welcome.
