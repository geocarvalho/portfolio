---
layout: post
title: "From CRAM to FASTQ: A Practical Guide for Bioinformaticians"
date: 2025-07-02
categories: bioinformatics
tags: [bioinformatics, cram, fastq, nextflow, workflow]
---

In the world of genomics, data often comes in CRAM format to save disk space. But what if you need the original FASTQ files for downstream processing like re-alignment or variant calling with a different tool?

In this post, I'll walk you through the steps to convert a CRAM file into high-quality, paired-end FASTQ files using Dockerized bioinformatics tools like samtools, fastp, and bbmap. You'll also learn how to validate references, check file integrity, and repair read pairs.

## Prerequisites

- At least moderate familiarity with the command line.
- Docker is installed and running on your machine.
- Access to the corresponding FASTA reference file (critical for decompression).
- A CRAM file aligned to hg38 (no-alt or UCSC version depending on the lab standard).

## Step 1: Check Reference Compatibility

Before working with a CRAM, make sure your reference genome matches the one used for compression. This is critical, or decompression will fail or give incorrect results.

```bash
docker run -v $PWD:$PWD quay.io/biocontainers/samtools:1.22--h96c455f_0 samtools samples \
  -h -f $PWD/hg38.fa $PWD/sample.sorted.noalt.hg38.cram
```

Pro Tip: If you see a `.` in the last column of the output, your reference does not match.

## Step 2: Validate Your CRAM

Ensure the file isn't truncated or corrupted:

```bash
docker run -v $PWD:$PWD quay.io/biocontainers/samtools:1.22--h96c455f_0 samtools quickcheck \
  $PWD/sample.sorted.noalt.hg38.cram \
  && echo OK
```

## Step 3: Convert to FASTQ

Now convert the name-sorted BAM to paired-end FASTQ files:

```bash
docker run -v $PWD:$PWD quay.io/biocontainers/samtools:1.22--h96c455f_0 samtools fastq \
  -@ 36 \
  --reference $PWD/ucsc/hg38.fa \
  -1 $PWD/sample.R1.fastq.gz \
  -2 $PWD/sample.R2.fastq.gz \
  $PWD/sample.sorted.noalt.hg38.cram
```

## Step 4: Repair Paired-end Reads

Trimming or filtering can sometimes desynchronize paired-end reads. Fix that with bbmap's `repair.sh`:

```bash
docker run -v $PWD:$PWD quay.io/biocontainers/bbmap:39.26--he5f24ec_0 repair.sh \
    in1=$PWD/sample.R1.fastq.gz \
    in2=$PWD/sample.R2.fastq.gz \
    out1=$PWD/sample_R1.fixed.fastq.gz \
    out2=$PWD/sample_R2.fixed.fastq.gz \
    outs=$PWD/sample_singletons.fastq.gz repair
```

## Step 5: Check FASTQ with Fastp

Use fastp to trim adapters, filter low-quality reads, and generate QC reports.

```bash
docker run -v $PWD:$PWD quay.io/biocontainers/fastp:1.0.0--heae3180_0 fastp \
  -w 36 \
  -i $PWD/sample_R1.fixed.fastq.gz \
  -I $PWD/sample_R2.fixed.fastq.gz \
  -o $PWD/results/sample_R1.fastp.fastq.gz \
  -O $PWD/results/sample_R2.fastp.fastq.gz \
  -j $PWD/results/sample.fastp.json \
  -h $PWD/results/sample.fastp.html \
  --detect_adapter_for_pe
```

## Final Output

After these steps, you will have:

- Quality control reports (`.json`, `.html`) from fastp.
- Singleton reads (optional): `sample_singletons.fastq.gz`
- Cleaned and validated FASTQ files: `sample_R1.fixed.fastq.gz` and `sample_R2.fixed.fastq.gz`

## Troubleshooting Tips

- Repairing reads is important if your pipeline expects synchronized pairs (e.g., for STAR, BWA, etc).
- If Fastp removes too many reads, inspect the quality thresholds or adapter settings.
- If `samtools fastq` fails, double-check that the BAM is name-sorted.
- For example, in one case, a company used `hg38_no_alt.fasta` that was not available, and I could not find it anywhere on the internet. It worked just fine with the UCSC `hg38.fasta`.
- Try to use the exact reference genome used during alignment, but in case you have a derived reference, one possible solution is to use the base reference used to create the derived reference.

## Conclusion

CRAM files are efficient for storage, but getting back to raw reads for reanalysis requires care. With this step-by-step pipeline, you can reliably convert and clean CRAMs into FASTQ using reproducible Docker-based commands.

## Extra: Nextflow Simple Pipeline

Here's a Nextflow pipeline that automates the entire workflow:

1. `check_cram_reference`
2. `validate_cram`
3. `sort_cram_to_bam`
4. `bam_to_fastq`
5. `fastp_cleanup`
6. `repair_reads`

### Directory Structure (suggested)

```
project/
├── main.nf
├── nextflow.config
├── data/
│   ├── input.cram
│   └── reference.fa
```

### main.nf

```groovy
nextflow.enable.dsl=2

params.cram        = "data/input.cram"
params.reference   = "data/reference.fa"
params.sample_id   = "sample"
params.outdir      = "results"

workflow {
    Channel
        .fromPath(params.cram)
        .set { cram_ch }

    check_cram_reference(cram_ch, file(params.reference))
        | validate_cram
        | sort_cram_to_bam
        | bam_to_fastq
        | fastp_cleanup
        | repair_reads
}

process check_cram_reference {
    tag "$sample_id"

    input:
    path cram
    path reference

    output:
    path "cram.ok", emit: ok

    script:
    sample_id = cram.getBaseName().tokenize('.')[0]
    """
    mkdir -p refcheck
    docker run -v \$PWD:\$PWD quay.io/biocontainers/samtools:1.22--h96c455f_0 \\
        samtools samples -h -f $reference $cram > refcheck/${sample_id}_ref_check.txt
    touch cram.ok
    """
}

process validate_cram {
    tag "$sample_id"

    input:
    path cram_ok

    output:
    path "${params.cram}", emit: cram

    script:
    """
    docker run -v \$PWD:\$PWD quay.io/biocontainers/samtools:1.22--h96c455f_0 \\
        samtools quickcheck ${params.cram} && echo OK
    cp ${params.cram} .
    """
}

process sort_cram_to_bam {
    tag "$sample_id"

    input:
    path cram
    path reference from file(params.reference)

    output:
    path "sorted.bam", emit: bam

    script:
    """
    docker run -v \$PWD:\$PWD quay.io/biocontainers/samtools:1.22--h96c455f_0 \\
        samtools sort -@ 4 -O bam -n --reference $reference \\
        -o sorted.bam $cram
    """
}

process bam_to_fastq {
    tag "$sample_id"

    input:
    path bam
    path reference from file(params.reference)

    output:
    tuple path("R1.fastq.gz"), path("R2.fastq.gz")

    script:
    """
    docker run -v \$PWD:\$PWD quay.io/biocontainers/samtools:1.22--h96c455f_0 \\
        samtools fastq -@ 4 --reference $reference \\
        -1 R1.fastq.gz -2 R2.fastq.gz $bam
    """
}

process fastp_cleanup {
    tag "$sample_id"

    input:
    tuple path(r1), path(r2)

    output:
    tuple path("clean_R1.fastq.gz"), path("clean_R2.fastq.gz")

    script:
    """
    mkdir -p ${params.outdir}/fastp
    docker run -v \$PWD:\$PWD quay.io/biocontainers/fastp:1.0.0--heae3180_0 \\
        fastp -w 4 -i $r1 -I $r2 \\
        -o clean_R1.fastq.gz -O clean_R2.fastq.gz \\
        -j ${params.outdir}/fastp/fastp.json \\
        -h ${params.outdir}/fastp/fastp.html \\
        --detect_adapter_for_pe
    """
}

process repair_reads {
    tag "$sample_id"

    input:
    tuple path(r1), path(r2)

    output:
    tuple path("fixed_R1.fastq.gz"), path("fixed_R2.fastq.gz"), path("singletons.fastq.gz")

    script:
    """
    docker run -v \$PWD:\$PWD quay.io/biocontainers/bbmap:39.26--he5f24ec_0 \\
        repair.sh \\
        in1=$r1 in2=$r2 \\
        out1=fixed_R1.fastq.gz \\
        out2=fixed_R2.fastq.gz \\
        outs=singletons.fastq.gz \\
        repair
    """
}
```

### nextflow.config

```groovy
params {
  cram      = "data/input.cram"
  reference = "data/reference.fa"
  sample_id = "sample"
  outdir    = "results"
}

process.container = ''
process.executor = 'local'
process.errorStrategy = 'terminate'
```

### Run the Pipeline

```bash
nextflow run main.nf --cram data/input.cram --reference data/reference.fa --sample_id SAMPLE01
```

---

*Published originally on [Medium](https://medium.com/@geocarvalho/from-cram-to-fastq-a-practical-guide-for-bioinformaticians-d24c82fde5eb).*
