# Amplicon processing

ONT FASTQ inputs and intermediate reads are tracked with Git LFS. After cloning,
run `git lfs install` and `git lfs pull` to retrieve their contents.

This folder stores platform-specific amplicon processing workflows and their run outputs for the *Ophiocordyceps nutans* project.

## Subfolders

| Path | Purpose |
| --- | --- |
| `Illumina/` | Multiplex paired-end amplicon workflow inputs and outputs. |
| `ONT/` | Nanopore amplicon workflow script, manifests, and run outputs. |

## Illumina workflow (`Illumina/`)

### Important files

- `multiplex_marker_pipeline.sh`  
  Multiplex marker pipeline (`fastp` → `cutadapt` → `vsearch` → consensus extraction → `blastn`/`taxonkit`) for ITS/SSU/LSU/TEF amplicons.
- `Julia_samples.csv`  
  Sample metadata and primer definitions used as pipeline input.
- `Julia_data/*.fastq.gz`  
  Input paired-end reads for samples `203_A`, `204_B`, `205_A`, `206_A`.

### Key processing steps (from script)

1. Read sample names from the CSV input.
2. Run quality control (`fastp`; optional `multiqc`).
3. Trim marker primers with `cutadapt`.
4. Merge read pairs and quality-filter merged reads (`vsearch`, `seqkit`).
5. Dereplicate and select top-abundance per-sample consensuses.
6. Assign taxonomy with `blastn` and format lineages with `taxonkit`.
7. Export summary statistics/tables.

### Outputs present in this repository

- `Julia_samples/qc/`: per-sample `fastp` HTML/JSON, QC FASTQs, and a MultiQC report.
- `Julia_samples/taxonomy/`: BLAST results/top hits, centroids, and lineage/taxid tables.
- `Julia_samples/final_results/`: `barcode_summary.tsv`, `consensus_sequences.tsv`, `pipeline_stats.tsv`.

## ONT workflow (`ONT/`)

### Important files

- `ONT_Amplicon_Pipeline_museomics.sh`  
  ONT demultiplexing/primer-trimming pipeline with `cutadapt`, `seqkit`, NGSpeciesID consensus generation, and BLAST taxonomy assignment.
- `All_input.csv`  
  Sample/marker/tag/primer manifest used by the ONT script.
- `all_barcodes.fastq`  
  Input ONT reads referenced by the workflow.

### Key processing steps (from script)

1. Read sample tags/primers from `All_input.csv`.
2. Demultiplex by linked tags (`cutadapt`) in forward and reverse orientations.
3. Trim marker primers and length/quality filter reads.
4. Run NGSpeciesID (+ Medaka) per sample-marker, then filter consensus sequences by supporting-read threshold.
5. Build global consensus tables.
6. Run `blastn` + `taxonkit` taxonomy mapping and write final summary tables.

### Outputs present in this repository

- `Onutans_203/` and `ONT_main_070726/`: saved pipeline runs with manifests and staged outputs. `ONT_main_070726/` is the larger multiplex run; `Onutans_203/` contains an earlier sample-focused run.
- `*/stage1_tagdemux/`: demultiplexed FASTQs and `demux_summary.tsv`.
- `*/samples/`: per-sample/per-marker intermediate files and NGSpeciesID outputs.
- `*/final_results/`: consolidated outputs such as `all_consensuses_annotated.tsv`, `sample_marker_status.tsv`, `selected_consensus.tsv`, BLAST/taxonomy summaries, and coverage/read-count reports.

The ONT workflow uses `All_input.csv` for sample, marker, tag, and primer metadata
and writes run-specific manifests, logs, demultiplexed reads, per-sample/per-marker
NGSpeciesID results, and final consensus tables.
