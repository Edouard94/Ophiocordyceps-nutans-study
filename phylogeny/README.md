# Phylogeny

This folder contains scripts and reference resources used to build multigene phylogenetic datasets for the *Ophiocordyceps nutans* complex.

## Subfolders

| Path | Purpose |
| --- | --- |
| `misc scripts/` | Amplicon extraction, alignment QC, and platform-comparison scripts. |
| `single genes phylogeny/` | Multigene pipeline script and single-gene/reference dataset outputs. |
| `final phylogeny/` | Final marker alignments, concatenated files, IQ-TREE outputs, and PDF tree. |
| `references for phylogeny/` | Representative reference list, retrieval script, downloaded reference FASTAs, and NCBI XML cache/audits. |

## `misc scripts/`

- `csv_to_amplicons.py`  
  Extracts per-marker FASTA files from a combined CSV table (expects marker, phylo_id, sequence fields).
- `qc_alignments.py`  
  Filters alignment records by N proportion, mismatch-to-consensus proportion, and minimum informative sites.
- `compare_platforms.py`  
  Compares sequence identity across platforms per sample/marker and exports CSV matrices and optional heatmaps.
## `single genes phylogeny/`

- `multigene_phylo.sh`  
  End-to-end phylogeny pipeline: amplicon extraction, optional consensus generation, MAFFT alignment, optional trimming/QC, AMAS concatenation and partitioning, and IQ-TREE inference.

This folder also contains single-gene amplicons, alignments, phylogeny outputs,
and a strict combined-reference dataset used for marker-level analyses.

### Expected inputs

- Combined sequence CSV (e.g., from `docs/All_seq_data.xlsx` exported to CSV).
- Reference FASTA sets from `references for phylogeny/references_folder/` or the single-gene reference datasets.
- External tools noted in scripts (`mafft`, `trimal`, `AMAS.py`, `iqtree`, etc.).

### Typical outputs produced by the scripts

- Per-gene amplicons and alignments (`*_amplicons.fasta`, `*_aligned.fasta`, `*_trimmed.fasta`).
- QC logs (`qc_*.tsv`) and reversal/platform comparison logs.
- Concatenated NEXUS alignment + partition file from AMAS.
- IQ-TREE outputs (`*.treefile`, `*.contree`, logs).

## `references for phylogeny/`

- `references_file.csv`: representative species/voucher/accession table used as retrieval input.
- `fetch_refs_from_reps_csv.py`: retrieves references from NCBI (or validates with `--parse-only`), caches XML records, audits missing/rejected records, and writes per-gene FASTAs.
- `references_folder/`: current retrieved datasets (`ITS_refs.fasta`, `SSU_refs.fasta`, `LSU_refs.fasta`, `TEF_refs.fasta`, `RPB1_refs.fasta`, `RPB2_refs.fasta`), `xml_cache/`, and audit tables in `audits/`.

## `final phylogeny/`

The committed final analysis contains nucleotide alignments for nSSU, nLSU, TEF,
RPB1, and RPB2, concatenated NEXUS/PHYLIP files, IQ-TREE model and tree outputs,
and `Ophio_nutans.pdf`. The `scripts/` subfolder contains helper scripts for
GenBank extraction and alignment concatenation.
