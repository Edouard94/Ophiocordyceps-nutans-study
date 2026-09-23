# Phylogeny

This folder contains scripts and reference resources used to build multigene phylogenetic datasets for the *Ophiocordyceps nutans* complex.

## Subfolders

| Path | Purpose |
| --- | --- |
| `amplicon & phylogenetic scripts/` | Data extraction, QC, platform comparison, and multigene tree pipeline scripts. |
| `references for phylogeny/` | Representative reference list, retrieval script, downloaded reference FASTAs, and NCBI XML cache/audits. |

## `amplicon & phylogenetic scripts/`

- `csv_to_amplicons.py`  
  Extracts per-marker FASTA files from a combined CSV table (expects marker, phylo_id, sequence fields).
- `qc_alignments.py`  
  Filters alignment records by N proportion, mismatch-to-consensus proportion, and minimum informative sites.
- `compare_platforms.py`  
  Compares sequence identity across platforms per sample/marker and exports CSV matrices and optional heatmaps.
- `multigene_phylo.sh`  
  End-to-end phylogeny pipeline: amplicon extraction, optional consensus generation, MAFFT alignment, optional trim/QC, AMAS concatenation/partitioning, IQ-TREE inference, and platform-comparison outputs.

### Expected inputs

- Combined sequence CSV (e.g., from `docs/All_seq_data.xlsx` exported to CSV).
- Reference FASTA sets from a references directory (`references_last/` or fallback `references/` in script logic).
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
