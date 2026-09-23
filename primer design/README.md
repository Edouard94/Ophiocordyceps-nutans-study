# Primer design

This folder contains the primer-design workflow and generated assets used for the *Ophiocordyceps nutans* study.

## Main workflow script

- `Primer_Design_Ophio_cordyceps_genes.sh`
  - Downloads fungal (*Ophiocordyceps*/*Cordyceps*) and Hemiptera gene records from NCBI.
  - Converts/filter records, removes duplicates (`seqkit`, `cd-hit-est`), aligns (`mafft`), trims (`trimal`), and generates entropy plots (`plotcon`).
  - Includes optional blocks for MBC-prime, ecoPrimers, and PrimerProspector analyses.

### Expected inputs

- NCBI access via `esearch`/`efetch`.
- External tools referenced in the script (OBITools, seqkit, cd-hit-est, MAFFT, trimAl, EMBOSS, ecoPrimers, PrimerProspector).

### Outputs represented in this folder

- `alignments/`: raw and trimmed per-gene alignments (`*_c99.aln`, `*_c99_trimmed.aln.fasta`).
- `entropy_plots/` and `entropy_plots1_non-trimmed/`: per-gene entropy PNGs.
- `PrimerProspector_base_frequencies_museum_samples/`: per-gene base-frequency summaries.
- `Primer_Blast_Museomics/`: Primer-BLAST PDF result reports.
- `Primers_museomics_opti.xlsx`: curated/optimized primer table.
