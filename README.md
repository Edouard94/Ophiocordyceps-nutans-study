# Disentangling the *Ophiocordyceps nutans* Complex Through Museomics

This repository contains sequencing datasets, amplicon-processing pipelines, primer-design resources, and phylogenetic workflows used to study the *Ophiocordyceps nutans* species complex from museum material.

## Repository structure

| Folder | Purpose |
| --- | --- |
| `amplicon processing/` | ONT and Illumina multiplex amplicon pipelines, run manifests, and generated consensus/taxonomy outputs. |
| `docs/` | Study metadata tables (sample sheets, extraction QC, sequence tracking). |
| `phylogeny/` | Scripts for reference-sequence retrieval and multigene phylogeny preparation/inference. |
| `primer design/` | Primer design workflow script plus alignments, entropy plots, and Primer-BLAST/PrimerProspector outputs. |
| `raw data/` | Placeholder folder for additional raw inputs not currently tracked in Git. |
| `results/` | Top-level destination folders for project figures and summary tables. |

## Folder-level documentation

- [`amplicon processing/README.md`](amplicon%20processing/README.md)
- [`docs/README.md`](docs/README.md)
- [`phylogeny/README.md`](phylogeny/README.md)
- [`primer design/README.md`](primer%20design/README.md)
- [`raw data/README.md`](raw%20data/README.md)
- [`results/README.md`](results/README.md)

## Citation

Please use the citation metadata in `CITATION.cff`.

## License

- **Software and code** (for example, `.py` and `.sh` files) are licensed under the MIT License. See [LICENSE](LICENSE).
- **Data, documentation, and other non-software repository content** are licensed under CC BY 4.0. See [LICENSE-DATA](LICENSE-DATA).
