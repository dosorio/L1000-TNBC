# Revision 3 analyses

Code for the analyses in the revision-3 response to reviewers. Run everything with:

```bash
bash run_all.sh
```

`run_all.sh` downloads the public data (`00_download.sh`), runs the scripts below in order, and writes all outputs to `results/`, with figures in `figures/`. A full run takes about two minutes after the downloads.

Requirements: R (>= 4.3) with `ccdata`, `preprocessCore`, `fgsea`, `readxl`, `ggplot2`, `ggrepel` and `patchwork`, plus metaLINCS 0.9.0 (`R CMD INSTALL metaLINCS` after `00_download.sh` clones it).

| Script | Reviewer point | Output |
|---|---|---|
| `01_lincs_metadata.R` | 4 | LINCS-L1000 Phase I (GSE92742) and `ccdata` catalogue counts |
| `02_signatures.R` | 1, 3 | retriever (published and exact-matching), naive average and PRL signatures; audit of the published profiles |
| `03_benchmark.R` | 1, 2 | Accuracy against GDSC1 sensitivity in 21 held-out ER-/HER2- breast cancer lines; consistency-filter analysis |
| `04_stability.R` | 1, 3 | Leave-one-cell-line-out robustness |
| `05_public_sensitivity.R` | 2 | GDSC1 and HMS LINCS data for QL-XII-47, GSK-690693, ipatasertib and Torin2 |
| `06_viability.R` | 2 | Bliss reanalysis of Supplementary Table 6 (Figure 5F) |
| `07_threshold.R` | 3 | Sensitivity to the retriever correlation threshold |
| `08_deg_direction.R` | 5 | Direction of the 205 differentially expressed genes |
| `09_enrichment.R` | 5 | GO Biological Process and Hallmark GSEA, leading-edge genes |
| `10_figures.R` | 1, 5 | `Fig_R3_benchmark.png` (Supplementary Fig. S4), `Fig_R3_GO_BP_heatmap.png` (S6) |
| `11_public_sensitivity_figure.R` | 2 | `Fig_R3_public_sensitivity.png` (Supplementary Fig. S5) |

`functions.R` holds the shared code: data loading, retriever with exact metadata matching (identical to the corrected package), the PRL implementation (Iorio et al., PNAS 2010), and the metaLINCS wrapper.

Data sources: LINCS-L1000 signatures from Bioconductor `ccdata`; the single-cell TNBC signature, viability data and HMS LINCS dataset 20367 from this repository, and the previously published retriever profiles from its history (commit ad5332e); GDSC1 release 8.5; DepMap 22Q2 cell line annotation; HMS LINCS dataset 20344; Enrichr GO Biological Process 2023 and MSigDB Hallmark 2020 gene sets; GEO GSE92742 metadata.
