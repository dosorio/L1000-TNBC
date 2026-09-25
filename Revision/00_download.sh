#!/usr/bin/env bash
# Downloads the public data used by the revision-3 analyses into ./public and metaLINCS 0.9.0 into ./metaLINCS.
set -euo pipefail
cd "$(dirname "$0")"
[ -d metaLINCS ] || git clone --depth 1 https://github.com/bigomics/metaLINCS.git
mkdir -p public && cd public
get() { [ -s "$1" ] || curl -fL --retry 5 -C - -o "$1" "$2"; }
G=https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5
get GDSC1_fitted_dose_response_27Oct23.xlsx $G/GDSC1_fitted_dose_response_27Oct23.xlsx
get sample_info.csv https://ndownloader.figshare.com/files/35020903          # DepMap 22Q2
H=https://lincs.hms.harvard.edu/_static/db/prod-20200624
get Screen20344_DrugSensitivity2.xlsx $H/Screen20344_DrugSensitivity2.xlsx   # HMS LINCS 20344
get GO_BP_2023.gmt "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023"
get Hallmark_2020.gmt "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
E=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl
get GSE92742_Broad_LINCS_pert_info.txt.gz $E/GSE92742_Broad_LINCS_pert_info.txt.gz
get GSE92742_Broad_LINCS_sig_info.txt.gz $E/GSE92742_Broad_LINCS_sig_info.txt.gz
