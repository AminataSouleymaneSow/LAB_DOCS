#!/bin/bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate GRN_env
cd /mnt/nfs/WORKSP/aminata.sow/LAB_DOCS/GRN_SCENIC/MS

echo "===== STEP 1: GRN ====="
pyscenic grn input/expression.loom input/allTFs_hg38.txt -o output/adjacencies.tsv --num_workers 20
echo "===== GRN FINISHED ====="

echo "===== STEP 2: CTX ====="
pyscenic ctx output/adjacencies.tsv homo_sapiens/*.feather --annotations_fname homo_sapiens/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname input/expression.loom --mask_dropouts --min_genes 5 --output output/regulons.csv --num_workers 20
echo "===== CTX FINISHED ====="

echo "===== STEP 3: AUCell ====="
pyscenic aucell input/expression.loom output/regulons.csv --output output/scenic_output.loom --num_workers 20
echo "===== AUCell FINISHED ====="

echo "===== FULL PIPELINE FINISHED ====="
