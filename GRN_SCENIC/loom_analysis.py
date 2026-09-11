import sys
import os
import pandas as pd
import loompy as lp
import numpy as np
import json
import base64
import zlib
from scipy.stats import rankdata

# ── Get condition from command line ────────────────────────────────────────
if len(sys.argv) < 2:
    print("Usage: python loom_analysis.py <condition>")
    print("Example: python loom_analysis.py Low")
    sys.exit(1)

condition = sys.argv[1]
print(f"Running SCENIC analysis for: {condition}")

# ── Paths (auto-built from condition) ─────────────────────────────────────
base_path   = f"/mnt/nfs/WORKSP/aminata.sow/LAB_DOCS/GRN_SCENIC/{condition}"
loom_path   = f"{base_path}/output/scenic_output.loom"
annot_path  = f"{base_path}/input/annotation.txt"
output_path = f"{base_path}/output/"

# ── Check files exist ──────────────────────────────────────────────────────
if not os.path.exists(loom_path):
    print(f"ERROR: loom file not found: {loom_path}")
    sys.exit(1)
if not os.path.exists(annot_path):
    print(f"ERROR: annotation file not found: {annot_path}")
    sys.exit(1)

# ── Load loom ──────────────────────────────────────────────────────────────
print("Loading loom file...")
lf = lp.connect(loom_path, mode='r+', validate=False)

# ── AUC matrix ─────────────────────────────────────────────────────────────
print("Extracting AUC matrix...")
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
print(f"AUC matrix shape: {auc_mtx.shape}")

# ── Load & align annotation ────────────────────────────────────────────────
print("Loading annotation...")
annot  = pd.read_csv(annot_path, sep='\t', header=0, index_col=0)
annot2 = annot.reindex(lf.ca.CellID)
print(f"Annotation aligned: {all(annot2.index == lf.ca.CellID)}")

# ── Update loom metadata ───────────────────────────────────────────────────
print("Updating loom metadata...")
meta = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))
meta["annotations"] = [
    {
        "name": "Clusters",
        "values": list(set(annot2['Clusters'].astype(str)))
    }
]
lf.attrs['MetaData'] = base64.b64encode(
    zlib.compress(json.dumps(meta).encode('ascii'))
).decode('ascii')
lf.ca.Clusters = np.array(annot2['Clusters'].values)

# ── Cell annotation table ──────────────────────────────────────────────────
cellAnnot = pd.DataFrame(lf.ca.Clusters, index=lf.ca.CellID)
cellAnnot.columns = ['Clusters']

# ── Close loom ─────────────────────────────────────────────────────────────
lf.close()
print("Loom file closed!")

# ── Regulon Specificity Scores (numpy fix) ─────────────────────────────────
print("Computing Regulon Specificity Scores...")

def rss_fix(aucs, labels):
    celltypes = sorted(labels.unique())
    regulons  = aucs.columns
    rss = pd.DataFrame(index=celltypes, columns=regulons, dtype=float)
    for ct in celltypes:
        mask = labels == ct
        for reg in regulons:
            vals = aucs[reg].values
            r = rankdata(vals) / len(vals)
            rss.loc[ct, reg] = r[mask].mean()
    return rss

rss_cellType = rss_fix(auc_mtx, cellAnnot['Clusters'])
print("RSS done!")

# ── Export ─────────────────────────────────────────────────────────────────
auc_mtx.to_csv(output_path + "auc_mtx.csv")
rss_cellType.to_csv(output_path + "rss_cellType.csv")
print("Done! Files saved:")
print(f"  -> {output_path}auc_mtx.csv")
print(f"  -> {output_path}rss_cellType.csv")
