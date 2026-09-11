import loompy
import scipy.io
import pandas as pd
import numpy as np
print("Creating LOOM file in the input folder...")
matrix_file   = "input/matrix.mtx"
genes_file    = "input/genes.tsv"
barcodes_file = "input/barcodes.tsv"
output_loom   = "input/expression.loom"
matrix = scipy.io.mmread(matrix_file).tocsc()
genes = pd.read_csv(genes_file, header=None)[0].values
cells = pd.read_csv(barcodes_file, header=None)[0].values
assert matrix.shape[0] == len(genes), "Mismatch: genes vs matrix rows"
assert matrix.shape[1] == len(cells), "Mismatch: cells vs matrix columns"
_, idx = np.unique(genes, return_index=True)
genes = genes[idx]
matrix = matrix[idx, :]
loompy.create(output_loom, matrix, row_attrs={"Gene": genes}, col_attrs={"CellID": cells})
print("LOOM file created successfully!")
print("Matrix shape:", matrix.shape)
