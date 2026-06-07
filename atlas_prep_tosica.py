# imports
# ----------------------------------------------------------------------------------------------------------------------
import utilTOSICA.preprocessing as uTp
import scanpy as sc
import numpy as np
from pathlib import Path
from scipy.sparse import csr_matrix
import urllib.request
import random
# ----------------------------------------------------------------------------------------------------------------------

def rand_split_adata_obs(obs_idx_list: list, chunk_size: int) -> list:
    shuffled_index = np.random.permutation(obs_idx_list)
    return [shuffled_index[i:i + chunk_size] for i in range(0, len(shuffled_index), chunk_size)]


if __name__ == "__main__":
    out_dir = Path.cwd() / "AtlasTOSICA_out"
    out_dir.mkdir(exist_ok=True, parents=True)

    # retrieve h5ad GBM atlas (cell x gene)
    # https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c; retrieved 02.06.2026
    url = "https://datasets.cellxgene.cziscience.com/b07b7eb7-1613-4d5e-b7d1-1b19aada8534.h5ad"
    tag = "GBM_atlas"
    retrive_as_path = out_dir / f"{tag}.h5ad"
    urllib.request.urlretrieve(url, retrive_as_path)
    print(f"Dataset < {url} > has been retrieved as {retrive_as_path}")

    # import anndata
    adata_gbm = sc.read(retrive_as_path)
    # unique tags only!
    adata_gbm.obs_names_make_unique()
    adata_gbm.var_names_make_unique()

    # dim reduction of adata
    sc.pp.highly_variable_genes(adata_gbm,
                                n_top_genes=3000,
                                flavor="seurat_v3_paper",
                                subset=True)

    obs_new_celltags = [f"{t3} | {t1}" for t3, t1 in zip(list(adata_gbm.obs["annotation_level_3"]), list(adata_gbm.obs["annotation_level_1"]))]
    adata_gbm.obs["curated_cell_tag"] = obs_new_celltags

    dict_keep = uTp.reduce_cellclasses_by_cluster(adata_gbm, "curated_cell_tag", n_count_max=5000)
    list_keep = []
    for list_obs_class in dict_keep.values():
        list_keep.extend(list_obs_class)
    adata_gbm_red = adata_gbm.copy()[list_keep, :]
    # make X dense
    if type(adata_gbm_red.X) != np.ndarray:
        adata_gbm_red.X = csr_matrix.todense(adata_gbm_red.X.copy())
    adata_gbm_red.X = np.nan_to_num(adata_gbm_red.X.copy(), nan=0.0)
    adata_gbm_red.write(out_dir / "gbm_atlas_genes3k_class5k.h5")

    # split big adata_gbm into smaller chunks
    list_shuffled_idx = rand_split_adata_obs(list(adata_gbm.obs.index), 100000)
    for i, list_subindex in enumerate(list_shuffled_idx):
        adata_gbm_red = adata_gbm.copy()[list_subindex, :]
        # make X dense
        if type(adata_gbm_red.X) != np.ndarray:
            adata_gbm_red.X = csr_matrix.todense(adata_gbm_red.X.copy())
        adata_gbm_red.X = np.nan_to_num(adata_gbm_red.X.copy(), nan=0.0)
        adata_gbm_red.write(out_dir / f"slice_{i}__gbm_atlas_genes3k.h5")