import os
import psutil

process = psutil.Process(os.getpid())
base_memory_usage = process.memory_info().rss

from datetime import datetime
start_time = datetime.now()
##################################################
import pandas as pd
import anndata
import sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
from matplotlib import rcParams
# scRNA-seq reference data
scsim = "/storage/group/qul12/default/retrofit/Mouse/CerebellumPuck_sc_count.tsv"
scsimlabel = "/storage/group/qul12/default/retrofit/Mouse/CerebellumPuck_sc_label.tsv"
adata = anndata.read_csv(scsim, delimiter='\t')
obs = pd.read_csv(scsimlabel, delimiter='\t')
adata.obs = obs
adata_ref = adata
# prepare scRNA-seq anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        #batch_key='bio_celltype',
                        # cell type, covariate used for constructing signatures
                        labels_key='bio_celltype'#,
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        #categorical_covariate_keys=['Method']
                       )
# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()
mod.train(max_epochs=2500, use_gpu=True)
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)
# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
######load ST sim data
sim = pd.read_csv("/storage/group/qul12/default/retrofit/Mouse/CerebellumPuck_counts.csv",
                   index_col=0,
                   delimiter=',')
adata_vis = sc.AnnData(sim.T)
location =  pd.read_csv("/storage/group/qul12/default/retrofit/Mouse/Cerebellum_coords.csv",
                        delimiter=',')
adata_vis.obs = location
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model 
cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=1,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()
mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )
# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)
proportion = adata_vis.obsm['means_cell_abundance_w_sf']
proportion = proportion.div(proportion.sum(axis=1), axis=0)
for i in proportion.columns:
    proportion.columns=proportion.columns.str.replace(i, i.split("_")[4])
proportion.index=adata_vis.obs['barcode']

#################################################
memory_usage = process.memory_info().rss
loop_memory_usage = memory_usage - base_memory_usage
print(loop_memory_usage/1024**2)

end_time = datetime.now()
print((end_time - start_time).total_seconds()/60,"min")

