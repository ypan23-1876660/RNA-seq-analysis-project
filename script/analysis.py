   
# ## RNA-seq analysis
# Paper: Genetic and genomic analyses of Drosophila melanogaster models of chromatin modification disorders

   
## Setup
import numpy as np
import scipy 
from scipy import stats
import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols
import seaborn as sns
import csv
import os
import sys
import argparse
from time import time
import pandas as pd

import altair as alt
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

pip install pydeseq2

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# # Filtering data 
# Loading files: raw gene counts, description.xlsx
from google.colab import files
uploaded = files.upload()


#Read in data
counts = pd.read_csv("GSE213763_read_counts_forGEO.csv", index_col = 0).iloc[:, 1:].T 
counts.index = counts.index.str.replace('_sorted.bam', '')
description = pd.read_excel("description.xlsx", index_col=0)#.index.name=None#.rename(columns = {'sample_name':'Geneid'})#.reset_index(drop=True)#.set_index('Geneid')
description.index.name=None
description['disease'] = np.where(description.gene.str.contains('osa|brm|Snr1'),'SSRIDD','CdLS')
description.iloc[36:42, -1] = 'Control'

 
# Convert to log2
counts_log2 = np.log2(counts + 0.5)
counts_log2_values = counts_log2.values.flatten()  # flatten the matrix into a 1D array

 
# Histogram for counts
'''
plt.hist(counts_log2_values, bins=50)
plt.title("Histogram of Log2 RNAseq counts")
plt.xlabel("Counts")
plt.ylabel("Frequency")
plt.show()
'''

 
#first calculate the mean expression for each gene (remember we did add a small count to zeroes)
mean = np.mean(counts_log2,axis = 0)


#next calculate the variability within each gene
sqrt_std = np.sqrt(counts_log2.std())

 
plt.title("Mean-Variance Plot")
plt.xlabel("Mean Log2(counts)")
plt.ylabel("sqrt(Std.Dev.)") 

#fit polynomial models up to degree 6
model1 = np.poly1d(np.polyfit(mean, sqrt_std, 1))
model2 = np.poly1d(np.polyfit(mean, sqrt_std, 2))
model3 = np.poly1d(np.polyfit(mean, sqrt_std, 3))
model4 = np.poly1d(np.polyfit(mean, sqrt_std, 4))
model5 = np.poly1d(np.polyfit(mean, sqrt_std, 5))
model6 = np.poly1d(np.polyfit(mean, sqrt_std, 6))

#create scatterplot
polyline = np.linspace(-1, 15,50 )
plt.scatter(mean, sqrt_std, alpha=0.4,s = 1)
plt.xticks(np.arange(0, 20, 1))


#add fitted polynomial lines to scatterplot 
plt.plot(polyline, model1(polyline), color='green')
plt.plot(polyline, model2(polyline), color='black')
plt.plot(polyline, model3(polyline), color='purple')
plt.plot(polyline, model4(polyline), color='blue')
plt.plot(polyline, model5(polyline), color='red')
plt.plot(polyline, model6(polyline), color='orange')
plt.show()


 
# Set drop off to be at 2 or 3, selected 2
keep_mean = mean > 2

# Filter counts_log2 based on selected threshold log2 counts of 2
counts_log2_filtered = counts_log2.loc[:,keep_mean]
counts_filtered = counts.loc[:,keep_mean]

 
# Histogram for filtered counts 
counts_log2_filtered_values = counts_log2_filtered.values.flatten()  # flatten the matrix into a 1D array
'''
plt.hist(counts_log2_filtered_values, bins=50)
plt.title("Histogram of Log2 RNAseq filtered counts")
plt.xlabel("Counts")
plt.ylabel("Frequency")
plt.show()
'''

 
# Checking filtered scatterplot
mean_filtered = np.mean(counts_log2_filtered,axis = 0)
sqrt_std_filtered = np.sqrt(counts_log2_filtered.std())
plt.scatter(mean_filtered, sqrt_std_filtered, alpha=0.4,s = 1)
plt.title("Mean-Variance Plot")
plt.xlabel("Mean Log2(counts)")
plt.ylabel("sqrt(Std.Dev.)") 
plt.show()

# # DESEQ analysis
counts_filtered_w_description = counts_filtered.merge(description, left_index=True, right_index=True)
counts_filtered_w_description.index.name=None

 
# filtered counts file for clds vs control and ssridd vs control 
ssridd_control = counts_filtered_w_description.iloc[np.where(counts_filtered_w_description.disease.str.contains('SSRIDD|Control'))].iloc[:,2:15192]
cdls_control = counts_filtered_w_description.iloc[np.where(counts_filtered_w_description.disease.str.contains('CdLS|Control'))].iloc[:,2:15192]

 

# description file for cdls vs control and ssridd vs control
ssridd_control_descripton = description.iloc[np.where(description.disease.str.contains('SSRIDD|Control'))]
cdls_control_descripton = description.iloc[np.where(description.disease.str.contains('CdLS|Control'))]
ssridd_control_descripton.name=None
cdls_control_descripton.name=None

 
def DESeq2Wrapper(counts,meta,comparison):
#Start by initializing the data, using the DeseqDataSet() function. The function requires 3 inputs: 
# the CountsData data frame (containing the gene expression data in counts format), 
# the meta data frame (containing the metadata associated with each sample),
# and the name of the column in meta that contains the group information for which comparison you want to make. To compare normal vs. metastatic, use "Stage" 
  dds = DeseqDataSet(counts = counts, clinical = meta, design_factors=comparison)

  dds.fit_size_factors()
  scaling_deseq = list(dds.obsm["size_factors"])
  dds.fit_genewise_dispersions()

  dds.deseq2()

  dds.calculate_cooks()
  if dds.refit_cooks:
    # Replace outlier counts
    dds.refit()

  CountsNormal = pd.DataFrame(dds.layers['normed_counts'], 
                                     index=dds.obs.index, 
                                     columns=dds.var.index)

  deseq_stats = DeseqStats(dds,alpha = 0.05)
  deseq_stats.summary()
  return CountsNormal, deseq_stats.results_df

 
# deseq analys cdls vs control and ssridd vs control
count_normal_ssridd_control, deseq_results_ssridd_control = DESeq2Wrapper(ssridd_control, ssridd_control_descripton, "disease")
count_normal_cdls_control, deseq_results_cdls_control = DESeq2Wrapper(cdls_control, cdls_control_descripton, "disease")

 
# Saving ssridd vs control results
count_normal_ssridd_control.to_csv('counts_normal_ssridd_control.csv')
deseq_results_ssridd_control.to_csv('deseq_result2_ssridd_control.csv')

# Saving cdls vs control results
count_normal_cdls_control.to_csv('counts_normal_cdls_control.csv')
deseq_results_cdls_control.to_csv('deseq_result2_cdls_control.csv')

 
# Loading files: raw gene counts, description.xlsx
from google.colab import files
uploaded = files.upload()

 
# Read in deseq results and normalized counts 
# ssridd
ssridd_deseq_results = pd.read_csv('deseq_result2_ssridd_control.csv', index_col = 0)
ssridd_counts = pd.read_csv('counts_normal_ssridd_control.csv', index_col = 0)
ssridd_deseq_results.index.name = None
ssridd_counts_w_description = ssridd_counts.merge(ssridd_control_descripton, left_index=True, right_index=True, how='outer').reset_index().set_index('sample_title')
ssridd_counts_w_description.index.name=None


#clds
cdls_deseq_results = pd.read_csv('deseq_result2_cdls_control.csv', index_col = 0)
cdls_counts = pd.read_csv('counts_normal_cdls_control.csv', index_col = 0)
cdls_deseq_results.index.name = None
cdls_counts_w_description = cdls_counts.merge(cdls_control_descripton, left_index=True, right_index=True, how='outer').reset_index().set_index('sample_title')
cdls_counts_w_description.index.name=None

 
# Find the most significant genes from cdls deseq results 
significant_genes = cdls_deseq_results[(cdls_deseq_results['padj'] < 0.05)]
significant_genes_ordered = significant_genes.sort_values(by='log2FoldChange', key=abs, ascending = False)
significant_genes_ordered = significant_genes_ordered.assign(rank=range(len(significant_genes_ordered)))

 
plt.scatter(significant_genes_ordered['rank'], abs(significant_genes_ordered['log2FoldChange']))
plt.title("Elbow plot for fold changes")
plt.xlabel("Rank")
plt.ylabel("Absolute log2 Fold change")
plt.show()

 
# Filter and identify gene of interest according to padj and log2fold change threshold
goi = cdls_deseq_results[(cdls_deseq_results['padj'] < 0.05) & (abs(cdls_deseq_results['log2FoldChange']) >= 2)]
goi = goi.sort_values(by = 'log2FoldChange', key=abs, ascending = False)
goi_25 = goi.iloc[:25,]

# Get filtered counts for gene of interst 
'''
goi_counts2 = cdls_counts_w_description.iloc[:,2:15193].loc[:,goi_25.index]
'''
counts_filtered_w_description2 = counts_filtered_w_description.reset_index().set_index('sample_title')
counts_filtered_w_description2.index.name=None
#goi_counts = counts_filtered_w_description2.iloc[:,2:15193].loc[:,goi_25.index]
goi_counts = cdls_control.merge(cdls_control_descripton, left_index=True, right_index=True).reset_index().set_index('sample_title')
goi_counts.index.name=None
goi_counts = goi_counts.iloc[:,2:15193].loc[:,goi_25.index]

 
cdls_control_descripton = cdls_control_descripton.reset_index().set_index('sample_title')
cdls_control_descripton.index.name=None

   
## Visualization 1: Clustermap for CdLS vs Control
goi_counts_log2 = np.log2(goi_counts + 0.5)
goi_counts_log2_trip_mean = goi_counts_log2.iloc[18:,:].mean()
goi_counts_log2_diff = -(goi_counts_log2 - goi_counts_log2_trip_mean)

 
diseasetype = set(cdls_control_descripton.disease)
disease_color = dict(zip(diseasetype,sns.color_palette()))
disease = cdls_control_descripton.loc[:,'disease']
row_colors1 = disease.map(disease_color)


 
# Comparison between each disease with control 
clustermap = sns.clustermap(goi_counts_log2_diff, row_colors=row_colors1,cmap = 'coolwarm', row_cluster=False)
#clustermap.figure.savefig("clustermap_cdls_control.png")

   
## Visualization 2: Volcano plot for most significant fold change genes
cdls_counts = goi_counts
cdls_goi_deseq_result = cdls_deseq_results.loc[cdls_counts.columns]
ssridd_goi_deseq_result = ssridd_deseq_results.loc[cdls_counts.columns]

 
cdls_goi_deseq_result['-log10_padj'] = -np.log10(cdls_goi_deseq_result['padj'])
ssridd_goi_deseq_result['-log10_padj'] = -np.log10(ssridd_goi_deseq_result['padj'])

 
cdls_goi_deseq_result['padj_significant'] = np.where(cdls_goi_deseq_result['padj'] < 0.05, True, False)
ssridd_goi_deseq_result['padj_significant'] = np.where(ssridd_goi_deseq_result['padj'] < 0.05, True, False)

 
cdls_goi_deseq_result = cdls_goi_deseq_result.add_suffix('_cdls')
ssridd_goi_deseq_result = ssridd_goi_deseq_result.add_suffix('_ssridd')

 
combined = pd.concat([cdls_goi_deseq_result, ssridd_goi_deseq_result], axis=1).dropna()


 
DESeq2Results_AllPairwise = combined
interval = alt.selection_interval(encodings = ['x','y'])

DESeq2Results_AllPairwise_volcano = alt.Chart(DESeq2Results_AllPairwise).mark_point().encode(
    x = "log2FoldChange_cdls",
    y = "-log10_padj_cdls",
    color = alt.condition(interval, "padj_significant_ssridd", alt.value('gray')),

).properties(
    selection = interval
)

DESeq2Results_AllPairwise_volcano | DESeq2Results_AllPairwise_volcano.encode(x = "log2FoldChange_ssridd", y = "-log10_padj_ssridd")

   
## Visualization 3: PCA plot
counts_w_description = counts_filtered.merge(description, left_index=True, right_index=True, how='outer').reset_index().set_index('sample_title')
counts_w_description.index.name=None
diseasetypes = set(counts_w_description.disease)

 
# Convert the data and transform data for PCA 
from   sklearn.decomposition import PCA
from   sklearn.manifold import TSNE
from   sklearn.preprocessing import StandardScaler

 
goi_counts = counts_filtered.loc[:,cdls_counts.columns.to_list()]

 
goi_counts_scaled = StandardScaler().fit_transform(goi_counts)
goi_pca = PCA().fit(goi_counts_scaled)
goi_pca_data = goi_pca.transform(goi_counts_scaled)

 
# Plotting PCA 
fig = plt.figure(figsize = (8,5))
ax = fig.add_subplot(111)

for disease in diseasetypes:
  ax.scatter(goi_pca_data[counts_filtered_w_description['disease'] == disease, 0], 
             goi_pca_data[counts_filtered_w_description['disease'] == disease, 1], 
             alpha=0.75, marker=".", s=50, label=disease)
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.legend(loc="best")

 
# Genes most strongly associated with PC1
pc1_loadings = goi_pca.components_[0, :].T
pc1_weighted_loadings = pc1_loadings* goi_pca.explained_variance_[0]

gene_names = counts_w_description.iloc[:,2:15194].loc[:,goi.index].columns.values

n_top_genes=1
top_pc1_gene_idxs = np.argsort(np.abs(pc1_weighted_loadings))[-n_top_genes:]
for gene_idx in top_pc1_gene_idxs:
  print(gene_names[gene_idx])
  ax.annotate(text = gene_names[gene_idx],
              xy = [0, 0],
              xytext = [goi_pca.components_.T[gene_idx, 0] * goi_pca.explained_variance_[0], 
                        goi_pca.components_.T[gene_idx, 1] * goi_pca.explained_variance_[1]],
              arrowprops=dict(arrowstyle='<-',linewidth=1, shrinkA=0.9)) 


 
# Genes most strongly associated with PC1 and PC2 
pc2_loadings = goi_pca.components_[1, :].T
pc2_weighted_loadings = pc2_loadings* goi_pca.explained_variance_[2]


gene_names = counts_w_description.iloc[:,2:15194].loc[:,goi.index].columns.values

n_top_genes=1
top_pc2_gene_idxs = np.argsort(np.abs(pc2_weighted_loadings))[-n_top_genes:]

for gene_idx in top_pc2_gene_idxs:
  print(gene_names[gene_idx])
  ax.annotate(text = gene_names[gene_idx],
              xy = [0, 0],
              xytext = [goi_pca.components_.T[gene_idx, 0] * goi_pca.explained_variance_[0], 
                        goi_pca.components_.T[gene_idx, 1] * goi_pca.explained_variance_[1]],
              arrowprops=dict(arrowstyle='<-',linewidth=1, shrinkA=0.9)) 
  
fig.figure.savefig("PCA_cdls_control.png")
fig

   
## Most significant genes
goi.sort_values(by='padj')
ssridd_deseq_results.iloc['padj'] < 0.05
ssridd_deseq_results.loc[goi_25.index]
ssridd_counts_w_description.loc[:,cdls_top_10]

 



