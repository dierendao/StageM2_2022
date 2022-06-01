#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# import modules

import sys, getopt 
from pathlib import Path 
import pandas as pd 
import stream as st 
import pickle 
import matplotlib.pyplot as plt 
import numpy as np 
import re as re 
import random

random.seed(1993)


# ## Single Cell Trajectory REconstruction, Exploration And Mapping of omics data (STREAM)

# In[2]:


st.__version__


# In[3]:


st.set_figure_params(dpi=80,style='white',figsize=[5.4,4.8], rc={'image.cmap': 'viridis'})


# In[4]:


# Parameters and paths to files

epgAlpha = 0.01    
epgMu = 0.05     
epgLambda = 0.05 
inputFile = "/home/ndao/Images/prepareSTREAM/All_scaleDataForStream.tsv"   
conditionCells ="/home/ndao/Images/prepareSTREAM/All_seuratClusters.tsv"     
conditionColors ="/home/ndao/Images/prepareSTREAM/All_colorClusters.tsv"
outputDir = "/home/ndao/Images/streamAnalysis"


# In[5]:


adata=st.read(file_name = inputFile, workdir = outputDir)   
st.add_cell_labels(adata, file_name=conditionCells)   
st.add_cell_colors(adata, file_name=conditionColors)


# In[6]:


# Dimmension of the data

adata.obs.shape


# In[7]:


st.cal_qc(adata, assay='rna')


# In[8]:


#select top principal components 
st.select_top_principal_components(adata, first_pc=True, n_pc=15)


# # Dimension reduction

# In[9]:


st.dimension_reduction(adata, method="mlle", feature="top_pcs", n_components=4, n_jobs=12, n_neighbors=50)


# In[23]:


adata.obs.columns
adata.var_names


# In[10]:


st.plot_dimension_reduction(adata,color=['label','Gata1','n_genes'],
                            n_components=2,show_graph=False,show_text=False)


# #Trajectory inference

# In[11]:


st.seed_elastic_principal_graph(adata, n_clusters=10)


# In[27]:


st.plot_dimension_reduction(adata,color=['label','Gata1','n_genes'],n_components=2,show_graph=True,show_text=False)
st.plot_branches(adata,show_text=True)


# In[12]:


st.elastic_principal_graph(adata=adata, epg_alpha=epgAlpha, epg_mu=epgMu, epg_lambda=epgLambda)   


# In[13]:


st.plot_dimension_reduction(adata,color=['label','Gata1','n_genes'],n_components=2,show_graph=True,show_text=False)
st.plot_branches(adata,show_text=False)


#  # Adjusting trajectories (optional)

# In[14]:


###Extend leaf branch to reach further cells 
st.extend_elastic_principal_graph(adata, epg_ext_mode='WeigthedCentroid',epg_ext_par=0.8)
st.plot_dimension_reduction(adata,color=['label'],n_components=2,show_graph=True,show_text=True)
st.plot_branches(adata,show_text=True)


# # Trajectory visualization

# #flat tree

# In[15]:


st.plot_flat_tree(adata,color=['label','branch_id_alias','S3_pseudotime'],
                  dist_scale=0.5,show_graph=True,show_text=True)


# #stream plot at single cell level

# In[16]:


st.plot_stream_sc(adata,root='S3',color=['label','Ctsg'],
                  dist_scale=0.3,show_graph=True,show_text=True)


# #stream plots

# In[18]:


st.plot_stream(adata,root='S3',color=['label','Ctsg'])


# In[19]:


#Save plot as png
st.plot_stream(adata, root="S3", color=['S3_pseudotime'], log_scale=True, factor_zoomin=200, save_fig=True, fig_path=outputDir, fig_format="png")


# # Marker genes detection

# 1) detect marker genes for each leaf branch

# In[46]:


st.detect_leaf_markers(adata,cutoff_zscore=1.0,cutoff_pvalue=0.01,
                       root='S3',n_jobs=4)


# In[47]:


adata.uns['leaf_markers_all'].head()


# In[48]:


adata.uns['leaf_markers'][('S1','S2')].head()


# 2) detect transition genes for each branch

# In[50]:


st.detect_transition_markers(adata,cutoff_spearman=0.4,cutoff_logfc=0.25,
                             root='S3',n_jobs=4)


# In[52]:


adata.uns['transition_markers'][('S3','S1')].head()


# In[53]:


st.plot_transition_markers(adata,fig_size=(10,5))


# 3) detect marker genes that are differentially expressed between pairs of branches

# In[54]:


st.detect_de_markers(adata,cutoff_zscore=1,cutoff_logfc=0.25,
                     root='S3',n_jobs=4)


# In[56]:


adata.uns['de_markers_greater'][(('S1', 'S2'), ('S1', 'S0'))].head()


# In[57]:


adata.uns['de_markers_less'][(('S1', 'S2'), ('S1', 'S0'))].head()


# In[58]:


st.plot_de_markers(adata)


# In[66]:


#adata.obs.to_csv("".join([outputDir,"stream_metadata2.csv"]))


# In[65]:


filepath = Path('/home/ndao/Images/streamAnalysis/stream_metadata.csv')
adata.obs.to_csv(filepath)

