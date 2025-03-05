#!/usr/bin/env python
# coding: utf-8

# In[24]:


# read vcf data
import sys
import os
import re
import argparse
import pandas as pd
import numpy as np  
import gzip
#import pysam
import io
from struct_lmm import StructLMM


# In[25]:


HAP_P  = "dataPrepScripts/PCoutput_local_ancestry_pcs.csv"
PHE_P = "testPlink/synthetic_small_v1.pheno1"
SUBSET = "Phenotype(binary)"
SNV_P = "dataPrepScripts/SNVoutput.txt"


# In[26]:


HAPS = pd.read_csv(HAP_P, header=1, index_col=0)
PHE = pd.read_csv(PHE_P, sep="\t")
SNV = pd.read_csv(SNV_P, sep="\t")


# In[27]:


y = PHE[SUBSET].to_numpy()

#%%
y=y-np.mean(y)
y=y/np.std(y)


#%%
def standardize_columns(data):

    # Calculate mean for each column
    means = data.mean(axis=0)
    
    # Calculate standard deviation for each column
    stds = data.std(axis=0)
    
    # Handle potential zero standard deviation (avoid division by zero)
    stds[stds == 0] = 1.0
    
    # Standardize the data
    return((data - means) / stds)



# In[28]:


E = HAPS.iloc[:, 1:].to_numpy()
#%%
E=standardize_columns(E)
E=E/np.sqrt(E.shape[1])
# In[29]:


M = np.concatenate([np.ones((E.shape[0], 1)), E], axis=1)


# In[30]:


g = SNV.iloc[:, 1].to_numpy()


#%%
def ensure_minor_allele_coding(genotype_vector):
    """
    Ensures genotype vector has minor allele coded as 1/2 by flipping if needed.
    Flips coding if there are more 2s than 0s.
    """
    import numpy as np
    geno = np.array(genotype_vector)
    
    # Count 0s and 2s
    zeros = np.sum(geno == 0)
    twos = np.sum(geno == 2)
    
    # Flip if there are more 2s than 0s
    if twos > zeros:
        flipped = np.zeros_like(geno)
        flipped[geno == 0] = 2
        flipped[geno == 1] = 1
        flipped[geno == 2] = 0
        return flipped
    else:
        return geno

#%%
g=ensure_minor_allele_coding(g)

# In[31]:

y, M, E, g, y.shape, E.shape, M.shape, g.shape



# %%time
lmm = StructLMM(y, M, E)
lmm.fit()
lmm.score_2dof_inter(g)

