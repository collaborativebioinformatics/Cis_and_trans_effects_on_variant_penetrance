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
import pysam
import io
from struct_lmm import StructLMM



HAP_P  = "dataPrepScripts/outputFromPlink_local_ancestry_pcs.csv"
PHE_P = "testPlink/synthetic_small_v1.pheno1"
SUBSET = "Phenotype(liability)"
SNV_P = "dataPrepScripts/SNVoutputFromPlink.txt"
HAPS = pd.read_csv(HAP_P, header=1, index_col=0)
PHE = pd.read_csv(PHE_P, sep="\t")
SNV = pd.read_csv(SNV_P, sep="\t")
y = PHE[SUBSET].to_numpy()
E = HAPS.iloc[:, 1:].to_numpy()
M = np.concatenate([np.ones((E.shape[0], 1)), E], axis=1)
g = SNV.iloc[:, 1].to_numpy()
lmm = StructLMM(y, M, E)
lmm.fit()
lmm.score_2dof_inter(g)