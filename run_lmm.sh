#!/bin/bash

# Assign input arguments to variables
HAP_P=$1
PHE_P=$2
SUBSET=$3
SNV_P=$4

# Run Python script with provided arguments
python3 << EOF
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



# Assign variables
HAP_P = "$HAP_P"
PHE_P = "$PHE_P"
SUBSET = "$SUBSET"
SNV_P = "$SNV_P"

# Load data
HAPS = pd.read_csv(HAP_P, header=1, index_col=0)
PHE = pd.read_csv(PHE_P, sep="\t")
SNV = pd.read_csv(SNV_P, sep="\t")

# Prepare data
y = PHE[SUBSET].to_numpy()
E = HAPS.iloc[:, 1:].to_numpy()
M = np.concatenate([np.ones((E.shape[0], 1)), E], axis=1)
g = SNV.iloc[:, 1].to_numpy()

# Run analysis
lmm = StructLMM(y, M, E)
lmm.fit()
score = lmm.score_2dof_inter(g)
print(score)
EOF
