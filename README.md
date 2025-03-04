# Local Ancestry StructLMM for GxG Interactions

## Overview

This repository contains a methodology that adapts the StructLMM framework to study Gene-Gene (GxG) interactions instead of Gene-Environment (GxE) interactions. Our approach leverages local ancestry Principal Components (PCs) as a proxy for the "environment" in the StructLMM model.

## Methodology

We extend the StructLMM framework by Moore et al. (2018) to detect GxG interactions using the following inputs:

1. A query SNV that has a large effect size for a phenotype
2. A query region of interest to test for interactions (can be cis or trans)
3. The phenotype

The key innovation is using local ancestry Principal Components as the "environment" matrix in the StructLMM model. This allows us to capture interaction effects between the query SNV and genetic variants in the region of interest.

### Model

The adapted model is structured as:

$$\mathbf{y} = \mathbf{M}\boldsymbol{\alpha} + g\boldsymbol{\beta_0} + g\circ\boldsymbol{\beta_1} + \mathbf{e} + \boldsymbol{\varepsilon},$$

where:

$$\boldsymbol{\beta_0} \sim \mathcal{N}(\mathbf{0}, r_0 \cdot \rho\mathbf{I}), \quad \boldsymbol{\beta_1} \sim \mathcal{N}(\mathbf{0}, r_0(1-\rho)\mathbf{EE}^\top), \quad \mathbf{e} \sim \mathcal{N}(\mathbf{0}, r_1\mathbf{EE}^\top), \quad \text{and} \quad \boldsymbol{\varepsilon} \sim \mathcal{N}(\mathbf{0}, r_2\mathbf{I}).$$

- **y** is the phenotype vector
- **M** contains covariates
- **g** is the query SNV with large effect size
- **E** is a matrix representing local ancestry PCs from the query region
- $\rho \in [0, 1]$ dictates the relevance of the interaction effect

## Pipeline

### Step 1: Extract Local Ancestry PCs

We extract Principal Components from a specified genomic region to capture local ancestry patterns:

```bash
./local_ancestry_pc_extraction.sh input.vcf.gz chr1:1000000-2000000 output_prefix
```

The script:
1. Extracts variants from the region of interest
2. Performs quality control
3. Calculates PCs that explain up to 90% of variance
4. Outputs PC coordinates and variance explained

### Step 2: Run StructLMM with Local Ancestry PCs

We then apply the StructLMM framework using:
- The query SNV as the genetic variant
- Local ancestry PCs as the "environment"
- The phenotype of interest

## Requirements

- PLINK 2.0
- BCFtools
- Python 3.6+
  - numpy
  - pandas
  - scipy
  - [StructLMM](https://github.com/limix/struct-lmm)

## Usage

```bash
# Step 1: Extract local ancestry PCs
./local_ancestry_pc_extraction.sh input.vcf.gz chr1:1000000-2000000 output_prefix

# Step 2: Run GxG analysis
python run_gxg_analysis.py \
  --pcs output_prefix_local_ancestry_pcs.csv \
  --query-snv chr6_32609603_A_G_b38 or 6_32609603_A_G_b38 \
  --phenotype phenotype.csv \
  --covariates covariates.csv \
  --output results.csv
```

## Reference Implementation

This methodology builds upon the StructLMM framework. The original StructLMM implementation is available at: https://github.com/limix/struct-lmm

## Citation

If you use this method in your research, please cite:

Original StructLMM paper:
> Moore, R., Casale, F. P., Bonder, M. J., Horta, D., Franke, L., Barroso, I., & Stegle, O. (2018). A linear mixed-model approach to study multivariate gene–environment interactions. Nature Genetics, 50(7), 1167-1174.


