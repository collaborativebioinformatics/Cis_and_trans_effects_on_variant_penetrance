#!/bin/bash

# Local Ancestry Principal Component Extraction Pipeline
# Requires: PLINK, BCFtools, Python with numpy/pandas

# Argument parsing
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input_vcf> <region> <output_prefix>"
    echo "Example: $0 data.vcf.gz chr1:1000000-2000000 output_analysis"
    exit 1
fi

# Input parameters
INPUT_VCF=$1      # Input VCF file
REGION=$2         # Genomic region (e.g., chr1:1000000-2000000)
OUTPUT_PREFIX=$3  # Output file prefix

# Logging
LOG_FILE="${OUTPUT_PREFIX}_local_ancestry_extraction.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "Starting Local Ancestry PC Extraction Pipeline"
echo "Input VCF: $INPUT_VCF"
echo "Region: $REGION"
echo "Output Prefix: $OUTPUT_PREFIX"
echo "Date: $(date)"

# 1. Extract specific region using BCFtools
echo "Extracting region-specific variants..."
bcftools view -r "$REGION" "$INPUT_VCF" -O z -o "${OUTPUT_PREFIX}_region.vcf.gz"
bcftools index "${OUTPUT_PREFIX}_region.vcf.gz"

# 2. Convert VCF to PLINK binary format with enhanced variant ID setting
echo "Converting VCF to PLINK binary format..."
plink2 \
    --vcf "${OUTPUT_PREFIX}_region.vcf.gz" \
    --set-all-var-ids '@_#_$r_$a' \
    --max-alleles 2 \
    --geno 0.1 \
    --maf 0.05 \
    --hwe 1e-6 \
    --make-bed \
    --out "${OUTPUT_PREFIX}_qc"

# Check if binary files were created
if [ ! -f "${OUTPUT_PREFIX}_qc.bed" ] || [ ! -f "${OUTPUT_PREFIX}_qc.bim" ] || [ ! -f "${OUTPUT_PREFIX}_qc.fam" ]; then
    echo "Error: PLINK binary file conversion failed"
    exit 1
fi

# 3. Prune variants for PCA
echo "Pruning variants for PCA..."
plink2 \
    --bfile "${OUTPUT_PREFIX}_qc" \
    --indep-pairwise 50 5 0.2 \
    --out "${OUTPUT_PREFIX}_pruned"

# 4. Calculate Principal Components 
echo "Calculating Principal Components..."
plink2 \
    --bfile "${OUTPUT_PREFIX}_qc" \
    --extract "${OUTPUT_PREFIX}_pruned.prune.in" \
    --pca 10 \
    --out "${OUTPUT_PREFIX}_pca"

# 5. Process PCA results with dynamic PC selection
echo "Processing PCA results..."
python3 - << EOF
import numpy as np
import pandas as pd

# Read eigenvalues
eigenvalues = np.loadtxt("${OUTPUT_PREFIX}_pca.eigenval")

# Calculate relative variance proportion
total_variance = np.sum(eigenvalues)
relative_variance = eigenvalues / total_variance
cumulative_variance = np.cumsum(relative_variance)

# Select PCs up to 90% variance
n_pcs = np.argmax(cumulative_variance >= 0.9) + 1
n_pcs = min(n_pcs, 10)  # Ensure not more than 10 PCs

print(f"Number of PCs selected: {n_pcs}")
print("\nVariance Explained:")
for i in range(n_pcs):
    print(f"PC{i+1}: {relative_variance[i]*100:.2f}% (Cumulative: {cumulative_variance[i]*100:.2f}%)")

# Read eigenvectors
pca_df = pd.read_csv("${OUTPUT_PREFIX}_pca.eigenvec", 
                     delim_whitespace=True, 
                     names=['FID', 'IID'] + [f'PC{i}' for i in range(1, 11)])

# Select and save appropriate number of PCs
selected_pcs = pca_df[['FID', 'IID'] + [f'PC{i}' for i in range(1, n_pcs+1)]]
selected_pcs.to_csv("${OUTPUT_PREFIX}_local_ancestry_pcs.csv", index=False)

# Save variance information
variance_df = pd.DataFrame({
    'PC': [f'PC{i+1}' for i in range(n_pcs)],
    'Eigenvalue': eigenvalues[:n_pcs],
    'Relative_Variance_Proportion': relative_variance[:n_pcs],
    'Cumulative_Variance_Proportion': cumulative_variance[:n_pcs]
})
variance_df.to_csv("${OUTPUT_PREFIX}_variance_explained.csv", index=False)
EOF

# 6. Cleanup intermediate files
echo "Cleaning up intermediate files..."
rm "${OUTPUT_PREFIX}_region.vcf.gz"
rm "${OUTPUT_PREFIX}_qc.log"
rm "${OUTPUT_PREFIX}_pruned"*

echo "Local Ancestry PC Extraction Complete!"
echo "Results in: ${OUTPUT_PREFIX}_local_ancestry_pcs.csv"
echo "Variance explained in: ${OUTPUT_PREFIX}_variance_explained.csv"
