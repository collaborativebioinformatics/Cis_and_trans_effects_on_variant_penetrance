#!/bin/bash

# Local Ancestry Principal Component Extraction Pipeline for PLINK format data
# Requires: PLINK2, Python with numpy/pandas

# Argument parsing
if [ $# -lt 3 ]; then
    echo "Usage: $0 <plink_files_directory> <region> <output_prefix>"
    echo "Example: $0 /path/to/plink_files chr1:1000000-2000000 output_analysis"
    exit 1
fi

# Input parameters
PLINK_DIR=$1      # Directory containing PLINK files (.bed, .bim, .fam)
REGION=$2         # Genomic region (e.g., chr1:1000000-2000000)
OUTPUT_PREFIX=$3  # Output file prefix

# Logging
LOG_FILE="${OUTPUT_PREFIX}_local_ancestry_extraction.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "Starting Local Ancestry PC Extraction Pipeline"
echo "PLINK files directory: $PLINK_DIR"
echo "Region: $REGION"
echo "Output Prefix: $OUTPUT_PREFIX"
echo "Date: $(date)"

# Parse the region to get chromosome
CHROM=$(echo $REGION | cut -d':' -f1)
RANGE=$(echo $REGION | cut -d':' -f2)
START=$(echo $RANGE | cut -d'-' -f1)
END=$(echo $RANGE | cut -d'-' -f2)

echo "Parsed region - Chromosome: $CHROM, Start: $START, End: $END"

# Find the appropriate chromosome file in the directory
# Try common naming patterns
PLINK_PREFIX=""
PATTERNS=("$CHROM" "${CHROM#chr}" "chr${CHROM#chr}")

for pattern in "${PATTERNS[@]}"; do

    FOUND=$(find $PLINK_DIR -name "*.bed" | grep -E "(^|[^0-9a-zA-Z])${pattern}([^0-9a-zA-Z]|$)" | head -1)
    if [ -n "$FOUND" ]; then
        PLINK_PREFIX="${FOUND%.bed}"
        echo "Found PLINK files with precise matching, prefix: $PLINK_PREFIX"
        break
    fi

done

if [ -z "$PLINK_PREFIX" ]; then
    echo "Error: Could not find PLINK files for chromosome $CHROM in directory $PLINK_DIR"
    exit 1
fi

# 1. Extract specific region using PLINK
echo "Extracting region-specific variants..."
plink2 \
    --bfile "$PLINK_PREFIX" \
    --chr ${CHROM#chr} \
    --from-bp $START \
    --to-bp $END \
    --make-bed \
    --out "${OUTPUT_PREFIX}_region"

# Check if binary files were created
if [ ! -f "${OUTPUT_PREFIX}_region.bed" ] || [ ! -f "${OUTPUT_PREFIX}_region.bim" ] || [ ! -f "${OUTPUT_PREFIX}_region.fam" ]; then
    echo "Error: Region extraction failed"
    exit 1
fi

# 2. Apply QC filters
echo "Applying quality control filters..."
plink2 \
    --bfile "${OUTPUT_PREFIX}_region" \
    --set-all-var-ids '@_#_$r_$a' \
    --max-alleles 2 \
    --geno 0.1 \
    --maf 0.05 \
    --hwe 1e-6 \
    --make-bed \
    --out "${OUTPUT_PREFIX}_qc"

# Check if binary files were created
if [ ! -f "${OUTPUT_PREFIX}_qc.bed" ] || [ ! -f "${OUTPUT_PREFIX}_qc.bim" ] || [ ! -f "${OUTPUT_PREFIX}_qc.fam" ]; then
    echo "Error: PLINK QC filtering failed"
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
# 6. Cleanup intermediate files
echo "Cleaning up intermediate files..."
# Remove region files
rm -f "${OUTPUT_PREFIX}_region"*

# Remove QC files
rm -f "${OUTPUT_PREFIX}_qc.bed"
rm -f "${OUTPUT_PREFIX}_qc.bim"
rm -f "${OUTPUT_PREFIX}_qc.fam"
rm -f "${OUTPUT_PREFIX}_qc.log"

# Remove pruning files
rm -f "${OUTPUT_PREFIX}_pruned"*

# Keep only these important output files:
# - ${OUTPUT_PREFIX}_local_ancestry_pcs.csv
# - ${OUTPUT_PREFIX}_variance_explained.csv
# - ${OUTPUT_PREFIX}_pca.eigenval
# - ${OUTPUT_PREFIX}_pca.eigenvec
# - ${OUTPUT_PREFIX}_pca.log

# Remove log file if requested
# rm -f "${OUTPUT_PREFIX}_local_ancestry_extraction.log"


echo "Local Ancestry PC Extraction Complete!"
echo "Results in: ${OUTPUT_PREFIX}_local_ancestry_pcs.csv"
echo "Variance explained in: ${OUTPUT_PREFIX}_variance_explained.csv"
