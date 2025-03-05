#!/usr/bin/env python3

import argparse
import sys
import os
import re
import pandas as pd
import numpy as np
from struct_lmm import StructLMM

def standardize_columns(data):
    """Standardize each column of a matrix."""
    means = data.mean(axis=0)
    stds = data.std(axis=0)
    stds[stds == 0] = 1.0
    return (data - means) / stds

def ensure_minor_allele_coding(genotype_vector):
    """Ensure genotype vector has minor allele coded as 1/2."""
    geno = np.array(genotype_vector)
    zeros = np.sum(geno == 0)
    twos = np.sum(geno == 2)
    
    if twos > zeros:
        print(f"Flipping SNV coding (more 2s than 0s: {twos} vs {zeros})")
        flipped = np.zeros_like(geno)
        flipped[geno == 0] = 2
        flipped[geno == 1] = 1
        flipped[geno == 2] = 0
        return flipped
    else:
        print(f"SNV coding is correct (more 0s than 2s: {zeros} vs {twos})")
        return geno

def main():
    parser = argparse.ArgumentParser(description="Run StructLMM analysis for GxG interactions using local ancestry PCs")
    
    # Required arguments
    parser.add_argument("--pcs", required=True, help="File containing local ancestry PCs")
    parser.add_argument("--snv", required=True, help="File containing query SNV genotypes")
    parser.add_argument("--phenotype", required=True, help="File containing phenotype data")
    
    # Optional arguments
    parser.add_argument("--phenotype-column", help="Name of the phenotype column (default: first column with 'phenotype' in name)")
    parser.add_argument("--covariates", help="File containing covariates (optional)")
    parser.add_argument("--output", help="File to save results (default: stdout)")
    parser.add_argument("--sep", default=None, help="Separator for input files (default: auto-detect)")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")
    
    args = parser.parse_args()
    
    # Read input files
    try:
        # Read PCs file
        print(f"Reading local ancestry PCs from {args.pcs}")
        # Handle special case with two header rows, where second starts with #
        with open(args.pcs, 'r') as f:
            first_lines = [next(f) for _ in range(3) if f]
        
        if len(first_lines) >= 2 and first_lines[1].startswith('#'):
            print("Detected PC file with comment line as second row")
            # Skip the second line (comment line)
            pcs_df = pd.read_csv(args.pcs, header=0, comment='#', index_col=1, sep=args.sep, engine='python')
        else:
            try:
                pcs_df = pd.read_csv(args.pcs, header=0, index_col=0, sep=args.sep, engine='python')
            except:
                # Try with header in second row
                pcs_df = pd.read_csv(args.pcs, header=1, index_col=0, sep=args.sep, engine='python')
        
        # Read phenotype file
        print(f"Reading phenotype data from {args.phenotype}")
        pheno_df = pd.read_csv(args.phenotype, sep=args.sep, engine='python')
        
        # Read SNV file
        print(f"Reading SNV data from {args.snv}")
        snv_df = pd.read_csv(args.snv, sep=args.sep, engine='python')
        
        # Read covariates if provided
        if args.covariates:
            print(f"Reading covariates from {args.covariates}")
            cov_df = pd.read_csv(args.covariates, sep=args.sep, engine='python')
    except Exception as e:
        print(f"Error reading input files: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Determine phenotype column
    if args.phenotype_column:
        phenotype_column = args.phenotype_column
    else:
        # Find first column with 'phenotype' in name (case insensitive)
        phenotype_columns = [col for col in pheno_df.columns if 'phenotype' in col.lower()]
        if phenotype_columns:
            phenotype_column = phenotype_columns[0]
            print(f"Using phenotype column: {phenotype_column}")
        else:
            # If no column with 'phenotype', use the second column (assuming first is ID)
            if len(pheno_df.columns) > 1:
                phenotype_column = pheno_df.columns[1]
                print(f"No column with 'phenotype' found. Using column: {phenotype_column}")
            else:
                print("Error: Cannot determine phenotype column", file=sys.stderr)
                sys.exit(1)
    
    # Check if phenotype column exists
    if phenotype_column not in pheno_df.columns:
        print(f"Error: Phenotype column '{phenotype_column}' not found in phenotype file", file=sys.stderr)
        sys.exit(1)
    
    # Extract and standardize phenotype
    y = pheno_df[phenotype_column].to_numpy()
    orig_y_mean = np.mean(y)
    orig_y_std = np.std(y)
    print(f"Phenotype '{phenotype_column}' mean: {orig_y_mean:.3f}, std: {orig_y_std:.3f}")
    
    # Standardize phenotype
    y = (y - orig_y_mean) / orig_y_std
    
    # Process and extract PC data (environment matrix)
    if args.verbose:
        print(f"PC data shape before processing: {pcs_df.shape}")
        print(f"PC columns: {', '.join(pcs_df.columns)}")
    
    # Extract PC columns only (skip ID columns if they exist)
    pc_cols = [col for col in pcs_df.columns if col.startswith('PC')]
    
    if not pc_cols:
        print("Warning: No columns starting with 'PC' found. Using all numeric columns.")
        pc_cols = pcs_df.select_dtypes(include=[np.number]).columns
    
    if len(pc_cols) < len(pcs_df.columns):
        print(f"Using these PC columns: {', '.join(pc_cols)}")
        pcs_df = pcs_df[pc_cols]
    
    E = pcs_df.to_numpy()
    if args.verbose:
        print(f"E shape after extraction: {E.shape}")
    
    # Standardize PC matrix
    E = standardize_columns(E)
    E = E / np.sqrt(E.shape[1])
    print(f"Using {E.shape[1]} PCs from local ancestry")
    
    # Create covariate matrix
    if args.covariates:
        # Extract covariates (assuming first column is ID)
        cov_data = cov_df.iloc[:, 1:].to_numpy()
        print(f"Using {cov_data.shape[1]} covariates")
        M = np.concatenate([np.ones((y.shape[0], 1)), cov_data], axis=1)
    else:
        # Use intercept only if no covariates
        print("No covariates provided, using intercept only")
        M = np.ones((y.shape[0], 1))
    
    # Process SNV data
    # Assuming first column is ID and second column is genotype
    if len(snv_df.columns) < 2:
        print("Error: SNV file must have at least 2 columns (ID and genotype)", file=sys.stderr)
        sys.exit(1)
    
    g = snv_df.iloc[:, 1].to_numpy()
    
    # Ensure minor allele coding
    g = ensure_minor_allele_coding(g)
    
    # Check for size consistency
    if not (len(y) == E.shape[0] == M.shape[0] == len(g)):
        print(f"Error: Inconsistent data sizes: y: {len(y)}, E: {E.shape[0]}, M: {M.shape[0]}, g: {len(g)}", file=sys.stderr)
        sys.exit(1)
    
    # Print dimensions
    print(f"Data dimensions: samples: {len(y)}, PCs: {E.shape[1]}, covariates: {M.shape[1]-1}")
    
    # Run StructLMM analysis
    print("\nRunning StructLMM analysis...")
    try:
        lmm = StructLMM(y, M, E)
        lmm.fit()
        
        # In the basic installation, score_2dof_inter only returns p-value
        pv = lmm.score_2dof_inter(g)
        
        # Try to get alternate and null likelihoods if available
        try:
            # For basic install, this might not be available
            null_lml = lmm.null_lml
            alt_lml = lmm.alt_lml
            if null_lml is not None and alt_lml is not None:
                var_explained = 1 - np.exp((null_lml - alt_lml) / y.shape[0])
            else:
                var_explained = None
        except:
            var_explained = None
        
        # Prepare results
        results = {
            "p_value": pv,
            "n_samples": y.shape[0],
            "n_pcs": E.shape[1],
            "phenotype": phenotype_column,
            "phenotype_mean": orig_y_mean,
            "phenotype_std": orig_y_std
        }
        
        # Add variance explained if available
        if var_explained is not None:
            results["variance_explained"] = var_explained
        
        # Create results DataFrame
        results_df = pd.DataFrame([results])
        
        # Print or save results
        print("\nResults:")
        print(f"P-value: {pv:.8e}")
        if var_explained is not None:
            print(f"Variance explained: {var_explained:.4f}")
        
        if args.output:
            print(f"Saving results to {args.output}")
            results_df.to_csv(args.output, index=False)
        else:
            print("\nFull results:")
            print(results_df.to_string(index=False))
        
    except Exception as e:
        print(f"Error in StructLMM analysis: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
