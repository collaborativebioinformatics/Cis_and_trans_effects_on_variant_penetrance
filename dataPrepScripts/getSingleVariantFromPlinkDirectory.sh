#!/bin/bash

# Script to extract a single variant from PLINK files and output as 0,1,2 coding
# Requires: PLINK2

# Argument parsing
if [ $# -lt 3 ]; then
    echo "Usage: $0 <plink_files_directory> <variant_id> <output_file>"
    echo "Example: $0 /path/to/plink_files rs12345 variant_genotypes.txt"
    echo "Alternative with chr:pos: $0 /path/to/plink_files chr1:12345 variant_genotypes.txt"
    echo "Alternative with chr:pos:ref:alt: $0 /path/to/plink_files chr1:12345:A:G variant_genotypes.txt"
    exit 1
fi

# Input parameters
PLINK_DIR=$1      # Directory containing PLINK files (.bed, .bim, .fam)
VARIANT=$2        # Variant ID (rsID, chr:pos, or chr:pos:ref:alt)
OUTPUT_FILE=$3    # Output file name

echo "Starting Variant Extraction"
echo "PLINK files directory: $PLINK_DIR"
echo "Variant ID: $VARIANT"
echo "Output file: $OUTPUT_FILE"

# Parse the variant identifier to determine chromosome
if [[ $VARIANT == rs* ]]; then
    # For rsIDs, we need to search all BIM files to find the chromosome
    echo "Searching for rsID in BIM files..."
    
    CHROM=""
    FOUND_BIM=""
    for bim in "$PLINK_DIR"/*.bim; do
        if grep -q "$VARIANT" "$bim"; then
            # Extract chromosome from the first column of matching line
            CHROM=$(grep "$VARIANT" "$bim" | head -1 | awk '{print $1}')
            FOUND_BIM="$bim"
            echo "Found variant in chromosome $CHROM (file: $FOUND_BIM)"
            break
        fi
    done
    
    if [ -z "$CHROM" ]; then
        echo "Error: Variant $VARIANT not found in any BIM file"
        exit 1
    fi
    
    # Check if it's multi-allelic by counting occurrences in the BIM file
    ALLELE_COUNT=$(grep "$VARIANT" "$FOUND_BIM" | wc -l)
    if [ "$ALLELE_COUNT" -gt 1 ]; then
        echo "Error: Variant $VARIANT is multi-allelic with $ALLELE_COUNT alleles"
        echo "Showing all alleles:"
        grep "$VARIANT" "$FOUND_BIM"
        echo "Please specify a specific allele using chr:pos:ref:alt format"
        exit 1
    else
        VARIANT_SPEC="--snp $VARIANT"
    fi
else
    # For chr:pos format, extract chromosome
    CHROM=$(echo $VARIANT | cut -d':' -f1)
    
    # Handle chr:pos or chr:pos:ref:alt format
    if [[ $VARIANT == *:*:*:* ]]; then
        # chr:pos:ref:alt format - specific alleles
        POS=$(echo $VARIANT | cut -d':' -f2)
        REF=$(echo $VARIANT | cut -d':' -f3)
        ALT=$(echo $VARIANT | cut -d':' -f4)
        
        echo "Looking for specific alleles: REF=$REF, ALT=$ALT at position $POS"
        
        # Find the chromosome file to check if this specific allele combination exists
        PATTERNS=("$CHROM" "${CHROM#chr}" "chr${CHROM#chr}")
        FOUND_BIM=""

        for pattern in "${PATTERNS[@]}"; do
            # Try direct pattern matching first
            if ls "$PLINK_DIR/$pattern.bim" 2>/dev/null; then
                FOUND_BIM="$PLINK_DIR/$pattern.bim"
                break
            fi
            
            # Find all matching files and use regex for precise matching
            FOUND_FILES=$(find "$PLINK_DIR" -name "*.bim" | grep -E "(^|[^0-9a-zA-Z])${pattern}([^0-9a-zA-Z]|$)")
            
            # Count the number of matching files
            FOUND_COUNT=$(echo "$FOUND_FILES" | grep -v "^$" | wc -l)
            
            if [ "$FOUND_COUNT" -eq 1 ]; then
                # Exactly one match found
                FOUND_BIM="$FOUND_FILES"
                echo "Found a single BIM file $FOUND_BIM"
                break
            elif [ "$FOUND_COUNT" -gt 1 ]; then
                # Multiple matches found - show error and exit
                echo "Error: Found multiple matching files for chromosome pattern $pattern:"
                echo "$FOUND_FILES"
                exit 1
            fi
        done
        
        # Check if a BIM file was found, if not show error and list available files
        if [ -z "$FOUND_BIM" ]; then
            echo "Error: Could not find BIM file for chromosome $CHROM in directory $PLINK_DIR"
            echo "Available .bim files in directory:"
            find "$PLINK_DIR" -name "*.bim" | sort
            exit 1
        fi

        # Set PLINK_PREFIX based on the found BIM file
        PLINK_PREFIX="${FOUND_BIM%.bim}"
        echo "Using PLINK prefix: $PLINK_PREFIX"

        if [ -n "$FOUND_BIM" ]; then
            # Check if the exact allele combination exists
            MATCHES=$(awk -v pos="$POS" -v ref="$REF" -v alt="$ALT" '$4 == pos && (($5 == ref && $6 == alt) || ($5 == alt && $6 == ref)) {print}' "$FOUND_BIM")
            
            if [ -z "$MATCHES" ]; then
                echo "Warning: Specified alleles $REF/$ALT not found at position $POS"
                echo "Available variants at this position:"
                awk -v pos="$POS" '$4 == pos {print}' "$FOUND_BIM"
                
                echo "Error: The requested variant with specific alleles does not exist. Please choose from available variants."
                exit 1
            else
                # Found the specific allele combination
                echo "Found variant with requested alleles:"
                echo "$MATCHES"
                
                # Extract variant ID to use with PLINK
                VARIANT_ID=$(echo "$MATCHES" | head -1 | awk '{print $2}')
                VARIANT_SPEC="--snp $VARIANT_ID"
                echo "Selected variant: $VARIANT_ID"
            fi
        else
            echo "Error: Could not find chromosome file to verify alleles"
            exit 1
        fi
    else
        # chr:pos format - check for multiple variants
        POS=$(echo $VARIANT | cut -d':' -f2)
        
        # FOUND_BIM and PLINK_PREFIX were already set in the previous code block
        
        if [ -n "$FOUND_BIM" ]; then
            # Check if position has multiple variants
            MATCHES=$(awk -v pos="$POS" '$4 == pos {print}' "$FOUND_BIM")
            MATCH_COUNT=$(echo "$MATCHES" | grep -v "^$" | wc -l)
            
            if [ "$MATCH_COUNT" -eq 0 ]; then
                echo "Error: No variants found at position $POS"
                exit 1
            elif [ "$MATCH_COUNT" -gt 1 ]; then
                echo "Error: Position $POS has $MATCH_COUNT variants:"
                echo "$MATCHES"
                echo "Please specify a specific variant using chr:pos:ref:alt format"
                exit 1
            else
                # Single variant at this position
                VARIANT_SPEC="--chr ${CHROM#chr} --from-bp $POS --to-bp $POS"
            fi
        else
            echo "Error: Could not find chromosome file to check for variants at position $POS"
            exit 1
        fi
    fi
    
    echo "Using chromosome $CHROM and position $POS"
fi

# Create temporary directory
TEMP_DIR=$(mktemp -d)
TEMP_PREFIX="$TEMP_DIR/variant"

# Extract the variant using PLINK2 with recode A
echo "Extracting variant using: plink2 --bfile $PLINK_PREFIX $VARIANT_SPEC --recode A --out $TEMP_PREFIX"
plink2 \
    --bfile "$PLINK_PREFIX" \
    $VARIANT_SPEC \
    --recode A \
    --out "$TEMP_PREFIX"

# Check if extraction was successful
if [ ! -f "${TEMP_PREFIX}.raw" ]; then
    echo "Error: Failed to extract variant"
    echo "Check if the variant exists and if PLINK command is valid"
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Count the number of variants in the raw file
VARIANT_COLUMN_COUNT=$(head -1 "${TEMP_PREFIX}.raw" | wc -w)
VARIANT_COUNT=$(( $VARIANT_COLUMN_COUNT - 6 ))  # 6 standard columns before variants

if [ "$VARIANT_COUNT" -eq 0 ]; then
    echo "Error: No variants were extracted. The variant may not exist in the PLINK files."
    rm -rf "$TEMP_DIR"
    exit 1
fi

# Check if multiple variants were extracted (shouldn't happen with our modifications)
if [ "$VARIANT_COUNT" -gt 1 ]; then
    echo "Error: Multiple variants ($VARIANT_COUNT) were extracted."
    echo "This should not happen - please specify a single variant."
    rm -rf "$TEMP_DIR"
    exit 1
fi

echo "Successfully extracted variant"

# Process the single variant
echo "Processing variant data..."
awk 'BEGIN {FS=" "; OFS="\t"; print "IID\tGENOTYPE"} NR > 1 {print $2, $7}' "${TEMP_PREFIX}.raw" > "${OUTPUT_FILE}"

# Count genotypes for reporting
COUNTS=$(awk 'NR>1 {count[$2]++} END {for (val in count) printf("%s:%s ", val, count[val])}' "${OUTPUT_FILE}")
TOTAL=$(awk 'NR>1 {count++} END {print count}' "${OUTPUT_FILE}")

echo "Final genotype counts: $COUNTS"
echo "Extracted genotypes for $TOTAL samples"
echo "Coding scheme: PLINK default (0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate)"

# Clean up temporary files
rm -rf "$TEMP_DIR"

echo "Variant extraction complete! Results saved to $OUTPUT_FILE"
