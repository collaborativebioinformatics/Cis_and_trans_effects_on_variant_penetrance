#!/bin/bash

# Script to extract a single variant from PLINK files and output as 0,1,2 coding
# When only position is specified, 0 is assigned to the major allele
# Requires: PLINK2, Python with pandas

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
    for bim in $PLINK_DIR/*.bim; do
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
        echo "Warning: Variant $VARIANT is multi-allelic with $ALLELE_COUNT alleles"
        echo "Showing all alleles:"
        grep "$VARIANT" "$FOUND_BIM"
        
        # Get first instance position, ref and alt alleles
        VARIANT_INFO=$(grep "$VARIANT" "$FOUND_BIM" | head -1)
        POS=$(echo "$VARIANT_INFO" | awk '{print $4}')
        REF=$(echo "$VARIANT_INFO" | awk '{print $5}')
        ALT=$(echo "$VARIANT_INFO" | awk '{print $6}')
        
        echo "Using first instance: POS=$POS, REF=$REF, ALT=$ALT"
        VARIANT_SPEC="--snp $VARIANT"
    else
        VARIANT_SPEC="--snp $VARIANT"
    fi
    
    # For rsIDs, we'll use PLINK's allele coding (ref=0)
    CODE_MAJOR_AS_ZERO=false
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
        
        # When specific alleles are provided, use PLINK's default coding
        CODE_MAJOR_AS_ZERO=false
    else
        # chr:pos format - potentially multi-allelic
        POS=$(echo $VARIANT | cut -d':' -f2)
        
        # Find the chromosome file to check for variants at this position
        PATTERNS=("$CHROM" "${CHROM#chr}" "chr${CHROM#chr}")
        FOUND_BIM=""
        
        for pattern in "${PATTERNS[@]}"; do
            if ls $PLINK_DIR/$pattern.bim 2>/dev/null; then
                FOUND_BIM="$PLINK_DIR/$pattern.bim"
                break
            fi
            
            TEMP=$(find $PLINK_DIR -name "*.bim" | grep -i "$pattern" | head -1)
            if [ -n "$TEMP" ]; then
                FOUND_BIM="$TEMP"
                break
            fi
        done
        
        if [ -n "$FOUND_BIM" ]; then
            # Check if position has multiple variants
            MATCHES=$(awk -v pos="$POS" '$4 == pos {print}' "$FOUND_BIM")
            MATCH_COUNT=$(echo "$MATCHES" | grep -v "^$" | wc -l)
            
            if [ "$MATCH_COUNT" -eq 0 ]; then
                echo "Error: No variants found at position $POS"
                exit 1
            elif [ "$MATCH_COUNT" -gt 1 ]; then
                echo "Warning: Position $POS has $MATCH_COUNT variants:"
                echo "$MATCHES"
                
                # Create a temporary file to extract all variants at this position
                TEMP_EXTRACT_FILE=$(mktemp)
                echo "$MATCHES" | awk '{print $2}' > "$TEMP_EXTRACT_FILE"
                
                echo "Extracting all variants at this position..."
                VARIANT_SPEC="--extract $TEMP_EXTRACT_FILE"
                MULTI_ALLELIC=true
            else
                # Single variant at this position
                VARIANT_SPEC="--chr ${CHROM#chr} --from-bp $POS --to-bp $POS"
                MULTI_ALLELIC=false
            fi
        else
            # Cannot check, proceed with position only but warn user
            echo "Warning: Could not find chromosome file to check for variants at position $POS"
            echo "Proceeding with position only - results may be unexpected if the position doesn't exist"
            VARIANT_SPEC="--chr ${CHROM#chr} --from-bp $POS --to-bp $POS"
            MULTI_ALLELIC=false
        fi
        
        # For position-only queries, we want major allele as 0
        CODE_MAJOR_AS_ZERO=true
    fi
    
    echo "Using chromosome $CHROM and position $POS"
fi

# Find the appropriate chromosome file in the directory - IMPROVED MATCHING
PLINK_PREFIX=""
PATTERNS=("$CHROM" "${CHROM#chr}" "chr${CHROM#chr}")

# Get all .bed files in directory
BED_FILES=$(find "$PLINK_DIR" -name "*.bed")

# First, look for exact file matches
for pattern in "${PATTERNS[@]}"; do
    if [ -f "$PLINK_DIR/$pattern.bed" ]; then
        PLINK_PREFIX="$PLINK_DIR/$pattern"
        echo "Found exact PLINK files with prefix: $PLINK_PREFIX"
        break
    fi
done

# If no exact match, search for files that end with the target chromosome
if [ -z "$PLINK_PREFIX" ]; then
    echo "Looking for files containing chromosome pattern..."
    
    # For each pattern, try to find files that end with the pattern
    for pattern in "${PATTERNS[@]}"; do
        # Look specifically for files that end with the chromosome pattern
        # but avoid matching chr10 when looking for chr1
        for bed_file in $BED_FILES; do
            # Get basename without extension
            base_name=$(basename "$bed_file" .bed)
            
            # Check if filename ends with exact chromosome (e.g., ends with "chr1" not within "chr10")
            if [[ "$base_name" == *"_${pattern}" || "$base_name" == *"_${pattern#chr}" ]]; then
                PLINK_PREFIX="${bed_file%.bed}"
                echo "Found prefix-ending match: $PLINK_PREFIX"
                break 2
            elif [[ "$base_name" == *"${pattern}" && "$base_name" != *"${pattern}0"* && "$base_name" != *"${pattern}1"* && "$base_name" != *"${pattern}2"* ]]; then
                # Match "something_chr1" or "something_chr1_something" but not "something_chr10"
                PLINK_PREFIX="${bed_file%.bed}"
                echo "Found chromosome-specific match: $PLINK_PREFIX"
                break 2
            fi
        done
    done
fi

# As a last resort, look for any file containing the pattern
if [ -z "$PLINK_PREFIX" ]; then
    for bed_file in $BED_FILES; do
        base_name=$(basename "$bed_file")
        # For each chromosome pattern
        for pattern in "${PATTERNS[@]}"; do
            # Find files containing the exact chromosome string
            if [[ "$base_name" == *"${pattern}"* ]]; then
                # Exclude matches to chr10, chr11 etc when looking for chr1
                # by checking that the next character isn't a digit
                next_char="${base_name:${#pattern}:1}"
                if [[ ! "$next_char" =~ [0-9] ]]; then
                    PLINK_PREFIX="${bed_file%.bed}"
                    echo "Found pattern match: $PLINK_PREFIX (matched $pattern)"
                    break 2
                fi
            fi
        done
    done
fi

if [ -z "$PLINK_PREFIX" ]; then
    echo "Error: Could not find PLINK files for chromosome $CHROM in directory $PLINK_DIR"
    echo "Available .bed files in directory:"
    find "$PLINK_DIR" -name "*.bed" | sort
    exit 1
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

echo "Successfully extracted $VARIANT_COUNT variant(s)"

# Handle multi-allelic sites by combining or selecting a single variant
if [ "$VARIANT_COUNT" -gt 1 ] || [ "$MULTI_ALLELIC" = "true" ]; then
    echo "Processing multi-allelic site with $VARIANT_COUNT variants"
    
    # Print individual allele statistics
    echo "Individual allele statistics:"
    echo "---------------------------------------------"
    HEADER=$(head -1 "${TEMP_PREFIX}.raw")
    
    # For each variant column, print statistics and check if flipping is needed
    FLIP_COLUMNS=""
    
    for i in $(seq 1 $VARIANT_COUNT); do
        COL_INDEX=$(( i + 6 ))  # Add 6 to get actual column index
        
        # Extract variant name from header
        VARIANT_NAME=$(echo "$HEADER" | awk -v col=$COL_INDEX '{print $col}')
        
        # Get counts for this variant column (0, 1, 2)
        COUNTS_0=$(awk -v col=$COL_INDEX 'NR>1 && $col==0 {count++} END {print count+0}' "${TEMP_PREFIX}.raw")
        COUNTS_1=$(awk -v col=$COL_INDEX 'NR>1 && $col==1 {count++} END {print count+0}' "${TEMP_PREFIX}.raw")
        COUNTS_2=$(awk -v col=$COL_INDEX 'NR>1 && $col==2 {count++} END {print count+0}' "${TEMP_PREFIX}.raw")
        TOTAL=$(awk 'NR>1 {count++} END {print count}' "${TEMP_PREFIX}.raw")
        
        # Determine if we need to flip this column (if 2 is more common than 0)
        if [ $COUNTS_2 -gt $COUNTS_0 ]; then
            FLIP_COLUMNS="$FLIP_COLUMNS $COL_INDEX"
            echo "Variant $i ($VARIANT_NAME): 0:$COUNTS_0 1:$COUNTS_1 2:$COUNTS_2 (Total: $TOTAL) - FLIPPED CODING"
        else
            echo "Variant $i ($VARIANT_NAME): 0:$COUNTS_0 1:$COUNTS_1 2:$COUNTS_2 (Total: $TOTAL)"
        fi
    done
    echo "---------------------------------------------"
    
    # Use a Python script for more reliable handling
    PYTHON_SCRIPT=$(mktemp)
    cat > "$PYTHON_SCRIPT" << 'EOF'
#!/usr/bin/env python3
import sys
import os
from collections import Counter

# Get input arguments
raw_file = sys.argv[1]
output_file = sys.argv[2]
flip_columns = sys.argv[3].split()
flip_columns = [int(col) for col in flip_columns if col]

# Read the raw file
with open(raw_file, 'r') as f:
    lines = [line.strip() for line in f]

# Parse header
header = lines[0].split()
variant_start = 6  # Index where variant columns start (0-based)
variant_count = len(header) - variant_start

print(f"Processing {variant_count} variants with Python")
if flip_columns:
    print(f"Columns to flip: {flip_columns}")

# Function to flip a value (0->2, 2->0, 1->1)
def flip_value(val):
    if val == 0:
        return 2
    elif val == 2:
        return 0
    return val

# Combined genotype counts
combined_counts = Counter()

# Process data rows
results = []
results.append("IID\tGENOTYPE\n")  # Header line

# Process each data row
for line_idx in range(1, len(lines)):
    parts = lines[line_idx].split()
    
    if len(parts) <= variant_start:
        continue  # Skip rows with insufficient columns
    
    iid = parts[1]  # IID is column 1
    
    # Process all variant columns with flipping as needed
    values = []
    for i in range(variant_count):
        col_idx = i + variant_start
        if col_idx < len(parts):
            try:
                val = int(float(parts[col_idx]))
                # Apply flipping if needed
                if col_idx in flip_columns:
                    val = flip_value(val)
                values.append(val)
            except (ValueError, IndexError):
                pass
    
    # Find maximum value after flipping
    max_val = max(values) if values else 0
    
    # Add to results
    results.append(f"{iid}\t{max_val}\n")
    combined_counts[max_val] += 1

# Write the output file
with open(output_file, 'w') as f:
    for line in results:
        f.write(line)

# Print final statistics
print("\nFinal combined genotype counts after flipping:")
total = sum(combined_counts.values())
for val in sorted(combined_counts.keys()):
    percent = (combined_counts[val] / total) * 100 if total > 0 else 0
    print(f"  {val}: {combined_counts[val]} ({percent:.2f}%)")

# Calculate MAF
if total > 0:
    minor_allele_count = combined_counts[1] + (2 * combined_counts[2])
    total_alleles = 2 * total
    maf = minor_allele_count / total_alleles if total_alleles > 0 else 0
    print(f"Minor allele frequency (combined): {maf:.3f}")

print(f"Processed {total} samples")
EOF

    # Run the Python script for combining variants with proper flipping
    echo "Combining variants with Python script..."
    python3 "$PYTHON_SCRIPT" "${TEMP_PREFIX}.raw" "${OUTPUT_FILE}" "$FLIP_COLUMNS"
    
    # If Python failed, fall back to a simpler approach
    if [ $? -ne 0 ]; then
        echo "Warning: Python processing failed, falling back to simple AWK processing"
        awk 'BEGIN {FS=" "; OFS="\t"; print "IID\tGENOTYPE"} NR > 1 {print $2, $7}' "${TEMP_PREFIX}.raw" > "${OUTPUT_FILE}"
    fi
    
    # Clean up
    rm -f "$PYTHON_SCRIPT" "$TEMP_EXTRACT_FILE"
else
    # Single variant - just extract the genotype column
    echo "Processing single variant data..."
    awk 'BEGIN {FS=" "; OFS="\t"; print "IID\tGENOTYPE"} NR > 1 {print $2, $7}' "${TEMP_PREFIX}.raw" > "${OUTPUT_FILE}"
    
    # Count genotypes for reporting
    COUNTS=$(awk 'NR>1 {count[$2]++} END {for (val in count) printf("%s:%s ", val, count[val])}' "${OUTPUT_FILE}")
    TOTAL=$(awk 'NR>1 {count++} END {print count}' "${OUTPUT_FILE}")
    
    echo "Final genotype counts: $COUNTS"
    echo "Extracted genotypes for $TOTAL samples"
    
    # Note about coding scheme
    if [ "$CODE_MAJOR_AS_ZERO" = "true" ]; then
        echo "Warning: Cannot guarantee major allele is coded as 0 in this mode"
        echo "Coding scheme: PLINK default (0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate)"
    else
        echo "Coding scheme: PLINK default (0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate)"
    fi
fi

# Clean up temporary files
rm -rf "$TEMP_DIR"

echo "Variant extraction complete! Results saved to $OUTPUT_FILE"
