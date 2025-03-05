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
        
        # We need to be specific with alleles
        VARIANT_SPEC="--chr ${CHROM#chr} --from-bp $POS --to-bp $POS --snps-only"
        
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
            
            if [ "$MATCH_COUNT" -gt 1 ]; then
                echo "Found multi-allelic site with $MATCH_COUNT variants at position $POS:"
                echo "$MATCHES"
                
                # Create a temp file to extract all variants at this position
                TEMP_EXTRACT_FILE=$(mktemp)
                echo "$MATCHES" | awk '{print $2}' > "$TEMP_EXTRACT_FILE"
                
                echo "Combining all minor alleles into one variant extraction..."
                # Extract all variants and combine them
                VARIANT_SPEC="--extract $TEMP_EXTRACT_FILE"
                
                # Store flag to combine multi-allelic variants
                MULTI_ALLELIC=true
            else
                # Single variant at this position
                VARIANT_SPEC="--chr ${CHROM#chr} --from-bp $POS --to-bp $POS"
                MULTI_ALLELIC=false
            fi
        else
            # Cannot check, proceed with position only
            echo "Warning: Could not check for multi-allelic variants - proceeding with position only"
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

# Extract the variant using PLINK2 with recode ad (additive coding)
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

echo "Processing variant data using raw output..."

# If multi-allelic, first report on individual alleles
if [ "$MULTI_ALLELIC" = "true" ]; then
    echo "Individual allele statistics before combining:"
    echo "---------------------------------------------"
    
    # Use a loop to extract columns and print stats for each variant
    COLUMN_COUNT=$(head -1 "${TEMP_PREFIX}.raw" | wc -w)
    VARIANT_COLUMNS=$(expr $COLUMN_COUNT - 6) # Subtract 6 header columns
    
    # Get the header line to extract variant names
    HEADER=$(head -1 "${TEMP_PREFIX}.raw")
    
    # Debug: show column headers
    if [ "$DEBUG_MODE" = "true" ]; then
        echo "DEBUG: Header has $COLUMN_COUNT columns, $VARIANT_COLUMNS are variants"
        echo "DEBUG: Header: $HEADER"
    fi
    
    # Process each variant column
    for i in $(seq 1 $VARIANT_COLUMNS); do
        COL_INDEX=$(expr $i + 6) # Add 6 to get actual column index
        
        # Extract variant name from header
        VARIANT_NAME=$(echo "$HEADER" | awk -v col=$COL_INDEX '{print $col}')
        
        # Get counts for this variant column
        COUNTS=$(awk -v col=$COL_INDEX 'NR>1 {count[$col]++} END {for (val in count) printf("%s:%s ", val, count[val])}' "${TEMP_PREFIX}.raw")
        TOTAL=$(awk 'NR>1 {count++} END {print count}' "${TEMP_PREFIX}.raw")
        
        echo "Variant $i ($VARIANT_NAME): $COUNTS (Total: $TOTAL samples)"
        
        # Debug: Dump first few values from this column
        if [ "$DEBUG_MODE" = "true" ]; then
            echo "DEBUG: First 5 values of column $COL_INDEX:"
            awk -v col=$COL_INDEX 'NR>1 && NR<=6 {print $col}' "${TEMP_PREFIX}.raw"
        fi
    done
    
    echo "---------------------------------------------"
    echo "Combining all variants using maximum dosage approach..."
fi

# Process the .raw file with awk to extract IID and genotype
# For multi-allelic sites, we need to combine all variant columns
# Process the .raw file with awk to extract IID and genotype
# For multi-allelic sites, we need to combine all variant columns
AWK_SCRIPT='
BEGIN { 
    FS=" "; # PLINK raw files use space as separator, not tab
    print "IID\tGENOTYPE" 
}
NR > 1 {
    # If we have 7 or more columns, extract IID and genotype data
    if (NF >= 7) {
        # For multi-allelic sites, we need to check all genotype columns (7 and beyond)
        # and use the maximum value across all alternate alleles
        # Important: convert to numeric value for proper comparison
        max_value = 0 + $7;  # Force numeric context
        
        # Debug line
        #printf("Sample: %s, First value: %s, ", $2, max_value);
        
        for (i=8; i<=NF; i++) {
            curr_value = 0 + $i;  # Force numeric context
            if (curr_value > max_value) {
                max_value = curr_value;
                #printf("New max at col %d: %s, ", i, max_value);
            }
        }
        #printf("\n");
        
        printf("%s\t%d\n", $2, max_value);
    }
}
'

# If multi-allelic, use special processing to combine all alternate alleles
if [ "$MULTI_ALLELIC" = "true" ]; then
    # Create a temporary AWK script file for better control
    TMP_AWK_SCRIPT=$(mktemp)
    
    cat > "$TMP_AWK_SCRIPT" << 'EOF'
BEGIN { 
    FS=" "; # PLINK raw files use space as separator, not tab
    print "IID\tGENOTYPE"; 
}
NR > 1 {
    # If we have 7 or more columns, extract IID and genotype data
    if (NF >= 7) {
        # For multi-allelic sites, we need to check all genotype columns (7 and beyond)
        # and use the maximum value across all alternate alleles
        # Important: convert to numeric value for proper comparison
        max_value = 0;
        for (i=7; i<=NF; i++) {
            curr_value = $i + 0;  # Force numeric context
            if (curr_value > max_value) {
                max_value = curr_value;
            }
        }
        
        printf("%s\t%d\n", $2, max_value);
    }
}
EOF
    
    # Use the awk script file for processing
    awk -f "$TMP_AWK_SCRIPT" "${TEMP_PREFIX}.raw" > "${OUTPUT_FILE}"
    
    # Clean up temporary files
    rm -f "$TMP_AWK_SCRIPT" "$TEMP_EXTRACT_FILE"
    
    # Debug output - dump the first 10 lines of the output file
    if [ "$DEBUG_MODE" = "true" ]; then
        echo "DEBUG: First 10 lines of final output file:"
        head -10 "${OUTPUT_FILE}"
    fi
else
    # For single allelic sites, just extract the first variant column
    awk 'BEGIN {FS=" "; OFS="\t"; print "IID\tGENOTYPE"} 
         NR > 1 {print $2, $7}' "${TEMP_PREFIX}.raw" > "${OUTPUT_FILE}"
fi

# Count genotypes for reporting
COUNTS=$(awk 'NR>1 {count[$2]++} END {for (val in count) printf("%s:%s ", val, count[val])}' "${OUTPUT_FILE}")
TOTAL=$(awk 'NR>1 {count++} END {print count}' "${OUTPUT_FILE}")

echo "Final combined genotype counts: $COUNTS"
echo "Extracted genotypes for $TOTAL samples"

# Additional statistics for combined data
ZEROS=$(awk 'NR>1 && $2==0 {count++} END {print count+0}' "${OUTPUT_FILE}")
ONES=$(awk 'NR>1 && $2==1 {count++} END {print count+0}' "${OUTPUT_FILE}")
TWOS=$(awk 'NR>1 && $2==2 {count++} END {print count+0}' "${OUTPUT_FILE}")

# Calculate frequencies
if [ $TOTAL -gt 0 ]; then
    ZERO_FREQ=$(awk -v z=$ZEROS -v t=$TOTAL 'BEGIN {printf "%.2f%%", (z/t)*100}')
    ONE_FREQ=$(awk -v o=$ONES -v t=$TOTAL 'BEGIN {printf "%.2f%%", (o/t)*100}')
    TWO_FREQ=$(awk -v w=$TWOS -v t=$TOTAL 'BEGIN {printf "%.2f%%", (w/t)*100}')
    
    echo "Genotype distribution: 0=$ZEROS ($ZERO_FREQ), 1=$ONES ($ONE_FREQ), 2=$TWOS ($TWO_FREQ)"
    
    # Calculate allele frequency (for combined variant)
    MAF=$(awk -v o=$ONES -v w=$TWOS -v t=$TOTAL 'BEGIN {printf "%.3f", (o+2*w)/(2*t)}')
    echo "Minor allele frequency (combined): $MAF"
fi

# Note about coding scheme
if [ "$MULTI_ALLELIC" = "true" ]; then
    echo "Note: Values have been properly recoded to follow minor allele coding:"
    echo "Coding scheme: 0 = homozygous major allele, 1 = heterozygous, 2 = homozygous minor allele"
elif [ "$CODE_MAJOR_AS_ZERO" = "true" ]; then
    echo "Warning: Without Python processing, we cannot guarantee major allele is coded as 0"
    echo "Coding scheme: PLINK default (0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate)"
else
    echo "Coding scheme: PLINK default (0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate)"
fi

# Clean up temporary files
rm -rf "$TEMP_DIR"

echo "Variant extraction complete! Results saved to $OUTPUT_FILE"
