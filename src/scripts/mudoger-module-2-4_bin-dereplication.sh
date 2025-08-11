#!/bin/bash

# 4 BIN REDUNDANCY REMOVAL - CORRECTED VERSION
echo '------- START MODULE 2-4 BIN REDUNDANCY REMOVAL'

# Set error handling (allow unset variables for compatibility)
set -eo pipefail

# Get arguments
lib_folder=$1
bin_count=0

# Extract sample name
lib_name=$(basename "$lib_folder")

echo "Processing sample: $lib_name"
echo "Working directory: $lib_folder"

# Create unique_bins directory (FIXED: missing in original)
mkdir -p "$lib_folder"/prokaryotes/binning/unique_bins

echo "Created unique bins directory: $lib_folder/prokaryotes/binning/unique_bins"

# Verify refinement directories exist
if [[ ! -d "$lib_folder/prokaryotes/binning/refinement-bac" ]]; then
    echo "ERROR: Bacterial refinement directory not found: $lib_folder/prokaryotes/binning/refinement-bac"
    exit 1
fi

if [[ ! -d "$lib_folder/prokaryotes/binning/refinement-arc" ]]; then
    echo "ERROR: Archaeal refinement directory not found: $lib_folder/prokaryotes/binning/refinement-arc"
    exit 1
fi

echo "✅ Refinement directories found"

# Calculate MD5 checksums for bacterial refined bins
echo "Calculating MD5 checksums for bacterial refined bins..."
if ls "$lib_folder"/prokaryotes/binning/refinement-bac/metawrap*bins/*fa >/dev/null 2>&1; then
    md5sum "$lib_folder"/prokaryotes/binning/refinement-bac/metawrap*bins/*fa > "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum
    bac_count=$(ls "$lib_folder"/prokaryotes/binning/refinement-bac/metawrap*bins/*fa | wc -l)
    echo "Found $bac_count bacterial refined bins"
else
    echo "No bacterial refined bins found"
    touch "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum
    bac_count=0
fi

# Calculate MD5 checksums for archaeal refined bins
echo "Calculating MD5 checksums for archaeal refined bins..."
if ls "$lib_folder"/prokaryotes/binning/refinement-arc/metawrap*bins/*fa >/dev/null 2>&1; then
    md5sum "$lib_folder"/prokaryotes/binning/refinement-arc/metawrap*bins/*fa >> "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum
    arc_count=$(ls "$lib_folder"/prokaryotes/binning/refinement-arc/metawrap*bins/*fa | wc -l)
    echo "Found $arc_count archaeal refined bins"
else
    echo "No archaeal refined bins found"
    arc_count=0
fi

total_refined_bins=$((bac_count + arc_count))
echo "Total refined bins: $total_refined_bins"

# Check if any bins were found
if [[ $total_refined_bins -eq 0 ]]; then
    echo "ERROR: No refined bins found for dereplication"
    exit 1
fi

# Verify md5_sum file was created and has content
if [[ ! -f "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum ]]; then
    echo "ERROR: MD5 sum file was not created"
    exit 1
fi

if [[ ! -s "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum ]]; then
    echo "ERROR: MD5 sum file is empty"
    exit 1
fi

echo "MD5 checksums calculated for $(wc -l < "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum) bins"

# Extract unique MD5 hashes
cat "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum | cut -f1 -d' ' | sort | uniq > "$lib_folder"/prokaryotes/binning/unique_bins/md5_unique

unique_hash_count=$(wc -l < "$lib_folder"/prokaryotes/binning/unique_bins/md5_unique)
duplicate_count=$((total_refined_bins - unique_hash_count))

echo "Found $unique_hash_count unique bins (removed $duplicate_count duplicates)"

# Copy unique bins with standardized names
echo "Copying unique bins to final location..."

while read l; do
    # Extract path and trim leading/trailing whitespace
    bininit="$(grep -F "$l" "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum | head -n1 | cut -f2- -d' ' | xargs)"
    binafter="$lib_folder"/prokaryotes/binning/unique_bins/"$lib_name"-bin."$bin_count".fa

    # Verify source file exists before copying
    if [[ -f "$bininit" ]]; then
        if cp "$bininit" "$binafter"; then
            echo " Copied bin $bin_count: $(basename "$bininit") -> $(basename "$binafter")"
            bin_count=$((bin_count + 1))
        else
            echo " ERROR: Failed to copy $bininit to $binafter"
        fi
    else
        echo " WARNING: Source bin file not found: $bininit"
    fi
done < "$lib_folder"/prokaryotes/binning/unique_bins/md5_unique

# Clean up temporary files
rm -f "$lib_folder"/prokaryotes/binning/unique_bins/md5_unique
rm -f "$lib_folder"/prokaryotes/binning/unique_bins/md5_sum

# Verify final results
final_bin_count=$(find "$lib_folder"/prokaryotes/binning/unique_bins -name "*.fa" 2>/dev/null | wc -l)

echo ""
echo "=== BIN DEREPLICATION SUMMARY ==="
echo "Original refined bins: $total_refined_bins"
echo "Unique bins identified: $unique_hash_count"
echo "Duplicates removed: $duplicate_count"
echo "Final unique bins created: $final_bin_count"

if [[ $final_bin_count -gt 0 ]]; then
    echo "✅ Bin dereplication completed successfully!"
    echo "Unique bins saved in: $lib_folder/prokaryotes/binning/unique_bins/"
else
    echo "❌ ERROR: No final bins were created during dereplication"
    exit 1
fi

echo ""
echo "Bin redundancy removal completed successfully"