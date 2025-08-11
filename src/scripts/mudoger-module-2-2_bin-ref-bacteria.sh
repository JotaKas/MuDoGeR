#!/bin/bash

########## PROKARYOTIC REFINEMENT FOR BACTERIA (50,10) WITHOUT CONCOCT ###################
echo '------- START MODULE 2-2 BIN REFINEMENT FOR BACTERIA (NO CONCOCT)'

# Set strict error handling
set -eo pipefail

# Logging functions
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}

log_success() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SUCCESS] ✅ $*"
}

log_warning() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARNING] ⚠️  $*"
}

# Load conda environment
conda activate mudoger_env
config_path="$(which config.sh)"
database="${config_path/config/database}"
source "$config_path"
conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env"
source "$config_path"
source "$database"

# Arguments declaration
output_folder="$1"              # output folder to be created inside master output folder
cores="$2"                      # number of cores
assembly="$3"                   # assembly file
maxbin2_bins="$4"              # maxbin2 bins directory
metabat2_bins="$5"             # metabat2 bins directory
memory="$6"                     # memory allocation

log_info "Starting bacterial bin refinement without CONCOCT"
log_info "Output folder: $output_folder"
log_info "Cores: $cores"
log_info "Assembly: $assembly"
log_info "MaxBin2 bins: $maxbin2_bins"
log_info "MetaBAT2 bins: $metabat2_bins"
log_info "Memory: ${memory}G"

# Validate input directories
if [[ ! -d "$maxbin2_bins" ]]; then
    log_error "MaxBin2 bins directory not found: $maxbin2_bins"
    exit 1
fi

if [[ ! -d "$metabat2_bins" ]]; then
    log_error "MetaBAT2 bins directory not found: $metabat2_bins"
    exit 1
fi

# Check if bins exist
maxbin2_count=$(find "$maxbin2_bins" -name "*.fa" 2>/dev/null | wc -l)
metabat2_count=$(find "$metabat2_bins" -name "*.fa" 2>/dev/null | wc -l)

log_info "MaxBin2 bins available: $maxbin2_count"
log_info "MetaBAT2 bins available: $metabat2_count"

if [[ $((maxbin2_count + metabat2_count)) -eq 0 ]]; then
    log_error "No bins found in input directories"
    exit 1
fi

# Create output directory
mkdir -p "$output_folder"

# Run only once during database installation configuration
CHECKM_DB="$DATABASES_LOCATION/checkm"
echo "${CHECKM_DB}" | checkm data setRoot "${CHECKM_DB}"

log_info "Running metaWRAP bin refinement for bacteria (completeness ≥50%, contamination ≤10%)"

# Run bin refinement with only two bin sets (no CONCOCT)
if metawrap bin_refinement -o "$output_folder" -t "$cores" \
    -A "$metabat2_bins" -B "$maxbin2_bins" \
    -c 50 -x 10 -m "$memory"; then

    log_success "Bacterial bin refinement completed successfully"

    # Verify output
    stats_file="$output_folder/metawrap_50_10_bins.stats"
    refined_bins_dir="$output_folder/metawrap_50_10_bins"

    if [[ -f "$stats_file" ]]; then
        log_success "Stats file created: $stats_file"

        # Count refined bins
        if [[ -d "$refined_bins_dir" ]]; then
            refined_count=$(find "$refined_bins_dir" -name "*.fa" 2>/dev/null | wc -l)
            log_success "Bacterial refinement produced $refined_count refined bins"

            if [[ $refined_count -eq 0 ]]; then
                log_warning "No refined bins met the quality criteria (≥50% completeness, ≤10% contamination)"
            fi
        else
            log_warning "Refined bins directory not found: $refined_bins_dir"
        fi
    else
        log_error "Stats file not created, refinement may have failed"
        exit 1
    fi

else
    log_error "metaWRAP bin refinement failed"
    exit 1
fi

conda deactivate

log_success "Bacterial bin refinement module completed successfully"