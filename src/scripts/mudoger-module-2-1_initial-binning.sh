#!/bin/bash

##########  INITIAL PROKARYOTIC BINNING WITHOUT CONCOCT  ###################
echo '------- START MODULE 2-1 INITIAL BINNING (NO CONCOCT)'

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
source "$config_path"

conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env"

# Arguments declaration
assembly="$1"                     # assembly fasta file
forward_library="$2"              # forward library path
reverse_library="$3"              # reverse library path
output_folder="$4"                # output folder to be created inside master output folder
num_cores="$5"                    # number of threads
memory="$6"                       # memory in GB

log_info "Starting initial binning without CONCOCT"
log_info "Assembly: $assembly"
log_info "Forward reads: $forward_library"
log_info "Reverse reads: $reverse_library"
log_info "Output folder: $output_folder"
log_info "Cores: $num_cores"
log_info "Memory: ${memory}G"

# Validate input files
if [[ ! -f "$assembly" ]]; then
    log_error "Assembly file not found: $assembly"
    exit 1
fi

if [[ ! -f "$forward_library" ]]; then
    log_error "Forward reads file not found: $forward_library"
    exit 1
fi

if [[ ! -f "$reverse_library" ]]; then
    log_error "Reverse reads file not found: $reverse_library"
    exit 1
fi

# Create output directory
mkdir -p "$output_folder"

# Run metaWRAP binning with only MetaBAT2 and MaxBin2
log_info "Running metaWRAP binning with MetaBAT2 and MaxBin2 (CONCOCT excluded)"

if metawrap binning -o "$output_folder" -t "$num_cores" -a "$assembly" \
    --metabat2 --maxbin2 "$forward_library" "$reverse_library"; then

    log_success "Initial binning completed successfully"

    # Verify outputs
    local metabat2_bins="$output_folder/metabat2_bins"
    local maxbin2_bins="$output_folder/maxbin2_bins"

    if [[ -d "$metabat2_bins" ]]; then
        local metabat2_count=$(find "$metabat2_bins" -name "*.fa" 2>/dev/null | wc -l)
        log_success "MetaBAT2 produced $metabat2_count bins"
    else
        log_error "MetaBAT2 output directory not found"
        exit 1
    fi

    if [[ -d "$maxbin2_bins" ]]; then
        local maxbin2_count=$(find "$maxbin2_bins" -name "*.fa" 2>/dev/null | wc -l)
        log_success "MaxBin2 produced $maxbin2_count bins"
    else
        log_error "MaxBin2 output directory not found"
        exit 1
    fi

    local total_bins=$((metabat2_count + maxbin2_count))
    log_success "Total bins produced: $total_bins"

    if [[ $total_bins -eq 0 ]]; then
        log_error "No bins were produced by any binner"
        exit 1
    fi

else
    log_error "metaWRAP binning failed"
    exit 1
fi

conda deactivate

log_success "Initial binning module completed successfully"