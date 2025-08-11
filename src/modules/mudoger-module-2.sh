#!/bin/bash
# PROKARYOTES MODULES FILESYSTEM STRUCTURE
#   BINNING
#     INITIAL-BINNING
#     REFINEMENT-BACTERIA
#     REFINEMENT-ARCHAEA
#     BIN-FINAL-SET (REDUNDANCY REMOVAL)
#   METRICS
#     CHECKM
#     GTDBtk
#     PROKKA
#     STATS (GENOMIC)
#       N50, NUM_NUCLEOTIDE, NUM_CONTIGS, ATCG and more...

# Set strict error handling
set -eo pipefail

# Global variables for logging
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_FILE=""
PROGRESS_STEPS=8
CURRENT_STEP=0

# Logging functions
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*" | tee -a "${LOG_FILE:-/dev/null}"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" | tee -a "${LOG_FILE:-/dev/null}" >&2
}

log_success() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SUCCESS] ✅ $*" | tee -a "${LOG_FILE:-/dev/null}"
}

log_warning() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARNING] ⚠️  $*" | tee -a "${LOG_FILE:-/dev/null}"
}

progress_step() {
    CURRENT_STEP=$((CURRENT_STEP + 1))
    local progress=$((CURRENT_STEP * 100 / PROGRESS_STEPS))
    log_info "[$CURRENT_STEP/$PROGRESS_STEPS] ($progress%) $*"
}

# Verification functions
verify_file_exists() {
    local file_path="$1"
    local description="$2"

    if [[ -f "$file_path" ]]; then
        log_success "$description exists: $file_path"
        return 0
    else
        log_error "$description not found: $file_path"
        return 1
    fi
}

verify_directory_exists() {
    local dir_path="$1"
    local description="$2"

    if [[ -d "$dir_path" ]]; then
        log_success "$description exists: $dir_path"
        return 0
    else
        log_error "$description not found: $dir_path"
        return 1
    fi
}

verify_directory_not_empty() {
    local dir_path="$1"
    local description="$2"

    if [[ -d "$dir_path" && "$(find "$dir_path" -maxdepth 1 -name "*.fa" -o -name "*.fasta" | wc -l)" -gt 0 ]]; then
        local count=$(find "$dir_path" -maxdepth 1 -name "*.fa" -o -name "*.fasta" | wc -l)
        log_success "$description contains $count bins: $dir_path"
        return 0
    else
        log_error "$description is empty or doesn't exist: $dir_path"
        return 1
    fi
}

verify_binning_output() {
    local output_dir="$1"
    local required_bins=("metabat2_bins" "maxbin2_bins")

    log_info "Verifying binning outputs in $output_dir"

    local success=true
    for bin_dir in "${required_bins[@]}"; do
        if ! verify_directory_not_empty "$output_dir/$bin_dir" "$bin_dir directory"; then
            success=false
        fi
    done

    # Check if at least one binner produced results
    local total_bins=0
    for bin_dir in "${required_bins[@]}"; do
        if [[ -d "$output_dir/$bin_dir" ]]; then
            local count=$(find "$output_dir/$bin_dir" -maxdepth 1 -name "*.fa" -o -name "*.fasta" 2>/dev/null | wc -l)
            total_bins=$((total_bins + count))
        fi
    done

    if [[ $total_bins -gt 0 ]]; then
        log_success "Total bins found: $total_bins"
        return 0
    else
        log_error "No bins found in any binner output"
        return 1
    fi
}

verify_refinement_output() {
    local refinement_dir="$1"
    local stats_file="$2"
    local description="$3"

    log_info "Verifying $description output"

    if verify_file_exists "$stats_file" "$description stats file" &&
       verify_directory_exists "$refinement_dir" "$description directory"; then

        # Check if refined bins directory has content
        local refined_bins_dir=$(find "$refinement_dir" -name "*bins" -type d | head -1)
        if [[ -d "$refined_bins_dir" ]]; then
            local bin_count=$(find "$refined_bins_dir" -name "*.fa" -o -name "*.fasta" 2>/dev/null | wc -l)
            log_success "$description produced $bin_count refined bins"
            return 0
        fi
    fi

    log_error "$description verification failed"
    return 1
}

help_message() {
    cat << EOF

Usage: bash -i MuDoGeR/src/modules/mudoger-module-2.sh -1 reads_1.fastq -2 reads_2.fastq -o output_dir -a assembly.fasta

Options:
    -1 STR          forward fastq reads path
    -2 STR          reverse fastq reads path
    -a STR          assembly path
    -o STR          output directory path
    -m INT          given RAM in GB (default=50G)
    -t INT          number of threads/cores (default=1)
    -h --help       print this message

EOF
}

# Step verification functions
verify_initial_binning() {
    local binning_dir="$1"
    verify_binning_output "$binning_dir"
}

verify_bacterial_refinement() {
    local refinement_dir="$1"
    local stats_file="$refinement_dir/metawrap_50_10_bins.stats"
    verify_refinement_output "$refinement_dir" "$stats_file" "Bacterial refinement"
}

verify_archaeal_refinement() {
    local refinement_dir="$1"
    local stats_file="$refinement_dir/metawrap_40_30_bins.stats"
    verify_refinement_output "$refinement_dir" "$stats_file" "Archaeal refinement"
}

verify_dereplication() {
    local unique_bins_dir="$1"
    verify_directory_not_empty "$unique_bins_dir" "Unique bins directory"
}

verify_taxonomy() {
    local taxonomy_file="$1"
    verify_file_exists "$taxonomy_file" "GTDBtk taxonomy results"
}

verify_checkm() {
    local checkm_file="$1"
    verify_file_exists "$checkm_file" "CheckM results"
}

verify_prokka() {
    local prokka_dir="$1"
    verify_directory_exists "$prokka_dir" "Prokka annotation directory"
}

verify_metrics() {
    local metrics_file="$1"
    verify_file_exists "$metrics_file" "Genome metrics"
}

# Step execution functions
execute_initial_binning() {
    local assembly="$1"
    local forward_library="$2"
    local reverse_library="$3"
    local binning_dir="$4"
    local cores="$5"
    local memory="$6"

    mkdir -p "$binning_dir"

    # Modified script call without CONCOCT
    bash -i "$MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-1_initial-binning.sh" \
        "$assembly" \
        "$forward_library" \
        "$reverse_library" \
        "$binning_dir" \
        "$cores" \
        "$memory"
}

execute_bacterial_refinement() {
    refinement_dir="$1"
    cores="$2"
    assembly="$3"
    binning_dir="$4"
    memory="$5"

    # Only use MetaBAT2 and MaxBin2 bins for refinement
    bash -i "$MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-2_bin-ref-bacteria.sh" \
        "$refinement_dir" \
        "$cores" \
        "$assembly" \
        "$binning_dir/maxbin2_bins" \
        "$binning_dir/metabat2_bins" \
        "$memory"
}

execute_archaeal_refinement() {
    refinement_dir="$1"
    cores="$2"
    assembly="$3"
    binning_dir="$4"
    memory="$5"

    # Only use MetaBAT2 and MaxBin2 bins for refinement
    bash -i "$MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-3_bin-ref-archea.sh" \
        "$refinement_dir" \
        "$cores" \
        "$assembly" \
        "$binning_dir/maxbin2_bins" \
        "$binning_dir/metabat2_bins" \
        "$memory"
}

execute_dereplication() {
    local output_folder="$1"

    bash -i "$MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-4_bin-dereplication.sh" "$output_folder"
}

execute_taxonomy() {
    local prokaryotes_dir="$1"
    local cores="$2"

    bash -i "$MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-5_bin-taxonomy.sh" "$prokaryotes_dir" "$cores"
}

execute_checkm() {
    local prokaryotes_dir="$1"
    local cores="$2"

    bash -i "$MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-6_bin-QC.sh" "$prokaryotes_dir" "$cores"
}

execute_prokka() {
    local prokaryotes_dir="$1"
    local cores="$2"

    bash -i "$MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-7_prokka.sh" "$prokaryotes_dir" "$cores"
}

execute_metrics() {
    local prokaryotes_dir="$1"
    local cores="$2"

    bash -i "$MUDOGER_CONDA_ENVIRONMENT_PATH/bin/mudoger-module-2-8_genomics-metrics.sh" "$prokaryotes_dir" "$cores"
}

# Main execution function
main() {
    # Default parameters
    local libname_folder=$(pwd)
    local memory=50
    local cores=1
    local forward_library=""
    local reverse_library=""
    local assembly=""

    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -1) forward_library="$2"; shift 2;;
            -2) reverse_library="$2"; shift 2;;
            -a) assembly="$2"; shift 2;;
            -o) libname_folder="$2"; shift 2;;
            -t) cores="$2"; shift 2;;
            -m) memory="$2"; shift 2;;
            -h|--help) help_message; exit 0;;
            *) log_error "Unknown option: $1"; help_message; exit 1;;
        esac
    done

    # Validate required parameters
    if [[ -z "$forward_library" || -z "$reverse_library" || -z "$assembly" ]]; then
        log_error "Missing required parameters"
        help_message
        exit 1
    fi

    # Set up logging
    mkdir -p "$libname_folder/logs"
    LOG_FILE="$libname_folder/logs/prokaryotes_module.log"

    log_info "=== MuDoGeR Prokaryotes Module Started ==="
    log_info "Run ID: $(date '+%Y%m%d_%H%M%S')"
    log_info "Parameters:"
    log_info "  Forward reads: $forward_library"
    log_info "  Reverse reads: $reverse_library"
    log_info "  Assembly: $assembly"
    log_info "  Output directory: $libname_folder"
    log_info "  Cores: $cores"
    log_info "  Memory: ${memory}G"
    log_info "  Log file: $LOG_FILE"

    # Activate conda environment
    conda activate mudoger_env
    config_path="$(which config.sh)"
    source "$config_path"

    # Define paths
    local binning_dir="$libname_folder/prokaryotes/binning/initial-binning"
    local refinement_bac_dir="$libname_folder/prokaryotes/binning/refinement-bac"
    local refinement_arc_dir="$libname_folder/prokaryotes/binning/refinement-arc"
    local unique_bins_dir="$libname_folder/prokaryotes/binning/unique_bins"
    local metrics_dir="$libname_folder/prokaryotes/metrics"

    # Ensure metrics directory exists
    mkdir -p "$metrics_dir"

    # Execute pipeline steps with verification
    local success=true

    # Step 1: Initial Binning (without CONCOCT)
    progress_step "Starting Initial binning (MetaBAT2 + MaxBin2)"

    if verify_initial_binning "$binning_dir"; then
        log_success "Initial binning (MetaBAT2 + MaxBin2) already completed successfully - skipping"
    else
        log_info "Running Initial binning (MetaBAT2 + MaxBin2)"

        if execute_initial_binning "$assembly" "$forward_library" "$reverse_library" "$binning_dir" "$cores" "$memory"; then
            if verify_initial_binning "$binning_dir"; then
                log_success "Initial binning (MetaBAT2 + MaxBin2) completed successfully"
            else
                log_error "Initial binning execution completed but verification failed"
                success=false
            fi
        else
            log_error "Initial binning execution failed"
            success=false
        fi
    fi

    # Step 2: Bacterial refinement
    if [[ "$success" == true ]]; then
        progress_step "Starting Bacterial refinement"

        if verify_bacterial_refinement "$refinement_bac_dir"; then
            log_success "Bacterial refinement already completed successfully - skipping"
        else
            log_info "Running Bacterial refinement"

            if execute_bacterial_refinement "$refinement_bac_dir" "$cores" "$assembly" "$binning_dir" "$memory"; then
                if verify_bacterial_refinement "$refinement_bac_dir"; then
                    log_success "Bacterial refinement completed successfully"
                else
                    log_error "Bacterial refinement execution completed but verification failed"
                    success=false
                fi
            else
                log_error "Bacterial refinement execution failed"
                success=false
            fi
        fi
    fi

    # Step 3: Archaeal refinement
    if [[ "$success" == true ]]; then
        progress_step "Starting Archaeal refinement"

        if verify_archaeal_refinement "$refinement_arc_dir"; then
            log_success "Archaeal refinement already completed successfully - skipping"
        else
            log_info "Running Archaeal refinement"

            if execute_archaeal_refinement "$refinement_arc_dir" "$cores" "$assembly" "$binning_dir" "$memory"; then
                if verify_archaeal_refinement "$refinement_arc_dir"; then
                    log_success "Archaeal refinement completed successfully"
                else
                    log_error "Archaeal refinement execution completed but verification failed"
                    success=false
                fi
            else
                log_error "Archaeal refinement execution failed"
                success=false
            fi
        fi
    fi

    # Step 4: Bin dereplication
    if [[ "$success" == true ]]; then
        progress_step "Starting Bin dereplication"

        if verify_dereplication "$unique_bins_dir"; then
            log_success "Bin dereplication already completed successfully - skipping"
        else
            log_info "Running Bin dereplication"

            if execute_dereplication "$libname_folder"; then
                if verify_dereplication "$unique_bins_dir"; then
                    log_success "Bin dereplication completed successfully"
                else
                    log_error "Bin dereplication execution completed but verification failed"
                    success=false
                fi
            else
                log_error "Bin dereplication execution failed"
                success=false
            fi
        fi
    fi

    # Step 5: GTDBtk taxonomy assignment
    if [[ "$success" == true ]]; then
        progress_step "Starting GTDBtk taxonomy assignment"
        local taxonomy_file="$metrics_dir/GTDBtk_taxonomy/gtdbtk_result.tsv"

        if verify_taxonomy "$taxonomy_file"; then
            log_success "GTDBtk taxonomy assignment already completed successfully - skipping"
        else
            log_info "Running GTDBtk taxonomy assignment"

            if execute_taxonomy "$libname_folder/prokaryotes" "$cores"; then
                if verify_taxonomy "$taxonomy_file"; then
                    log_success "GTDBtk taxonomy assignment completed successfully"
                else
                    log_warning "GTDBtk taxonomy assignment execution completed but verification failed"
                fi
            else
                log_warning "GTDBtk taxonomy assignment execution failed, continuing with other steps"
            fi
        fi
    fi

    # Step 6: CheckM quality control
    if [[ "$success" == true ]]; then
        progress_step "Starting CheckM quality control"
        local checkm_file="$metrics_dir/checkm_qc/outputcheckm.tsv"

        if verify_checkm "$checkm_file"; then
            log_success "CheckM quality control already completed successfully - skipping"
        else
            log_info "Running CheckM quality control"

            if execute_checkm "$libname_folder/prokaryotes" "$cores"; then
                if verify_checkm "$checkm_file"; then
                    log_success "CheckM quality control completed successfully"
                else
                    log_warning "CheckM quality control execution completed but verification failed"
                fi
            else
                log_warning "CheckM quality control execution failed, continuing with other steps"
            fi
        fi
    fi

    # Step 7: Prokka annotation
    if [[ "$success" == true ]]; then
        progress_step "Starting Prokka annotation"
        local prokka_dir="$metrics_dir/prokka"

        if verify_prokka "$prokka_dir"; then
            log_success "Prokka annotation already completed successfully - skipping"
        else
            log_info "Running Prokka annotation"

            if execute_prokka "$libname_folder/prokaryotes" "$cores"; then
                if verify_prokka "$prokka_dir"; then
                    log_success "Prokka annotation completed successfully"
                else
                    log_warning "Prokka annotation execution completed but verification failed"
                fi
            else
                log_warning "Prokka annotation execution failed, continuing with other steps"
            fi
        fi
    fi

    # Step 8: Genome metrics
    if [[ "$success" == true ]]; then
        progress_step "Starting Genome metrics calculation"
        local metrics_file="$metrics_dir/genome_statistics/prok_genomes_stats.tsv"

        if verify_metrics "$metrics_file"; then
            log_success "Genome metrics calculation already completed successfully - skipping"
        else
            log_info "Running Genome metrics calculation"

            if execute_metrics "$libname_folder/prokaryotes" "$cores"; then
                if verify_metrics "$metrics_file"; then
                    log_success "Genome metrics calculation completed successfully"
                else
                    log_warning "Genome metrics calculation execution completed but verification failed"
                fi
            else
                log_warning "Genome metrics calculation execution failed"
            fi
        fi
    fi

    # Final summary
    if [[ "$success" == true ]]; then
        log_success "=== MuDoGeR Prokaryotes Module Completed Successfully ==="
        log_info "Results can be found in: $libname_folder/prokaryotes"
        log_info "Unique bins: $unique_bins_dir"
        log_info "Metrics: $metrics_dir"
        log_info "Logs: $LOG_FILE"
    else
        log_error "=== MuDoGeR Prokaryotes Module Failed ==="
        log_error "Check the log file for details: $LOG_FILE"
        exit 1
    fi
}

# Run main function
main "$@"
