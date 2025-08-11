#!/bin/bash

# Enhanced Module 2-5: Bin Taxonomy Assignment with improved logging and verification
# Usage: mudoger-module-2-5_bin-taxonomy.sh <prokaryotes_folder> <cores>

echo '=================================================='
echo '      MuDoGeR Module 2-5: Bin Taxonomy Assignment'
echo '=================================================='

# Initialize logging
LOG_DIR="$1/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/module-2-5_bin-taxonomy.log"
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

# Function for logging with timestamps
log_message() {
    echo "[$TIMESTAMP] $1" | tee -a "$LOG_FILE"
}

# Function to check if previous run was successful
check_previous_success() {
    local master_folder="$1"
    local expected_files=(
        "$master_folder/metrics/GTDBtk_taxonomy/gtdbtk_result.tsv"
        "$master_folder/metrics/GTDBtk_taxonomy/gtdbtk.bac120.summary.tsv"
        "$master_folder/metrics/GTDBtk_taxonomy/gtdbtk.ar53.summary.tsv"
    )

    # Check if at least the main result file exists and is not empty
    if [[ -f "${expected_files[0]}" && -s "${expected_files[0]}" ]]; then
        local line_count=$(wc -l < "${expected_files[0]}")
        if [[ $line_count -gt 1 ]]; then  # More than just header
            return 0  # Success
        fi
    fi
    return 1  # Not successful
}

# Function to verify input requirements
verify_inputs() {
    local master_folder="$1"
    local cores="$2"

    log_message "Verifying input requirements..."

    # Check if unique_bins directory exists and has .fa files
    if [[ ! -d "$master_folder/binning/unique_bins" ]]; then
        log_message "ERROR: unique_bins directory not found at $master_folder/binning/unique_bins"
        return 1
    fi

    local bin_count=$(find "$master_folder/binning/unique_bins" -name "*.fa" | wc -l)
    if [[ $bin_count -eq 0 ]]; then
        log_message "ERROR: No .fa files found in unique_bins directory"
        return 1
    fi

    log_message "Found $bin_count bin files to process"

    # Check cores parameter
    if ! [[ "$cores" =~ ^[0-9]+$ ]] || [[ $cores -lt 1 ]]; then
        log_message "ERROR: Invalid cores parameter: $cores"
        return 1
    fi

    log_message "Using $cores CPU cores"
    return 0
}

# Function to setup conda environment
setup_environment() {
    log_message "Setting up conda environment..."

    conda activate mudoger_env
    if [[ $? -ne 0 ]]; then
        log_message "ERROR: Failed to activate mudoger_env"
        return 1
    fi

    config_path="$(which config.sh)"
    if [[ ! -f "$config_path" ]]; then
        log_message "ERROR: config.sh not found"
        return 1
    fi

    database="${config_path/config/database}"
    source "$config_path"
    source "$database"

    # Activate GTDBtk environment
    conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH/gtdbtk_env"
    if [[ $? -ne 0 ]]; then
        log_message "ERROR: Failed to activate gtdbtk_env"
        return 1
    fi

    # Setup GTDBTK database path
    GTDBTK_DATA_PATH=$(realpath "$DATABASES_LOCATION"/gtdbtk/release*)
    if [[ ! -d "$GTDBTK_DATA_PATH" ]]; then
        log_message "ERROR: GTDBTK database not found at $DATABASES_LOCATION/gtdbtk/"
        return 1
    fi

    export GTDBTK_DATA_PATH="$GTDBTK_DATA_PATH"
    log_message "GTDBTK database path set to: $GTDBTK_DATA_PATH"

    return 0
}

# Function to run GTDBtk classification
run_gtdbtk() {
    local master_folder="$1"
    local cores="$2"

    log_message "Starting GTDBtk classification workflow..."

    local output_dir="$master_folder/metrics/GTDBtk_taxonomy"
    local genome_dir="$master_folder/binning/unique_bins"

    mkdir -p "$output_dir"

    # Run GTDBtk classify workflow
    log_message "Running: gtdbtk classify_wf --skip_ani_screen --extension fa --cpus $cores --genome_dir $genome_dir --out_dir $output_dir"

    gtdbtk classify_wf --skip_ani_screen --extension fa --cpus "$cores" --genome_dir "$genome_dir" --out_dir "$output_dir" 2>&1 | tee -a "$LOG_FILE"

    local gtdbtk_exit_code=${PIPESTATUS[0]}
    if [[ $gtdbtk_exit_code -ne 0 ]]; then
        log_message "ERROR: GTDBtk failed with exit code $gtdbtk_exit_code"
        return 1
    fi

    log_message "GTDBtk classification completed successfully"
    return 0
}

# Function to create merged results file
create_merged_results() {
    local master_folder="$1"
    local output_dir="$master_folder/metrics/GTDBtk_taxonomy"

    log_message "Creating merged results file..."

    # Check if summary files exist
    local bac_summary="$output_dir/gtdbtk.bac120.summary.tsv"
    local arc_summary="$output_dir/gtdbtk.ar53.summary.tsv"
    local merged_results="$output_dir/gtdbtk_result.tsv"

    if [[ -f "$bac_summary" ]] || [[ -f "$arc_summary" ]]; then
        # Use awk to merge files, handling headers properly
        awk 'FNR==1 && NR!=1 {next;}{print}' "$output_dir"/gtdbtk.*summ*.tsv > "$merged_results"

        if [[ -s "$merged_results" ]]; then
            local result_count=$(tail -n +2 "$merged_results" | wc -l)
            log_message "Merged results file created with $result_count classified genomes"
            return 0
        else
            log_message "ERROR: Merged results file is empty"
            return 1
        fi
    else
        log_message "ERROR: No GTDBtk summary files found"
        return 1
    fi
}

# Function to verify outputs
verify_outputs() {
    local master_folder="$1"
    local output_dir="$master_folder/metrics/GTDBtk_taxonomy"

    log_message "Verifying outputs..."

    # Check main result file
    local merged_results="$output_dir/gtdbtk_result.tsv"
    if [[ ! -f "$merged_results" ]] || [[ ! -s "$merged_results" ]]; then
        log_message "ERROR: Main results file missing or empty"
        return 1
    fi

    # Count results
    local total_bins=$(find "$master_folder/binning/unique_bins" -name "*.fa" | wc -l)
    local classified_bins=$(tail -n +2 "$merged_results" | wc -l)

    log_message "Summary: $classified_bins out of $total_bins bins were classified"

    if [[ $classified_bins -eq 0 ]]; then
        log_message "WARNING: No bins were successfully classified"
        return 1
    fi

    # Check for expected output files
    local expected_files=(
        "$output_dir/classify/gtdbtk.bac120.summary.tsv"
        "$output_dir/classify/gtdbtk.ar53.summary.tsv"
        "$output_dir/identify/gtdbtk.bac120.markers_summary.tsv"
        "$output_dir/identify/gtdbtk.ar53.markers_summary.tsv"
    )

    local found_files=0
    for file in "${expected_files[@]}"; do
        if [[ -f "$file" ]]; then
            ((found_files++))
            log_message "Found: $(basename "$file")"
        fi
    done

    log_message "Found $found_files expected output files"

    return 0
}

# Main execution
main() {
    local master_folder="$1"
    local cores="$2"

    log_message "Starting Module 2-5: Bin Taxonomy Assignment"
    log_message "Input folder: $master_folder"
    log_message "CPU cores: $cores"

    # Check if already completed
    if check_previous_success "$master_folder"; then
        log_message "SUCCESS: Bin taxonomy assignment already completed successfully"
        log_message "Results available at: $master_folder/metrics/GTDBtk_taxonomy/gtdbtk_result.tsv"

        # Display summary
        local classified_count=$(tail -n +2 "$master_folder/metrics/GTDBtk_taxonomy/gtdbtk_result.tsv" | wc -l)
        log_message "Previously classified $classified_count genomes"

        echo "-> Bin taxonomy assignment is done. Please check: $master_folder/metrics/GTDBtk_taxonomy"
        return 0
    fi

    log_message "Previous run not found or incomplete, starting fresh analysis..."

    # Verify inputs
    if ! verify_inputs "$master_folder" "$cores"; then
        log_message "ERROR: Input verification failed"
        return 1
    fi

    # Setup environment
    if ! setup_environment; then
        log_message "ERROR: Environment setup failed"
        return 1
    fi

    # Run GTDBtk
    if ! run_gtdbtk "$master_folder" "$cores"; then
        log_message "ERROR: GTDBtk classification failed"
        conda deactivate
        return 1
    fi

    # Create merged results
    if ! create_merged_results "$master_folder"; then
        log_message "ERROR: Failed to create merged results"
        conda deactivate
        return 1
    fi

    # Verify outputs
    if ! verify_outputs "$master_folder"; then
        log_message "ERROR: Output verification failed"
        conda deactivate
        return 1
    fi

    conda deactivate

    log_message "SUCCESS: Module 2-5 completed successfully"
    log_message "Results saved to: $master_folder/metrics/GTDBtk_taxonomy/"
    echo "-> Bin taxonomy assignment completed. Results: $master_folder/metrics/GTDBtk_taxonomy"

    return 0
}

# Script entry point
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <prokaryotes_folder> <cores>"
    echo "Example: $0 /path/to/sample/prokaryotes 8"
    exit 1
fi

main "$1" "$2"
exit_code=$?

if [[ $exit_code -eq 0 ]]; then
    echo "Module 2-5 completed successfully!"
else
    echo "Module 2-5 failed. Check log file: $1/logs/module-2-5_bin-taxonomy.log"
fi

exit $exit_code