#!/bin/bash

# Enhanced Module 2-6: Bin Quality Control with improved logging and verification
# Usage: mudoger-module-2-6_bin-QC.sh <prokaryotes_folder> <cores>

echo '=================================================='
echo '      MuDoGeR Module 2-6: Bin Quality Control'
echo '=================================================='

# Initialize logging
LOG_DIR="$1/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/module-2-6_bin-QC.log"
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

# Function for logging with timestamps
log_message() {
    echo "[$TIMESTAMP] $1" | tee -a "$LOG_FILE"
}

# Function to check if previous run was successful
check_previous_success() {
    local prokaryotes_folder="$1"
    local output_file="$prokaryotes_folder/metrics/checkm_qc/outputcheckm.tsv"

    # Check if output file exists and is not empty
    if [[ -f "$output_file" && -s "$output_file" ]]; then
        local line_count=$(wc -l < "$output_file")
        if [[ $line_count -gt 1 ]]; then  # More than just header
            # Check if the file has proper CheckM format (should have multiple columns)
            local column_count=$(head -n 1 "$output_file" | tr -cd '\t' | wc -c)
            if [[ $column_count -ge 10 ]]; then  # CheckM output should have many columns
                return 0  # Success
            fi
        fi
    fi
    return 1  # Not successful
}

# Function to verify input requirements
verify_inputs() {
    local prokaryotes_folder="$1"
    local cores="$2"

    log_message "Verifying input requirements..."

    # Check if unique_bins directory exists and has .fa files
    local bins_dir="$prokaryotes_folder/binning/unique_bins"
    if [[ ! -d "$bins_dir" ]]; then
        log_message "ERROR: unique_bins directory not found at $bins_dir"
        return 1
    fi

    local bin_count=$(find "$bins_dir" -name "*.fa" | wc -l)
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

    # Activate metawrap environment (which includes CheckM)
    conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH/metawrap_env"
    if [[ $? -ne 0 ]]; then
        log_message "ERROR: Failed to activate metawrap_env"
        return 1
    fi

    # Setup CheckM database
    local checkm_db="$DATABASES_LOCATION/checkm"
    if [[ ! -d "$checkm_db" ]]; then
        log_message "ERROR: CheckM database not found at $checkm_db"
        return 1
    fi

    log_message "Setting CheckM database root to: $checkm_db"
    echo "$checkm_db" | checkm data setRoot "$checkm_db" 2>&1 | tee -a "$LOG_FILE"

    if [[ ${PIPESTATUS[1]} -ne 0 ]]; then
        log_message "ERROR: Failed to set CheckM database root"
        return 1
    fi

    return 0
}

# Function to run CheckM quality assessment
run_checkm() {
    local prokaryotes_folder="$1"
    local cores="$2"

    log_message "Starting CheckM quality assessment..."

    local output_folder="$prokaryotes_folder/metrics/checkm_qc"
    local input_bins_folder="$prokaryotes_folder/binning/unique_bins"
    local extension="fa"

    # Create output folder
    mkdir -p "$output_folder"

    # Run CheckM lineage workflow
    log_message "Running: checkm lineage_wf -t $cores --reduced_tree --tab_table -x $extension -f $output_folder/outputcheckm.tsv $input_bins_folder $output_folder"

    checkm lineage_wf -t "$cores" --reduced_tree --tab_table -x "$extension" -f "$output_folder/outputcheckm.tsv" "$input_bins_folder" "$output_folder" 2>&1 | tee -a "$LOG_FILE"

    local checkm_exit_code=${PIPESTATUS[0]}
    if [[ $checkm_exit_code -ne 0 ]]; then
        log_message "ERROR: CheckM failed with exit code $checkm_exit_code"
        return 1
    fi

    log_message "CheckM quality assessment completed successfully"
    return 0
}

# Function to verify outputs and provide quality summary
verify_outputs() {
    local prokaryotes_folder="$1"
    local output_file="$prokaryotes_folder/metrics/checkm_qc/outputcheckm.tsv"

    log_message "Verifying outputs..."

    # Check if output file exists and has content
    if [[ ! -f "$output_file" ]] || [[ ! -s "$output_file" ]]; then
        log_message "ERROR: CheckM output file missing or empty"
        return 1
    fi

    # Count processed bins
    local total_bins=$(find "$prokaryotes_folder/binning/unique_bins" -name "*.fa" | wc -l)
    local processed_bins=$(tail -n +2 "$output_file" | wc -l)

    log_message "Summary: $processed_bins out of $total_bins bins were processed"

    if [[ $processed_bins -eq 0 ]]; then
        log_message "ERROR: No bins were successfully processed"
        return 1
    fi

    # Generate quality statistics
    log_message "Generating quality statistics..."

    # Extract completeness and contamination columns (assuming standard CheckM format)
    # CheckM output typically has: Bin Id, Marker lineage, # genomes, # markers, # marker sets, 0, 1, 2, 3, 4, 5+, Completeness, Contamination, Strain heterogeneity
    if command -v awk >/dev/null 2>&1; then
        local stats_output=$(awk -F'\t' '
        NR > 1 {
            comp = $(NF-2)  # Completeness is typically 3rd from last column
            cont = $(NF-1)  # Contamination is typically 2nd from last column
            if (comp != "" && cont != "") {
                total++
                if (comp >= 90 && cont <= 5) high_quality++
                else if (comp >= 50 && cont <= 10) medium_quality++
                else low_quality++
                total_comp += comp
                total_cont += cont
            }
        }
        END {
            if (total > 0) {
                printf "High quality (≥90%% complete, ≤5%% contamination): %d\n", high_quality
                printf "Medium quality (≥50%% complete, ≤10%% contamination): %d\n", medium_quality
                printf "Low quality: %d\n", low_quality
                printf "Average completeness: %.2f%%\n", total_comp/total
                printf "Average contamination: %.2f%%\n", total_cont/total
            }
        }' "$output_file")

        log_message "Quality distribution:"
        echo "$stats_output" | while read line; do
            log_message "  $line"
        done
    fi

    # Check for additional CheckM output files
    local checkm_output_dir="$prokaryotes_folder/metrics/checkm_qc"
    local expected_files=(
        "$checkm_output_dir/bins"
        "$checkm_output_dir/lineage.ms"
        "$checkm_output_dir/storage"
    )

    local found_files=0
    for item in "${expected_files[@]}"; do
        if [[ -e "$item" ]]; then
            ((found_files++))
            log_message "Found: $(basename "$item")"
        fi
    done

    log_message "Found $found_files expected CheckM output items"

    return 0
}

# Function to cleanup temporary files
cleanup_temp_files() {
    local prokaryotes_folder="$1"
    local checkm_output_dir="$prokaryotes_folder/metrics/checkm_qc"

    log_message "Cleaning up temporary files..."

    # CheckM creates some temporary files that can be safely removed
    if [[ -d "$checkm_output_dir" ]]; then
        find "$checkm_output_dir" -name "*.tmp" -delete 2>/dev/null || true
        find "$checkm_output_dir" -name "*.log" -size +100M -delete 2>/dev/null || true  # Remove very large log files
    fi
}

# Main execution
main() {
    local prokaryotes_folder="$1"
    local cores="$2"

    log_message "Starting Module 2-6: Bin Quality Control"
    log_message "Input folder: $prokaryotes_folder"
    log_message "CPU cores: $cores"

    # Check if already completed
    if check_previous_success "$prokaryotes_folder"; then
        log_message "SUCCESS: Bin quality control already completed successfully"
        log_message "Results available at: $prokaryotes_folder/metrics/checkm_qc/outputcheckm.tsv"

        # Display summary of previous results
        local processed_count=$(tail -n +2 "$prokaryotes_folder/metrics/checkm_qc/outputcheckm.tsv" | wc -l)
        log_message "Previously processed $processed_count bins"

        echo "-> MAGs quality control is done. Please check: $prokaryotes_folder/metrics/checkm_qc/outputcheckm.tsv"
        return 0
    fi

    log_message "Previous run not found or incomplete, starting fresh analysis..."

    # Verify inputs
    if ! verify_inputs "$prokaryotes_folder" "$cores"; then
        log_message "ERROR: Input verification failed"
        return 1
    fi

    # Setup environment
    if ! setup_environment; then
        log_message "ERROR: Environment setup failed"
        return 1
    fi

    # Run CheckM
    if ! run_checkm "$prokaryotes_folder" "$cores"; then
        log_message "ERROR: CheckM quality control failed"
        conda deactivate
        return 1
    fi

    # Verify outputs
    if ! verify_outputs "$prokaryotes_folder"; then
        log_message "ERROR: Output verification failed"
        conda deactivate
        return 1
    fi

    # Cleanup temporary files
    cleanup_temp_files "$prokaryotes_folder"

    conda deactivate

    log_message "SUCCESS: Module 2-6 completed successfully"
    log_message "Results saved to: $prokaryotes_folder/metrics/checkm_qc/"
    echo "-> MAGs quality control completed. Results: $prokaryotes_folder/metrics/checkm_qc/outputcheckm.tsv"

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
    echo "Module 2-6 completed successfully!"
else
    echo "Module 2-6 failed. Check log file: $1/logs/module-2-6_bin-QC.log"
fi

exit $exit_code