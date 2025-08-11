#!/bin/bash

# Enhanced Module 2-7: Prokka Annotation with improved logging and verification
# Usage: mudoger-module-2-7_prokka.sh <prokaryotes_folder> <cores>

echo '=================================================='
echo '      MuDoGeR Module 2-7: Prokka Annotation'
echo '=================================================='

# Initialize logging
LOG_DIR="$1/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/module-2-7_prokka.log"
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

# Function for logging with timestamps
log_message() {
    echo "[$TIMESTAMP] $1" | tee -a "$LOG_FILE"
}

# Function to check if previous run was successful
check_previous_success() {
    local prokaryotes_folder="$1"
    local prokka_dir="$prokaryotes_folder/metrics/prokka"

    # Check if prokka directory exists
    if [[ ! -d "$prokka_dir" ]]; then
        return 1
    fi

    # Count bins and check if all have been processed
    local total_bins=$(find "$prokaryotes_folder/binning/unique_bins" -name "*.fa" | wc -l)
    local processed_bins=0

    # Check each bin for successful prokka output
    for bin_file in "$prokaryotes_folder/binning/unique_bins"/*.fa; do
        if [[ -f "$bin_file" ]]; then
            local bin_name=$(basename "$bin_file" .fa)
            local bin_prokka_dir="$prokka_dir/$bin_name"

            # Check for key prokka output files
            if [[ -f "$bin_prokka_dir/PROKKA_"*".tsv" ]] && \
               [[ -f "$bin_prokka_dir/PROKKA_"*".faa" ]] && \
               [[ -f "$bin_prokka_dir/PROKKA_"*".gff" ]]; then
                # Verify that the TSV file has content (more than just header)
                local tsv_file=$(find "$bin_prokka_dir" -name "PROKKA_*.tsv" | head -1)
                if [[ -s "$tsv_file" ]]; then
                    local line_count=$(wc -l < "$tsv_file")
                    if [[ $line_count -gt 1 ]]; then
                        ((processed_bins++))
                    fi
                fi
            fi
        fi
    done

    # Consider successful if all bins are processed
    if [[ $processed_bins -eq $total_bins ]] && [[ $total_bins -gt 0 ]]; then
        return 0
    fi

    return 1
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

    local bin_count=$(find "$bins_dir" -name "*.fa" -type f | wc -l)
    if [[ $bin_count -eq 0 ]]; then
        log_message "ERROR: No .fa files found in unique_bins directory"
        return 1
    fi

    log_message "Found $bin_count bin files to process"

    # Verify that bin files are valid FASTA format
    local invalid_files=0
    for bin_file in "$bins_dir"/*.fa; do
        if [[ -f "$bin_file" ]]; then
            if ! head -1 "$bin_file" | grep -q "^>"; then
                log_message "WARNING: $bin_file does not appear to be valid FASTA format"
                ((invalid_files++))
            fi
        fi
    done

    if [[ $invalid_files -gt 0 ]]; then
        log_message "WARNING: Found $invalid_files files with invalid FASTA format"
    fi

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

    # Activate prokka environment
    conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH/prokka_env"
    if [[ $? -ne 0 ]]; then
        log_message "ERROR: Failed to activate prokka_env"
        return 1
    fi

    # Verify prokka is available
    if ! command -v prokka >/dev/null 2>&1; then
        log_message "ERROR: prokka command not found in PATH"
        return 1
    fi

    local prokka_version=$(prokka --version 2>&1 | head -1)
    log_message "Using prokka: $prokka_version"

    return 0
}

# Function to run prokka annotation
run_prokka() {
    local prokaryotes_folder="$1"
    local cores="$2"

    log_message "Starting Prokka annotation..."

    local input_bins_folder="$prokaryotes_folder/binning/unique_bins"
    local output_results="$prokaryotes_folder/metrics/prokka"

    mkdir -p "$output_results"

    local total_bins=0
    local successful_bins=0
    local failed_bins=0

    # Process each bin
    for bin_file in "$input_bins_folder"/*.fa; do
        if [[ -f "$bin_file" ]]; then
            ((total_bins++))
            local bin_name=$(basename "$bin_file" .fa)
            local bin_output_dir="$output_results/$bin_name"

            log_message "Processing bin $total_bins: $bin_name"

            # Check if this bin was already processed successfully
            if [[ -f "$bin_output_dir/PROKKA_"*".tsv" ]] && \
               [[ -f "$bin_output_dir/PROKKA_"*".faa" ]] && \
               [[ -f "$bin_output_dir/PROKKA_"*".gff" ]]; then
                local tsv_file=$(find "$bin_output_dir" -name "PROKKA_*.tsv" | head -1)
                if [[ -s "$tsv_file" ]]; then
                    local line_count=$(wc -l < "$tsv_file")
                    if [[ $line_count -gt 1 ]]; then
                        log_message "  Bin $bin_name already processed successfully, skipping"
                        ((successful_bins++))
                        continue
                    fi
                fi
            fi

            # Run prokka for this bin
            log_message "  Running prokka for $bin_name..."

            # Remove existing output directory to ensure clean run
            if [[ -d "$bin_output_dir" ]]; then
                rm -rf "$bin_output_dir"
            fi

            prokka "$bin_file" --cpus "$cores" --outdir "$bin_output_dir" --prefix "PROKKA_$bin_name" --metagenome --quiet 2>&1 | tee -a "$LOG_FILE"

            local prokka_exit_code=${PIPESTATUS[0]}
            if [[ $prokka_exit_code -eq 0 ]]; then
                # Verify output files were created
                if [[ -f "$bin_output_dir/PROKKA_${bin_name}.tsv" ]] && \
                   [[ -s "$bin_output_dir/PROKKA_${bin_name}.tsv" ]]; then
                    local gene_count=$(tail -n +2 "$bin_output_dir/PROKKA_${bin_name}.tsv" | wc -l)
                    log_message "  SUCCESS: $bin_name annotated with $gene_count genes"
                    ((successful_bins++))
                else
                    log_message "  ERROR: $bin_name - prokka completed but output files missing or empty"
                    ((failed_bins++))
                fi
            else
                log_message "  ERROR: $bin_name - prokka failed with exit code $prokka_exit_code"
                ((failed_bins++))
            fi
        fi
    done

    log_message "Prokka annotation summary:"
    log_message "  Total bins: $total_bins"
    log_message "  Successful: $successful_bins"
    log_message "  Failed: $failed_bins"

    if [[ $failed_bins -gt 0 ]]; then
        log_message "WARNING: Some bins failed annotation"
        if [[ $successful_bins -eq 0 ]]; then
            log_message "ERROR: No bins were successfully annotated"
            return 1
        fi
    fi

    log_message "Prokka annotation completed"
    return 0
}

# Function to verify outputs and generate summary
verify_outputs() {
    local prokaryotes_folder="$1"
    local prokka_dir="$prokaryotes_folder/metrics/prokka"

    log_message "Verifying outputs and generating summary..."

    if [[ ! -d "$prokka_dir" ]]; then
        log_message "ERROR: Prokka output directory not found"
        return 1
    fi

    local total_bins=$(find "$prokaryotes_folder/binning/unique_bins" -name "*.fa" | wc -l)
    local processed_bins=0
    local total_genes=0
    local total_known_genes=0
    local total_hypothetical_genes=0

    # Check each bin's prokka output
    for bin_file in "$prokaryotes_folder/binning/unique_bins"/*.fa; do
        if [[ -f "$bin_file" ]]; then
            local bin_name=$(basename "$bin_file" .fa)
            local bin_prokka_dir="$prokka_dir/$bin_name"

            if [[ -d "$bin_prokka_dir" ]]; then
                local tsv_file=$(find "$bin_prokka_dir" -name "PROKKA_*.tsv" | head -1)
                if [[ -f "$tsv_file" && -s "$tsv_file" ]]; then
                    ((processed_bins++))

                    # Count genes
                    local bin_genes=$(tail -n +2 "$tsv_file" | wc -l)
                    local bin_known=$(tail -n +2 "$tsv_file" | grep -v "hypothetical" | wc -l)
                    local bin_hypothetical=$(tail -n +2 "$tsv_file" | grep "hypothetical" | wc -l)

                    ((total_genes += bin_genes))
                    ((total_known_genes += bin_known))
                    ((total_hypothetical_genes += bin_hypothetical))

                    log_message "  $bin_name: $bin_genes genes ($bin_known known, $bin_hypothetical hypothetical)"
                fi
            fi
        fi
    done

    log_message "Annotation summary:"
    log_message "  Processed bins: $processed_bins/$total_bins"
    log_message "  Total genes: $total_genes"
    log_message "  Known genes: $total_known_genes"
    log_message "  Hypothetical genes: $total_hypothetical_genes"

    if [[ $total_genes -gt 0 ]]; then
        local known_percentage=$(awk "BEGIN {printf \"%.1f\", ($total_known_genes/$total_genes)*100}")
        log_message "  Known genes percentage: $known_percentage%"
    fi

    if [[ $processed_bins -eq 0 ]]; then
        log_message "ERROR: No bins were successfully processed"
        return 1
    fi

    # Check for expected file types
    local expected_extensions=("tsv" "gff" "faa" "ffn" "fna")
    for ext in "${expected_extensions[@]}"; do
        local count=$(find "$prokka_dir" -name "*.$ext" | wc -l)
        log_message "  Found $count .$ext files"
    done

    return 0
}

# Function to cleanup temporary files
cleanup_temp_files() {
    local prokaryotes_folder="$1"
    local prokka_dir="$prokaryotes_folder/metrics/prokka"

    log_message "Cleaning up temporary files..."

    if [[ -d "$prokka_dir" ]]; then
        # Remove prokka temporary files and logs
        find "$prokka_dir" -name "*.tmp" -delete 2>/dev/null || true
        find "$prokka_dir" -name "*.log" -delete 2>/dev/null || true
        find "$prokka_dir" -name "*.err" -delete 2>/dev/null || true
    fi
}

# Main execution
main() {
    local prokaryotes_folder="$1"
    local cores="$2"

    log_message "Starting Module 2-7: Prokka Annotation"
    log_message "Input folder: $prokaryotes_folder"
    log_message "CPU cores: $cores"

    # Check if already completed
    if check_previous_success "$prokaryotes_folder"; then
        log_message "SUCCESS: Prokka annotation already completed successfully"
        log_message "Results available at: $prokaryotes_folder/metrics/prokka"

        # Display summary of previous results
        local total_bins=$(find "$prokaryotes_folder/binning/unique_bins" -name "*.fa" | wc -l)
        local processed_bins=0

        for bin_file in "$prokaryotes_folder/binning/unique_bins"/*.fa; do
            if [[ -f "$bin_file" ]]; then
                local bin_name=$(basename "$bin_file" .fa)
                local tsv_file=$(find "$prokaryotes_folder/metrics/prokka/$bin_name" -name "PROKKA_*.tsv" 2>/dev/null | head -1)
                if [[ -f "$tsv_file" && -s "$tsv_file" ]]; then
                    ((processed_bins++))
                fi
            fi
        done

        log_message "Previously processed $processed_bins/$total_bins bins"

        echo "-> Rapid prokka annotation is done. Please check: $prokaryotes_folder/metrics/prokka"
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

    # Run prokka
    if ! run_prokka "$prokaryotes_folder" "$cores"; then
        log_message "ERROR: Prokka annotation failed"
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

    log_message "SUCCESS: Module 2-7 completed successfully"
    log_message "Results saved to: $prokaryotes_folder/metrics/prokka/"
    echo "-> Prokka annotation completed. Results: $prokaryotes_folder/metrics/prokka"

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
    echo "Module 2-7 completed successfully!"
else
    echo "Module 2-7 failed. Check log file: $1/logs/module-2-7_prokka.log"
fi

exit $exit_code