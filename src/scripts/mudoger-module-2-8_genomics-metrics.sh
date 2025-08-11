#!/bin/bash

# Enhanced Module 2-8: Genomics Metrics with improved logging and verification
# Usage: mudoger-module-2-8_genomics-metrics.sh <prokaryotes_folder> <cores>

echo '=================================================='
echo '      MuDoGeR Module 2-8: Genomics Metrics'
echo '=================================================='

# Initialize logging
LOG_DIR="$1/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/module-2-8_genomics-metrics.log"
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

# Function for logging with timestamps
log_message() {
    echo "[$TIMESTAMP] $1" | tee -a "$LOG_FILE"
}

# Function to check if previous run was successful
check_previous_success() {
    local prokaryotes_folder="$1"
    local expected_files=(
        "$prokaryotes_folder/metrics/genome_statistics/genome_metrics.tsv"
        "$prokaryotes_folder/metrics/genome_statistics/prok_genomes_stats.tsv"
        "$prokaryotes_folder/metrics/genome_statistics/bbtools.tsv"
        "$prokaryotes_folder/final_outputs/allbins_metrics_summary.tsv"
        "$prokaryotes_folder/final_outputs/mags_results_summary.tsv"
    )

    # Check if all expected files exist and are not empty
    for file in "${expected_files[@]}"; do
        if [[ ! -f "$file" ]] || [[ ! -s "$file" ]]; then
            return 1
        fi

        # Check if files have more than just headers
        local line_count=$(wc -l < "$file")
        if [[ $line_count -le 1 ]]; then
            return 1
        fi
    done

    # Check if final_outputs directory structure exists
    local final_dirs=(
        "$prokaryotes_folder/final_outputs/all_bins_seq"
        "$prokaryotes_folder/final_outputs/only_mags_seq"
        "$prokaryotes_folder/final_outputs/bins_metrics_summary"
        "$prokaryotes_folder/final_outputs/bins_genes_prokka_summary"
    )

    for dir in "${final_dirs[@]}"; do
        if [[ ! -d "$dir" ]]; then
            return 1
        fi
    done

    return 0
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

    # Check for required input files from previous modules
    local required_files=(
        "$prokaryotes_folder/metrics/GTDBtk_taxonomy/gtdbtk_result.tsv"
        "$prokaryotes_folder/metrics/checkm_qc/outputcheckm.tsv"
    )

    for file in "${required_files[@]}"; do
        if [[ ! -f "$file" ]] || [[ ! -s "$file" ]]; then
            log_message "ERROR: Required input file missing or empty: $file"
            return 1
        fi
    done

    # Check prokka directory
    local prokka_dir="$prokaryotes_folder/metrics/prokka"
    if [[ ! -d "$prokka_dir" ]]; then
        log_message "ERROR: Prokka directory not found at $prokka_dir"
        return 1
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

    return 0
}

# Function to calculate genome statistics
calculate_genome_stats() {
    local prokaryotes_folder="$1"

    log_message "Calculating genome statistics..."

    local output_path="$prokaryotes_folder/metrics/genome_statistics"
    local output_file="$output_path/prok_genomes_stats.tsv"

    mkdir -p "$output_path"

    # Create header
    echo -e "genome_name\tgenome_size\tnumber_of_scaffolds\tlargest_scaffold_size\tN50\tN90" > "$output_file"

    local processed_genomes=0
    local failed_genomes=0

    for genome in "$prokaryotes_folder/binning/unique_bins"/*.fa; do
        if [[ -f "$genome" ]]; then
            local genome_name=$(basename "$genome")
            log_message "Processing genome statistics for: $genome_name"

            # Create temporary directory for calculations
            local temp_path="$prokaryotes_folder/metrics/temp_$$"
            mkdir -p "$temp_path"

            # Calculate contig lengths and statistics
            if calculate_single_genome_stats "$genome" "$temp_path" "$output_file"; then
                ((processed_genomes++))
                log_message "  Successfully processed $genome_name"
            else
                ((failed_genomes++))
                log_message "  ERROR: Failed to process $genome_name"
            fi

            # Clean up temporary files
            rm -rf "$temp_path"
        fi
    done

    log_message "Genome statistics summary: $processed_genomes processed, $failed_genomes failed"

    if [[ $processed_genomes -eq 0 ]]; then
        log_message "ERROR: No genomes were successfully processed"
        return 1
    fi

    return 0
}

# Function to calculate statistics for a single genome
calculate_single_genome_stats() {
    local genome="$1"
    local temp_path="$2"
    local output_file="$3"

    local genome_name=$(basename "$genome")

    # Extract contig lengths, ordered from large to small
    cat "$genome" | awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' | \
    sed '/^>/ d' | awk '{ print length($0) }' | sort -gr > "$temp_path/contig_lengths.txt"

    if [[ ! -s "$temp_path/contig_lengths.txt" ]]; then
        return 1
    fi

    # Calculate basic statistics
    local Y=$(cat "$temp_path/contig_lengths.txt" | wc -l)  # number of contigs
    local X=$(paste -sd+ "$temp_path/contig_lengths.txt" | bc)  # sum of contig lengths

    if [[ -z "$X" ]] || [[ $X -eq 0 ]]; then
        return 1
    fi

    # Calculate cumulative lengths
    awk 'BEGIN {sum=0} {sum= sum+$0; print sum}' "$temp_path/contig_lengths.txt" > "$temp_path/contig_lengths_cum.txt"

    # Calculate cumulative percentages
    awk -v var=$X '{print $0/var}' "$temp_path/contig_lengths_cum.txt" > "$temp_path/cum_perc.txt"

    # Join results
    paste "$temp_path/contig_lengths.txt" "$temp_path/cum_perc.txt" > "$temp_path/matrix.txt"

    # Calculate N50 and N90
    local N50=$(awk '$2 >= 0.50' "$temp_path/matrix.txt" | head -1 | awk '{ print $1}')
    local N90=$(awk '$2 >= 0.90' "$temp_path/matrix.txt" | head -1 | awk '{ print $1}')
    local large_contig=$(head -1 "$temp_path/contig_lengths.txt")

    # Set defaults if calculations failed
    if [[ -z "$N50" ]]; then N50=0; fi
    if [[ -z "$N90" ]]; then N90=0; fi
    if [[ -z "$large_contig" ]]; then large_contig=0; fi

    # Write results
    echo -e "$genome_name\t$X\t$Y\t$large_contig\t$N50\t$N90" >> "$output_file"

    return 0
}

# Function to run BBTools statistics
run_bbtools() {
    local prokaryotes_folder="$1"

    log_message "Running BBTools statistics..."

    # Activate bbtools environment
    conda activate "$MUDOGER_DEPENDENCIES_ENVS_PATH/bbtools_env"
    if [[ $? -ne 0 ]]; then
        log_message "ERROR: Failed to activate bbtools_env"
        return 1
    fi

    local output_path="$prokaryotes_folder/metrics/genome_statistics"
    local bbtools_output="$output_path/bbtools.tsv"

    # Change to bins directory and run BBTools
    cd "$prokaryotes_folder/binning/unique_bins/"

    if [[ $(ls *.fa 2>/dev/null | wc -l) -eq 0 ]]; then
        log_message "ERROR: No .fa files found in unique_bins directory"
        cd - >/dev/null
        return 1
    fi

    log_message "Running statswrapper.sh on *.fa files..."
    "$MUDOGER_DEPENDENCIES_ENVS_PATH/bbtools_env/bin/statswrapper.sh" *.fa > "$bbtools_output" 2>&1

    local bbtools_exit_code=$?
    cd - >/dev/null

    if [[ $bbtools_exit_code -ne 0 ]]; then
        log_message "ERROR: BBTools failed with exit code $bbtools_exit_code"
        return 1
    fi

    if [[ ! -s "$bbtools_output" ]]; then
        log_message "ERROR: BBTools output file is empty"
        return 1
    fi

    local stat_count=$(tail -n +2 "$bbtools_output" | wc -l)
    log_message "BBTools statistics completed for $stat_count sequences"

    return 0
}

# Function to create comprehensive genome metrics
create_comprehensive_metrics() {
    local prokaryotes_folder="$1"

    log_message "Creating comprehensive genome metrics..."

    local output_path="$prokaryotes_folder/metrics/genome_statistics"
    local comprehensive_file="$output_path/genome_metrics.tsv"

    # Create header
    echo -e "OTU\tcompleteness\tcontamination\tstr.heterogeneity\ttaxonomy\tgenome_size\t#scaffolds\tlargest_scaff\tN50\tN90\tprokka_known\tprokka_unknown" > "$comprehensive_file"

    local processed_bins=0
    local total_bins=$(find "$prokaryotes_folder/binning/unique_bins" -name "*.fa" | wc -l)

    # Process each bin
    for bin_path in "$prokaryotes_folder/binning/unique_bins"/*; do
        if [[ -f "$bin_path" ]]; then
            local bin=$(basename "$bin_path" | sed "s/.fa//g")

            log_message "Creating comprehensive metrics for: $bin"

            # Extract data from various sources
            local tax=$(grep "$bin" "$prokaryotes_folder/metrics/GTDBtk_taxonomy/gtdbtk_result.tsv" | cut -f2 | head -1)
            local qual=$(grep "$bin" "$prokaryotes_folder/metrics/checkm_qc/outputcheckm.tsv" | cut -f12,13,14 | head -1)
            local metrics=$(grep "$bin" "$output_path/prok_genomes_stats.tsv" | cut -f2-6 | head -1)

            # Count prokka genes
            local prokka_files=("$prokaryotes_folder/metrics/prokka/$bin"/PROKKA*.tsv)
            local prokka_known=0
            local prokka_unknown=0

            if [[ -f "${prokka_files[0]}" ]]; then
                prokka_known=$(tail -n +2 "${prokka_files[0]}" | grep -v "hypothetical" | wc -l)
                prokka_unknown=$(tail -n +2 "${prokka_files[0]}" | grep "hypothetical" | wc -l)
            fi

            # Handle missing data
            if [[ -z "$tax" ]]; then tax="Unclassified"; fi
            if [[ -z "$qual" ]]; then qual="0\t0\t0"; fi
            if [[ -z "$metrics" ]]; then metrics="0\t0\t0\t0\t0"; fi

            # Write comprehensive metrics
            echo -e "$bin\t$qual\t$tax\t$metrics\t$prokka_known\t$prokka_unknown" >> "$comprehensive_file"
            ((processed_bins++))
        fi
    done

    log_message "Comprehensive metrics created for $processed_bins/$total_bins bins"

    if [[ $processed_bins -eq 0 ]]; then
        log_message "ERROR: No comprehensive metrics were created"
        return 1
    fi

    return 0
}

# Function to filter MAGs and create final outputs
create_final_outputs() {
    local prokaryotes_folder="$1"

    log_message "Creating final outputs and filtering MAGs..."

    local output_path="$prokaryotes_folder/metrics/genome_statistics"
    local genome_metrics="$output_path/genome_metrics.tsv"

    # Filter high-quality MAGs (Completeness - 5*Contamination >= 50)
    log_message "Filtering MAGs based on quality score (Completeness - 5*Contamination >= 50)..."

    awk 'BEGIN {FS="\t"; OFS="\t"} NR==1 {print} NR>1 && ($2 - (5*$3)) >= 50' "$genome_metrics" > "$prokaryotes_folder/MAGS_results.tsv"

    local mag_count=$(tail -n +2 "$prokaryotes_folder/MAGS_results.tsv" | wc -l)
    log_message "Identified $mag_count high-quality MAGs"

    # Create final output directories
    log_message "Creating final output directory structure..."

    local final_outputs="$prokaryotes_folder/final_outputs"
    mkdir -p "$final_outputs/all_bins_seq"
    mkdir -p "$final_outputs/only_mags_seq"
    mkdir -p "$final_outputs/bins_metrics_summary"
    mkdir -p "$final_outputs/bins_genes_prokka_summary"

    # Copy all bin sequences
    log_message "Copying all bin sequences..."
    cp "$prokaryotes_folder/binning/unique_bins"/*.fa "$final_outputs/all_bins_seq/" 2>/dev/null || {
        log_message "ERROR: Failed to copy bin sequences"
        return 1
    }

    # Copy taxonomy and quality summaries
    log_message "Copying taxonomy and quality summaries..."
    cp "$prokaryotes_folder/metrics/GTDBtk_taxonomy/gtdbtk_result.tsv" "$final_outputs/bins_metrics_summary/taxa_bins_gtdbtk_summary.tsv"
    cp "$prokaryotes_folder/metrics/checkm_qc/outputcheckm.tsv" "$final_outputs/bins_metrics_summary/qual_bins_checkm_summary.tsv"

    # Copy gene annotations
    log_message "Copying gene annotations..."
    local copied_annotations=0
    for bin_path in "$prokaryotes_folder/binning/unique_bins"/*.fa; do
        if [[ -f "$bin_path" ]]; then
            local bin=$(basename "$bin_path" | sed "s/.fa//g")
            local prokka_tsv=("$prokaryotes_folder/metrics/prokka/$bin"/PROKKA*.tsv)

            if [[ -f "${prokka_tsv[0]}" ]]; then
                cp "${prokka_tsv[0]}" "$final_outputs/bins_genes_prokka_summary/${bin}_genes_prokka.tsv"
                ((copied_annotations++))
            fi
        fi
    done

    log_message "Copied $copied_annotations gene annotation files"

    # Copy complete summary
    cp "$genome_metrics" "$final_outputs/allbins_metrics_summary.tsv"

    # Move MAGs summary
    mv "$prokaryotes_folder/MAGS_results.tsv" "$final_outputs/mags_results_summary.tsv"

    # Copy only MAG sequences
    log_message "Copying MAG sequences..."
    local copied_mags=0
    if [[ -s "$final_outputs/mags_results_summary.tsv" ]]; then
        cut -f1 "$final_outputs/mags_results_summary.tsv" | tail -n+2 > "$final_outputs/tmp_mag_list"

        while IFS= read -r mag; do
            if [[ -f "$prokaryotes_folder/binning/unique_bins/$mag.fa" ]]; then
                cp "$prokaryotes_folder/binning/unique_bins/$mag.fa" "$final_outputs/only_mags_seq/"
                ((copied_mags++))
            fi
        done < "$final_outputs/tmp_mag_list"

        rm -f "$final_outputs/tmp_mag_list"
    fi

    log_message "Copied $copied_mags MAG sequences"

    return 0
}

# Function to verify final outputs
verify_final_outputs() {
    local prokaryotes_folder="$1"

    log_message "Verifying final outputs..."

    local final_outputs="$prokaryotes_folder/final_outputs"
    local expected_files=(
        "$final_outputs/allbins_metrics_summary.tsv"
        "$final_outputs/mags_results_summary.tsv"
        "$final_outputs/bins_metrics_summary/taxa_bins_gtdbtk_summary.tsv"
        "$final_outputs/bins_metrics_summary/qual_bins_checkm_summary.tsv"
    )

    local expected_dirs=(
        "$final_outputs/all_bins_seq"
        "$final_outputs/only_mags_seq"
        "$final_outputs/bins_genes_prokka_summary"
    )

    # Check files
    local missing_files=0
    for file in "${expected_files[@]}"; do
        if [[ ! -f "$file" ]] || [[ ! -s "$file" ]]; then
            log_message "ERROR: Missing or empty file: $file"
            ((missing_files++))
        fi
    done

    # Check directories
    local missing_dirs=0
    for dir in "${expected_dirs[@]}"; do
        if [[ ! -d "$dir" ]]; then
            log_message "ERROR: Missing directory: $dir"
            ((missing_dirs++))
        fi
    done

    # Count contents
    local all_bins_count=$(find "$final_outputs/all_bins_seq" -name "*.fa" | wc -l)
    local mags_count=$(find "$final_outputs/only_mags_seq" -name "*.fa" | wc -l)
    local annotation_count=$(find "$final_outputs/bins_genes_prokka_summary" -name "*_genes_prokka.tsv" | wc -l)

    log_message "Final output verification:"
    log_message "  All bins: $all_bins_count sequences"
    log_message "  MAGs: $mags_count sequences"
    log_message "  Gene annotations: $annotation_count files"
    log_message "  Missing files: $missing_files"
    log_message "  Missing directories: $missing_dirs"

    # Generate final summary
    if [[ -f "$final_outputs/mags_results_summary.tsv" ]]; then
        local mag_summary_count=$(tail -n +2 "$final_outputs/mags_results_summary.tsv" | wc -l)
        if [[ $mag_summary_count -gt 0 ]]; then
            log_message "Quality summary of MAGs:"

            # Calculate quality statistics for MAGs
            local avg_completeness=$(awk -F'\t' 'NR>1 {sum+=$2; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}' "$final_outputs/mags_results_summary.tsv")
            local avg_contamination=$(awk -F'\t' 'NR>1 {sum+=$3; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}' "$final_outputs/mags_results_summary.tsv")

            log_message "  Average completeness: $avg_completeness%"
            log_message "  Average contamination: $avg_contamination%"

            # Count MAGs by quality levels
            local high_quality=$(awk -F'\t' 'NR>1 && $2>=90 && $3<=5' "$final_outputs/mags_results_summary.tsv" | wc -l)
            local medium_quality=$(awk -F'\t' 'NR>1 && $2>=50 && $2<90 && $3<=10' "$final_outputs/mags_results_summary.tsv" | wc -l)

            log_message "  High quality MAGs (≥90% complete, ≤5% contamination): $high_quality"
            log_message "  Medium quality MAGs (≥50% complete, ≤10% contamination): $medium_quality"
        fi
    fi

    if [[ $missing_files -gt 0 ]] || [[ $missing_dirs -gt 0 ]]; then
        log_message "ERROR: Final output verification failed"
        return 1
    fi

    return 0
}

# Function to cleanup temporary files
cleanup_temp_files() {
    local prokaryotes_folder="$1"

    log_message "Cleaning up temporary files..."

    # Remove any remaining temporary directories
    find "$prokaryotes_folder/metrics" -name "temp_*" -type d -exec rm -rf {} + 2>/dev/null || true

    # Remove temporary files
    find "$prokaryotes_folder" -name "tmp_*" -delete 2>/dev/null || true
}

# Main execution function
main() {
    local prokaryotes_folder="$1"
    local cores="$2"

    log_message "Starting Module 2-8: Genomics Metrics"
    log_message "Input folder: $prokaryotes_folder"
    log_message "CPU cores: $cores"

    # Check if already completed
    if check_previous_success "$prokaryotes_folder"; then
        log_message "SUCCESS: Genomics metrics already completed successfully"
        log_message "Results available at: $prokaryotes_folder/final_outputs/"

        # Display summary of previous results
        local all_bins=$(find "$prokaryotes_folder/final_outputs/all_bins_seq" -name "*.fa" | wc -l)
        local mags=$(find "$prokaryotes_folder/final_outputs/only_mags_seq" -name "*.fa" | wc -l)

        log_message "Previously processed: $all_bins total bins, $mags MAGs"

        echo "-> MAGs statistics is done. Please check: $prokaryotes_folder/metrics/genome_statistics/prok_genomes_stats.tsv"
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

    # Calculate genome statistics
    if ! calculate_genome_stats "$prokaryotes_folder"; then
        log_message "ERROR: Genome statistics calculation failed"
        return 1
    fi

    # Run BBTools
    if ! run_bbtools "$prokaryotes_folder"; then
        log_message "ERROR: BBTools analysis failed"
        return 1
    fi

    # Create comprehensive metrics
    if ! create_comprehensive_metrics "$prokaryotes_folder"; then
        log_message "ERROR: Comprehensive metrics creation failed"
        return 1
    fi

    # Create final outputs
    if ! create_final_outputs "$prokaryotes_folder"; then
        log_message "ERROR: Final outputs creation failed"
        return 1
    fi

    # Verify final outputs
    if ! verify_final_outputs "$prokaryotes_folder"; then
        log_message "ERROR: Final output verification failed"
        return 1
    fi

    # Cleanup temporary files
    cleanup_temp_files "$prokaryotes_folder"

    # Deactivate conda environment
    conda deactivate

    log_message "SUCCESS: Module 2-8 completed successfully"
    log_message "Results saved to: $prokaryotes_folder/final_outputs/"
    echo "-> MAGs statistics completed. Results: $prokaryotes_folder/final_outputs/"

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
    echo "Module 2-8 completed successfully!"
else
    echo "Module 2-8 failed. Check log file: $1/logs/module-2-8_genomics-metrics.log"
fi

exit $exit_code