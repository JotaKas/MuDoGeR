#!/bin/bash


input_uvigs="$1"
output_checkv="$2"
cores="$3"


source activate conda-vcheck
export CHECKVDB=/data/msb/tools/vcheck/checkv-db-v0.6


checkv end_to_end  "$input_uvigs" "$output_checkv" -t "$cores"
