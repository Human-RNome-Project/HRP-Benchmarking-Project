#!/usr/bin/env bash
set -euo pipefail

METADATA="unified_metadata_ont_0430.tsv"

# Ensure we're executing in wdl_pipelines directory
cd "$(dirname "$0")"

# Clean up OUTPUT_JSON format in the generator script to remove the _2 suffix
sed -i 's/OUTPUT_JSON="tRNA\/inputs_${SAMPLE_ID}_2.json"/OUTPUT_JSON="tRNA\/inputs_${SAMPLE_ID}.json"/' generate_tRNA_inputs.sh

for i in {1..9}; do
    bc="bc${i}"
    
    # Process both HRP_A_020 and HRP_A_021 sample prefixes
    for prefix in "HRP_A_020" "HRP_A_021"; do
        # Extract the correct sample_id matching the prefix and barcode from metadata
        SAMPLE_ID=$(awk -F'\t' -v pfx="$prefix" -v barcode="$bc" '$3 ~ "^" pfx && $61 == barcode {print $3; exit}' "$METADATA")
        
        if [[ -n "$SAMPLE_ID" ]]; then
            if [[ "$prefix" == "HRP_A_020" ]]; then
                POD5_DIR="/global/cfs/cdirs/m5243/analysis/ont/all_runs/raw/NOVOA_E/ONT_tRNA_native_vs_gentegra/bc_${i}"
            else
                POD5_DIR="/global/cfs/cdirs/m5243/analysis/ont/all_runs/raw/NOVOA_E/ONT_tRNA_native_vs_IVT/bc_${i}"
            fi
            echo "Processing ${SAMPLE_ID} (Barcode: ${bc}) -> ${POD5_DIR}"
            
            # Update POD5_DIR in generate script
            sed -i "s|^POD5_DIR=.*|POD5_DIR=\"${POD5_DIR}\"|" generate_tRNA_inputs.sh
            
            # Generate FAST JSON
            sed -i "s|^SAMPLE_ID=.*|SAMPLE_ID=\"${SAMPLE_ID}_FAST\"|" generate_tRNA_inputs.sh
            sed -i "s|^BASECALLING_TYPE=.*|BASECALLING_TYPE=\"fast\"|" generate_tRNA_inputs.sh
            ./generate_tRNA_inputs.sh
            
            # Generate SUP JSON
            sed -i "s|^SAMPLE_ID=.*|SAMPLE_ID=\"${SAMPLE_ID}_SUP\"|" generate_tRNA_inputs.sh
            sed -i "s|^BASECALLING_TYPE=.*|BASECALLING_TYPE=\"sup\"|" generate_tRNA_inputs.sh
            ./generate_tRNA_inputs.sh
        fi
    done
done

echo "All JSONs generated successfully."
