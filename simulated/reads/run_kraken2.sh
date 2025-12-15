#!/bin/bash

# Require a database path
if [ -z "$1" ]; then
    echo "Usage: $0 <kraken2_database_path>"
    exit 1
fi

db="$1"

# Loop over all prefixes detected from *.R1.fq
for r1 in *.R1.fq; do
    prefix="${r1%.R1.fq}"
    r2="${prefix}.R2.fq"

    # Skip if R2 file missing
    if [ ! -f "$r2" ]; then
        echo "Warning: Missing pair for $r1, skipping..."
        continue
    fi

    echo "Running Kraken2 on prefix: $prefix"

    kraken2 \
        --db "$db" \
        --output "${prefix}.output" \
        --report "${prefix}.kreport" \
        --paired "$r1" "$r2"
done
