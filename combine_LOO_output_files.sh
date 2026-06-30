#!/bin/bash

LOO_dir=$1

output_file="${LOO_dir}/combined_results.csv"

# Find all output files holdout_*.csv
csv_files=($(find "$LOO_dir" -maxdepth 1 -type f -name "holdout_*.csv" | sort))

if [ ${#csv_files[@]} -eq 0 ]; then
    echo "No holdout*.csv files found in $LOO_dir"
    exit 1
fi

if [ ! -f "$output_file" ]; then
    # Write header from the first file
    head -n 1 "${csv_files[0]}" > "$output_file"
    # Append all CSVs, skipping header lines
    for f in "${csv_files[@]}"; do
        tail -n +2 "$f" >> "$output_file"
    done
else
    # combined_results.csv exists, append only new rows (skip header)
    for f in "${csv_files[@]}"; do
        tail -n +2 "$f" >> "$output_file"
    done
fi

# Remove the old csv files
for f in "${csv_files[@]}"; do
    rm "$f"
done

