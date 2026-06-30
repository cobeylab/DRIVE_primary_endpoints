#!/bin/bash

fit_specs_dir=$1 # Specification directory (which model is being fitted)
loglik_unit=$2 # LOO unit (1 = individual/strain/year, 2 = individual)
mem_per_cpu=${3:-8000}
partition=${4:-cobey}

fit_id=$(basename "$fit_specs_dir")_LOO

if [[ "$partition" != *"cobey"* ]] && [[ "$partition" != *"cobey-hm"* ]]; then
  time_limit="5:00:00"
else
    time_limit="48:00:00"
fi

# Determine number of LOO units (array size) from processed_data.csv using cut/awk
csv="bayesian_fits_real_data/processed_data.csv"

# Helper to get column index by name from CSV header
get_col_idx() {
  awk -F',' -v col="$1" 'NR==1{
    for(i=1;i<=NF;i++){
      gsub(/^[ \t]+|[ \t]+$/, "", $i)
      if($i==col){print i; exit}
    }
  }' "$csv"
}

ind_idx=$(get_col_idx "individual")
str_idx=$(get_col_idx "strain")
year_idx=$(get_col_idx "year")

# Automatically determine the number of LOO units based on processed_data.csv
# to create a job array.
# This will be greater than the actual number because of downstream processing,
#  (e.g. if looking only at day 30), but the Rscript will handle array indices 
# in excess of actual units
if [ "$loglik_unit" -eq 1 ]; then
  # Unique combinations of individual,strain,year
  num_indices=$(cut $csv -d "," -f $ind_idx,$str_idx,$year_idx | tail -n+2 | sort -u | wc -l)
else
  # Unique individuals
  num_indices=$(cut $csv -d "," -f $ind_idx | tail -n+2 | sort -u | wc -l)
fi

if [ "$loglik_unit" -eq 1 ]; then
  LOO_dir="LOO_Hpre_group"
else
  LOO_dir="LOO_individual"
fi

# Determine indices to run (exclude already completed ones and those with .tmp files)
LOO_dir_path="${fit_specs_dir}/${LOO_dir}"
combined_csv="${LOO_dir_path}/combined_results.csv"

if [ -f "$combined_csv" ]; then
  # Get already run indices (assume first column)
  already_run=$(cut -d, -f1 "$combined_csv" | tail -n+2 | sort -u)
  # Generate all indices
  all_indices=$(seq 1 $num_indices)
  # Exclude already run indices (ensure both are sorted numerically)
  indices_to_run=$(comm -23 <(echo "$all_indices" | sort -n) <(echo "$already_run" | sort -n))
else
  indices_to_run=$(seq 1 $num_indices)
fi

# Exclude indices for which holdout_[index].csv.tmp or holdout_[index].csv exists in LOO_dir
indices_to_run=$(for idx in $indices_to_run; do
  if [ ! -f "${LOO_dir_path}/holdout_${idx}.csv.tmp" ] && [ ! -f "${LOO_dir_path}/holdout_${idx}.csv" ]; then
    echo $idx
  fi
done)

# If partition is neither cobey nor cobey-hm, keep only the first 1000 indices
if [[ "$partition" != *"cobey"* ]] && [[ "$partition" != *"cobey-hm"* ]]; then
  indices_to_run=$(echo "$indices_to_run" | head -n 1000)
fi

if [ -z "$indices_to_run" ]; then
  echo "No indices to run."
  exit 0
fi

# Generate a unique suffix for this run (timestamp)
run_suffix=$(date +%Y%m%d_%H%M%S)

# Write indices to a uniquely named file in out_err_files/
indices_file="out_err_files/${fit_id}_indices_to_run_${run_suffix}.txt"
echo "$indices_to_run" > "$indices_file"
num_indices_to_run=$(wc -l < "$indices_file")

# Convert to SLURM array range
array_range="1-${num_indices_to_run}"

# Decide which Rscript command to use
if [[ "$fit_specs_dir" == *"null_model"* ]]; then
  # If using the null model, a simplified R script takes only the holdout index as an argument.
  echo "Warning: Fitting to the null model. loglik_unit is set to 1 by default, regardless of the provided argument."
  rscript_cmd='Rscript fit_null_model_real_data.R "$holdout_index"'
else
  # If using the main class of models
  rscript_cmd='Rscript fit_bayesian_model_real_data.R "'"$fit_specs_dir"'" "'"$loglik_unit"'" "$holdout_index"'
fi

# Set qos based on partition
if [[ "$partition" == *"cobey"* ]] || [[ "$partition" == *"cobey-hm"* ]]; then
  qos="cobey"
else
  qos="$partition"
fi

sbatch --job-name=${fit_id}_LOO --output=out_err_files/${fit_id}_LOO_rep%A_%a.out --error=out_err_files/${fit_id}_LOO_rep%A_%a.err --account=pi-cobey --qos=${qos} --time=${time_limit} --partition=${partition} --nodes=1 --ntasks-per-node=4 --mem-per-cpu=${mem_per_cpu} --array=${array_range} << EOF
#!/bin/bash
module load R/4.4.1

# Read the actual holdout index from the uniquely named file in out_err_files/
indices_file="${indices_file}"
holdout_index=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "\$indices_file")

eval $rscript_cmd
EOF