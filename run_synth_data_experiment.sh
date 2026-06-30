#!/bin/bash

experiment_dir=$1
detailed_replicate=$2
n_reps=$3
mem_per_cpu=${4:-4000}
partition=${5:-cobey}

experiment_id=$(basename "$experiment_dir")

if [[ "$partition" != *"cobey"* ]] && [[ "$partition" != *"cobey-hm"* ]]; then
    time_limit="5:00:00"
else
    time_limit="48:00:00"
fi

sbatch --job-name=${experiment_id}_array --output=out_err_files/${experiment_id}_rep%A_%a.out --error=out_err_files/${experiment_id}_rep%A_%a.err --account=pi-cobey --qos=cobey --time=${time_limit} --partition=${partition} --nodes=1 --ntasks-per-node=4 --mem-per-cpu=${mem_per_cpu} --array=1-${n_reps} << EOF
#!/bin/bash
module load R/4.4.1

Rscript fit_bayesian_model_synthetic_data.R $experiment_dir $detailed_replicate
EOF