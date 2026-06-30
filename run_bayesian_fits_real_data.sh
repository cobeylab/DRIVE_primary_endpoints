#!/bin/bash

fit_specs_dir=$1
mem_per_cpu=${2:-4000}
partition=${3:-cobey}

fit_id=$(basename "$fit_specs_dir")

if [[ "$partition" != *"cobey"* ]] && [[ "$partition" != *"cobey-hm"* ]]; then
    time_limit="5:00:00"
else
    time_limit="48:00:00"
fi

sbatch --job-name=${fit_id}_array --output=out_err_files/${fit_id}_rep%A_%a.out --error=out_err_files/${fit_id}_rep%A_%a.err --account=pi-cobey --qos=cobey --time=${time_limit} --partition=${partition} --nodes=1 --ntasks-per-node=4 --mem-per-cpu=${mem_per_cpu} << EOF
#!/bin/bash
module load R/4.4.1

Rscript fit_bayesian_model_real_data.R $fit_specs_dir
EOF