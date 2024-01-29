#!/bin/bash
# execute file using: path/to/log_make_products.sh GALAXY

source activate /arc/home/thbrown/.conda/envs/gistenv
galaxy=$1

sleep 1s

# use psrecord to log the mem and cpu usage over time
# any value gt than -1 will run scalene line by line profiling
psrecord "/arc/projects/mauve/products/scripts/make_products.sh ${galaxy} 1" \
--log /arc/projects/mauve/resource_usage/${galaxy}_$(date +"%FT%H%M")_psrecord.txt \
--plot /arc/projects/mauve/resource_usage/${galaxy}_$(date +"%FT%H%M")_psrecord.png --include-children --interval 1