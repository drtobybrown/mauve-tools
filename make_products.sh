#!/bin/bash
# execute file using: path/to/make_products.sh GALAXY

# run even if a command fails with set -e
particular_script()
{
    false
}

# exit if a command fails
set -e

source activate /arc/home/thbrown/.conda/envs/gistenv
# mamba activate gistenv

sleep 1s

# get the session id
sessionid=$(env | grep -i 'skaha_sessionid=' | awk -F= '{print $2}')

galaxy=$1

# if log flag passed run scalene
scalene_log=${2:-1}

sleep 1s

# send a mail notification the job has started
python /arc/home/thbrown/bin/sendmail.py --sessionid=$sessionid --jobname=$galaxy --subject='job started'

# gist config file
configFile="/arc/projects/mauve/products/configFiles/${galaxy}/${galaxy}_MAUVE_MasterConfig_v5.4.2_sn40.yaml"

# gist default file
defaultDir="/arc/projects/mauve/products/configFiles/${galaxy}/${galaxy}.defaultDir"

# set the output directory
gistMain=/arc/projects/mauve/bin/src/gistpipeline/gistPipeline/MainPipeline.py
programPath=/arc/home/thbrown/mauve/dev/gist-geckos/gistPipeline
logPath=/arc/projects/mauve/resource_usage

# run GIST, recording memory usage if flag set
if [ "$2" -gt "-1" ]
then
    echo "running gist with scalene"
    sleep 2s
        
    scalene --cli --reduced-profile --profile-all --profile-only ppxf vorbin \
    --outfile $logPath/${galaxy}_$(date +"%FT%H%M")_scalene.txt \
    $gistMain --config=$configFile --default-dir=$defaultDir
    
else
    gistPipeline --config=$configFile --default-dir=$defaultDir
fi

sleep 3s 

# convert the output fits cubes to hdf5
# which command requires that fits2idia exists before executing the bash script
fits2idia_script="/arc/projects/mauve/products/scripts/fits2hdf5.py"
which fits2idia && python $fits2idia_script --config=$configFile --default-dir=$defaultDir

# make the moment maps
python /arc/projects/mauve/products/scripts/make_moments.py --galaxy=$galaxy --config=$configFile --default-dir=$defaultDir

sleep 1s

# move the Gist generated pipeline files to a subfolder
source $defaultDir
resultsDir="${outputDir%/}/${galaxy}"
pipleineFileDir="${resultsDir}/pipeline_files"
[ -d $pipleineFileDir ] || mkdir $pipleineFileDir # make dir only if it doesn't exist

for file in $(cat /arc/projects/mauve/products/scripts/pipeline_files.txt); do mv $resultsDir/*"$file" $pipleineFileDir; done

# send a mail notification
particular_script || python /arc/home/thbrown/bin/sendmail.py --sessionid=$sessionid --jobname=$galaxy --subject='job complete'

