#!/bin/bash -f

#$ -N Klimakverna 
#$ -l h_rt=24:00:00
#$ -l h_rss=3G,mem_free=3G
#$ -S /bin/bash
#$ -q all.q
#$ -M kajalhb@met.no
#$ -m bae
#$ -o /lustre/storeC-ext/users/klimakverna/development/jobs/OUT_$JOB_NAME.$JOB_ID
#$ -e /lustre/storeC-ext/users/klimakverna/development/jobs/ERR_$JOB_NAME.$JOB_ID
#$ -R y
echo "*** Starting ***"
source /modules/centos7/conda/prod_04_2021/etc/profile.d/conda.sh
conda activate /lustre/storeC-ext/users/klimakverna/development/conda/klimakverna
echo "*** Running snakemake ***"
snakemake --cores=3 --rerun-incomplete --snakefile /lustre/storeC-ext/users/klimakverna/development/Klimakverna-Pilot1/workflow/Snakefile --configfile /lustre/storeC-ext/users/klimakverna/development/Klimakverna-Pilot1/config/testing_config.yaml
echo "*** Done ***"