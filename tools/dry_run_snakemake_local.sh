#!/bin/bash
# source /modules/rhel8/conda/install/etc/profile.d/conda.sh
# conda activate /lustre/storeC-ext/users/klimakverna/development/conda/klimakverna
snakemake -n --snakefile /lustre/storeC-ext/users/klimakverna/development/Klimakverna-Pilot1/workflow/Snakefile --configfile /lustre/storeC-ext/users/klimakverna/development/Klimakverna-Pilot1/config/testcase_8_config.yaml
