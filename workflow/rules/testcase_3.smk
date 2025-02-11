import os
import shutil
import pandas as pd
import glob
import json

config = config["testcase_3"]
home = '../'

# Extract the list from the "file_pattern" key
INPUT_BASE = config["input_base"]
OUTPUT_DIR = home + config["output_dir"]
FILE_PATTERNS = config["file_pattern"]
CUTOUTS = config["cutouts"]
VAR_RANGES = config["var_ranges"]
VARIABLES = list(config["file_pattern"].keys())

ALL_FILES = []
for var, pattern in FILE_PATTERNS.items():
    for file_path in glob.glob(INPUT_BASE + '/' + pattern):
        infile = os.path.relpath(file_path, INPUT_BASE)
        infile_base = os.path.splitext(infile)[0]  # Discard the extension
        #print(infile_base)
        ALL_FILES.append(infile_base)  # Append the file name without extension

print('INPUT_BASE:', INPUT_BASE)
print('OUTPUT_DIR:', OUTPUT_DIR)
print('VARS:', VARIABLES)
print('RANGES:', VAR_RANGES)
print('CUTOUTS:', CUTOUTS)
print('IN-FILES:', len(ALL_FILES))

rule variable_ranges:
    input:
        "%s/{sample}.nc" % INPUT_BASE
    output:
        "%s/{sample}.txt" % OUTPUT_DIR
    params:
        config=config,
        variables=VARIABLES,
        var_ranges=VAR_RANGES,
        cutouts=CUTOUTS
    script:
        "../scripts/variable_ranges.py"
    #shell:
        #    """
        #    echo {input} > {output}
        #    """