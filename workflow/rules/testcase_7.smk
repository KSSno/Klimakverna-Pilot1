
import os
import shutil
import pandas as pd
import glob
import json

config =config["testcase_7"]
ATTRIBUTES = config["attributes"]
VARIABLES = config["variables_pr"]
INDEX = config["index_tasmax20ge"]
MEAN = config["30_year_mean"]

UPDATED_NC_OUT = config["updated_nc"]
INPUT_DIR = config["input_nc"]
INPUT_BASE = config["input_base"]
with open(INPUT_DIR, "r") as f:
    data = json.load(f)

# Extract the list from the "file_pattern" key
file_pattern_list = data["file_pattern"]
all_file_paths = []

for pattern_type,pattern  in file_pattern_list.items():
    matching_files = glob.glob(pattern, recursive=True)

    for file_path in matching_files:
        struct = os.path.relpath(file_path, INPUT_BASE)
        file_base, _ = os.path.splitext(struct)  # Discard the extension
        all_file_paths.append(file_base)  # Append the file name without extension
        

rule metadata_filling:
    input:
        lambda wildcards: glob.glob(f"{INPUT_BASE}{wildcards.filename}.nc4") or glob.glob(f"{INPUT_BASE}{wildcards.filename}.nc"),
        ATTRIBUTES,
        VARIABLES,
        INDEX,
        MEAN,
        INPUT_DIR
    output:
        nc_output=os.path.join(UPDATED_NC_OUT, "{filename}.nc4"),
    script:
        "../scripts/metadata_filling.py"
