
import os
import shutil
import pandas as pd
import glob
import json

testcase_config = config["testcase_7"]
ATTRIBUTES = testcase_config["attributes"]
VARIABLES = testcase_config["variables_pr"]
INDEX = testcase_config["index_tasmax20ge"]
MEAN = testcase_config["30_year_mean"]

UPDATED_NC_OUT = config["updated_nc"]
INPUT_DIR = config["input_nc"]
INPUT_BASE = config["input_base"]
TEMPLATE_VAR = config["template_variables"]
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
        OBJECTS,
        INPUT_DIR,
        TEMPLATE_VAR
    output:
        nc_output=os.path.join(UPDATED_NC_OUT, "{filename}.nc4"),
    script:
        "../scripts/metadata_filling.py"
