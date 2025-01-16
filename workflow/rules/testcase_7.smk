
import os
import shutil
import pandas as pd
import glob
import json

config =config["testcase_7"]
ATTRIBUTES = config["attributes"]
VARIABLES = config["variables_pr"]
UPDATED_NC_OUT = config["updated_nc"]
INPUT_DIR = config["input_nc"]
INPUT_BASE = config["input_base"]
with open(INPUT_DIR, "r") as f:
    data = json.load(f)

# Extract the list from the "file_pattern" key
file_pattern_list = data["file_pattern"]

all_file_paths = []

for pattern in file_pattern_list:
    matching_files = glob.glob(pattern, recursive=True)
    for file_path in matching_files:
        struct = os.path.relpath(file_path, INPUT_BASE)
        file_base, _ = os.path.splitext(struct)  # Discard the extension
        all_file_paths.append(file_base)  # Append the file name without extension
        
print(all_file_paths)

rule metadata_filling:
    input:
        f"{INPUT_BASE}{{filename}}.nc4",
        ATTRIBUTES,
        VARIABLES

    output:
        nc_output=os.path.join(UPDATED_NC_OUT, "{filename}.nc4"),
    script:
        "../scripts/metadata_filling.py"
