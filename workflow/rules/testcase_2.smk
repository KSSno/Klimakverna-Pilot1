import os
import shutil
import pandas as pd
import glob
import json

config =config["testcase_2"]

NC_OUT_DIR = config["nc_output"]
INPUT_DIR = config["input_json"]
INPUT_BASE = config["input_base"]

with open(INPUT_DIR, "r") as f:
    data = json.load(f)

# Extract the list from the "file_pattern" key
# Extract the list from the "file_pattern" key
file_pattern_list = data["file_pattern"]

all_file_paths = []

for pattern in file_pattern_list:
    matching_files = glob.glob(pattern, recursive=True)

    for file_path in matching_files:
        struct = os.path.relpath(file_path, INPUT_BASE)
        file_base, _ = os.path.splitext(struct)  # Discard the extension
        all_file_paths.append(file_base)  # Append the file name without extension
rule validate_and_move:
    input:
        f"{INPUT_BASE}{{filename}}.nc4"
    output:
        nc_output=os.path.join(NC_OUT_DIR, "{filename}.nc4")
    
    script:
        "../scripts/validate.py"