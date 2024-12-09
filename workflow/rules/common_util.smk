import os
import shutil
import pandas as pd
import glob
import json

config =config["testcase_2"]
MMD_OUT_DIR = config["mmd_output"]
NC_OUT_DIR = config["nc_output"]

INPUT_DIR = config["input_json"]
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
