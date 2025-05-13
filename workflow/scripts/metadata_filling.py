from netCDF4 import Dataset
import json
from shutil import copyfile
from datetime import date
import fnmatch
import os

modified = False


def update_metadata(input_file, output_file, global_attributes, attributes_variable_pr,attributes_variable_tasmax20ge,variable_diff_norheatwave, file_pattern_list):
    # Copy input file to the output location before editing
    print(input_file)
    filename, _ = os.path.splitext(input_file)  # Discard the extension
    print(filename)
    parts = filename.split("_")
    if len(parts)>= 4:
        scenario = parts[3]

    copyfile(input_file, output_file)
    add_common_global_attribute( output_file, global_attributes)
    for  type, pattern  in file_pattern_list.items():
        if (fnmatch.fnmatch(input_file, pattern)):
            if type=="pr" :
                add_attributes_variables(output_file, type, attributes_variable_pr,scenario)
    
            if type=="tasmax20ge":
                add_attributes_variables(output_file,type , attributes_variable_tasmax20ge)

            if type=="norheatwave":
                add_attributes_variables(output_file,type , variable_diff_norheatwave)  
    global modified
    if modified:
        with Dataset(output_file, "r+") as nc:
            today = date.today().isoformat()
            nc.setncattr("date_metadata_modified", today)
            print(f"Updated date {today}")
            modified = False


def add_common_global_attribute(output_file,  global_attributes):
    global modified 

    # Open the copied file in 'r+' mode to edit in place
    with Dataset(output_file, "r+") as nc:
        # Check and update metadata attributes
        for key, value in global_attributes.items():
            if key not in nc.ncattrs():
                nc.setncattr(key, value) 
                modified=True
                print(f"Updated global metadata {key}:{value}")   
        

def add_attributes_variables(output_file,variable_name, variable_attributes,scenario):
    global modified 

    with Dataset(output_file, "r+") as nc:
    # Check and update metadata attributes
        for key, value in variable_attributes.get("global", {}).items():
            if key not in nc.ncattrs():
                if '{scenario}' in value:
                    value = value.format(scenario=scenario)
                nc.setncattr(key, value) 
                modified=True
                print(f"Updated global attribute {variable_name} {key}:{value}")
            if variable_name in nc.variables:
                pr_var = nc.variables[variable_name]
                for key, value in variable_attributes.get(variable_name, {}).items():
                    if key not in pr_var.ncattrs():
                        pr_var.setncattr(key, value)
                        modified=True
                        print(f"Updated variables {variable_name} {key}:{value}")    


# def print_global_attributes(nc_file):
#     # Open the NetCDF file in read mode
#     with Dataset(nc_file, "r") as nc:
#         # Iterate through all global attributes
#         print("Global Attributes in the NetCDF file:")
#         for attr_name in nc.ncattrs():
#             attr_value = getattr(nc, attr_name)
#             print(f"{attr_name}: {attr_value}")


if __name__ == "__main__":
    try:
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
        attributes = snakemake.input[1]
        variable_pr = snakemake.input[2]
        variable_tasmax20ge = snakemake.input[3]
        variable_diff_norheatwave = snakemake.input[4]
        pattern_list = snakemake.input[5]


       
    except NameError:
        # Default values for standalone testing
        input_file = "example.nc"
        output_file = "example_updated.nc"
        attributes = "attributes.json"
        variable_pr="variables.json"
        variable_diff_norheatwave="variable_diff_norheatwave.json"
        pattern_list= "pattern_list.json"
        
    # Load attributes from JSON
    with open(attributes, "r") as f:
        global_attributes = json.load(f)
    with open(variable_pr, "r") as f:
        attributes_variable_pr = json.load(f)
    with open(variable_tasmax20ge, "r") as f:
        attributes_variable_tasmax20ge = json.load(f)
    with open(variable_diff_norheatwave, "r") as f:
        attributes_variable_diff_norheatwave = json.load(f)    
    with open(pattern_list, "r") as f:
        data = json.load(f)
        
    # Extract the list from the "file_pattern" key
    file_pattern_list = data["file_pattern"]
    
    update_metadata(input_file, output_file, global_attributes, attributes_variable_pr,attributes_variable_tasmax20ge,attributes_variable_diff_norheatwave, file_pattern_list)

