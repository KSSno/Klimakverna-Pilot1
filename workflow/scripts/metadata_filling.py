from netCDF4 import Dataset
import json
from shutil import copyfile
from datetime import date
import fnmatch
import os

modified = False


def update_metadata(
    input_file,
    output_file,
    global_attributes,
    objects_collection,
    file_pattern_list,
    scenario_list
):

    # Copy input file to the output location before editing
    filename = os.path.basename(input_file)  # Discard the extension   
    copyfile(input_file, output_file)

    add_common_global_attribute(output_file, global_attributes)
    for object_type, pattern  in file_pattern_list.items():
        if (fnmatch.fnmatch(input_file, pattern)):
            object = objects_collection.get(object_type)
            
            # Deriving scenario from filename for different type of inputs(pr, tasmax20ge, norheatwave)
            if object:
                parts = filename.split("_")
                index = object["scenario_index"]
                derived_scenario = parts[index] if len(parts) > index else None

                if derived_scenario:
                    scenario = scenario_list.get(derived_scenario)
                    if scenario:
                        print(f"Scenario: {scenario}")
                    else:
                        print(f"Warning: Scenario not found in scenario_list") 

                add_attributes_variables(
                        output_file,
                        object_type,
                        object["attributes"],
                        scenario
                    )

    global modified
    if modified:
        with Dataset(output_file, "r+") as nc:
            today = date.today().isoformat()
            nc.setncattr("date_metadata_modified", today)
            print(f"Updated date {today}")
            modified = False


def add_common_global_attribute(output_file, global_attributes):
    global modified

    # Open the copied file in 'r+' mode to edit in place
    with Dataset(output_file, "r+") as nc:
        # Check and update metadata attributes
        for key, value in global_attributes.items():
            if key not in nc.ncattrs():
                nc.setncattr(key, value) 
                modified=True
                print(f"Updated global metadata {key}:{value}")   
        
#Adding variable attributes for different types of inputs(pr, tasmax20ge,norheatwave)
def add_attributes_variables(output_file,variable_name, variable_attributes, scenario= None):
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
                var = nc.variables[variable_name]
                for key, value in variable_attributes.get("variable", {}).items():
                    if key not in var.ncattrs():
                        var.setncattr(key, value)
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
        objects = snakemake.input[2]
        pattern_list = snakemake.input[3]
        template_variables = snakemake.input[4]

    except NameError:
        # Default values for standalone testing
        input_file = "example.nc"
        output_file = "example_updated.nc"
        attributes = "attributes.json"
        objects="objects.json"
        pattern_list= "pattern_list.json"
        
    # Load attributes from JSON
    with open(attributes, "r") as f:
        global_attributes = json.load(f)
    with open(objects, "r") as f:
        objects_collection = json.load(f)
       
    with open(pattern_list, "r") as f:
        data = json.load(f)
    with open(template_variables, "r") as f:
        templates = json.load(f)    
    # Extract the list from the "file_pattern" key
    file_pattern_list = data["file_pattern"]
    scenario_list = templates["scenario"]

    update_metadata(
        input_file,
        output_file,
        global_attributes,
        objects_collection,
        file_pattern_list,
        scenario_list
    )
