from netCDF4 import Dataset
import json
from shutil import copyfile
from datetime import date

def update_metadata(input_file, output_file, global_attributes, variables):
    # Copy input file to the output location before editing
    copyfile(input_file, output_file)
    modified=False

    # Open the copied file in 'r+' mode to edit in place
    with Dataset(output_file, "r+") as nc:
        # Check and update metadata attributes
        for key, value in global_attributes.items():
            if key not in nc.ncattrs():
                nc.setncattr(key, value) 
                modified=True
                print(f"Updated metadata {key}:{value}")
        # if "pr" in nc.variables:
        #     pr_var = nc.variables['pr']
        #     for key, value in variables.items():
        #         if key not in pr_var.ncattrs():
        #             pr_var.setncattr(key, value)
        #             modified=True
        #             print(f"Updated variables {key}:{value}")    
        if modified:
            today = date.today().isoformat()
            nc.setncattr("date_metadata_modified", today)
            print(f"Updated date {today}")

def print_global_attributes(nc_file):
    # Open the NetCDF file in read mode
    with Dataset(nc_file, "r") as nc:
        # Iterate through all global attributes
        print("Global Attributes in the NetCDF file:")
        for attr_name in nc.ncattrs():
            attr_value = getattr(nc, attr_name)
            print(f"{attr_name}: {attr_value}")


if __name__ == "__main__":
    try:
        input_file = snakemake.input[0]
        output_file = snakemake.output[0]
        attributes = snakemake.input[1]
        variables = snakemake.input[2]
    except NameError:
        # Default values for standalone testing
        input_file = "example.nc"
        output_file = "example_updated.nc"
        attributes = "attributes.json"
        variables="variables.json"
    # Load global attributes from JSON
    with open(attributes, "r") as f:
        global_attributes = json.load(f)
    with open(variables, "r") as f:
        variables = json.load(f)
    
    # print_global_attributes(input_file)
    update_metadata(input_file, output_file, global_attributes, variables)
    # print_global_attributes(output_file)
