import netCDF4 
from pydantic import BaseModel, Field, ValidationError , model_validator
from dateutil import parser
from shutil import copyfile
import os


class NetCDFMetadata(BaseModel):
    title: str = Field(..., description="Title of the dataset")
    institution: str = Field(..., description="Institution providing the dataset")
    bc_method: str = Field(default='Empirical quantile mapping (EQM) or EQM + 3DBC', description="Source of the data")
    bc_method_id:str = Field(default='NVE-EQM or MET-3DBC')
    bc_period:str = Field(default='1985-2014')
    dataset_production_status:str = Field(default='Complete')
    time_coverage_end :str = Field()
    time_coverage_start:str = Field()
    time_coverage_resolution:str = Field()
    history:str = Field()
    date_created:str = Field()
    date_modified:str = Field()
    observation:str =Field()
    observation_id:str = Field()
    product_name :str = Field()


class NetCDFVariablePr(BaseModel):
    FillValue: str = Field()
    long_name: str = Field()
    norwegian_name: str = Field()
    description: str = Field()
    norwegian_description: str = Field()
    units: str = Field()
    coordinates: str = Field()
    grid_mapping: str = Field()
    cell_methods: str = Field()


class NetCDFValidation(BaseModel):
    global_attrs: NetCDFMetadata
    variable_pr: NetCDFVariablePr

@model_validator(mode="after")
def check_datetime_iso(self) -> "NetCDFMetadata":
    try:
        dt = parser.isoparse(self.time_coverage_end)
    except ValueError:
        raise ValueError(f"{self.time_coverage_end} not in ISO format(YYYY-MM-DDTHH:MM:SSZ)")
    except Exception as e:
        raise e


def validate_nc_file(file_path,ouput_dir):

    dataset = netCDF4.Dataset(file_path, mode='r')

    # Extract all global attributes dynamically
    global_attrs = {attr: getattr(dataset, attr) for attr in dataset.ncattrs()}


    # Extract attributes for variable 'pr' 
    var_name = "pr"
    
    if var_name in dataset.variables:
        variable = dataset.variables[var_name]
        pr_attributes = {}
        for attr in variable.ncattrs():
            new_attr_name = attr.lstrip("_")  # Remove leading underscore if present
            pr_attributes[new_attr_name] = getattr(variable, attr)
        # Store attributes in a dictionary
        
    else:
        print(f"Variable '{var_name}' not found in the NetCDF file.")

    dataset.close()

    # Validate using Pydantic
    try:
        validated_data = NetCDFValidation(
            global_attrs=NetCDFMetadata(**global_attrs),
            variable_pr = NetCDFVariablePr(**pr_attributes)
        )
        destination = os.path.join(ouput_dir, input_file)
        copyfile(input_file, ouput_dir)

        print("Validation successful:", validated_data)
    except ValidationError as e:
        print("Validation failed:", e)

if __name__ == "__main__":

    try:
        input_file = snakemake.input[0]
        output_dir = snakemake.output[0]

    except NameError:
            # Default values for standalone testing
            input_file = "example.nc"
            
    # Usage
    validate_nc_file(input_file,output_dir)
