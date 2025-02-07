import netCDF4 
from pydantic import BaseModel, Field, ValidationError , model_validator
from dateutil import parser



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



class NetCDFVariable(BaseModel):
    name: str = Field(..., description="Variable name")
    dimensions: list[str] = Field(..., description="Variable dimensions")
    dtype: str = Field(..., description="Variable data type")
class NetCDFValidation(BaseModel):
    global_attrs: NetCDFMetadata
    # variables: list[NetCDFVariable]

@model_validator(mode="after")
def check_datetime_iso(self) -> "NetCDFMetadata":
    try:
        dt = parser.isoparse(self.time_coverage_end)
    except ValueError:
        raise ValueError(f"{self.time_coverage_end} not in ISO format(YYYY-MM-DDTHH:MM:SSZ)")
    except Exception as e:
        raise e


def validate_nc_file(file_path):

    dataset = netCDF4.Dataset(file_path, mode='r')

    # Extract all global attributes dynamically
    global_attrs = {attr: getattr(dataset, attr) for attr in dataset.ncattrs()}

    dataset.close()

    # Validate using Pydantic
    try:
        validated_data = NetCDFValidation(
            global_attrs=NetCDFMetadata(**global_attrs)
        )
        print("Validation successful:", validated_data)
    except ValidationError as e:
        print("Validation failed:", e)

        
        # print("\n==== Variables ====")
        # variables = [
        #     {
        #         "name": var_name,
        #         "dimensions": list(var.dimensions),
        #         "dtype": str(var.dtype),
        #         "shape": var.shape,  # Optional: Include shape
        #         "attributes": {attr: getattr(var, attr) for attr in var.ncattrs()}  # Extract variable-specific attributes
        #     }
        #     for var_name, var in nc.variables.items()
        #     ]
    
        # try:
        #     validated_data = NetCDFValidation(
        #         global_attrs=NetCDFMetadata(**global_attrs),
        #         variables=[NetCDFVariable(**var) for var in variables]
        #     )
        #     print("Validation successful:", validated_data)
        # except ValidationError as e:
        #     print("Validation failed:", e)

# Usage
validate_nc_file("/home/shamlym/workspace/klima-kverna/nc/cnrm-r1i1p1f2-hclim_ssp370_3dbc-eqm-sn2018v2005_rawbc_norway_1km_tasmax20ge_annual_merged.nc")
