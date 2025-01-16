from netCDF4 import Dataset

# Open the NetCDF file in read mode
file_path = "/home/shamlym/workspace/klima-kverna/nc/cnrm-r1i1p1f2-hclim_ssp370_3dbc-eqm-sn2018v2005_rawbc_norway_1km_tasmax20ge_annual_merged.nc"
with Dataset(file_path, "r") as nc:
    print("==== Global Attributes ====")
    for attr in nc.ncattrs():
        print(f"{attr}: {getattr(nc, attr)}")
    
    print("\n==== Dimensions ====")
    for dim_name, dimension in nc.dimensions.items():
        print(f"Dimension: {dim_name}, Size: {len(dimension)}")
    
    print("\n==== Variables ====")
    for var_name, variable in nc.variables.items():
        print(f"Variable: {var_name}")
        print(f"  Dimensions: {variable.dimensions}")
        print(f"  Shape: {variable.shape}")
        print(f"  Attributes:")
        for attr in variable.ncattrs():
            print(f"    {attr}: {getattr(variable, attr)}")
        print(f"  Data: {variable[:]}")  # Be careful with large data!