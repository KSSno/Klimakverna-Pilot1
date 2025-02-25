import json
import os
import glob
import xarray as xr
import numpy as np


def select_subset(da, cutouts):
    if cutouts["method"] == "lonlatbox":
        # Select the data within the specified bounds
        lon = slice(cutouts["xmin"], cutouts["xmax"])
        lat = slice(cutouts["ymin"], cutouts["ymax"])
        return da.sel(lon=lon, lat=lat)
    return da


def check_variable_minmax(nc_file, var, var_ranges, cutouts, outfp):
    with xr.open_dataset(nc_file) as da:
        da = select_subset(da, cutouts)
        min_value = var_ranges[f'{var}_min']
        max_value = var_ranges[f'{var}_max']

        print(var, 'minmax:', min_value, max_value)
        # Find indices where 'tas' values are below and above the thresholds
        outside_range = np.argwhere((da[var].values < min_value) | (da[var].values > max_value))

        if len(outside_range) == 0:
            print(f'# {nc_file}: OK: values are within range ({min_value}, {max_value})', file=outfp)
        else:
            print(f'# {nc_file}: WARN: values out of range:', file=outfp)
            for t, x, y in outside_range:
                print('\t%s\t%s\t%g\t%s\t%.2f\t%.2f\t%s' %
                      (var,
                      f'>{max_value}' if da[var].values[t][x][y] > max_value else f'<{min_value}',
                      da[var].values[t][x][y],
                      da['time'].values[t],
                      da['lat'].values[x],
                      da['lon'].values[y],
                      nc_file), file=outfp)
            da.to_netcdf(output_file + '.nc') # Save the area outside


if __name__ == "__main__":
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    variables = snakemake.params.variables
    var_ranges = snakemake.params.var_ranges
    cutouts = snakemake.params.cutouts

    var = os.path.basename(input_file).split('_', 1)[0]
    #print('IN:', input_file)
    #print('OUT', output_file)
    if var in variables:
        with open(output_file, 'w') as outfp:
            #outfp.write(input_file + '\n')
            check_variable_minmax(input_file, var, var_ranges, cutouts, outfp)
