import arrow
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
from osgeo import gdal
import rasterio
from shapely.geometry import mapping
# import numpy as np


def save_areal_mean_for_region_to_netcdf(
    config: dict, indicator_id: str, scenario: str, period: str, region: str, inputs: list[str], output: str
):
    """
    Saves mean and percentiles of indicator for one scenario, period and region as netcdf.

    Args:
        config: configuration from KAPy that contains information from all configuration files
        indicator_id: which parameter in the input file to calculate ensemble statistics for
        scenario: ssp370 for CMIP6 and rcp26 and rcp45 for CMIP5 inputs
        period: historical, near future or far future
        region: id for region in shapefile defined in configuration file to cut out of input file and calculate ensemble statistics for
        inputs: original datasets with 30 year means
        output: path to output directory including netcdf file name

    Returns:
        Nothing, but saves the areal mean and percentiles to netcdf for the input scenario, period and region, as well as a spatial plot
    """
    indicator_name = config["indicators"][indicator_id]["variables"]
    region_datasets = {}
    polygons_file = gpd.read_file(config["region"][region]["shapefile"])
    this_region = polygons_file[polygons_file.Region == int(region)].geometry
    polygon = gpd.GeoSeries(data=this_region, crs=polygons_file.crs)

    for input_file in inputs:
        dataset = xr.open_dataset(input_file)
        projection = gdal.Open(input_file).GetProjection()
        dataset = dataset.rio.write_crs(projection)
        polygon = polygon.to_crs(projection)

        mask = rasterio.features.geometry_mask(
            [mapping(polygon.geometry.iloc[0])],
            out_shape=(len(dataset.Yc), len(dataset.Xc)),
            transform=dataset.rio.transform(),
            invert=True,
        )
        # all_touched=False by default, means "If False, only pixels whose center is within the polygon or that are selected by Bresenham’s line algorithm will be burned in"
        mask = xr.DataArray(mask, dims=("Yc", "Xc"))
        clipped_ds = dataset.where(mask, drop=True)
        region_datasets[input_file] = clipped_ds

    if period == "hist":
        period_name = "ref_period"
    else:
        period_name = period
    datasets = [region_datasets[filename] for filename in region_datasets if period_name in filename]

    if not datasets:  # doesn't work for daily datasets, have to extract time period
        for idx in range(1, len(config["periods"]) + 1):
            if period == config["periods"][f"{idx}"]["short_name"]:
                period_config = config["periods"][f"{idx}"]
        year = arrow.get(str(period_config["start"]))
        end = arrow.get(str(period_config["end"]))

        while year <= end:
            datasets.extend(
                [
                    region_datasets[filename]
                    for filename in region_datasets
                    if f"daily_{year.format('YYYY')}" in filename
                ]
            )
            year = year.shift(years=1)

        dataset_total = xr.combine_nested(datasets, concat_dim="time")
        dataarray_areal_mean = dataset_total.mean(dim=["Yc", "Xc"])
        dataarray_areal_mean = dataarray_areal_mean.assign({"scenario": [scenario]})
        dataarray_areal_mean = dataarray_areal_mean.drop("projection_utm")
        dataarray_areal_mean.to_netcdf(output)

    else:
        dataset_total = xr.combine_nested(datasets, concat_dim="realization")
        # Take mean over all bias adjustments
        dataset_mean = dataset_total[indicator_name].mean(dim=["realization"])

        # To save region dataset for ncview
        dataset_mean.to_netcdf(
            f"{output.split(indicator_id)[0]}/{indicator_id}_{scenario}_{period}_region_{region}_area.nc"
        )

        # Take areal mean
        dataarray_areal_mean = dataset_mean.mean(dim=["Yc", "Xc"])
        dataarray_areal_mean = dataarray_areal_mean.rename("indicator_mean")
        dataarray_areal_mean = dataarray_areal_mean.rename({"time": "periodID"})
        dataarray_areal_mean = dataarray_areal_mean.to_dataset()
        dataarray_areal_mean = dataarray_areal_mean.assign({"scenario": [scenario]})
        dataarray_areal_mean = dataarray_areal_mean.drop("projection_utm")
        dataarray_areal_mean.to_netcdf(output)

        limit = [-15, 15]
        plt.figure()
        dataset_mean.plot(robust=True, vmin=limit[0], vmax=limit[-1])
        plt.savefig(f"{output.split(indicator_id)[0]}/{indicator_id}_{scenario}_{period}_region_{region}_pr_mean.png")

    # plt.figure()
    # dataset_total.sel(Xc=list(np.arange(195500, 205500, 1000))).sel(Yc=list(np.arange(6889500, 6840500, -1000))).plot.scatter(x="Xc", y="Yc", hue=indicator_name, cmap="viridis")
    # plt.savefig(f"{output.split(indicator_id)[0]}/{indicator_id}_{scenario}_{period}_region_{region}_test.png")
