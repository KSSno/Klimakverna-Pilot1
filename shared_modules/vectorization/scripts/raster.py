import logging
import shutil
from pathlib import Path

import numpy as np
import rasterio
from scipy import ndimage

from .models import Smoothing

logger = logging.getLogger(__name__)


def get_raster_stats(dataset: Path, band: int = 1) -> dict:
    with rasterio.open(dataset, "r") as src:
        # Read the data
        data = src.read(band)

        # Get the stats
        stats = {
            "min": data.min(),
            "max": data.max(),
            "mean": data.mean(),
            "std": data.std(),
            "nodata": src.nodata,
        }

        return stats


def add_nodata_values(
    in_ds: Path,
    out_ds: Path,
    nodata: list[int | float] = [],
    band: int = 1,
):
    """
    Add nodata values to a raster dataset.
    """

    logger.debug(f"Adding nodata values to {in_ds}")
    logger.debug(f"Output dataset: {out_ds}")
    logger.debug(f"Nodata values: {nodata}")
    logger.debug(f"Band: {band}")

    # Copy the input dataset to the output dataset
    shutil.copy(in_ds, out_ds)

    if len(nodata) == 0:
        logger.warning(
            f"No nodata values provided. No changes will be made to {out_ds}"
        )
        return out_ds

    with rasterio.open(out_ds, "r+") as src:
        if src.nodata is None:
            logger.warning(
                f"Dataset {in_ds} does not have a nodata value. Setting nodata to {nodata[0]}"
            )
            src.nodata = nodata[0]

        # Read the data
        data = src.read(band)

        # Create a new mask with nodata values
        nodata_mask = np.isin(data, nodata)

        # Replace values in data with nodata values whren mask is True
        data[nodata_mask] = src.nodata
        # Write the modified data back to the dataset
        src.write(data, band)

        return out_ds


def uniform_filter(
    in_ds: Path,
    out_ds: Path,
    size: int = 3,
    band: int = 1,
):
    """
    Apply an average filter to a raster dataset.
    """

    logger.debug(f"Applying uniform filter to {in_ds}")
    logger.debug(f"Output dataset: {out_ds}")
    logger.debug(f"Band: {band}")
    logger.debug(f"Filter size: {size}")

    with rasterio.open(in_ds) as src:
        data = src.read(band)

        # Apply the average filter
        filtered_data = ndimage.uniform_filter(data, size=size)

        # Write the filtered data to the output dataset
        with rasterio.open(out_ds, "w", **src.meta) as dst:
            dst.write(filtered_data, band)

        return out_ds


def median_filter(
    in_ds: Path,
    out_ds: Path,
    size: int = 3,
    band: int = 1,
):
    """
    Apply a median filter to a raster dataset.
    """

    logger.debug(f"Applying median filter to {in_ds}")
    logger.debug(f"Output dataset: {out_ds}")
    logger.debug(f"Band: {band}")
    logger.debug(f"Filter size: {size}")

    with rasterio.open(in_ds) as src:
        data = src.read(band)

        # Apply the median filter
        filtered_data = ndimage.median_filter(data, size=size)

        # Write the filtered data to the output dataset
        with rasterio.open(out_ds, "w", **src.meta) as dst:
            dst.write(filtered_data, band)

        return out_ds


def gaussian_filter(
    in_ds: Path,
    out_ds: Path,
    sigma: float = 1.0,
    radius: int | None = None,
    band: int = 1,
):
    """
    Apply a Gaussian filter to a raster dataset.
    """

    logger.debug(f"Applying Gaussian filter to {in_ds}")
    logger.debug(f"Output dataset: {out_ds}")
    logger.debug(f"Sigma: {sigma}")
    logger.debug(f"Band: {band}")

    with rasterio.open(in_ds) as src:
        data = src.read(band)

        # Apply the Gaussian filter
        filtered_data = ndimage.gaussian_filter(data, sigma=sigma, radius=radius)

        # Write the filtered data to the output dataset
        with rasterio.open(out_ds, "w", **src.meta) as dst:
            dst.write(filtered_data, band)

        return out_ds


def smooth_raster(in_ds: Path, out_ds: Path, smoothing_config: Smoothing):
    if smoothing_config.method == "uniform":
        logger.info(
            f"Smoothing dataset with uniform filter, size: {smoothing_config.size}"
        )
        uniform_filter(
            in_ds=in_ds,
            out_ds=out_ds,
            size=smoothing_config.size,
        )
    elif smoothing_config.method == "gaussian":
        logger.info(
            f"Smoothing dataset with Gaussian filter, sigma: {smoothing_config.sigma}"
        )
        gaussian_filter(
            in_ds=in_ds,
            out_ds=out_ds,
            sigma=smoothing_config.sigma,
            radius=smoothing_config.radius,
        )
    elif smoothing_config.method == "median":
        logger.info(
            f"Smoothing dataset with median filter, size: {smoothing_config.size}"
        )
        median_filter(
            in_ds=in_ds,
            out_ds=out_ds,
            size=smoothing_config.size,
        )

    return out_ds


if __name__ == "__main__":
    # Example usage of the add_nodata_values function
    in_ds = Path("data/test.tif")
    out_ds = Path("data/test_out.tif")
    nodata = [0]

    add_nodata_values(in_ds, out_ds, nodata)

    # Example usage of the average_filter function
    in_ds = out_ds
    out_ds = Path("data/test_filtered.tif")
    size = 5
    uniform_filter(in_ds, out_ds, size)

    # Example usage of the gaussian_filter function
    in_ds = out_ds
    out_ds = Path("data/test_gaussian.tif")
    sigma = 1.0
    gaussian_filter(in_ds, out_ds, sigma)

    for k, v in get_raster_stats(in_ds).items():
        print(f"{k}: {v}")
