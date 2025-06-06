import logging
import subprocess
from pathlib import Path
from subprocess import CalledProcessError
from typing import Dict

from .config import MAX_FIELD_NAME, MIN_FIELD_NAME

logger = logging.getLogger(__name__)


def subprocess_runner(command: list):
    """
    A simple wrapper around subprocess.run to execute a command and handle errors.
    """
    logger.debug(f"Running command: {' '.join(command)}")

    try:
        result = subprocess.run(command, check=True, capture_output=True)
        logger.debug(f"Command succeeded with output: {result.stdout.decode().strip()}")
    except CalledProcessError as e:
        logger.error(f"Command failed with error: {e.stderr.decode()}")
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred: {str(e)}")
        raise
    logger.debug(f"Command completed with return code: {result.returncode}")

    return result.returncode


def gdal_calc(
    input_rasters: Dict[str, Path],
    output_raster: Path,
    expression: str,
    nodata_value: float | None = None,
    overwrite: bool = True,
    quiet: bool = False,
):
    command = ["gdal_calc"]

    if nodata_value is not None:
        command += [f"--NoDataValue={str(nodata_value)}"]
    if quiet:
        command += ["--quiet"]
    if overwrite:
        command += ["--overwrite"]

    # Add input rasters as named variables
    for name, raster in input_rasters.items():
        command += [f"-{name}", f"{str(raster)}"]

    # Add output raster
    command += [f"--outfile={str(output_raster)}"]

    # Add expression
    command += [f'--calc="{expression}"']

    return subprocess_runner(command)


def gdaldem_TPI(
    input_raster: Path,
    output_raster: Path,
    band: int | None = None,
    output_format: str | None = None,
    compute_edges: bool = False,
    quiet: bool = False,
):
    command = ["gdaldem", "tpi"]

    command.append(str(input_raster))
    command.append(str(output_raster))

    if band:
        command += ["-b", str(band)]
    if output_format:
        command += ["-of", output_format]
    if compute_edges:
        command += ["-compute_edges"]
    if quiet:
        command += ["-q"]

    return subprocess_runner(command)


def gdal_translate(
    input_raster: Path,
    output_raster: Path,
    intput_format: str = None,
    output_format: str = None,
    bands: list = None,
    resampling: str = None,
    outsize: list = None,
    quiet: bool = False,
    srs: str = None,
):
    command = ["gdal_translate"]

    if bands:
        for band in bands:
            command += ["-b", str(band)]

    if srs:
        command += ["-a_srs", srs]

    if resampling:
        command += ["-r", resampling]
    if intput_format:
        command += ["-if", intput_format]
    if output_format:
        command += ["-of", output_format]

    if outsize:
        command += ["-outsize", str(outsize[0]), str(outsize[1])]

    if quiet:
        command += ["-q"]

    command.append(str(input_raster))
    command.append(str(output_raster))

    return subprocess_runner(command)


def gdal_contour(
    input_raster: Path,
    output_vector: Path,
    output_format: str = None,
    output_layer: str = None,
    contour_interval: float = None,
    contour_offset: float = None,
    contour_base: float = None,
    contour_levels: list = None,
    input_band: int = None,
    polygonize: bool = False,
    elevation_field_name: str = "value",
    elevation_max_field_name: str = MAX_FIELD_NAME,
    elevation_min_field_name: str = MIN_FIELD_NAME,
    quiet: bool = False,
):
    command = ["gdal_contour"]

    if input_band:
        command += ["-b", str(input_band)]

    if polygonize:
        command += [
            "-p",
            "-amax",
            elevation_max_field_name,
            "-amin",
            elevation_min_field_name,
        ]
    else:
        command += ["-a", elevation_field_name]

    # Contour options
    if contour_interval:
        command += ["-i", str(contour_interval)]
    if contour_offset:
        command += ["-off", str(contour_offset)]
    if contour_base:
        command += ["-e", str(contour_base)]
    if contour_levels:
        command += ["-fl", " ".join(map(str, contour_levels))]

    if output_format:
        command += ["-f", output_format]
    if output_layer:
        command += ["-nln", output_layer]

    if quiet:
        command += ["-q"]

    # Input and output files
    command.append(str(input_raster))
    command.append(str(output_vector))

    return subprocess_runner(command)


def ogr2ogr(
    src_dataset: Path,
    dst_dataset: Path,
    layer_name: str = None,
    input_format: str = None,
    output_format: str = None,
    overwrite: bool = True,
    clip_ds: Path = None,
    clip_layer: str = None,
    make_valid: bool = False,
    explode: bool = False,
):
    command = ["ogr2ogr"]

    if overwrite:
        command += ["-overwrite"]

    if make_valid:
        command += ["-makevalid"]

    if explode:
        command += ["-explodecollections"]

    if clip_ds:
        command += ["-clipsrc", str(clip_ds)]
    if clip_layer:
        command += ["-clipsrclayer", clip_layer]

    if input_format:
        command += ["-if", input_format]
    if output_format:
        command += ["-of", output_format]

    command.append(str(dst_dataset))
    command.append(str(src_dataset))

    if layer_name:
        command.append(layer_name)

    return subprocess_runner(command)


def gdal_fillnodata(
    input_raster: Path,
    output_raster: Path,
    max_distance: int = None,
    smooth_iterations: int = None,
):
    command = ["gdal_fillnodata"]

    if max_distance:
        command += ["-md", str(max_distance)]
    if smooth_iterations:
        command += ["-si", str(smooth_iterations)]

    command += [str(input_raster)]
    command += [str(output_raster)]

    return subprocess_runner(command)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)

    gdal_calc(
        input_rasters={"A": "test.tif"},
        output_raster="output2.tif",
        expression="A",
        nodata_value=0,
    )
