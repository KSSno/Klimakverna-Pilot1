import logging
import shutil
from importlib import resources as imp_resources
from pathlib import Path
from tempfile import TemporaryDirectory

import pyogrio

from . import resources
from .config import DEBUG_FOLDER, DefultContourConfig
from .gdal import (
    gdal_calc,
    gdal_contour,
    gdal_fillnodata,
    gdal_translate,
    gdaldem_TPI,
    ogr2ogr,
)
from .models import ContourGeneratorConfig, DynamicSmoothing
from .raster import (
    add_nodata_values,
    gaussian_filter,
    get_raster_stats,
    smooth_raster,
)
from .vector import (
    color_layer,
    eliminate_small_polygons,
    load_vector_layer,
    replace_with_nodata,
    save_vector_layer,
)

logger = logging.getLogger(__name__)


def get_feature_data():
    with imp_resources.path(resources, "feature_data.gpkg") as path:
        layers = pyogrio.list_layers(path)

        return path, layers


def list_feature_layers():
    path, layers = get_feature_data()
    return [layer[0] for layer in layers]


def get_feature_layer(layer_name: str):
    path, layers = get_feature_data()

    for layer in layers:
        if layer[0].upper() == layer_name.upper():
            return {"dataset": path, "layer": layer[0]}

    # If the layer is not found, raise an error
    err_msg = f"Feature layer {layer_name} not found in package. Available layers: {[layer[0] for layer in layers]}"
    logger.error(err_msg)
    raise ValueError(err_msg)


def get_contour_definition(
    name: str, contour_dict: dict = DefultContourConfig.contour_levels
):
    """
    Get contour levels and colormap from the contour level dictionary.
    """
    # Check if the name exists in the dictionary. Case insensitive.
    for key in contour_dict.keys():
        if key.upper() == name.upper():
            return contour_dict[key]

    err_msg = f"Contour definition '{name}' not found in the contour level dictionary. Available definitions: {list(contour_dict.keys())}"
    # Log the error message
    logger.error(err_msg)
    # Raise a ValueError with the error message
    raise ValueError(err_msg)


def smooth_raster_dynamic(
    in_ds: Path,
    out_ds: Path,
    smoothing_config: DynamicSmoothing,
    tempdir_params: dict = {},
):
    """
    Implemented using this method:
    https://www.tandfonline.com/doi/full/10.1080/23729333.2017.1300998

    Some details from this repository:
    https://github.com/MathiasGroebe/Smooth-Contours

    """
    with TemporaryDirectory(**tempdir_params) as tmpdir:
        # Dynamic smoothing files
        tmp_tpi_ds = Path(tmpdir) / "tmp_tpi.tif"
        tmp_tpi_reclass_ds = Path(tmpdir) / "tmp_tpi_reclass.tif"
        tmp_tpi_smooth_ds = Path(tmpdir) / "tmp_tpi_smooth.tif"
        tmp_tpi_smooth_reclass_ds = Path(tmpdir) / "tmp_tpi_smooth_reclass.tif"
        tmp_tpi_norm_ds = Path(tmpdir) / "tmp_tpi_norm.tif"
        tmp_smooth_low_ds = Path(tmpdir) / "smooth_low_ds.tif"
        tmp_smooth_high_ds = Path(tmpdir) / "smooth_high_ds.tif"

        NODATA_VALUE = -9999

        # Create raster with low smoothing
        smooth_raster(
            in_ds=in_ds,
            out_ds=tmp_smooth_low_ds,
            smoothing_config=smoothing_config.low,
        )
        # Create raster with high smoothing
        smooth_raster(
            in_ds=in_ds,
            out_ds=tmp_smooth_high_ds,
            smoothing_config=smoothing_config.high,
        )
        # Create TPI
        gdaldem_TPI(
            input_raster=in_ds,
            output_raster=tmp_tpi_ds,
        )
        # Reclassify TPI
        gdal_calc(
            input_rasters={"A": tmp_tpi_ds},
            output_raster=tmp_tpi_reclass_ds,
            expression="((-1)*A*(A<0))+(A*(A>=0))",
            nodata_value=NODATA_VALUE,
        )
        # Smooth TPI
        gaussian_filter(
            in_ds=tmp_tpi_reclass_ds,
            out_ds=tmp_tpi_smooth_ds,
            sigma=smoothing_config.low.sigma,
            radius=smoothing_config.low.radius,
        )
        # Reclassify smooth TPI
        gdal_calc(
            input_rasters={"A": tmp_tpi_smooth_ds},
            output_raster=tmp_tpi_smooth_reclass_ds,
            expression="(0*A*(A<0))+(A*(A>=0))",
            nodata_value=NODATA_VALUE,
        )
        # Normalize TPI
        gdal_calc(
            input_rasters={"A": tmp_tpi_smooth_reclass_ds},
            output_raster=tmp_tpi_norm_ds,
            expression=f"A / {get_raster_stats(tmp_tpi_smooth_reclass_ds)['max']}",
            nodata_value=NODATA_VALUE,
        )

        # Calculate the final smoothed dataset
        gdal_calc(
            input_rasters={
                "A": tmp_tpi_norm_ds,
                "B": tmp_smooth_low_ds,
                "C": tmp_smooth_high_ds,
            },
            output_raster=out_ds,
            expression="A*B+(1-A)*C",
            nodata_value=NODATA_VALUE,
        )

    return out_ds


def contour_generator(config: ContourGeneratorConfig, debug: bool = False):
    logger.info(f"Processing: {config.output_layer}")

    # Write model config to debug log
    logger.debug(f"ContourGeneratorConfig: {config.model_dump_json(indent=4)}")

    tempdir_params = {}

    if debug:
        # Create a temporary directory in the debug folder
        tempdir_params = {
            "dir": DEBUG_FOLDER,
            "delete": False,
            "prefix": config.output_layer,
        }

    with TemporaryDirectory(**tempdir_params) as tmpdir:
        # Input dataset
        in_ds = f"{config.format}:{str(config.dataset)}:{config.subdataset}"

        # Temporary files
        tmp_src_ds = Path(tmpdir) / "src_ds.tif"
        tmp_nodata_ds = Path(tmpdir) / "nodata_ds.tif"
        tmp_fill_ds = Path(tmpdir) / "fill_ds.tif"
        tmp_resample_ds = Path(tmpdir) / "resampled.tif"
        tmp_smooth_ds = Path(tmpdir) / "smooth_ds.tif"
        tmp_gpkg_contour = Path(tmpdir) / "contour.gpkg"
        tmp_gpkg_valid = Path(tmpdir) / "contour_valid.gpkg"
        tmp_gpkg_clip = Path(tmpdir) / "clip.gpkg"

        # Extract subdataset from the NetCDF file
        logger.info(f"Extracting subdataset {config.subdataset} from {config.dataset}")
        gdal_translate(
            input_raster=in_ds,
            output_raster=tmp_src_ds,
            intput_format=config.format,
            output_format="GTiff",
            bands=[config.band],
            srs=config.srs,
        )

        # Set nodata value to zero if specified
        if config.add_nodata_values:
            logger.info(
                f"Setting pixel values to nodata values: {config.add_nodata_values}"
            )
            add_nodata_values(
                in_ds=tmp_src_ds,
                out_ds=tmp_nodata_ds,
                nodata=config.add_nodata_values,
            )
        else:
            tmp_nodata_ds = tmp_src_ds

        # Fill nodata values
        logger.info("Filling nodata values")
        gdal_fillnodata(
            input_raster=tmp_nodata_ds,
            output_raster=tmp_fill_ds,
        )

        # Resample the dataset
        logger.info(
            f"Resampling dataset with method: {config.resample_method}, factor: {config.resample_factor}"
        )
        gdal_translate(
            input_raster=tmp_fill_ds,
            output_raster=tmp_resample_ds,
            intput_format="GTiff",
            output_format="GTiff",
            resampling=config.resample_method,
            outsize=[f"{config.resample_factor * 100}%", "0"],
        )

        # Smooth the dataset

        if config.dynamic_smoothing:
            logger.info("Use dynamic smoothing")
            smooth_raster_dynamic(
                in_ds=tmp_resample_ds,
                out_ds=tmp_smooth_ds,
                smoothing_config=config.dynamic_smoothing,
                tempdir_params=tempdir_params,
            )

        elif config.smoothing:
            logger.info("Use static smoothing")
            smooth_raster(
                in_ds=tmp_resample_ds,
                out_ds=tmp_smooth_ds,
                smoothing_config=config.smoothing,
            )
        else:
            logger.info("No smoothing applied")
            tmp_smooth_ds = tmp_resample_ds

        # Create contours
        logger.info(
            f"Creating contours with {len(config.contour_levels.levels)} levels"
        )
        gdal_contour(
            input_raster=tmp_smooth_ds,
            output_vector=tmp_gpkg_contour,
            output_format=config.output_format,
            output_layer=config.output_layer,
            contour_levels=config.contour_levels.levels,
            input_band=1,
            polygonize=True,
        )
        # Validate contours
        logger.info("Exploding and validating contours")
        ogr2ogr(
            src_dataset=tmp_gpkg_contour,
            dst_dataset=tmp_gpkg_valid,
            layer_name=config.output_layer,
            input_format=config.output_format,
            output_format=config.output_format,
            overwrite=config.overwrite,
            explode=True,
            make_valid=True,
        )

        # Clip layer
        logger.info(f"Clipping layer to {config.clip_layer.layer}")
        ogr2ogr(
            src_dataset=tmp_gpkg_valid,
            dst_dataset=tmp_gpkg_clip,
            layer_name=config.output_layer,
            input_format=config.output_format,
            output_format=config.output_format,
            overwrite=config.overwrite,
            clip_ds=config.clip_layer.dataset,
            clip_layer=config.clip_layer.layer,
        )

        gdf = load_vector_layer(
            dataset=tmp_gpkg_clip,
            layer=config.output_layer,
        )

        logger.info(f"Adding color and labels to contours in {config.output_layer}")
        gdf = color_layer(
            gdf=gdf,
            levels=config.contour_levels.levels,
            colormap=config.contour_levels.colormap,
        )
        if debug:
            gdf.to_file(
                Path(tmpdir) / "colors.gpkg",
                driver="GPKG",
                layer=config.output_layer,
            )

        if config.eliminate_polygon_threshold:
            logger.info(
                f"Eliminating small polygons with threshold: {config.eliminate_polygon_threshold}"
            )
            gdf = eliminate_small_polygons(
                gdf=gdf,
                threshold=config.eliminate_polygon_threshold,
            )
            if debug:
                gdf.to_file(
                    Path(tmpdir) / "eliminated.gpkg",
                    driver="GPKG",
                    layer=config.output_layer,
                )

        if config.nodata_feature_layer:
            logger.info(
                f"Replacing features with nodata features from {config.nodata_feature_layer.layer}"
            )
            gdf = replace_with_nodata(
                gdf=gdf,
                nodata_gdf=load_vector_layer(
                    dataset=config.nodata_feature_layer.dataset,
                    layer=config.nodata_feature_layer.layer,
                ),
            )
            if debug:
                gdf.to_file(
                    Path(tmpdir) / "nodata.gpkg",
                    driver="GPKG",
                    layer=config.output_layer,
                )

        logger.info(
            f"Saving vector layer {config.output_layer} to {config.output_dataset}"
        )
        save_vector_layer(
            gdf=gdf,
            output_dataset=config.output_dataset,
            output_layer=config.output_layer,
        )

        if debug:
            shutil.copy(
                config.output_dataset,
                Path(tmpdir) / "result.gpkg",
            )

        logger.info(f"Completed processing dataset: {config.output_layer}")


if __name__ == "__main__":
    pass
