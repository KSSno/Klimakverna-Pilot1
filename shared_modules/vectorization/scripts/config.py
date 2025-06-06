import logging
from pathlib import Path
from typing import Annotated, Any, Optional

from pydantic import BeforeValidator
from pydantic_settings import BaseSettings, SettingsConfigDict


def _folder_validator(v: Any) -> Optional[Path]:
    """
    Convert empty string, 0, False, or "." to None for LOG_FOLDER.
    This allows disabling file logging by setting KLIMAKART_LOG_FOLDER to an empty value.
    """
    if v in ["", 0, False, "0", "false", "False"]:
        return None
    return v


class Settings(BaseSettings):
    """
    Application settings using pydantic-settings.
    All settings can be overridden by environment variables with the KLIMAKART_ prefix.
    """

    # Configuration for environment variables
    model_config = SettingsConfigDict(env_prefix="KLIMAKART_")

    # Logging configuration
    LOG_LEVEL: int = logging.INFO
    LOG_FOLDER: Annotated[Path | None, BeforeValidator(_folder_validator)] = Path(
        "logs"
    )
    LOG_DAYS: int = 30  # Number of days to keep log files

    # Debug mode
    DEBUG_FOLDER: Annotated[Path | None, BeforeValidator(_folder_validator)] = Path(
        "debug"
    )

    # Infinite min and max values
    INF_MIN: int = -9999
    INF_MAX: int = 9999

    # Max and min field names
    MAX_FIELD_NAME: str = "max_value"
    MIN_FIELD_NAME: str = "min_value"

    # Nodata for features
    NODATA_COLOR: str = "#AAD3DF"
    NODATA_LABEL: str = "nodata"


# Create settings instance
settings = Settings()

# Export constants from settings for backward compatibility
LOG_LEVEL = settings.LOG_LEVEL
LOG_FOLDER = settings.LOG_FOLDER
LOG_DAYS = settings.LOG_DAYS
DEBUG_FOLDER = settings.DEBUG_FOLDER
INF_MIN = settings.INF_MIN
INF_MAX = settings.INF_MAX
MAX_FIELD_NAME = settings.MAX_FIELD_NAME
MIN_FIELD_NAME = settings.MIN_FIELD_NAME
NODATA_COLOR = settings.NODATA_COLOR
NODATA_LABEL = settings.NODATA_LABEL


class DefultContourConfig:
    """
    Default contour configuration for the contour levels and color maps.
    """

    srs = "EPSG:25833"
    data_folder = Path("data")
    output_dataset = data_folder / "output.gpkg"
    clip_layer = "N500_landflate"
    nodata_feature_layer = "N500_innsjo"
    resample_factor = 2
    resample_method = "cubicspline"

    """
    smoothing = {
        "method": "gaussian",
        "sigma": 1,
    }
    """

    dynamic_smoothing = {
        "low": {
            "method": "gaussian",
            "sigma": 1,
        },
        "high": {
            "method": "gaussian",
            "sigma": 3,
        },
    }

    # smoothing = "gaussian"
    eliminate_polygon_threshold = 1000000

    contour_levels = {
        "nullgradspassering_endring_9": {
            "levels": [INF_MIN, -40, -30, -20, -10, 10, 20, 30, 40, INF_MAX],
            "colormap": [
                "#1D6936",
                "#398439",
                "#72A472",
                "#ABC4AB",
                "#E5E5E5",
                "#D3B3E7",
                "#C282EA",
                "#B151ED",
                "#A020F0",
            ],
        },
        "nullgradspassering_endring_7": {
            "levels": [INF_MIN, -20, -10, -5, 5, 10, 20, INF_MAX],
            "colormap": [
                "#1D6936",
                # "#398439",
                "#72A472",
                "#ABC4AB",
                "#E5E5E5",
                "#D3B3E7",
                # "#C282EA",
                "#B151ED",
                "#A020F0",
            ],
        },
        "nullgradspassering_referanse_8": {
            "levels": [0, 25, 50, 75, 100, 125, 150, 175, 200],
            "colormap": [
                "#8C510A",
                "#BF812D",
                "#DFC27D",
                "#F6E8C3",
                # "#F5F5F5",
                "#C7EAE5",
                "#80CDC1",
                "#35978F",
                "#01665E",
            ],
        },
        "nullgradspassering_referanse_6": {
            "levels": [0, 50, 75, 100, 125, 150, 200],
            "colormap": [
                "#8C510A",
                "#BF812D",
                # "#DFC27D",
                "#F6E8C3",
                # "#F5F5F5",
                "#C7EAE5",
                # "#80CDC1",
                "#35978F",
                "#01665E",
            ],
        },
        "skisesong_referanse": {
            "levels": [0, 30, 60, 90, 120, 150, 180, 210, 240, INF_MAX],
            "colormap": [
                "#CCF57A",
                "#D9FFFF",
                "#B3FFFF",
                "#80ECFF",
                "#40CCFF",
                "#0099FF",
                "#0019FF",
                "#000099",
                "#00005B",
            ],
        },
        "skisesong_endring": {
            "levels": [INF_MIN, -90, -60, -30, 0, INF_MAX],
            "colormap": [
                "#8D7B68",
                "#A4907C",
                "#C8B6A6",
                "#F1DEC9",
                "#D3E8F7",
            ],
        },
        "kraftig_nedbor_referanse": {
            "levels": [0, 5, 10, 15, 20, 25, INF_MAX],
            "colormap": [
                "#CCAA66",
                "#FFECBF",
                "#B5DE9F",
                "#6ABBD8",
                "#3484C9",
                "#45508A",
            ],
        },
        "kraftig_nedbor_endring": {
            "levels": [INF_MIN, 0, 4, 8, 12, 16, INF_MAX],
            "colormap": [
                "#F6E8C3",
                "#C7EAE5",
                "#80CDC1",
                "#35978F",
                "#01665E",
                "#003C30",
            ],
        },
    }

    _ds_unique = [
        # ──────────────── NULLGRADSPASSERING ANNUAL ────────────────
        {
            "dataset": "nullgradspasseringer_referanseperiode_1991-2020_annual-mean.nc",
            "subdataset": "dzc",
            "band": 1,
            "contour_levels": "nullgradspassering_referanse_8",
            "output_layer": "nullgradspasseringer_referanseperiode_1991-2020_annual-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_annual",
            "band": 1,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2041-2070_annual-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_annual",
            "band": 2,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2071-2100_annual-mean",
        },
        # ──────────────── NULLGRADSPASSERING autumn ────────────────
        {
            "dataset": "nullgradspasseringer_referanseperiode_1991-2020_autumn-mean.nc",
            "subdataset": "dzc",
            "band": 1,
            "contour_levels": "nullgradspassering_referanse_8",
            "output_layer": "nullgradspasseringer_referanseperiode_1991-2020_autumn-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_autumn",
            "band": 1,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2041-2070_autumn-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_autumn",
            "band": 2,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2071-2100_autumn-mean",
        },
        # ──────────────── NULLGRADSPASSERING spring ────────────────
        {
            "dataset": "nullgradspasseringer_referanseperiode_1991-2020_spring-mean.nc",
            "subdataset": "dzc",
            "band": 1,
            "contour_levels": "nullgradspassering_referanse_8",
            "output_layer": "nullgradspasseringer_referanseperiode_1991-2020_spring-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_spring",
            "band": 1,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2041-2070_spring-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_spring",
            "band": 2,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2071-2100_spring-mean",
        },
        # ──────────────── NULLGRADSPASSERING summer ────────────────
        {
            "dataset": "nullgradspasseringer_referanseperiode_1991-2020_summer-mean.nc",
            "subdataset": "dzc",
            "band": 1,
            "contour_levels": "nullgradspassering_referanse_8",
            "output_layer": "nullgradspasseringer_referanseperiode_1991-2020_summer-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_summer",
            "band": 1,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2041-2070_summer-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_summer",
            "band": 2,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2071-2100_summer-mean",
        },
        # ──────────────── NULLGRADSPASSERING winter ────────────────
        {
            "dataset": "nullgradspasseringer_referanseperiode_1991-2020_winter-mean.nc",
            "subdataset": "dzc",
            "band": 1,
            "contour_levels": "nullgradspassering_referanse_8",
            "output_layer": "nullgradspasseringer_referanseperiode_1991-2020_winter-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_winter",
            "band": 1,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2041-2070_winter-mean",
        },
        {
            "dataset": "nullgradspasseringer_framtidsperioder.nc",
            "subdataset": "diff-dzc_winter",
            "band": 2,
            "contour_levels": "nullgradspassering_endring_9",
            "output_layer": "nullgradspasseringer_framtidsperiode_2071-2100_winter-mean",
        },
        # ─────────────────── SKI ───────────────────
        {
            "dataset": "skisesong_referanseperiode_1991-2020.nc",
            "subdataset": "swelangrenn",
            "band": 1,
            "contour_levels": "skisesong_referanse",
            "output_layer": "skisesong_referanseperiode_1991-2020",
            "nodata_feature_layer": nodata_feature_layer,
            "add_nodata_values": [0],
        },
        {
            "dataset": "skisesong_framtidsperioder_annual-mean.nc4",
            "subdataset": "diff-SWEge60",
            "band": 1,
            "contour_levels": "skisesong_endring",
            "output_layer": "skisesong_framtidsperiode_2041-2070",
            "nodata_feature_layer": nodata_feature_layer,
        },
        {
            "dataset": "skisesong_framtidsperioder_annual-mean.nc4",
            "subdataset": "diff-SWEge60",
            "band": 2,
            "contour_levels": "skisesong_endring",
            "output_layer": "skisesong_framtidsperiode_2071-2100",
            "nodata_feature_layer": nodata_feature_layer,
        },
        # ────────────────── KRAFTIG NEDBØR ──────────────────
        {
            "dataset": "kraftig_nedbor_referanseperiode_1991-2020_annual_mean.nc",
            "subdataset": "sdii",
            "band": 1,
            "contour_levels": "kraftig_nedbor_referanse",
            "output_layer": "kraftig_nedbor_referanseperiode_1991-2020",
        },
        {
            "dataset": "kraftig_nedbor_framtidsperioder_endringer_sdii.nc",
            "subdataset": "change-sdii_annual",
            "band": 1,
            "contour_levels": "kraftig_nedbor_endring",
            "output_layer": "kraftig_nedbor_framtidsperiode_2041-2070",
        },
        {
            "dataset": "kraftig_nedbor_framtidsperioder_endringer_sdii.nc",
            "subdataset": "change-sdii_annual",
            "band": 2,
            "contour_levels": "kraftig_nedbor_endring",
            "output_layer": "kraftig_nedbor_framtidsperiode_2071-2100",
        },
    ]

    _ds_base = {
        "resample_factor": resample_factor,
        "resample_method": resample_method,
        "clip_layer": clip_layer,
        "output_dataset": output_dataset,
        "srs": srs,
        # "smoothing": smoothing,
        "dynamic_smoothing": dynamic_smoothing,
        "eliminate_polygon_threshold": eliminate_polygon_threshold,
    }
    datasets = []
    for ds in _ds_unique:
        ds["dataset"] = data_folder / ds["dataset"]
        # Hack to fix ordering of the parameters in the JSON file
        datasets.append(ds | _ds_base | ds)


# Configuration functionality has been moved to CLI options in main.py
