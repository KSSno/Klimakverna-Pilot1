from pathlib import Path
from typing import Annotated, List, Literal

from pydantic import AfterValidator, BaseModel, BeforeValidator


def normalize_path(path: Path | str) -> Path:
    """
    Ensure path uses forward slashes for cross-platform compatibility.
    """
    # Create a new Path from the posix representation to ensure forward slashes
    return Path(Path(path).as_posix().replace("\\", "/"))


class ContourLevels(BaseModel):
    levels: List[int | float]
    colormap: List[str] | None = None


def validate_contour_levels_before(
    contour_levels: ContourLevels | List,
) -> ContourLevels:
    # If contour_levels is a list, convert it to ContourLevels
    if isinstance(contour_levels, list):
        return ContourLevels(levels=contour_levels)
    return contour_levels


def validate_contour_levels_after(contour_levels: ContourLevels) -> ContourLevels:
    # If colormap is None, create a default colormap with black color
    if contour_levels.colormap is None:
        # Create a list of black colors ("#000000") with length = levels - 1
        contour_levels.colormap = ["#000000"] * (len(contour_levels.levels) - 1)
    elif len(contour_levels.levels) != len(contour_levels.colormap) + 1:
        raise ValueError(
            f"The length of 'levels' ({len(contour_levels.levels)}) must be exactly one more than "
            f"the length of 'colormap' ({len(contour_levels.colormap)})"
        )
    return contour_levels


class FeatureLayer(BaseModel):
    dataset: Annotated[Path, AfterValidator(normalize_path)]
    layer: str
    format: str | None = None
    srs: str | None = None


class Smoothing(BaseModel):
    method: Literal["uniform", "gaussian", "median"]
    size: int = 5
    sigma: float = 1.0
    radius: int | None = None


def validate_smoothing(value: Smoothing | None):
    if isinstance(value, str):
        return Smoothing(method=value)
    return value


class DynamicSmoothing(BaseModel):
    low: Annotated[Smoothing, BeforeValidator(validate_smoothing)]
    high: Annotated[Smoothing, BeforeValidator(validate_smoothing)]


class ContourGeneratorConfig(BaseModel):
    dataset: Annotated[Path, AfterValidator(normalize_path)]
    subdataset: str
    format: str = "netcdf"
    band: int
    contour_levels: Annotated[
        ContourLevels | List,
        BeforeValidator(validate_contour_levels_before),
        AfterValidator(validate_contour_levels_after),
    ]
    resample_factor: float = 2
    resample_method: Literal[
        "nearest", "bilinear", "cubic", "cubicspline", "lanczos", "average", "mode"
    ] = "cubicspline"
    add_nodata_values: List[int | float] = []
    clip_layer: FeatureLayer
    nodata_feature_layer: FeatureLayer | None = None
    output_layer: str
    output_dataset: Annotated[Path, AfterValidator(normalize_path)]
    output_format: str = "GPKG"
    overwrite: bool = True
    srs: str | None = "EPSG:25833"
    smoothing: Annotated[Smoothing | None, BeforeValidator(validate_smoothing)] = None
    dynamic_smoothing: DynamicSmoothing | None = None
    eliminate_polygon_threshold: int | float | None = None


if __name__ == "__main__":
    # Example usage
    import json

    json_data_file = Path("tests/testdata/datasets.json")
    with json_data_file.open("r", encoding="utf-8") as f:
        datasets = json.load(f)
        print(f"Loaded datasets from {json_data_file}")
        print(datasets)

    for ds in datasets:
        print(normalize_path(ds["dataset"]))
