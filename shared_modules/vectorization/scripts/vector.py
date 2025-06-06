import logging
from pathlib import Path

import geopandas as gpd
import pandas as pd

from .config import (
    INF_MAX,
    INF_MIN,
    MAX_FIELD_NAME,
    MIN_FIELD_NAME,
    NODATA_COLOR,
    NODATA_LABEL,
)

logger = logging.getLogger(__name__)


def load_vector_layer(
    dataset: Path,
    layer: str,
) -> gpd.GeoDataFrame:
    """
    Load a vector layer from a dataset.
    """
    logger.debug(f"Loading vector layer {layer} from {dataset}")
    gdf = gpd.read_file(dataset, layer=layer)
    logger.debug(f"Loaded {len(gdf)} features from {layer}")

    return gdf


def save_vector_layer(
    gdf: gpd.GeoDataFrame,
    output_dataset: Path,
    output_layer: str,
    remove_cols: list = ["ID", "index"],
) -> None:
    """
    Save a vector layer to a dataset.
    """
    logger.debug(f"Saving vector layer {output_layer} to {output_dataset}")

    # Check if the output folder
    if not output_dataset.parent.exists():
        logger.debug(f"Creating output folder {output_dataset.parent}")
        output_dataset.parent.mkdir(parents=True, exist_ok=True)

    for col in remove_cols:
        if col in gdf.columns:
            gdf.drop(columns=col, inplace=True)
            logger.debug(f"Removed column {col} from {output_layer}")

    gdf.to_file(output_dataset, layer=output_layer, driver="GPKG", index=False)
    logger.debug(f"Saved {len(gdf)} features to {output_layer}")


def color_layer(
    gdf: gpd.GeoDataFrame,
    levels: list,
    colormap: list,
    max_field_name: str = MAX_FIELD_NAME,
    min_field_name: str = MIN_FIELD_NAME,
    inf_min: int = INF_MIN,
    inf_max: int = INF_MAX,
    explode: bool = True,
) -> gpd.GeoDataFrame:
    logger.debug(f"Levels: {levels}")
    logger.debug(f"Colormap: {colormap}")

    gdf.copy()
    gdf = gdf[gdf.geometry.notnull()]

    colors = []
    labels = []

    # iterate rows and set color based on levels
    for i, row in gdf.iterrows():
        row_colored = False
        for j in range(len(levels) - 1):
            if (
                row[min_field_name] >= levels[j]
                and row[max_field_name] <= levels[j + 1]
            ):
                colors.append(colormap[j])

                label_min = levels[j]
                label_max = levels[j + 1]
                if label_min == inf_min:
                    label = f"< {label_max}"
                elif label_max == inf_max:
                    label = f"> {label_min}"
                else:
                    label = f"{label_min} - {label_max}"

                labels.append(label)

                row_colored = True
                break

        if not row_colored:
            colors.append(None)
            labels.append(None)
            logger.warning(
                f"Row {i} not colored or labeled. Min: {row[min_field_name]}, Max: {row[max_field_name]}. Setting color and label to None."
            )

    gdf["color"] = colors
    gdf["label"] = labels

    # explode
    if explode:
        gdf = gdf.explode(column="geometry")
        logger.debug(f"Exploded {len(gdf)} features.")

    logger.debug(f"Colored {len(gdf)} features.")

    return gdf


def eliminate_small_polygons(
    gdf: gpd.GeoDataFrame,
    threshold: float,
    keep_islands: bool = True,
) -> gpd.GeoDataFrame:
    logger.debug(f"Threshold: {threshold}")

    # Get original columns
    original_columns = gdf.columns.tolist()

    logger.debug(f"Dropped {len(gdf.geometry.isna())} empty geometries")
    gdf = gdf.copy().dropna(subset=["geometry"])

    gdf = gdf.explode().reset_index()  # multi-part to single-part

    logger.debug(f"Original number of features: {len(gdf)}")

    merged = 0
    removed = 0

    gdf["area"] = gdf.geometry.area
    small_polygons_idx = list(gdf[gdf.area < threshold].sort_values("area").index)

    logger.debug(f"Number of small polygons: {len(small_polygons_idx)}")

    while small_polygons_idx:
        idx = small_polygons_idx.pop(0)
        if idx not in gdf.index:
            logger.debug(f"Polygon {idx} already removed. Skipping.")
            continue

        geom = gdf.loc[idx].geometry
        neighbors_touch = gdf.touches(geom)
        neighbors = gdf.loc[neighbors_touch]

        if idx in neighbors:
            neighbors = neighbors.drop(idx)

        if neighbors.empty:
            logger.debug(f"Polygon {idx} is an island. Keeping it as is.")
            continue

        largest_neighbor_idx = neighbors.geometry.area.idxmax()
        largest_neighbor = gdf.loc[largest_neighbor_idx]
        logger.debug(
            f"Largest neighbor of polygon {idx} is index: {largest_neighbor_idx}"
        )

        merged_polygon = gpd.GeoSeries(
            [
                largest_neighbor.geometry,
                geom,
            ]
        ).union_all()

        logger.debug(
            f"Merged polygon {idx} with {largest_neighbor_idx}. New area: {merged_polygon.area}"
        )

        gdf.at[largest_neighbor_idx, "geometry"] = merged_polygon
        gdf.at[largest_neighbor_idx, "area"] = merged_polygon.area

        gdf.drop(idx, inplace=True)
        removed += 1
        merged += 1

        small_polygons_idx = [
            i
            for i in small_polygons_idx
            if i in gdf.index and gdf.at[i, "area"] < threshold
        ]

    logger.info(
        f"Number of features removed: {removed} (merged: {merged}, dropped: {removed - merged})"
    )

    if not keep_islands:
        gdf = gdf[gdf.geometry.area >= threshold]
        logger.debug(f"Removed islands. Number of features: {len(gdf)}")

    gdf = gdf.reset_index(drop=True)
    gdf = gdf[original_columns]
    return gdf


def replace_with_nodata(
    gdf: gpd.GeoDataFrame,
    nodata_gdf: gpd.GeoDataFrame,
    nodata_color: str = NODATA_COLOR,
    nodata_label: str = NODATA_LABEL,
) -> gpd.GeoDataFrame:
    """
    Replace features in gdf with nodata_gdf where they overlap.
    """

    original_columns = gdf.columns.tolist()

    # Add color and label columns to nodata_gdf
    logger.debug(f"Nodata color: {nodata_color}")
    logger.debug(f"Nodata label: {nodata_label}")
    nodata_gdf["color"] = nodata_color
    nodata_gdf["label"] = nodata_label

    gdf = pd.concat([gdf.overlay(nodata_gdf, how="difference"), nodata_gdf])

    gdf = gdf.reset_index(drop=True)
    gdf = gdf[original_columns]
    return gdf
