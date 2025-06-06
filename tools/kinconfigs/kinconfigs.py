#!/usr/bin/python3

import pandas as pd
import json
from pathlib import Path

# Note: global metadata (metadata.json) will be made manually


# parse palette CSV file to extract palettes
def parse_palette_file(palette_file):
    """
    Columns of palette CSV:
        Collection:
            Refers to the collection of the variable,
            e.g. ReferenceIndices, ClimateStatistics
        Variable:
            (object level) variable name as used in other contexts
        Suffixes:
            Comma separated list of suffixes that will be used
            to create variable names. May be empty.
                Example: Variable "tas", Suffixes "winter,summer":
                variable names will be {tas_winter, tas_summer}
        Palette:
            Comma separated list of hex codes
        Intervals:
            Comma separated list of interval limits
        OpenLeft and OpenRight:
            Boolean indicators defining the type of interval
    """

    # cache palettes
    try:
        df_pal = pd.read_csv(palette_file, sep=",", header=0)
        df_pal = df_pal.set_index("Variable")
    except FileNotFoundError:
        print(f"Palettes file {palette_file} does not exist. Please create it.")
        return None
    except pd.errors.ParserError as e:
        print(f"Error parsing {palette_file}: {e}")
        return None
    else:
        return df_pal


# parse collection CSV files to extract metadata and palettes
def create_collection_config(input_dir, output_dir, collection, palettes):
    """
    Header rows of collection CSV files:
        Header row 1:
            Indicates type of attribute.
            One of {global, data variable, variable:<string>}.
        Header row 2:
            Index "Variabelnavn" plus pivoted attribute field names
    """

    # open and parse parameter CSV file for given collection
    filepath = input_dir / f"{collection}.csv"

    try:
        df_col = pd.read_csv(filepath, sep=",", header=[0, 1])
    except FileNotFoundError:
        print(f"Parameter file {filepath} not found.")
        return None
    except pd.errors.ParserError as e:
        print(f"Error parsing {filepath}: {e}")
        return None
    except ValueError as e:
        print(f"Error reading {filepath}: {e}")
        return None

    # sort by type of attribute
    df_col = df_col.T.sort_index().T

    # prepare config for all objects in this collection
    config = {}
    for _, row in df_col.iterrows():
        # extract and remove variable name
        varname = row["attr_typename"]["Variabelnavn"].strip()
        row = row.drop("attr_typename")

        config[varname] = object_config(collection, varname, row, palettes)

    # write to json file matching collection name
    config_filename = output_dir / f"{collection}.json"
    with open(config_filename, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2, ensure_ascii=False)


# return dict of all config for one object
def object_config(collection, varname, row, palettes):
    varconfig = {}

    # set metadata config if any
    metadata = object_config_metadata(collection, varname, row)
    if metadata is not None:
        varconfig["metadata"] = metadata

    # set palette config if any
    palettes, varnames = object_config_palettes(collection, varname, row, palettes)
    if palettes is not None and len(palettes) > 0:
        varconfig["palettes"] = palettes

    # set varnames config
    varconfig["varnames"] = varnames

    return varconfig


# define metadata config
def object_config_metadata(collection, varname, row):
    # get unique attribute types
    types = row.index.get_level_values(level=0).unique()

    # loop over each type (global, data variables, variable:<string>
    metadata_config = {}
    for attr_type in types:
        # pandas dataframe to dict
        dict_wnans = row[attr_type].to_dict()
        # omit nan fields
        dict_nonan = {
            k: dict_wnans[k] for k in dict_wnans if not pd.isna(dict_wnans[k])
        }
        # add to config if any fields were defined for this varname
        if len(dict_nonan) > 0:
            metadata_config[attr_type] = dict_nonan

    return metadata_config


# define palette config
def object_config_palettes(collection, varname, row, palettes):
    # skip if nothing defined for this collection
    df_pal_sub = palettes[palettes["Collection"] == collection]
    if len(df_pal_sub) == 0:
        return None, [varname]

    # skip if nothing defined for this variable
    if varname not in df_pal_sub.index:
        print(f"Did not find {collection} variable {varname} in palette CSV")
        return None, [varname]

    df_pal_rows = df_pal_sub.loc[[varname]]

    # loop over each sub definition
    palette_config = []
    all_varnames = []
    for _, row in df_pal_rows.iterrows():
        # deduce variable names for this palette config
        if pd.isna(row["Suffixes"]):
            varnames = [varname]
        else:
            suffixes = [x.strip() for x in row["Suffixes"].split(",")]
            varnames = [varname + "_" + v for v in suffixes]

        # save list of variable names for later
        all_varnames.extend(varnames)

        # are palettes and intervals both defined?
        nanpal = pd.isna(row["Palette"])
        nanint = pd.isna(row["Intervals"])
        if nanpal or nanint:
            print(f"Skipping incomplete definition for {collection} variable {varname}")
            continue

        # add row to config
        palette_config.append(
            {
                "colors": row["Palette"],
                "intervals": row["Intervals"],
                "openLowestInterval": row["OpenLeft"],
                "openHighestInterval": row["OpenRight"],
                "variables": varnames,
            }
        )

    return palette_config, all_varnames


# when run, create json for all defined collections
if __name__ == "__main__":
    # directories
    basedir = Path(__file__).resolve().parent
    input_dir = basedir / "input"
    output_dir = basedir / "output"
    palette_file = input_dir / "palettes.csv"

    # collections to use
    collections = [
        # "DailyReferenceData",
        "ReferenceIndices",
        "DailyTimeSeries",
        # "YearlyTimeSeries",
        # "30YearStatistics",
        "ClimateStatistics",
        # "RegionalStatistics"
    ]

    # parse palette file
    palettes = parse_palette_file(palette_file)
    if palettes is None:
        exit(1)

    # loop over all collections and create config files
    for collection in collections:
        create_collection_config(input_dir, output_dir, collection, palettes)
