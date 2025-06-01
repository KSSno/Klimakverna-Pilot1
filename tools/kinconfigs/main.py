#!/usr/bin/python3 

import pandas as pd
import json

# directories
basedir = 'tools/kinconfigs/'
input_dir = basedir + 'input/'
output_dir = basedir + 'output/'

# Note: global metadata (metadata.json) will be made manually

# Which collection files should be used
# Header row 1: Indicates type of attribute. One of {global, data variable, variable:<string>}.
# Header row 2: Index 'Variabelnavn' plus pivoted attribute field names
# Rows: One per variable
collections = [
    #'DailyReferenceData',
    'ReferenceIndices',
    'DailyTimeSeries',
    #'YearlyTimeSeries',
    #'30YearStatistics',
    'ClimateStatistics',
    #'RegionalStatistics'
]

# Path of palettes file (Note: specifically for collections ClimateStatistics and ReferenceIndices)
# Header: Variable, Type, Palette, Intervals, OpenLeft, OpenRight
# Variable: (object level) variable name as used in other contexts
# Type: One of {ref-annual, ref-seasonal, projection}. Will be used in combination with Variable to make list of actual varnames
#     Example: Variable 'tas', Type 'ref-seasonal': variable names will be {tas_winter, tas_spring, tas_summer, tas_autumn}
# Palette: Comma separated list of hex codes
# Intervals: Comma separated list of interval limits
# OpenLeft and OpenRight: Boolean indicators defining the type of interval
palette_file = input_dir + 'palettes.csv'

# First cache palettes
df_pal = pd.read_csv(palette_file, sep=',', header=0)
df_pal = df_pal.set_index('Variable')

# Setup a dict of palettes connected to collections 
palettes = {
    'ReferenceIndices': df_pal[ df_pal['Type'].isin( ['ref-annual','ref-seasonal'] ) ],
    'ClimateStatistics': df_pal[ df_pal['Type'] == 'projection' ]
}
# what Type means for the definition of variable palettes
pal_suffixes = {
    'projection': ['annual','winter','spring','summer','autumn'],
    'ref-annual': ['annual'],
    'ref-seasonal': ['winter','spring','summer','autumn']
}

# default suffixes when relevant
def_suffixes = {
    'YearlyTimeSeries': ['annual','winter','spring','summer','autumn'],
    '30YearStatistics': ['annual','winter','spring','summer','autumn'],
    'ReferenceIndices': ['annual','winter','spring','summer','autumn'],
    'ClimateStatistics': ['annual','winter','spring','summer','autumn'],
}


# create json file for given collection
def create_collection_config(collection):
    filepath = input_dir + collection + '.csv'
    df_col = pd.read_csv(filepath, sep=',', header=[0,1])
    # sort by type of attribute
    df_col = df_col.T.sort_index().T

    # prepare config for all objects in this collection
    config = {}
    for _, row in df_col.iterrows():
        # extract and remove variable name
        varname = row['attr_typename']['Variabelnavn'].strip()
        row = row.drop('attr_typename')

        config[varname] = object_config(collection, varname, row)

    # write to json file matching collection name
    config_filename = output_dir + collection + '.json'
    with open(config_filename, "w", encoding='utf-8') as f:
        json.dump(config, f, indent=2, ensure_ascii=False)


# return dict of all config for one object
def object_config(collection, varname, row):
    varconfig = {}

    # set metadata config if any
    metadata = object_config_metadata(collection, varname, row)
    if metadata is not None:
        varconfig['metadata'] = metadata

    # set palette config if any
    palettes = object_config_palettes(collection, varname, row)
    if palettes is not None:
        varconfig['palettes'] = palettes

    # set varnames config
    varconfig['varnames'] = object_config_varnames(collection, varname, row, palettes)

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
        dict_nonan = {k: dict_wnans[k] for k in dict_wnans if not pd.isna(dict_wnans[k])}
        # add to config if any fields were defined for this varname
        if len(dict_nonan) > 0:
            metadata_config[attr_type] = dict_nonan

    return metadata_config


# define palette config
def object_config_palettes(collection, varname, row):
    # skip if not relevant to this collection
    if collection not in palettes:
        return None

    df_pal_sub = palettes[collection]

    # skip if not defined for this variable
    if varname not in df_pal_sub.index:
        print('Did not find', collection, 'variable', varname, 'in palette CSV')
        return None

    df_pal_rows = df_pal_sub.loc[[varname]]

    # loop over each sub definition
    palette_config = {}
    for _, row in df_pal_rows.iterrows():
        # skip if this type of palette not recognized
        pal_type = row['Type']
        if pal_type not in pal_suffixes:
            print('Palette type', pal_type, 'unrecognized, found for variable', varname)
            continue

        # are palettes and intervals both defined?
        nanpal = pd.isna(row['Palette'])
        nanint = pd.isna(row['Intervals'])
        if nanpal or nanint:
            print('Skipping incomplete definition for', collection, 'variable', varname)
            continue

        # deduce variable names for this palette config
        suffixes = pal_suffixes[pal_type]
        varnames = [varname + '_' + v for v in suffixes]

        # add row to config
        for var in varnames:
            palette_config[var] = {
                'colors': row['Palette'],
                'intervals': row['Intervals'],
                'openLowestInterval': row['OpenLeft'] == 'TRUE',
                'openHighestInterval': row['OpenRight'] == 'TRUE',
                #'variables': varnames
            }

    return palette_config


# return list of variable names that may appear in files of this object
def object_config_varnames(collection, varname, row, palettes):
    # most collection do not have variable name variations
    if collection not in def_suffixes:
        return [varname]

    # what are the default suffixes for this collection?
    defaults = def_suffixes[collection]

    # if palettes defined for seasons then seasons exist
    # Note: is palettes is None, we don't really know!
    has_seasons = False
    if palettes is not None:
        for k in palettes.keys():
            if 'winter' in k:
                has_seasons = True

    # return full list if seasons, otherwise only _annual
    if has_seasons:
        return [varname + '_' + v for v in defaults]
    else:
        return [varname + '_annual']


# when run, create json for all defined collections
if __name__ == "__main__":
    for collection in collections:
        create_collection_config(collection)