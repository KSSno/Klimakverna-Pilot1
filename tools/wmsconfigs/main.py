#!/usr/bin/python3 

import pandas as pd
import json
import yaml

# setup style list for one variable in a collection
def style_config_from_variable(collection, varname, palettes):
    # one variable may have multiple palettes, e.g. for different seasons
    var_styles = []
    for palette in palettes:
        # parse colors into list
        if 'colors' not in palette:
            print(f'Skipping {varname} in {collection} - no colors defined for palette')
            return None

        colors = [x.strip() for x in palette['colors'].split(',')]

        # parse intervals into list
        if 'intervals' not in palette:
            print(f'Skipping {varname} in {collection} - no intervals defined for palette')
            return None

        intervals = [float(x) for x in palette['intervals'].split(',')]

        # open intervals are default
        open_lowest_interval = palette.get('openLowestInterval', True)
        open_highest_interval = palette.get('openHighestInterval', True)

        # add variable names
        variables = palette.get('variables', [])

        # make one style per variable
        for variable in variables:
            # define style config for given variables
            style_config = {
                'name': variable,
                'colors': colors,
                'intervals': intervals,
                'open_lowest_interval': open_lowest_interval,
                'open_highest_interval': open_highest_interval,
            }
            # append to config
            var_styles.append(style_config)

    return var_styles


# for one file of objects in a collection, extract the information needed for WMS config
def wms_config_from_file(input_dir, collection):
    # read object config from file
    filepath = input_dir + collection + '.json'
    with open(filepath, 'r', encoding='utf-8') as f:
        object_config = json.load(f)

    # loop over objects in the config
    elements = []
    for varname, config in object_config.items():
        # palettes defined at all?
        palettes = config.get('palettes', None)
        if palettes is None:
            print(f'Skipping {varname} in {collection} - no palettes defined')
            continue

        # TODO: deduce the pattern based on collection config and variable name

        # add palettes for each applicable variable name
        var_styles = style_config_from_variable(collection, varname, palettes)
        if var_styles is None:
            continue

        # define variable level config and append to elements
        elements.append({
            'pattern': '',
            'base_netcdf_directory': '',
            'module': 'mapgen.modules.generic_quicklook',
            'module_function': 'generic_quicklook',
            'mapfiles_path': '/mapfiles',
            'styles': var_styles,
        })

    return elements

#
def produce_wms_config(input_dir, output_dir, collections):
    # prepate config list for all objects
    wms_config = []

    # loop over collections and parse their object config files
    for collection in collections:
        collection_config = wms_config_from_file(input_dir, collection)
        if collection_config is not None:
            wms_config.extend(collection_config)

    # write to YAML file
    output_filename = output_dir + 'wms_config.yaml'
    with open(output_filename, 'w', encoding='utf-8') as f:
        yaml.dump(wms_config, f, default_flow_style=None, allow_unicode=True, explicit_start=True)


if __name__ == '__main__':
    # directories
    input_dir = 'tools/kinconfigs/output/'
    output_dir = 'tools/wmsconfigs/output/'

    # collections that should be published through WMS
    collections = [
        'ReferenceIndices',
        'ClimateStatistics',
    ]

    produce_wms_config(input_dir, output_dir, collections)