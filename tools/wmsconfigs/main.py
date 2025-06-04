#!/usr/bin/python3 

import pandas as pd
import json
import yaml

# yaml.dump(dict, filename, default_flow_style=False, allow_unicode=True)

# directories
input_dir = 'tools/kinconfigs/output/'
output_dir = 'tools/wmsconfigs/output/'

# collections that should be published through WMS
collections = [
    'ReferenceIndices',
    'ClimateStatistics',
]

# prepate config list for all objects
wms_config = []

# loop over collections
for collection in collections:
    # read collection config
    filepath = input_dir + collection + '.json'
    with open(filepath, 'r', encoding='utf-8') as f:
        object_config = json.load(f)
    
    print('Processing objects in collection', collection)

    # loop over objects in object config
    for varname, config in object_config.items():
        # palettes defined at all?
        palettes = config.get('palettes', None)
        if palettes is None:
            print(f'Skipping {varname} in {collection} - no palettes defined')
            continue

        # TODO: deduce the pattern based on collection config and variable name

        # any palettes defined for this variable
        var_styles = []
        for palette in palettes:
            # parse colors into list
            if 'colors' not in palette:
                print(f'Skipping {varname} in {collection} - no colors defined for palette')
                continue
            colors = [x.strip() for x in palette['colors'].split(',')]

            # parse intervals into list
            if 'intervals' not in palette:
                print(f'Skipping {varname} in {collection} - no intervals defined for palette')
                continue
            intervals = [float(x) for x in palette['intervals'].split(',')]

            # open intervals are default
            open_lowest_interval = palette.get('OpenLowestInterval', True)
            open_highest_interval = palette.get('OpenHighestInterval', True)

            # add variable names
            # TODO: drop these. Instead print once per suffix and use full variable name as name
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

        # define variable level config
        var_config = {
            'pattern': '',
            'base_netcdf_directory': '',
            'module': 'mapgen.modules.generic_quicklook',
            'module_function': 'generic_quicklook',
            'mapfiles_path': '/mapfiles',
            'styles': var_styles,
        }
        # append to WMS config
        wms_config.append(var_config)

# write to YAML file
output_filename = output_dir + 'wms_config.yaml'
with open(output_filename, 'w', encoding='utf-8') as f:
    yaml.dump(wms_config, f, default_flow_style=None, allow_unicode=True)