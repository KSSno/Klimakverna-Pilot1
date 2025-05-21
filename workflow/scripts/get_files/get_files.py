#
# Written by Tyge Løvset, NORCE Research, 2025 for Klimakverna.

import os
import fnmatch
import netCDF4 as nc4
import datetime
import json
import re

def check_common_keys(dict1, dict2):
    """
    Returns True only if all common keys are equal
    """
    common_keys = dict1.keys() & dict2.keys() # Find common keys using set intersection
    for key in common_keys:
        if dict1[key] != dict2[key]:
            return False
    return True


def get_files(config_path_or_dict, filters={}, base_path=None, return_groups=False):
    """
    Searches for files based on a flexible criteria. This function can be
    used in conjunction with a JSON "collection" configuration file.
    This config provides information about a specific collection of data files,
    such as EURO-CORDEX, CMIP or DailyTimeSeries files. All files in a collection
    should share a folder and filename structure. Se a full list of collections
    in the next header [TODO].

    Function Overview
    config        : A template dictionary for how to match files.
    base_path     : Root directory to start the search. If not specified, the one in config is used.
    filtering     : A dictionary with filtering options (e.g., file extensions, size, date).
    return_groups : Return both the list of files and a list of dicts with searchable path components.
    """
    if isinstance(config_path_or_dict, dict):
        config = config_path_or_dict
    else:
        with open(config_path_or_dict, 'r') as f:
            config = json.load(f)

    if base_path is None:
        if 'base_path' in config:
            base_path = config['base_path']
        else:
            return (None, None) if return_groups else None

    baselen = len(base_path)
    unmatched = set()
    outfiles = []
    outgroups = []

    re_folder = re.compile(config['folder_pattern'])
    re_file = re.compile(config['file_pattern'])

    start_year = None
    if 'start_year' in filters.keys() or 'end_year' in filters.keys():
        start_year = filters['start_year'] if 'start_year' in filters.keys() else -10000
        end_year = filters['end_year'] if 'end_year' in filters.keys() else 10000
    
    for root, dirs, files in os.walk(base_path):
        subdir = root[baselen:]

        for file in files:
            # Match the folder name
            match = re.search(re_folder, subdir)

            if not match:
                if not subdir in unmatched:
                    print('NOTE: Unmatched subfolder:', subdir)
                    unmatched.add(subdir)
                continue

            folder_group = match.groupdict()
            
            # Filter on the folder components
            found = True
            for key, val in filters.items():
                if key in folder_group.keys():
                    if not folder_group[key] in val:
                        found = False
                        break;
            if not found:
                continue

            # Match the file name againt RE patten
            match = re.search(re_file, file)
            if not match:
                print('no match: ', file)
                continue

            file_group = match.groupdict()

            # Match common keys from folder name against keys from file name
            if check_common_keys(folder_group, file_group) == False:
                print('mismatch: ', file)
                continue

            # Merge the two groups
            merged_groups = folder_group | file_group
            
            # Special handing of 'start_year', 'end_year' range filter keys:
            if start_year:
                g_start_year = int(merged_groups['start_year']) if 'start_year' in merged_groups.keys() else -10000
                g_end_year = int(merged_groups['end_year']) if 'end_year' in merged_groups.keys() else 10000

                if start_year > g_end_year or end_year < g_start_year:
                    #print('out of year range', file)
                    continue

            # Concat the subdir and the file name
            relative_path = os.path.join(subdir, file)
            
            # Add to the output
            outfiles.append(relative_path)
            if return_groups:
                outgroups.append(merged_groups)

    return (outfiles, outgroups) if return_groups else outfiles


if __name__ == '__main__':
    # ARC-44/AWI/MPI-M-MPI-ESM-LR/rcp85/r1i1p1/HIRHAM5/v2/day/tas/v20161011/
    # tas_ARC-44_MPI-M-MPI-ESM-LR_rcp85_r1i1p1_AWI-HIRHAM5_v2_day_20060101-20101231.nc

    euro_cordex_cmip5 = {
        "name": "EURO-CORDEX CMIP5",
        "folder_pattern": "(?P<domain>[^/]+)/(?P<institute>[^/]+)/(?P<model>[^/]+)/(?P<scenario>[^/]+)/(?P<ensemble>[^/]+)/(?P<rcm>[^/]+)/(?P<rcm_version>[^/]+)/(?P<frequency>[^/]+)/(?P<variable>[^/]+)/(?P<create_ver>[^/]+)",
        "file_pattern": "(?P<variable>[^_]+)_(?P<domain>[^_]+)_(?P<model>[^_]+)_(?P<scenario>[^_]+)_(?P<ensemble>[^_]+)_(?P<institute_rcm>[^_]+)_(?P<rcm_version>[^_]+)_(?P<frequency>[^_]+)_(?P<start_year>[0-9]{4})[0-9]*-(?P<end_year>[0-9]{4})[0-9]*\\.nc",
        "filter_keys": ["domain", "institute", "model", "scenario", "ensemble", "rcm", "rcm_version", "frequency", "variable", "create_ver", "start_year", "end_year"],
        "base_path": "/lustre/storeC-ext/users/kin2100/MET/cordex/output/"
    }

    filters = {
        'scenario': ['rcp45'],
        'variable': ['pr'],
        'start_year': 2010, # will also accept if only one of start_year and end_year exists
        'end_year': 2050,
    }

    #files = get_files(euro_cordex_cmip5, filters)
    #files = get_files(euro_cordex_cmip5, filters, '/lustre/storeC-ext/users/kin2100/MET/cordex/output/')
    files, groups = get_files('euro-cordex_cmip5.json', filters, return_groups=True)

    print('Print the first 10:')
    for i in range(min(len(files), 10)):
        print(files[i], ' : ', groups[i])

    print('Found:', len(files))
