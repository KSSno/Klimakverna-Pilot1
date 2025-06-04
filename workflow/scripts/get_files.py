''#
# Written by Tyge Løvset, NORCE Research, 2025 for Klimakverna.

import os
import sys
import json
import glob
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



def get_files(config_path_or_dict, filters={}, base_paths=None, return_groups=False, debug=False):
    """
    Searches for files based on a flexible criteria. This function can be
    used in conjunction with a JSON "collection" configuration file.
    This config provides information about a specific collection of data files,
    such as EURO-CORDEX, CMIP or DailyTimeSeries files. All files in a collection
    should share a folder and filename structure. Se a full list of collections
    in the next header [TODO].

    Function overview
    --------------------
    config        : A template dictionary for how to match files.
    filters       : A dictionary with filtering options (e.g., file extensions, size, date).
                    The available filter keys are specified in the config entry "filter_keys".
    base_paths    : Root directories to start the search. If not specified, the one in config is used.
    return_groups : Return both the list of files and a list of dicts with searchable path elements.

    Year range filtering
    --------------------
    The "start_year" and "end_year" filter keys are treated specially to filter on range.
    If only one of them is specified, it is treated as a half-open range.

    In addition, there is a special filter key "years", which expects a list of a mix
    of tuples (start, end) and single year entries. Both "years" and "start_year"/"end_year"
    will be used in the filtering when both are specified.
    """
    if type(config_path_or_dict) is dict:
        config = config_path_or_dict
    else:
        with open(config_path_or_dict, 'r') as f:
            config = json.load(f)

    if base_paths is None:
        if 'base_path' in config:
            base_paths = config['base_path']

    if type(base_paths) is str:
        base_paths = [base_paths]
    elif type(base_paths) is not list:
        print('Error: no "base_path" specified')
        return (None, None) if return_groups else None

    outfiles = []
    outgroups = []

    # Outer loop: each base path:
    for base_path in base_paths:
        if base_path[-1] != '/':
            base_path += '/'

        baselen = len(base_path)
        unmatched = set()

        re_folder = re.compile(config['folder_pattern'])
        re_file = re.compile(config['file_pattern'])

        flt_years = None
        if 'years' in filters.keys():
            flt_years = []  # allow both range and single year: "years": [[2020, 2030], 2050]
            for y in filters['years']:
                if type(y) in (list, tuple):
                    flt_years.append((int(y[0]), int(y[1])))
                else:
                    flt_years.append((int(y), int(y)))

        if 'start_year' in filters.keys() or 'end_year' in filters.keys():
            if not flt_years:
                flt_years = []
            yend = int(filters['end_year']) if 'end_year' in filters.keys() else 10000
            ystart = int(filters['start_year']) if 'start_year' in filters.keys() else -10000
            flt_years.append((int(ystart), int(yend)))

        # Main loop:
        for root, dirs, files in os.walk(base_path):
            subdir = root[baselen:]

            for file in files:
                # Match the folder name
                match = re.search(re_folder, subdir)

                if not match:
                    if not subdir in unmatched:
                        #print('Subfolder not matched:', os.path.join(base_path, subdir))
                        unmatched.add(subdir)
                    continue

                folder_group = match.groupdict()

                # Match the file name againt RE patten
                match = re.search(re_file, file)
                if not match:
                    if debug:
                        print('Filename not matched: %s' % file)
                    continue

                file_group = match.groupdict()

                # Match common keys from folder name against keys from file name
                if check_common_keys(folder_group, file_group) == False:
                    if True: # debug:
                        print('Folder/filename DRS mismatch:')
                        print('subdir:', subdir)
                        print('file:', file)
                        print('foldergroup:', folder_group)
                        print('filegroup:', file_group)
                    continue
                
                # Merge the two groups
                merged_groups = folder_group | file_group

                # Filter on the merged groups:
                found = True
                for key, val in filters.items():
                    if key in merged_groups.keys():
                        if not merged_groups[key] in val:
                            found = False
                            break;
                if not found:
                    if debug:
                        print('filtered group: %s/%s' % (subdir, file))
                        print('               ', merged_groups)
                    continue
                
                
                # Special handing of 'start_year', 'end_year' range filter keys:
                if flt_years and 'start_year' in merged_groups.keys():
                    ystart = int(merged_groups['start_year'])
                    yend = int(merged_groups['end_year']) if 'end_year' in merged_groups.keys() else ystart

                    in_range = False
                    for yr in flt_years:
                        if ystart >= yr[0] and yend <= yr[1]:
                            in_range = True
                            break
                    if not in_range:
                        continue # main loop

                # Concat the root and the file name
                out_path = os.path.join(root, file)
                
                # Add to the output
                outfiles.append(out_path)
                if return_groups:
                    outgroups.append(merged_groups)

    return (outfiles, outgroups) if return_groups else outfiles


def load_collections(collection_root):
    collections = {}
    for json_file in glob.glob(os.path.join(collection_root, '*.json')):
        with open(json_file, 'r') as f:
            print(json_file)
            config = json.load(f)
            collections[config['name']] = config
    return collections



if __name__ == '__main__':
    collections = load_collections('../../config/collections')

    filters = {
        "period": ["near_future_mean", "far_future_mean", "ref_period_mean"],
        "parameter": ["pr"],
        "scenario": ["rcp26", "rcp45"],
        "institution": "CLMcom-ETH",
        #'years': [[1980, 1990], [2020, 2030], 2090, [2095, 2100]],
    }

    if len(filters) > 0:
        print('Filters:')
        for k, v in filters.items():
            print(' ', k, ':', v)    

    coll = '30YearStatistics'
    if len(sys.argv) > 1:
        coll = sys.argv[1]

    print(collections[coll])
    files, groups = get_files(collections[coll], filters, return_groups=True)
    #files, groups = get_files('../../config/collections/cordex_cmip5.json', filters, return_groups=True)

    print('\nPrint the', len(files), 'found:')
    for i in range(len(files)):
        print(files[i])
        print('   ', groups[i])
