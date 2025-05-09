# Input from user to run testcase 8
# 30 year means and 24-hour time series

import argparse as ap
from pathlib import Path
import yaml
import csv
import json
import arrow

ssp370_models = ["cnrm_hclim", "cnrm_racmo", "ecearth_racmo", "ecearthveg_cclm", "ecearthveg_hclim", "miroc_icon", "mpi_hclim", "mpi_icon", "mpi_racmo", "noresm_hclim"]
rcp26_rcp45_models = ["cnrm_aladin", "ecearth_cclm", "ecearth_hirham", "ecearth_rca", "hadgem_rca", "hadgem_remo", "mpi_cclm", "mpi_remo", "noresm_rca", "noresm_remo"]
klimakverna_pilot1_path = Path("/lustre/storeC-ext/users/klimakverna/development/Klimakverna-Pilot1")
default_config_file = klimakverna_pilot1_path / "config/config.yaml"

with open(default_config_file, "r") as f:
    default_config = yaml.safe_load(f)
    default_config = default_config["testcase_8"]

def parse_config_tsv_file(file_path: str) -> dict:
    config = {}
    try:
        with open(file_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                for key, value in row.items():
                    if key in config:
                        config[key].append(value)
                    else:
                        config[key] = [value]

    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file {file_path} not found")

    return config

def write_config_to_tsv(config: dict, output_file: str) -> None:
    headers = list(config.keys())
    # Transpose the dictionary to rows
    rows = zip(*[config[key] for key in headers])

    with open(output_file, mode="w") as csvfile: # newline=""
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(rows)

parser = ap.ArgumentParser(description="Testcase 8: 30 year means and 24-hour time series in csv output format")
#this one 
parser.add_argument("-c", "--config-file", type=Path, help="Absolute path to configuration file, if not given use default settings; indicator=pr, region=1 and output location=")
# or these ones
parser.add_argument("-i", "--indicator", type=str, help="Indicator for calculation (pr or tas)", default="pr")
parser.add_argument("-s", "--scenarios", type=list, help="Scenario(s) to calculate for (rcp26, rcp45, ssp370, hist)", default=["rcp26", "rcp45"])
parser.add_argument("-p", "--periods", type=list, help="Period for calculation (For 30 year mean: 'hist' for 1991-2020, 'nf' for 2041-2070 or 'ff' for 2071-2100. For daily time series: 2041)", default=["nf"]) # --> calculation-type
parser.add_argument("-m", "--models", type=list, help="Model for calculation (rcp26/rcp45: cnrm_aladin, ssp370: ecearth_racmo and/or hist: cnrm_aladin (rcp26/rcp45) or ecearth_racmo (ssp370))", default=["cnrm_aladin"]) # --> CMIP version
parser.add_argument("-b", "--bias-method", type=str, help="Bias-adjustment method (EQM or 3DBC)", default="3DBC")
parser.add_argument("-ri", "--region-id", type=int, help="Region id in shapefile", default=1)

parser.add_argument("-r", "--region-shapefile", type=Path, help="Absolute path to shapefile for region to do calculation for", default="/lustre/storeC-ext/users/klimakverna/development/kaja/data/shapefile")
parser.add_argument("-o", "--output-location", type=Path, help="Absolute path to output location for csv file")

# defaults
# parser.add_argument("-t", "--calculation-type", type=str, help="Choose calculation to be 30 year means (30y) or 24-hour time series (24h). Default is 30 year means", default="30y")

args = parser.parse_args()
indicator = args.indicator 
periods = args.periods
scenarios = args.scenarios
models = args.models
bias_method = args.bias_method
region_id = args.region_id
region_shapefile = args.region_shapefile
output_location = args.output_location
calculation_type = "30y"

if args.config_file:
    try:
        with open(args.config_file, "r") as f:
            config = json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file {args.config_file} not found")
    except json.JSONDecodeError:
        raise ValueError(f"Configuration file {args.config_file} is not a valid JSON file")

    if "indicator" in config:
        indicator = config["indicator"]
        
    if "periods" in config:
        periods = config["periods"]

    if "scenarios" in config:
        scenarios = config["scenarios"]

    if "models" in config:
        scenarios = config["scenarios"]

    if "bias_method" in config:
        bias_method = config["bias_method"]

    if "region_id" in config:
        region_id = config["region_id"]

    if "region_shapefile" in config:
        region_shapefile = Path(config["region_shapefile"])

    if "output_location" in config:
        output_location = Path(config["output_location"])
           
    # if not indicator or not periods or not scenarios:
    #     raise ValueError("Configuration file must contain indicator, periods and scenarios")

# Validate input values
if indicator not in ["pr", "tas"]:
    raise ValueError(f"Indicator {indicator} not in ['pr', 'tas']")

if "rcp45" and "ssp370" in scenarios or "rcp26" and "ssp370" in scenarios:
    raise ValueError("Scenarios for different CMIP versions can't be mixed")
for scenario in scenarios:
    if scenario not in ["rcp26", "rcp45", "ssp370", "hist"]:
        raise ValueError(f"Scenario {scenario} not in ['rcp26', 'rcp45', 'ssp370', 'hist']")
    
for model in models:
    if "ssp370" in scenarios:
        cmip_version = 6
        if model not in ssp370_models:
            raise ValueError(f"Model {model} is invalid for scenario ssp370. It should be one of {ssp370_models}")
    elif "rcp26" in scenarios or "rcp45" in scenarios:
        cmip_version = 5
        if model not in rcp26_rcp45_models:
            raise ValueError(f"Model {model} is invalid for scenarios rcp26 and rcp45. It should be one of {rcp26_rcp45_models}")
    elif scenarios == ["hist"]: # sjekk denne, da må periods også være hist??
        if model not in rcp26_rcp45_models + ssp370_models:
            raise ValueError(f"Model {model} is invalid for scenario hist. It should be one of {ssp370_models + rcp26_rcp45_models}")
    else:
        raise ValueError(f"Model {model} is invalid. It should be one of {ssp370_models + rcp26_rcp45_models}")

if bias_method not in ["EQM", "3DBC"]:
    raise ValueError(f"Bias adjustment method {bias_method} not in ['EQM', '3DBC']")

if region_id > 7 or region_id < 1:
    raise ValueError(f"Region id {region_id} not in [1, 7]")

if not region_shapefile.exists():
    raise FileNotFoundError(f"Region shapefile {region_shapefile} not found")

if not output_location.exists():
    raise FileNotFoundError(f"Output location {output_location} not found")

# Modify default config according to inputs
# default_config_file contains path to csv and region shapefile, in addition to path to rest of config files (inputs, indicators, sceanrios, periods, region)
if indicator and indicator != "pr":
    # pr is the default indicator
    indicator_config = parse_config_tsv_file(f"{klimakverna_pilot1_path}/{default_config['configurationTables']['indicators']}")

    indicator_config["variables"] = [indicator]
    if indicator == "tas":
        indicator_config["id"] = ["101"]
        indicator_config["name"] = ["Annual mean temperature by period"]
        indicator_config["units"] = ["K"]

    write_config_to_tsv(indicator_config, "/lustre/storeC-ext/users/klimakverna/development/output/testcase_8/test/indicators.tsv")

if periods and periods != ["nf"]:
    # nf is the default period
    '''
    id	name	short_name	start	end
    1	Historical	hist	2018	2020
    2	Near future	nf	2041	2043
    3	Far future	ff	2071	2073
    '''
    periods_config = {"id": [], "name": [], "short_name": [], "start": [], "end": []}

    for period_id, period in enumerate(periods):
        periods_config["id"].append(f"{period_id+1}")
        
        if period == "hist":
            periods_config["short_name"].append(period)
            periods_config["name"].append("Historical")
            periods_config["start"].append("1991")
            periods_config["end"].append("2020")

        elif period == "nf":
            periods_config["short_name"].append(period)
            periods_config["name"].append("Near future")
            periods_config["start"].append("2041")
            periods_config["end"].append("2070")

        elif period == "ff":
            periods_config["short_name"].append(period)
            periods_config["name"].append("Far future")
            periods_config["start"].append("2071")
            periods_config["end"].append("2100")

        else:
            # '2043' or '2043-2045'
            calculation_type = "24h"
            if "-" in period:
                start, end = period.split("-")
                periods_config["start"].append(start)
                periods_config["end"].append(end)
                periods_config["name"].append(f"Years {period}")

            else:
                start = period
                end = period
                periods_config["name"].append(f"Year {period}")
                periods_config["start"].append(period)
                periods_config["end"].append(period)
            
            if arrow.get(end) <= arrow.get("2020") and arrow.get(start) >= arrow.get("1991"):
                name = "hist"
            elif arrow.get(end) <= arrow.get("2070") and arrow.get(start) >= arrow.get("2041"):
                name = "nf"
            elif arrow.get(end) <= arrow.get("2100") and arrow.get(start) >= arrow.get("2071"):
                name = "ff"
            else:
                raise ValueError(f"Period {period} is invalid. It should be in the range of 1991-2020, 2041-2070 or 2071-2100 for 30 year means. For daily time series, it should be a single year or a range of years within one period (e.g. '2043' or '2043-2045').")

            periods_config["short_name"].append(name)

    write_config_to_tsv(periods_config, "/lustre/storeC-ext/users/klimakverna/development/output/testcase_8/test/periods.tsv")

if scenarios and scenarios != ["rcp26", "rcp45"]:
    # rcp26 and rcp45 are the default scenarios
    '''
    id	description	scenarioStrings	hexcolour
    # For CMIP5 models
    #historical	Historical values	_hist_	66C2A5
    #rcp26	Low emissions scenario (RCP2.6)	_rcp26_	FC8D62
    #rcp45	Medium emissions scenario (RCP4.5)	_rcp45_	8DA0CB
    #
    # For CMIP6 models
    historical	Historical values	_hist_	66C2A5
    ssp370	2nd worst scenario (SSP370)	_ssp370_	FC8D62
    '''
    scenarios_config = {"id": [], "description": [], "scenarioStrings": [], "hexcolour": []}
    for scenario in scenarios:
        scenarios_config["id"].append(scenario)
        if scenario == "rcp26":
            scenarios_config["description"].append("Low emissions scenario (RCP2.6)")
            scenarios_config["scenarioStrings"].append("_rcp26_")
            scenarios_config["hexcolour"].append("FC8D62")
        elif scenario == "rcp45":
            scenarios_config["description"].append("Medium emissions scenario (RCP4.5)")
            scenarios_config["scenarioStrings"].append("_rcp45_")
            scenarios_config["hexcolour"].append("8DA0CB")
        elif scenario == "ssp370":
            scenarios_config["description"].append("2nd worst scenario (SSP370)")
            scenarios_config["scenarioStrings"].append("_ssp370_")
            scenarios_config["hexcolour"].append("FC8D62")
        elif scenario == "hist":
            scenarios_config["description"].append("Historical values")
            scenarios_config["scenarioStrings"].append("_hist_")
            scenarios_config["hexcolour"].append("66C2A5")
        else:
            raise ValueError(f"Scenario {scenario} is invalid. It should be one of ['rcp26', 'rcp45', 'ssp370', 'hist']")

    write_config_to_tsv(scenarios_config, "/lustre/storeC-ext/users/klimakverna/development/output/testcase_8/test/scenarios.tsv")

if region_id and region_id != 1 and region_shapefile:
    region_config = {"id": [region_id], "shapefile": [region_shapefile]}
    write_config_to_tsv(region_config, "/lustre/storeC-ext/users/klimakverna/development/output/testcase_8/test/region.tsv")


# inputs
    # with open(default_config["configurationTables"]["inputs"], "r") as f:
    #     input_config_file = csv.reader(f, delimiter="\t")

    #     for line_no, line in enumerate(input_config_file):
    #         #uncomment the correct line or create it if nout found
    #         print(line)

# change cconfig.yaml to use new config files