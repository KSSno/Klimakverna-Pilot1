# Input from user to run testcase 8
# 30 year means and 24-hour time series

import argparse as ap
from pathlib import Path
import yaml
import csv
import json
import arrow
import geopandas as gpd
import subprocess

model_info = {
    "cnrm_hclim": {"filename": "cnrm-r1i1p1f2-hclim", "cmip_version": 6},
    "cnrm_racmo": {"filename": "cnrm-r1i1p1f2-racmo", "cmip_version": 6},
    "ecearth_racmo": {"filename": "ecearth-r1i1p1f1-racmo", "cmip_version": 6},
    "ecearthveg_cclm": {"filename": "ecearthveg-r1i1p1f1-cclm", "cmip_version": 6},
    "ecearthveg_hclim": {"filename": "ecearthveg-r1i1p1f1-hclim", "cmip_version": 6},
    "miroc_icon": {"filename": "miroc-r1i1p1f1-icon", "cmip_version": 6},
    "mpi_hclim": {"filename": "mpi-r1i1p1f1-hclim", "cmip_version": 6},
    "mpi_icon": {"filename":"mpi-r1i1p1f1-icon", "cmip_version": 6},
    "mpi_racmo": {"filename": "mpi-r1i1p1f1-racmo", "cmip_version": 6},
    "noresm_hclim": {"filename": "noresm-r1i1p1f1-hclim", "cmip_version": 6},
    "cnrm_aladin": {"filename": "cnrm-r1i1p1-aladin", "cmip_version": 5},
    "ecearth_cclm": {"filename": "ecearth-r12i1p1-cclm", "cmip_version": 5},
    "ecearth_hirham": {"filename": "ecearth-r3i1p1-hirham", "cmip_version": 5},
    "ecearth_rca": {"filename": "ecearth-r12i1p1-rca", "cmip_version": 5},
    "hadgem_rca": {"filename": "hadgem-r1i1p1-rca", "cmip_version": 5},
    "hadgem_remo": {"filename": "hadgem-r1i1p1-remo", "cmip_version": 5},
    "mpi_cclm": {"filename": "mpi-r1i1p1-cclm", "cmip_version": 5},
    "mpi_remo": {"filename": "mpi-r2i1p1-remo", "cmip_version": 5},
    "noresm_rca": {"filename": "noresm-r1i1p1-rca", "cmip_version": 5},
    "noresm_remo": {"filename": "noresm-r1i1p1-remo", "cmip_version": 5}
}
ssp370_models = [model for model, info in model_info.items() if info["cmip_version"] == 6]
rcp26_rcp45_models = [model for model, info in model_info.items() if info["cmip_version"] == 5]

hist_period = {"start": "1991", "end": "2020", "name": "Historical"}
nf_period = {"start": "2041", "end": "2070", "name": "Near future"}
ff_period = {"start": "2071", "end": "2100", "name": "Far future"}

klimakverna_pilot1_path = Path("/lustre/storeC-ext/users/klimakverna/development/Klimakverna-Pilot1")
output_path = Path("/lustre/storeC-ext/users/klimakverna/development/output/testcase_8/model_results")
default_config_file = klimakverna_pilot1_path / "config/config.yaml"
path_eqm_24h = Path("/lustre/storeC-ext/users/kin2100/NVE/EQM")
path_dbc_24h = Path("/lustre/storeC-ext/users/kin2100/MET/3DBC/application")
path_30y_mean = Path("/lustre/storeC-ext/users/kin2100/MET/annmeans_bc/")

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

def get_input_path(model: str, bias_method: str, calculation_type: str, indicator: str) -> str:
    path = ""
    model_name = model_info[model]["filename"]
    cmip_version = model_info[model]["cmip_version"]

    if bias_method == "3DBC" and cmip_version == 5:
        if calculation_type == "24h":
            path = f"{path_dbc_24h}/{model_name}/{indicator}/[rrh][cci][pps][24t]*/*_????.nc4"
        else: # 30y, 30 year means
            path = f"{path_30y_mean}/[r,n,f]*_mean/{indicator}/30yrmean_[n,f,r][f,f,e][-,-,f]*_{model_name}_[r][c][p][24][65]*_3*.nc"

    elif bias_method == "EQM" and cmip_version == 5:
        if calculation_type == "24h":
            path = f"{path_eqm_24h}/{model_name}/{indicator}/*/*_????.nc4"
        else:
            path = f"{path_30y_mean}/[r,n,f]*_mean/{indicator}/30yrmean_[n,f,r][f,f,e][-,-,f]*_{model_name}_[r][c][p][24][65]*_e*.nc"

    elif bias_method == "3DBC" and cmip_version == 6:
        if calculation_type == "24h":
            path = f"{path_dbc_24h}/CMIP6/{model_name}/{indicator}/[sh][si][ps][3t]*/*_????.nc4"
        else:
            path = f"{path_30y_mean}/CMIP6/[r,n,f]*_mean/{indicator}/30yrmean_[n,f,r][f,f,e][-,-,f]*_{model_name}_ssp370_3*.nc"

    elif bias_method == "EQM" and cmip_version == 6:
        if calculation_type == "24h":
            path = f"{path_eqm_24h}/CMIP6/{model_name}/{indicator}/*/*_????.nc4"
        else:
            path = f"{path_30y_mean}/CMIP6/[r,n,f]*_mean/{indicator}/30yrmean_[n,f,r][f,f,e][-,-,f]*_{model_name}_ssp370_e*.nc"
    else:
        # Will never get here? Have already validated input
        raise ValueError(f"Combination of bias adjustment method {bias_method} and CMIP version {cmip_version} is not valid for model {model}")
    
    return path

def set_config_path(change: Path | bool, config_type: str):
    if change:
        default_config["configurationTables"][config_type] = str(change)
    else:
        default_config["configurationTables"][config_type] = f"{klimakverna_pilot1_path}/config/testcase_8/{config_type}.tsv"

calculation_type = "30y"

change_indicator = False
change_periods = False
change_scenarios = False
change_region_id = False
change_region_shapefile = False
change_output_location = False
change_inputs = False

default_indicator = "pr"
default_periods = ["nf"]
default_scenarios = ["rcp26", "rcp45"]
default_models = ["cnrm_aladin"]
default_bias_method = "3DBC"
default_region_id = 1
default_output_location = output_path
default_region_config = parse_config_tsv_file(f"{klimakverna_pilot1_path}/{default_config['configurationTables']['region']}")
default_region_shapefile_path = Path(default_region_config["shapefile"][0])

parser = ap.ArgumentParser(description=(
    "Testcase 8: 30 year means and 24-hour time series in csv output format. "
    "A configuration file or options set in command line overwrites default configuration. "
    f"Default configuration is indicator={default_indicator}, periods={default_periods}, scenarios={default_scenarios}, models={default_models}, "
    f"bias method={default_bias_method}, region id={default_region_id}, region shapefile={default_region_shapefile_path}, output location={default_output_location}"
    )
)      

parser.add_argument("-i", "--indicator", type=str, help="Indicator for calculation (pr or tas)", default=default_indicator)
parser.add_argument("-s", "--scenarios", nargs='+', type=str, help="Scenario(s) to calculate for (rcp26, rcp45, ssp370, hist)", default=default_scenarios)
parser.add_argument("-p", "--periods", nargs='+', type=str, help="Period for calculation (For 30 year mean: 'hist' for 1991-2020, 'nf' for 2041-2070 or 'ff' for 2071-2100. For daily time series: 2041 or 2041-2043)", default=default_periods) # --> calculation-type
parser.add_argument("-m", "--models", nargs='+', type=str, help="Model for calculation", default=default_models)
parser.add_argument("-b", "--bias-method", type=str, help="Bias-adjustment method (EQM or 3DBC)", default=default_bias_method)
parser.add_argument("-ri", "--region-id", type=int, help="Region id in shapefile", default=default_region_id)

parser.add_argument("-r", "--region-shapefile", type=Path, help="Absolute path to shapefile for region to do calculation for", default=default_region_shapefile_path)
parser.add_argument("-o", "--output-location", type=Path, help="Absolute path to output location for csv file and netcdf file", default= default_output_location)

parser.add_argument("-c", "--config-file", type=Path, help="Absolute path to json configuration file. See example /lustre/storeC-ext/users/klimakverna/development/Klimakverna-Pilot1/tools/example_config.json")

parser.add_argument("-n", "--dry-run", action='store_true',  help="If this option is set, a dry run of the snakemake pipeline will be performed where no calculations will be done")

args = parser.parse_args()
indicator = args.indicator 
periods = args.periods
scenarios = args.scenarios
models = args.models
bias_method = args.bias_method
region_id = args.region_id
region_shapefile = args.region_shapefile
output_location = args.output_location
dry_run = args.dry_run

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
        models = config["models"]

    if "bias_method" in config:
        bias_method = config["bias_method"]

    if "region_id" in config:
        region_id = config["region_id"]

    if "region_shapefile" in config:
        region_shapefile = Path(config["region_shapefile"])

    if "output_location" in config:
        output_location = Path(config["output_location"])

# Validate input values
if indicator not in ["pr", "tas"]:
    raise ValueError(f"Indicator {indicator} not in equal to 'pr' or 'tas'")

if "rcp45" in scenarios and "ssp370" in scenarios or "rcp26" in scenarios and "ssp370" in scenarios:
    raise ValueError("Scenarios for different CMIP versions can't be mixed")

for scenario in scenarios:
    if scenario not in ["rcp26", "rcp45", "ssp370", "hist"]:
        raise ValueError(f"Scenario {scenario} not in equal to 'rcp26', 'rcp45', 'ssp370' or 'hist'")

if "hist" in scenarios or "hist" in periods or ("nf" not in periods and "ff" not in periods):
    if "hist" in periods and "hist" not in scenarios:
        raise ValueError("If periods include 'hist', the scenarios also have to include 'hist'")
    
    hist_period_in_config = False
    for period in periods:
        if period == "hist":
            hist_period_in_config = True
        elif period == "nf" or period == "ff":
            continue
        elif "-" in period:
            start, end = period.split("-")
            if (arrow.get(start).is_between(arrow.get(hist_period["start"]), arrow.get(hist_period["end"]), "[]")
                or arrow.get(end).is_between(arrow.get(hist_period["start"]), arrow.get(hist_period["end"]), "[]")):
                hist_period_in_config = True
            
        elif arrow.get(period).is_between(arrow.get(hist_period["start"]), arrow.get(hist_period["end"]), "[]"):
            hist_period_in_config = True

    if hist_period_in_config and "hist" not in scenarios:
        raise ValueError("If periods include years in the 'hist' period, the scenarios also have to include 'hist'")

    if not hist_period_in_config and "hist" in scenarios:
        raise ValueError(f"If sceanrios include 'hist', the periods also have to include 'hist' or years between {hist_period['start']} and {hist_period['end']}")

for model in models:
    if "ssp370" in scenarios:
        if model not in ssp370_models:
            raise ValueError(f"Model {model} is invalid for scenario ssp370. Valid models for this sceanrio are {', '.join(f'{model}' for model in ssp370_models)}")
    elif "rcp26" in scenarios or "rcp45" in scenarios:
        if model not in rcp26_rcp45_models:
            raise ValueError(f"Model {model} is invalid for scenarios rcp26 and rcp45. Valid models for these sceanrios are {', '.join(f'{model}' for model in rcp26_rcp45_models)}")
    elif scenarios == ["hist"]:
        if model not in rcp26_rcp45_models + ssp370_models:
            raise ValueError(f"Model {model} is invalid for scenario hist. Valid models for this sceanrio are {', '.join(f'{model}' for model in ssp370_models + rcp26_rcp45_models)}")
    else:
        raise ValueError(f"Model {model} is invalid. Valid models are {', '.join(f'{model}' for model in ssp370_models + rcp26_rcp45_models)}")

if bias_method not in ["EQM", "3DBC"]:
    raise ValueError(f"Bias adjustment method {bias_method} not equal to 'EQM' or '3DBC'")

if not region_shapefile.exists():
    raise FileNotFoundError(f"Region shapefile {region_shapefile} not found")

if output_location:
    if not output_location.exists():
        raise FileNotFoundError(f"Output location {output_location} not found")
    change_output_location = output_location
else:
    output_location = Path(default_config["dirs"]["csv"])

# Sette default verdi i args eller lese det fra default_config??

# Modify default config according to inputs
# default_config_file contains paths to where the csv and the region netcdf fiels should be stored,
# in addition to paths to rest of config files (inputs, indicators, sceanrios, periods, region)
if indicator != default_indicator:
    indicator_config = parse_config_tsv_file(f"{klimakverna_pilot1_path}/{default_config['configurationTables']['indicators']}")

    indicator_config["variables"] = [indicator]
    if indicator == "tas":
        indicator_config["id"] = ["101"]
        indicator_config["name"] = ["Annual mean temperature by period"]
        indicator_config["units"] = ["K"]

    change_indicator = Path(f"{output_path}/config/indicators.tsv")
    write_config_to_tsv(indicator_config, f"{change_indicator}")

if periods != default_periods:
    periods_config = {"id": [], "name": [], "short_name": [], "start": [], "end": []}

    for period_id, period in enumerate(periods):
        periods_config["id"].append(f"{period_id+1}")
        
        if period == "hist":
            periods_config["short_name"].append(period)
            periods_config["name"].append(hist_period["name"])
            periods_config["start"].append(hist_period["start"])
            periods_config["end"].append(hist_period["end"])

        elif period == "nf":
            periods_config["short_name"].append(period)
            periods_config["name"].append(nf_period["name"])
            periods_config["start"].append(nf_period["start"])
            periods_config["end"].append(nf_period["end"])

        elif period == "ff":
            periods_config["short_name"].append(period)
            periods_config["name"].append(ff_period["name"])
            periods_config["start"].append(ff_period["start"])
            periods_config["end"].append(ff_period["end"])

        else:
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
            
            if arrow.get(end) >= arrow.get(hist_period["start"]) and arrow.get(start) <= arrow.get(hist_period["end"]):
                short_name = "hist"
            elif arrow.get(end) >= arrow.get(nf_period["start"]) and arrow.get(start) <= arrow.get(nf_period["end"]):
                short_name = "nf"
            elif arrow.get(end) >= arrow.get(ff_period["start"]) and arrow.get(start) <= arrow.get(ff_period["end"]):
                short_name = "ff"
            else:
                raise ValueError(f"Period {period} is invalid. It should be equal to nf, ff or hist for 30 year means. "
                                 "For daily time series, it should be a single year or a range of years within one period ('2043' or '2043-2045'). "
                                 f"hist: {hist_period['start']}-{hist_period['end']}, "
                                 f"nf: {nf_period['start']}-{nf_period['end']}, "
                                 f"ff: {ff_period['start']}-{ff_period['end']}")

            periods_config["short_name"].append(short_name)

    change_periods = Path(f"{output_path}/config/periods.tsv")
    write_config_to_tsv(periods_config, f"{change_periods}")

if scenarios != default_scenarios:
    scenarios_config = {"id": [], "description": [], "scenarioStrings": [], "hexcolour": []}
    # burde leses fra default config
    for scenario in scenarios:
        if scenario == "rcp26":
            scenarios_config["id"].append(scenario)
            scenarios_config["description"].append("Low emissions scenario (RCP2.6)")
            scenarios_config["scenarioStrings"].append("_rcp26_")
            scenarios_config["hexcolour"].append("FC8D62")
        elif scenario == "rcp45":
            scenarios_config["id"].append(scenario)
            scenarios_config["description"].append("Medium emissions scenario (RCP4.5)")
            scenarios_config["scenarioStrings"].append("_rcp45_")
            scenarios_config["hexcolour"].append("8DA0CB")
        elif scenario == "ssp370":
            scenarios_config["id"].append(scenario)
            scenarios_config["description"].append("2nd worst scenario (SSP370)")
            scenarios_config["scenarioStrings"].append("_ssp370_")
            scenarios_config["hexcolour"].append("FC8D62")
        elif scenario == "hist":
            scenarios_config["id"].append("historical")
            scenarios_config["description"].append("Historical values")
            scenarios_config["scenarioStrings"].append("_hist_")
            scenarios_config["hexcolour"].append("66C2A5")
        else:
            raise ValueError(f"Scenario {scenario} is invalid. Valid scenarios are 'rcp26', 'rcp45', 'ssp370' or 'hist'")

    change_scenarios = Path(f"{output_path}/config/scenarios.tsv")
    write_config_to_tsv(scenarios_config, f"{change_scenarios}")

if str(region_shapefile) == str(default_region_shapefile_path):
    if region_id != default_region_id:
        if region_id > 11 or region_id < 1:
            raise ValueError(f"Region id {region_id} should be between 1 and 11")
        
        region_config = {"id": [region_id], "shapefile": [region_shapefile]}
        change_region_id = Path(f"{output_path}/config/region.tsv")
        write_config_to_tsv(region_config, f"{change_region_id}")
    else:
        #  nothing? both region_id and region_shapefile are default values and in config already
        pass
else:
    # Check that region_id is valid for given shapefile
    df = gpd.read_file(region_shapefile)
    region_id_valid = False
    for col in df.columns:
        if region_id in df[col].values:
            region_id_valid = True
            break

    if not region_id_valid:
        raise ValueError(f"Region id {region_id} not found in shapefile {region_shapefile}. Please check the columns and values in the shapefile")
    
    region_config = {"id": [region_id], "shapefile": [region_shapefile]}
    change_region_shapefile = Path(f"{output_path}/config/region.tsv")
    write_config_to_tsv(region_config, f"{change_region_shapefile}")


# Update config.yaml to use new config files
default_config["configurationTables"]["seasons"] = f"{klimakverna_pilot1_path}/config/testcase_8/seasons.tsv"
set_config_path(change_indicator, "indicators")
set_config_path(change_periods, "periods")
set_config_path(change_scenarios, "scenarios")

if change_region_id:
    set_config_path(change_region_id, "region")
elif change_region_shapefile:
    set_config_path(change_region_shapefile, "region")

for model in models:
    if model != "cnrm_aladin" or bias_method != "3DBC" or indicator != "pr" or calculation_type != "30y":
        # input file name depends on models, bias_method, calcualtions_type, cmip_version (=scenario) and indicator 
        # periods are filtered afterwards
        input_config = {"id": [], "srcName": [], "varName": [], "path": [], "stemRegex": [], "internalVarName": [], "hasScenarios": [], "applyPreprocessor": []}
        stem_regex = "(.*).nc4"
        has_scenarios = "TRUE"
        apply_preprocessor = "FALSE"

        source_name = f"CMIP{model_info[model]['cmip_version']}"
        input_config["id"].append(f"CMIP{model_info[model]['cmip_version']}-{bias_method}-{calculation_type}-{model}-{indicator}")
        input_config["srcName"].append(source_name)
        input_config["varName"].append(indicator)
        input_config["path"].append(get_input_path(model, bias_method, calculation_type, indicator))
        input_config["stemRegex"].append(stem_regex)
        input_config["internalVarName"].append(indicator)
        input_config["hasScenarios"].append(has_scenarios)
        input_config["applyPreprocessor"].append(apply_preprocessor)

        change_inputs = Path(f"{output_path}/config/inputs.tsv")
        write_config_to_tsv(input_config, f"{change_inputs}")

    set_config_path(change_inputs, "inputs")

    if change_output_location:
        # Create if folder doesn't exist? think it works
        default_config["dirs"]["csv"] = str(change_output_location / f"{bias_method}/{calculation_type}/{model}")
        default_config["dirs"]["region"] = str(change_output_location / f"{bias_method}/{calculation_type}/{model}")
    else:
        default_config["dirs"]["csv"] = f"{output_path}/{bias_method}/{calculation_type}/{model}"
        default_config["dirs"]["region"] = f"{output_path}/{bias_method}/{calculation_type}/{model}"

    with open(default_config_file) as f:
        old_default_config = yaml.safe_load(f)

    old_default_config["testcase_8"] = default_config

    with open(klimakverna_pilot1_path / "config/testcase_8_config.yaml", "w") as f:
        yaml.dump(old_default_config, f)

    # Run calculation locally
    if dry_run:
        print(f"\n\nDry run for {model} with bias method {bias_method} and calculation type {calculation_type}")
        process = subprocess.Popen([f"{klimakverna_pilot1_path}/tools/dry_run_snakemake_local.sh"])
    else:
        print(f"\n\nRunning calculation for {model} with bias method {bias_method} and calculation type {calculation_type}")
        process = subprocess.Popen([f"{klimakverna_pilot1_path}/tools/run_snakemake_local.sh"])
    process.wait()

    # Run calculation in PPI queue
    # print(f"\n\nRunning calculation for {model} with bias method {bias_method} and calculation type {calculation_type}")
    # process = subprocess.Popen(["qsub", "-V", "-b", "n", "-cwd", f"{klimakverna_pilot1_path}/tools/run_snakemake_ppi_C.sh"])
    # process.wait()

# Add: possibility to choose both bias methods (save_areal mean supports  it)
# add possibility to dry run cmd option