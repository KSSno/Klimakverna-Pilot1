from tools import validate_config, get_output_files

config = validate_config(config["configuration"])
outputs = get_output_files(config)
output_base = config["output_base"]

# Cutting out data for region ------------------------------------------------------
if "region" in outputs and outputs['region']:
    rule save_areal_mean_for_region_to_netcdf:
        input:
            lambda wildcards: outputs['region'][os.path.join(outDirs["region"], f'{wildcards.indicator_id}_{wildcards.scenario}_{wildcards.period}_region_{wildcards.region}.nc')]
        output:
            os.path.join(output_base, '{indicator_id}_{scenario}_{period}_region_{region}.nc')
        run:
            KAPy.save_areal_mean_for_region_to_netcdf(config, wildcards.indicator_id, wildcards.scenario, wildcards.period, wildcards.region, input, str(output))

# Output as CSV for region timeseries or region 30 year means -----------------------                  
if "csv" in outputs and outputs["csv"]:
    rule timeseries_netcdf_to_csv:
        input:
            lambda wildcards: outputs['csv'][os.path.join(outDirs["csv"], f'{wildcards.indID}_region_{wildcards.region}.csv')]
        output:
            os.path.join(output_base, '{indicator_id}_region_{region}.csv')
        run:
            KAPy.timeseries_netcdf_to_csv(config, wildcards.indicator_id, wildcards.region, input, str(output))

# rule all --------------------------------------------------------------------------
rule all:
    input:
        outputs['all']
    default_target: True