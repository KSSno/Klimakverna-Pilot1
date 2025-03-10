config =config["testcase_hydro"]
modell = config["modell"]
input = config["input"]
nse = config["nse"]
output = config["output"]
print(output)
rule run_hbv_analysis:
    input:
        modell,
        input,
        nse

    output:
        hydro_output_var=os.path.join(output, "results/output.var"),
        hydro_output_nse=os.path.join(output, "results/nse.Rout"),

    shell:
        """
        tcsh -c '
        modell input;
        mv *var results/;
        R CMD BATCH nse
        '
        """