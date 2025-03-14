config =config["testcase_hydro"]
modell = config["modell"]
control_text = config["control_hbv"]
nse = config["nse"]
output = config["output"]
print(output,modell)
rule run_hbv_analysis:
    input:
        modell,
        control_text,
        nse

    output:
        hydro_output_var=os.path.join(output, "results/output.var"),
        hydro_output_nse=os.path.join(output, "results/nse.Rout"),

    shell:
        """
        {modell} {control_text};
        mv *var results/;
        R CMD BATCH {nse}
        """
