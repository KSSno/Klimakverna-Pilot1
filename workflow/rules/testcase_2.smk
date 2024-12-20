rule validate_and_move:
    input:
        f"{INPUT_BASE}{{filename}}.nc4",
    output:
        xml_output=os.path.join(MMD_OUT_DIR, "{filename}.xml"),
        nc_output=os.path.join(NC_OUT_DIR, "{filename}.nc4"),
    shell:
        """
        logdir="logs/{wildcards.filename}.txt"
        
        mkdir -p "$(dirname "$logdir")"

        # Get the directory path by removing the last component from the wildcards.filename
        base_dir=$(dirname "{wildcards.filename}")
        nc_dir="{NC_OUT_DIR}/$base_dir"
        mmd_dir="{MMD_OUT_DIR}/$base_dir"
        
        output=$(nc2mmd -i {input} -u "test" -o "$mmd_dir" 2>/dev/null)

        # Check if the output is not empty
        if [[ -n "$output" && "$output" =~ [^[:space:]] ]]; then
            echo "$output" > "$logdir"
        else
            # If no output, copy the input file to the output directory
            cp {input} "$nc_dir"
        fi

        """
