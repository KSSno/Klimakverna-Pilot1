rule print:
    output:
        "output/print_done.txt",
    shell:
        "echo 'here' > {output}"
