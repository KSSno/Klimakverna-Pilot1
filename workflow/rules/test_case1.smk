module kapy:
    snakefile:
        "/lustre/storeC-ext/users/klimakverna/development/KAPy/workflow/Snakefile"
    config:config['testcase_1']

use rule * from kapy as testcase1_*
