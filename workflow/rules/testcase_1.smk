module testcase_1:
    snakefile:
        "/lustre/storeC-ext/users/klimakverna/development/KAPy/workflow/Snakefile"
    config:config['testcase_1']

use rule * from testcase_1 as testcase_1_*
