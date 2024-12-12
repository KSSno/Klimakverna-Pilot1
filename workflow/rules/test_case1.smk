module kapy:
    snakefile:
        "/home/shamlym/workspace/klima-kverna/KAPy/workflow/Snakefile"
    config:config['testcase_1']

use rule * from kapy as testcase_1_*
