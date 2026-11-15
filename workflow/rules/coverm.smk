wildcard_constraints:
    domain="prok|euk",


rule coverm_genome_coverage:
    input:
        unpack(coverm_inputs),
    output:
        f"{RESULTS_DIR}/coverm/{{domain}}_genome_coverage.tsv",
    log:
        f"{RESULTS_DIR}/logs/coverm/{{sample}}/{{domain}}.log",
    conda:
        "../envs/coverm.yaml"
    threads: COVERM_THREADS
    params:
        input_mode=COVERM_INPUT_MODE,
        mapper=COVERM_MAPPER,
        methods=COVERM_METHODS,
        extension=COVERM_GENOME_FASTA_EXTENSION,
        extra=COVERM_EXTRA,
    script:
        "../scripts/coverm.py"
