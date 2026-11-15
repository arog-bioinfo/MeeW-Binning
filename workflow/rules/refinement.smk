if REFINEMENT_ENABLED:

    rule prepare_refinement_inputs:
        input:
            bins_dirs=refinement_binner_bins_inputs,
            assembly=lambda wc: sample_input(wc, "contigs"),
            tiara_tsv=f"{RESULTS_DIR}/domain/tiara_contigs.tsv",
        output:
            mappings=[
                f"{RESULTS_DIR}/refinement/das_tool/inputs/{tool}_prok_bins.tsv"
                for tool in SELECTED_BINNING_TOOLS
            ],
            assignments=f"{RESULTS_DIR}/refinement/bin_domain_assignments.tsv",
        log:
            f"{RESULTS_DIR}/logs/refinement/prepare_inputs/{{sample}}.log",
        params:
            tools=SELECTED_BINNING_TOOLS,
            fallback_domain=FALLBACK_DOMAIN,
        message:
            "Preparing sample-level DAS Tool mappings and domain assignments."
        script:
            "../scripts/prepare_refinement_inputs.py"

    rule refine_das_tool:
        input:
            assembly=lambda wc: sample_input(wc, "contigs"),
            contig2bin=rules.prepare_refinement_inputs.output.mappings,
            bin_domain_assignments=rules.prepare_refinement_inputs.output.assignments,
        output:
            bins=directory(f"{RESULTS_DIR}/refinement/das_tool/bins"),
            summary=f"{RESULTS_DIR}/refinement/das_tool/DASTool_summary.tsv",
            scores=f"{RESULTS_DIR}/refinement/das_tool/all_bins.eval.tsv",
            contig2bin=f"{RESULTS_DIR}/refinement/das_tool/DASTool_contig2bin.tsv",
            provenance=f"{RESULTS_DIR}/refinement/das_tool/input_provenance.tsv",
            run_status=f"{RESULTS_DIR}/refinement/das_tool/run_status.tsv",
            tool_log=f"{RESULTS_DIR}/refinement/das_tool/DAS_Tool.log",
        log:
            f"{RESULTS_DIR}/logs/refinement/das_tool/{{sample}}.log",
        benchmark:
            f"{RESULTS_DIR}/benchmarks/refinement/das_tool/{{sample}}.tsv"
        conda:
            "../envs/das_tool.yaml"
        threads: DASTOOL_THREADS
        resources:
            mem_mb=16000,
            disk_mb=50000,
            runtime=1440,
        params:
            labels=SELECTED_BINNING_TOOLS,
            score_threshold=DASTOOL_SCORE_THRESHOLD,
            search_engine="diamond",
            expected_version="1.1.7",
        message:
            "Refining prokaryotic bin inputs with DAS Tool 1.1.7."
        script:
            "../scripts/run_das_tool.py"
