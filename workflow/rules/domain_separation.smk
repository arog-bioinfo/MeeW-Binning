rule tiara_predict:
    input:
        contigs=lambda wc: sample_input(wc, "contigs"),
    output:
        tsv=f"{RESULTS_DIR}/domain/tiara_contigs.tsv",
    log:
        f"{RESULTS_DIR}/logs/domain/tiara/{{sample}}.log",
    conda:
        "../envs/tiara.yaml"
    threads: 4
    params:
        extra=config["domain_separation"]["tiara"].get("extra", ""),
    shell:
        """
        mkdir -p "$(dirname {output.tsv:q})"
        tiara -i {input.contigs:q} -o {output.tsv:q} --threads {threads} {params.extra} >{log:q} 2>&1
        """


if REFINEMENT_ENABLED:

    rule separate_bins_by_domain:
        input:
            prok_bins=rules.refine_das_tool.output.bins,
            source_bins=refinement_binner_bins_inputs,
            bin_domain_assignments=rules.prepare_refinement_inputs.output.assignments,
            tiara_tsv=f"{RESULTS_DIR}/domain/tiara_contigs.tsv",
            das_tool_status=rules.refine_das_tool.output.run_status,
        output:
            prok_bins=directory(f"{RESULTS_DIR}/domain/prok_bins"),
            euk_bins=directory(f"{RESULTS_DIR}/domain/euk_bins"),
            unknown_bins=directory(f"{RESULTS_DIR}/domain/unknown_bins"),
            summary=f"{RESULTS_DIR}/domain/domain_assignment.tsv",
            annotation_samples=f"{RESULTS_DIR}/domain/annotation_samples.tsv",
            rename_map=f"{RESULTS_DIR}/domain/bin_rename_map.tsv",
        log:
            f"{RESULTS_DIR}/logs/domain/separate_bins_by_domain/{{sample}}.log",
        params:
            tools=SELECTED_BINNING_TOOLS,
            euk_source_tool=EUK_SOURCE_TOOL,
            fallback_domain=FALLBACK_DOMAIN,
        message:
            "Materializing DAS Tool prokaryotic bins and source-binner eukaryotic and unknown bins."
        script:
            "../scripts/finalize_refined_domains.py"

else:

    rule separate_bins_by_domain:
        input:
            binner_output=lambda wc: selected_binner_tracked_output(wc.sample),
            bins_dir=lambda wc: (
                selected_binner_bins_dir(wc.sample)
                if SELECTED_BINNING_TOOLS[0] in {"comebin", "semibin2", "metabat2"}
                else []
            ),
            tiara_tsv=f"{RESULTS_DIR}/domain/tiara_contigs.tsv",
        output:
            prok_bins=directory(f"{RESULTS_DIR}/domain/prok_bins"),
            euk_bins=directory(f"{RESULTS_DIR}/domain/euk_bins"),
            unknown_bins=directory(f"{RESULTS_DIR}/domain/unknown_bins"),
            summary=f"{RESULTS_DIR}/domain/domain_assignment.tsv",
            annotation_samples=f"{RESULTS_DIR}/domain/annotation_samples.tsv",
            rename_map=f"{RESULTS_DIR}/domain/bin_rename_map.tsv",
        log:
            f"{RESULTS_DIR}/logs/domain/separate_bins_by_domain/{{sample}}.log",
        container:
            "docker://python:3.11-slim"
        params:
            bins_dir=lambda wc, input: (
                str(input.bins_dir)
                if input.bins_dir
                else selected_binner_bins_dir(wc.sample)
            ),
            fallback_domain=FALLBACK_DOMAIN,
        message:
            "Separating bins by Tiara domain calls using majority-bp rule."
        script:
            "../scripts/separate_bins_by_domain.py"
