from pathlib import Path


rule comebin:
    input:
        contigs=lambda wc: sample_input(wc, "contigs"),
        bam=lambda wc: sample_input(wc, "bam"),
        atomic_script=workflow.source_path("../scripts/run_comebin_atomic.sh"),
        checkm1_runner=CHECKM1_RUNNER,
        checkm1_helper=CHECKM1_HELPER,
    output:
        result_tsv=f"{RESULTS_DIR}/binning/comebin/comebin_res/comebin_res.tsv",
        bins=directory(f"{RESULTS_DIR}/binning/comebin/comebin_res/comebin_res_bins"),
    log:
        f"{RESULTS_DIR}/logs/binning/comebin/{{sample}}.log",
    conda:
        "../envs/comebin.yaml"
    threads: 8
    params:
        output_dir=lambda wc, output: str(Path(output.result_tsv).parents[1]),
        bam_dir=lambda wc, input: comebin_bam_dir(input.bam),
        checkm1_database=CHECKM1_DATABASE_PATH,
        checkm1_download_database=CHECKM1_DOWNLOAD_DATABASE,
        extra=config["binning"]["comebin"].get("extra", ""),
    shell:
        """
        bash {input.atomic_script:q} \
            --checkm-runner {input.checkm1_runner:q} \
            --checkm-helper {input.checkm1_helper:q} \
            --database {params.checkm1_database:q} \
            --download-database {params.checkm1_download_database:q} \
            --source-assembly {input.contigs:q} \
            --bam {input.bam:q} \
            --bam-dir {params.bam_dir:q} \
            --final-output {params.output_dir:q} \
            --threads {threads} \
            -- {params.extra} >{log:q} 2>&1
        """


rule metabinner_depth:
    input:
        bams=bam_inputs,
    output:
        depth=f"{RESULTS_DIR}/binning/metabinner/profiles/jgi_depth.txt",
    log:
        f"{RESULTS_DIR}/logs/binning/metabinner/{{sample}}.depth.log",
    conda:
        "../envs/metabinner.yaml"
    threads: 1
    params:
        output_dir=lambda wc, output: str(Path(output.depth).parent),
    shell:
        """
        mkdir -p {params.output_dir:q}
        if ! command -v jgi_summarize_bam_contig_depths >/dev/null 2>&1; then
            echo "ERROR: jgi_summarize_bam_contig_depths not found on PATH"
            exit 127
        fi
        jgi_summarize_bam_contig_depths --outputDepth {output.depth:q} --noIntraDepthVariance {input.bams:q} >{log:q} 2>&1
        """


rule metabinner_profiles:
    input:
        contigs=lambda wc: sample_input(wc, "contigs"),
        depth=rules.metabinner_depth.output.depth,
        profile_script=workflow.source_path("../scripts/prepare_metabinner_profiles.sh"),
    output:
        contigs=f"{RESULTS_DIR}/binning/metabinner/profiles/filtered_contigs.fa",
        coverage=f"{RESULTS_DIR}/binning/metabinner/profiles/coverage_profile.tsv",
        composition=f"{RESULTS_DIR}/binning/metabinner/profiles/tetranucleotide_profile.csv",
    log:
        f"{RESULTS_DIR}/logs/binning/metabinner/{{sample}}.profiles.log",
    conda:
        "../envs/metabinner.yaml"
    threads: 1
    params:
        output_dir=lambda wc, output: str(Path(output.contigs).parent),
        min_contig_length=METABINNER_CONFIG.get("min_contig_length", 1000),
    shell:
        """
        bash {input.profile_script:q} \
            --contigs {input.contigs:q} \
            --depth {input.depth:q} \
            --output-dir {params.output_dir:q} \
            --output-contigs {output.contigs:q} \
            --output-coverage {output.coverage:q} \
            --output-composition {output.composition:q} \
            --min-contig-length {params.min_contig_length} >{log:q} 2>&1
        """


rule metabinner:
    input:
        contigs=rules.metabinner_profiles.output.contigs,
        depth=rules.metabinner_profiles.output.coverage,
        kmer=rules.metabinner_profiles.output.composition,
        checkm1_runner=CHECKM1_RUNNER,
        checkm1_helper=CHECKM1_HELPER,
    output:
        result_tsv=f"{RESULTS_DIR}/binning/metabinner/upstream/metabinner_res/metabinner_result.tsv",
    log:
        f"{RESULTS_DIR}/logs/binning/metabinner/{{sample}}.log",
    conda:
        "../envs/metabinner.yaml"
    threads: 8
    params:
        output_root=lambda wc, output: str(Path(output.result_tsv).parents[1]),
        install_path=METABINNER_CONFIG.get("install_path", ""),
        scale=METABINNER_CONFIG.get("scale", "large"),
        checkm1_database=CHECKM1_DATABASE_PATH,
        checkm1_download_database=CHECKM1_DOWNLOAD_DATABASE,
    shell:
        """
        (
            configured_install={params.install_path:q}
            ambient_install="$(printenv METABINNER_HOME || true)"
            if [[ -z "$configured_install" && -n "$ambient_install" ]]; then
                configured_install="$ambient_install"
            fi

            if [[ -n "$configured_install" ]]; then
                if [[ -f "$configured_install/run_metabinner.sh" ]]; then
                    install_root="$(readlink -f "$configured_install")"
                elif [[ -f "$configured_install/bin/run_metabinner.sh" ]]; then
                    install_root="$(readlink -f "$configured_install/bin")"
                else
                    echo "ERROR: MetaBinner install_path/METABINNER_HOME must be a tool root, Bioconda bin directory, or Conda prefix"
                    exit 1
                fi
                runner="$install_root/run_metabinner.sh"
            else
                runner="$(command -v run_metabinner.sh || true)"
                if [[ -z "$runner" ]]; then
                    echo "ERROR: run_metabinner.sh not found on PATH"
                    exit 127
                fi
                runner="$(readlink -f "$runner")"
                install_root="$(dirname "$runner")"
            fi

            if [[ ! -f "$runner" || ! -d "$install_root/scripts" || ! -d "$install_root/auxiliary" ]]; then
                echo "ERROR: incomplete MetaBinner installation at $install_root"
                exit 1
            fi
            contigs="$(readlink -f {input.contigs:q})"
            coverage="$(readlink -f {input.depth:q})"
            composition="$(readlink -f {input.kmer:q})"

            output_root={params.output_root:q}
            output_parent="$(dirname "$output_root")"
            output_name="$(basename "$output_root")"
            mkdir -p "$output_parent"
            absolute_parent="$(readlink -f "$output_parent")"
            output_root="$absolute_parent/$output_name"
            result_tsv="$output_root/metabinner_res/metabinner_result.tsv"

            python {input.checkm1_runner:q} \
                --helper-source {input.checkm1_helper:q} \
                --database {params.checkm1_database:q} \
                --download-database {params.checkm1_download_database:q} \
                --tool metabinner \
                --tool-install-root "$install_root" \
                --protected-path "$contigs" \
                --protected-path "$coverage" \
                --protected-path "$composition" \
                --native-path "$contigs" \
                --native-path "$coverage" \
                --native-path "$composition" \
                --native-path "$output_root" \
                --cleanup-path "$output_root" \
                -- bash "$runner" \
                -a "$contigs" \
                -o "$output_root" \
                -d "$coverage" \
                -k "$composition" \
                -p "$install_root" \
                -t {threads} \
                -s {params.scale:q}

            if [[ ! -f "$result_tsv" ]]; then
                echo "ERROR: MetaBinner did not produce metabinner_result.tsv"
                exit 1
            fi
        ) >{log:q} 2>&1
        """


rule metabinner_bins:
    input:
        contigs=rules.metabinner_profiles.output.contigs,
        result=rules.metabinner.output.result_tsv,
    output:
        bins=directory(f"{RESULTS_DIR}/binning/metabinner/bins"),
    log:
        f"{RESULTS_DIR}/logs/binning/metabinner/{{sample}}.materialize.log",
    conda:
        "../envs/metabinner.yaml"
    threads: 1
    shell:
        """
        (
            if ! command -v gen_bins_from_tsv.py >/dev/null 2>&1; then
                echo "ERROR: gen_bins_from_tsv.py not found on PATH"
                exit 127
            fi
            rm -rf {output.bins:q}
            gen_bins_from_tsv.py \
                -f {input.contigs:q} \
                -r {input.result:q} \
                -o {output.bins:q}
            test -d {output.bins:q}
        ) >{log:q} 2>&1
        """


rule semibin2:
    input:
        contigs=lambda wc: sample_input(wc, "contigs"),
        bams=bam_inputs,
    output:
        bins_info=f"{RESULTS_DIR}/binning/semibin2/bins_info.tsv",
        bins=directory(f"{RESULTS_DIR}/binning/semibin2/output_bins"),
    log:
        f"{RESULTS_DIR}/logs/binning/semibin2/{{sample}}.log",
    conda:
        "../envs/semibin2.yaml"
    threads: 8
    params:
        output_dir=lambda wc, output: str(Path(output.bins_info).parent),
        environment=config["binning"]["semibin2"].get("environment", "global"),
        sequencing_type=config["binning"]["semibin2"].get(
            "sequencing_type", "short_read"
        ),
        extra=config["binning"]["semibin2"].get("extra", ""),
    script:
        "../scripts/semibin2.py"


rule metabat2_depth:
    input:
        bams=bam_inputs,
    output:
        depth=f"{RESULTS_DIR}/binning/metabat2/depth.txt",
    log:
        f"{RESULTS_DIR}/logs/binning/metabat2/{{sample}}.depth.log",
    conda:
        "../envs/metabat2.yaml"
    threads: 1
    params:
        output_dir=lambda wc, output: str(Path(output.depth).parent),
    shell:
        """
        mkdir -p {params.output_dir:q}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth:q} {input.bams:q} >{log:q} 2>&1
        """


rule metabat2:
    input:
        contigs=lambda wc: sample_input(wc, "contigs"),
        depth=rules.metabat2_depth.output.depth,
    output:
        bins=directory(f"{RESULTS_DIR}/binning/metabat2/bins"),
    log:
        f"{RESULTS_DIR}/logs/binning/metabat2/{{sample}}.log",
    conda:
        "../envs/metabat2.yaml"
    threads: 8
    params:
        prefix=lambda wc, output: str(Path(output.bins) / "bin"),
        min_contig_length=config["binning"]["metabat2"].get("min_contig_length", 1500),
        seed=config["binning"]["metabat2"].get("seed", 1),
        extra=config["binning"]["metabat2"].get("extra", ""),
    shell:
        """
        mkdir -p {output.bins:q}
        metabat2 -i {input.contigs:q} -a {input.depth:q} -o {params.prefix:q} -m {params.min_contig_length} -t {threads} --seed {params.seed} --saveCls {params.extra} >{log:q} 2>&1
        """
