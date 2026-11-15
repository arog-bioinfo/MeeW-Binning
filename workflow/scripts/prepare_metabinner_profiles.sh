#!/usr/bin/env bash

set -euo pipefail

usage() {
	echo "Usage: $0 --contigs PATH --depth PATH --output-dir PATH --output-contigs PATH --output-coverage PATH --output-composition PATH --min-contig-length N" >&2
}

contigs=""
depth=""
output_dir=""
output_contigs=""
output_coverage=""
output_composition=""
min_contig_length=""

while (($# > 0)); do
	case "$1" in
	--contigs)
		contigs="${2-}"
		shift 2
		;;
	--depth)
		depth="${2-}"
		shift 2
		;;
	--output-dir)
		output_dir="${2-}"
		shift 2
		;;
	--output-contigs)
		output_contigs="${2-}"
		shift 2
		;;
	--output-coverage)
		output_coverage="${2-}"
		shift 2
		;;
	--output-composition)
		output_composition="${2-}"
		shift 2
		;;
	--min-contig-length)
		min_contig_length="${2-}"
		shift 2
		;;
	*)
		echo "ERROR: unknown argument: $1" >&2
		usage
		exit 2
		;;
	esac
done

if [[ -z "$contigs" || -z "$depth" || -z "$output_dir" || -z "$output_contigs" || -z "$output_coverage" || -z "$output_composition" || -z "$min_contig_length" ]]; then
	usage
	exit 2
fi

if ! command -v Filter_tooshort.py >/dev/null 2>&1; then
	echo "ERROR: Filter_tooshort.py not found on PATH"
	exit 127
fi
if ! command -v gen_kmer.py >/dev/null 2>&1; then
	echo "ERROR: gen_kmer.py not found on PATH"
	exit 127
fi

mkdir -p "$output_dir"
profile_dir="$(readlink -f "$output_dir")"
source_contigs="$profile_dir/metabinner_input.fa"
generated_contigs="$profile_dir/metabinner_input_${min_contig_length}.fa"
length_file="$profile_dir/contig_length_filter${min_contig_length}.txt"
generated_kmer="$profile_dir/filtered_contigs_kmer_4_f0.csv"

rm -f \
	"$source_contigs" \
	"$generated_contigs" \
	"$length_file" \
	"$generated_kmer" \
	"$output_contigs" \
	"$output_coverage" \
	"$output_composition"

cp "$contigs" "$source_contigs"
Filter_tooshort.py "$source_contigs" "$min_contig_length"
mv "$generated_contigs" "$output_contigs"

awk -F '\t' -v OFS='\t' -v lengths="$length_file" '
    BEGIN {
        while ((getline < lengths) > 0) {
            keep[$1] = 1
        }
        close(lengths)
    }
    FNR == 1 || ($1 in keep) {
        printf "%s", $1
        for (field = 4; field <= NF; field++) {
            printf "%s%s", OFS, $field
        }
        printf "\n"
    }
' "$depth" >"$output_coverage"

gen_kmer.py "$output_contigs" 0 4
mv "$generated_kmer" "$output_composition"
rm -f "$source_contigs" "$length_file"
