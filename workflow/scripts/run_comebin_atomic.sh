#!/usr/bin/env bash

set -euo pipefail

usage() {
	echo "Usage: $0 --checkm-runner PATH --checkm-helper PATH --database PATH --download-database BOOL --source-assembly PATH --bam PATH --bam-dir PATH --final-output PATH --threads N [-- EXTRA...]" >&2
}

checkm_runner=""
checkm_helper=""
database=""
download_database=""
source_assembly=""
bam=""
bam_dir=""
final_output=""
threads=""
extra=()

while (($# > 0)); do
	case "$1" in
	--checkm-runner)
		checkm_runner="${2-}"
		shift 2
		;;
	--checkm-helper)
		checkm_helper="${2-}"
		shift 2
		;;
	--database)
		database="${2-}"
		shift 2
		;;
	--download-database)
		download_database="${2-}"
		shift 2
		;;
	--source-assembly)
		source_assembly="${2-}"
		shift 2
		;;
	--bam)
		bam="${2-}"
		shift 2
		;;
	--bam-dir)
		bam_dir="${2-}"
		shift 2
		;;
	--final-output)
		final_output="${2-}"
		shift 2
		;;
	--threads)
		threads="${2-}"
		shift 2
		;;
	--)
		shift
		extra=("$@")
		break
		;;
	*)
		echo "ERROR: unknown argument: $1" >&2
		usage
		exit 2
		;;
	esac
done

if [[ -z "$checkm_runner" || -z "$checkm_helper" || -z "$database" || -z "$download_database" || -z "$source_assembly" || -z "$bam" || -z "$bam_dir" || -z "$final_output" || -z "$threads" ]]; then
	usage
	exit 2
fi

if ! command -v run_comebin.sh >/dev/null 2>&1; then
	echo "ERROR: run_comebin.sh not found on PATH"
	exit 127
fi

final_parent="$(dirname "$final_output")"
final_name="$(basename "$final_output")"
mkdir -p "$final_parent"

workspace=""
cleanup() {
	status=$?
	trap - EXIT HUP INT TERM
	set +e
	if [[ -n "$workspace" && -d "$workspace" ]]; then
		rm -rf -- "$workspace"
		cleanup_status=$?
		if ((cleanup_status != 0)); then
			echo "ERROR: failed to clean temporary COMEBin workspace: $workspace"
			if ((status == 0)); then
				status=$cleanup_status
			fi
		fi
	fi
	exit "$status"
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

workspace="$(mktemp -d "$final_parent/.${final_name}.tmp.XXXXXX")"
assembly_dir="$workspace/assembly"
workspace_output="$workspace/output"
local_assembly="$assembly_dir/$(basename "$source_assembly")"
mkdir -p "$assembly_dir" "$workspace_output"

if ln -L -- "$source_assembly" "$local_assembly"; then
	echo "INFO: hard-linked assembly into temporary COMEBin workspace"
else
	echo "INFO: hard link unavailable; copying assembly into temporary COMEBin workspace"
	rm -f -- "$local_assembly"
	cp -L -- "$source_assembly" "$local_assembly"
fi
if [[ -L "$local_assembly" ]]; then
	echo "ERROR: workspace-local COMEBin assembly must not be a symlink"
	exit 1
fi

python "$checkm_runner" \
	--helper-source "$checkm_helper" \
	--database "$database" \
	--download-database "$download_database" \
	--protected-path "$source_assembly" \
	--protected-path "$bam" \
	--protected-path "$final_output" \
	--protected-path "$workspace" \
	-- run_comebin.sh \
	-a "$local_assembly" \
	-p "$bam_dir" \
	-t "$threads" \
	-o "$workspace_output" \
	"${extra[@]}"

workspace_result="$workspace_output/comebin_res/comebin_res.tsv"
workspace_bins="$workspace_output/comebin_res/comebin_res_bins"
if [[ ! -f "$workspace_result" ]]; then
	echo "ERROR: COMEBin did not produce comebin_res/comebin_res.tsv"
	exit 1
fi
if [[ ! -d "$workspace_bins" ]]; then
	echo "ERROR: COMEBin did not produce comebin_res/comebin_res_bins"
	exit 1
fi

rm -rf -- "$final_output"
mv -- "$workspace_output" "$final_output"
