#!/usr/bin/env bash
set -euo pipefail

# Benchmark persistent sequence insertion against a full genomeGenerate rebuild.
#
# Required:
#   BASE_GENOME_DIR=/path/to/prebuilt/index
#   INSERT_FASTA_FILES="/path/customer1.fa /path/customer2.fa"
#
# To run the full rebuild comparison, also provide:
#   BASE_FASTA_FILES="/path/reference.fa /path/ercc.fa"
#   SJDB_GTF_FILE=/path/combined.gtf # optional full rebuild GTF
#
# To benchmark inserted-sequence annotations, also provide:
#   INSERT_GTF_FILE=/path/inserted_sequences.gtf
#
# Output:
#   config.tsv, summary.tsv, *.time, *.command.txt, iostat snapshots, and
#   validation.tsv when both incremental and full-rebuild outputs are present.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"

star_bin="${STAR_BIN:-${repo_root}/source/STAR}"
time_bin="${TIME_BIN:-/usr/bin/time}"
out_root="${OUT_DIR:-/tmp/star-genome-insert-bench.$(date +%Y%m%d_%H%M%S).$$}"
threads="${THREADS:-$(nproc 2>/dev/null || echo 1)}"
base_genome_dir="${BASE_GENOME_DIR:-}"
insert_fasta_files_input="${INSERT_FASTA_FILES:-}"
base_fasta_files_input="${BASE_FASTA_FILES:-}"
run_full_rebuild="${RUN_FULL_REBUILD:-auto}"
sjdb_gtf_file="${SJDB_GTF_FILE:-}"
insert_gtf_file="${INSERT_GTF_FILE:-}"
sjdb_overhang="${SJDB_OVERHANG:-}"
genome_saindex_nbases="${GENOME_SAINDEX_NBASES:-14}"
genome_chrbin_nbits="${GENOME_CHRBIN_NBITS:-18}"
limit_genome_ram="${LIMIT_GENOME_GENERATE_RAM:-300000000000}"
extra_full_args="${EXTRA_FULL_REBUILD_ARGS:-}"
extra_insert_args="${EXTRA_GENOME_INSERT_ARGS:-}"
iostat_bin="${IOSTAT_BIN:-$(command -v iostat || true)}"
iostat_interval="${IOSTAT_INTERVAL:-1}"
iostat_samples="${IOSTAT_SAMPLES:-3}"

if [[ ! -x "${star_bin}" ]]; then
    echo "ERROR: STAR_BIN is not executable: ${star_bin}" >&2
    exit 1
fi

if [[ ! -x "${time_bin}" ]]; then
    echo "ERROR: TIME_BIN is not executable: ${time_bin}" >&2
    exit 1
fi

if [[ -z "${base_genome_dir}" || ! -d "${base_genome_dir}" ]]; then
    echo "ERROR: BASE_GENOME_DIR must point to an existing STAR genome directory" >&2
    exit 1
fi

if [[ -z "${insert_fasta_files_input}" ]]; then
    echo "ERROR: INSERT_FASTA_FILES is required" >&2
    exit 1
fi

read -r -a insert_fasta_files <<< "${insert_fasta_files_input}"
read -r -a base_fasta_files <<< "${base_fasta_files_input}"

for fasta_file in "${insert_fasta_files[@]}" "${base_fasta_files[@]}"; do
    if [[ -n "${fasta_file}" && ! -r "${fasta_file}" ]]; then
        echo "ERROR: FASTA file is not readable: ${fasta_file}" >&2
        exit 1
    fi
done

if [[ -n "${sjdb_gtf_file}" && ! -r "${sjdb_gtf_file}" ]]; then
    echo "ERROR: SJDB_GTF_FILE is not readable: ${sjdb_gtf_file}" >&2
    exit 1
fi

if [[ -n "${insert_gtf_file}" && ! -r "${insert_gtf_file}" ]]; then
    echo "ERROR: INSERT_GTF_FILE is not readable: ${insert_gtf_file}" >&2
    exit 1
fi

if [[ "${run_full_rebuild}" == "auto" ]]; then
    if [[ "${#base_fasta_files[@]}" -gt 0 && -n "${base_fasta_files[0]:-}" ]]; then
        run_full_rebuild=1
    else
        run_full_rebuild=0
    fi
fi

mkdir -p "${out_root}"
out_root="$(cd "${out_root}" && pwd)"

summary="${out_root}/summary.tsv"
validation="${out_root}/validation.tsv"
incremental_index="${out_root}/incremental_index"
full_rebuild_index="${out_root}/full_rebuild_index"

snapshot_iostat() {
    local label="$1"
    local out_file="${out_root}/${label}.iostat"
    if [[ -n "${iostat_bin}" && -x "${iostat_bin}" ]]; then
        "${iostat_bin}" -x "${iostat_interval}" "${iostat_samples}" > "${out_file}" 2>&1 || true
    else
        printf "iostat unavailable\n" > "${out_file}"
    fi
}

write_command() {
    local out_file="$1"
    shift
    printf '%q ' "$@" > "${out_file}"
    printf '\n' >> "${out_file}"
}

parse_time() {
    local time_file="$1"
    awk '
        /User time/ {user=$NF}
        /System time/ {sys=$NF}
        /Percent of CPU/ {cpu=$NF; gsub(/%/,"",cpu)}
        /Elapsed/ {wall=$NF}
        /Maximum resident/ {rss=$NF}
        /File system inputs/ {fsin=$NF}
        /File system outputs/ {fsout=$NF}
        END {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", wall, user, sys, cpu, rss, fsin, fsout}
    ' "${time_file}"
}

{
    printf 'key\tvalue\n'
    printf 'STAR_BIN\t%s\n' "${star_bin}"
    printf 'STAR_VERSION\t%s\n' "$("${star_bin}" --version)"
    printf 'THREADS\t%s\n' "${threads}"
    printf 'BASE_GENOME_DIR\t%s\n' "${base_genome_dir}"
    printf 'BASE_FASTA_FILES\t%s\n' "${base_fasta_files_input:-NA}"
    printf 'INSERT_FASTA_FILES\t%s\n' "${insert_fasta_files_input}"
    printf 'INSERT_GTF_FILE\t%s\n' "${insert_gtf_file:-NA}"
    printf 'RUN_FULL_REBUILD\t%s\n' "${run_full_rebuild}"
    printf 'SJDB_GTF_FILE\t%s\n' "${sjdb_gtf_file:-NA}"
    printf 'SJDB_OVERHANG\t%s\n' "${sjdb_overhang:-NA}"
    printf 'GENOME_SAINDEX_NBASES\t%s\n' "${genome_saindex_nbases}"
    printf 'GENOME_CHRBIN_NBITS\t%s\n' "${genome_chrbin_nbits}"
    printf 'LIMIT_GENOME_GENERATE_RAM\t%s\n' "${limit_genome_ram}"
} > "${out_root}/config.tsv"

printf 'scenario\twall_time\tuser_seconds\tsystem_seconds\tcpu_percent\tmax_rss_kb\tfs_inputs\tfs_outputs\tindex_dir\ttime_file\tcommand_file\n' > "${summary}"

insert_cmd=(
    "${star_bin}"
    --runMode genomeInsert
    --runThreadN "${threads}"
    --genomeDir "${base_genome_dir}"
    --genomeFastaFiles "${insert_fasta_files[@]}"
    --genomeInsertOutDir "${incremental_index}"
    --outFileNamePrefix "${out_root}/incremental_"
)
if [[ -n "${insert_gtf_file}" ]]; then
    insert_cmd+=(--sjdbGTFfile "${insert_gtf_file}")
fi
if [[ -n "${extra_insert_args}" ]]; then
    read -r -a extra_insert_args_array <<< "${extra_insert_args}"
    insert_cmd+=("${extra_insert_args_array[@]}")
fi

write_command "${out_root}/incremental.command.txt" "${insert_cmd[@]}"
snapshot_iostat "before_incremental"
"${time_bin}" -v -o "${out_root}/incremental.time" "${insert_cmd[@]}" > "${out_root}/incremental.stdout" 2>&1
snapshot_iostat "after_incremental"
read -r wall user sys cpu rss fsin fsout < <(parse_time "${out_root}/incremental.time")
printf 'genomeInsert\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${wall}" "${user}" "${sys}" "${cpu}" "${rss}" "${fsin}" "${fsout}" \
    "${incremental_index}" "${out_root}/incremental.time" "${out_root}/incremental.command.txt" >> "${summary}"

if [[ "${run_full_rebuild}" == "1" ]]; then
    full_cmd=(
        "${star_bin}"
        --runMode genomeGenerate
        --runThreadN "${threads}"
        --genomeDir "${full_rebuild_index}"
        --genomeFastaFiles "${base_fasta_files[@]}" "${insert_fasta_files[@]}"
        --limitGenomeGenerateRAM "${limit_genome_ram}"
        --genomeSAindexNbases "${genome_saindex_nbases}"
        --genomeChrBinNbits "${genome_chrbin_nbits}"
        --outFileNamePrefix "${out_root}/full_rebuild_"
    )
    if [[ -n "${sjdb_gtf_file}" ]]; then
        full_cmd+=(--sjdbGTFfile "${sjdb_gtf_file}")
    fi
    if [[ -n "${sjdb_overhang}" ]]; then
        full_cmd+=(--sjdbOverhang "${sjdb_overhang}")
    fi
    if [[ -n "${extra_full_args}" ]]; then
        read -r -a extra_full_args_array <<< "${extra_full_args}"
        full_cmd+=("${extra_full_args_array[@]}")
    fi

    write_command "${out_root}/full_rebuild.command.txt" "${full_cmd[@]}"
    snapshot_iostat "before_full_rebuild"
    "${time_bin}" -v -o "${out_root}/full_rebuild.time" "${full_cmd[@]}" > "${out_root}/full_rebuild.stdout" 2>&1
    snapshot_iostat "after_full_rebuild"
    read -r wall user sys cpu rss fsin fsout < <(parse_time "${out_root}/full_rebuild.time")
    printf 'full_rebuild\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${wall}" "${user}" "${sys}" "${cpu}" "${rss}" "${fsin}" "${fsout}" \
        "${full_rebuild_index}" "${out_root}/full_rebuild.time" "${out_root}/full_rebuild.command.txt" >> "${summary}"

    printf 'file\tstatus\n' > "${validation}"
    for file in Genome SAindex chrName.txt chrStart.txt chrLength.txt chrNameLength.txt sjdbInfo.txt sjdbList.out.tab sjdbList.fromGTF.out.tab exonInfo.tab exonGeTrInfo.tab geneInfo.tab transcriptInfo.tab; do
        if [[ -e "${incremental_index}/${file}" || -e "${full_rebuild_index}/${file}" ]]; then
            if cmp -s "${incremental_index}/${file}" "${full_rebuild_index}/${file}"; then
                printf '%s\tpass\n' "${file}" >> "${validation}"
            else
                printf '%s\tfail\n' "${file}" >> "${validation}"
            fi
        fi
    done
    if cmp -s "${incremental_index}/SA" "${full_rebuild_index}/SA"; then
        printf 'SA\tbyte-identical\n' >> "${validation}"
    else
        printf 'SA\tdifferent-equivalent-ordering\n' >> "${validation}"
    fi
fi

cat "${summary}"
