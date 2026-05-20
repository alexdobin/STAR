#!/usr/bin/env bash
set -euo pipefail

# Genome-generation benchmark and idempotence harness.
#
# Typical use after a genomeGenerate optimization:
#   make -C source STAR
#   BASELINE_STAR_BIN=/path/to/baseline/STAR THREADS="1 4 8 16" \
#       extras/tests/scripts/benchmarkGenomeGenerate.sh
#
# To compare against both a recent optimized baseline and a raw/base STAR build:
#   BASELINE_STAR_BIN=/path/to/optimized/STAR BASE_STAR_BIN=/path/to/base/STAR \
#       THREADS="64 128" extras/tests/scripts/benchmarkGenomeGenerate.sh
#
# For performance measurements, prefer a realistic FASTA:
#   BENCH_FASTA=/path/to/genome.fa LIMIT_GENOME_GENERATE_RAM=31000000000 \
#       THREADS="1 4 8 16 32" extras/tests/scripts/benchmarkGenomeGenerate.sh
#
# To benchmark a splice-junction-augmented genome index:
#   BENCH_FASTA=/path/to/genome.fa BENCH_GTF=/path/to/genes.gtf SJDB_OVERHANG=93 \
#       THREADS="32 64" extras/tests/scripts/benchmarkGenomeGenerate.sh
#
# To pass multiple FASTA files, use a space-separated list:
#   BENCH_FASTA_FILES="/path/to/genome.fa /path/to/spikeins.fa" ...
#
# The built-in FASTA is intentionally tiny. It is for fast byte-identity checks,
# not for meaningful wall-clock speedup measurements.
# Use SYNTHETIC_REPEAT=N to scale it up deterministically when BENCH_FASTA is
# not set.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"

star_bin="${STAR_BIN:-${repo_root}/source/STAR}"
baseline_star_bin="${BASELINE_STAR_BIN:-}"
base_star_bin="${BASE_STAR_BIN:-}"
out_root="${OUT_DIR:-/tmp/star-genome-generate-bench.$(date +%Y%m%d_%H%M%S).$$}"
threads_list="${THREADS:-1 4}"
runs_per_thread="${RUNS_PER_THREAD:-2}"
baseline_runs_per_thread="${BASELINE_RUNS_PER_THREAD:-${runs_per_thread}}"
base_runs_per_thread="${BASE_RUNS_PER_THREAD:-${baseline_runs_per_thread}}"
genome_saindex_nbases="${GENOME_SAINDEX_NBASES:-2}"
genome_chrbin_nbits="${GENOME_CHRBIN_NBITS:-4}"
limit_genome_ram="${LIMIT_GENOME_GENERATE_RAM:-20000}"
bench_gtf="${BENCH_GTF:-}"
sjdb_overhang="${SJDB_OVERHANG:-}"
time_bin="${TIME_BIN:-/usr/bin/time}"
synthetic_repeat="${SYNTHETIC_REPEAT:-1}"
iostat_bin="${IOSTAT_BIN:-$(command -v iostat || true)}"
iostat_interval="${IOSTAT_INTERVAL:-5}"
cpu_sample_interval="${CPU_SAMPLE_INTERVAL:-1}"

if [[ ! -x "${star_bin}" ]]; then
    echo "ERROR: STAR_BIN is not executable: ${star_bin}" >&2
    exit 1
fi

if [[ -n "${baseline_star_bin}" && ! -x "${baseline_star_bin}" ]]; then
    echo "ERROR: BASELINE_STAR_BIN is not executable: ${baseline_star_bin}" >&2
    exit 1
fi

if [[ -n "${base_star_bin}" && ! -x "${base_star_bin}" ]]; then
    echo "ERROR: BASE_STAR_BIN is not executable: ${base_star_bin}" >&2
    exit 1
fi

if [[ ! -x "${time_bin}" ]]; then
    echo "ERROR: TIME_BIN is not executable: ${time_bin}" >&2
    exit 1
fi

if [[ -n "${bench_gtf}" && ! -r "${bench_gtf}" ]]; then
    echo "ERROR: BENCH_GTF is not readable: ${bench_gtf}" >&2
    exit 1
fi

if [[ -n "${sjdb_overhang}" && -z "${bench_gtf}" ]]; then
    echo "ERROR: SJDB_OVERHANG requires BENCH_GTF" >&2
    exit 1
fi

mkdir -p "${out_root}"
out_root="$(cd "${out_root}" && pwd)"

fasta="${BENCH_FASTA:-${out_root}/synthetic.fa}"
bench_fasta_files="${BENCH_FASTA_FILES:-${BENCH_FASTA:-}}"
fasta_files=()

make_synthetic_fasta() {
    local out_fasta="$1"
    local repeat_i
    {
        echo ">chrA"
        for ((repeat_i=0; repeat_i<synthetic_repeat; repeat_i++)); do
            echo "ACGTTGCACCATGGTACGATCGTACGTTAGCTAGGCTAACCGTTAACGATGCTAGCTTACGATCGATGCGTACCATCGTTAACG"
            echo "TGCAACGGTACCTAGCATCGATCGTAGGCTAACGTACGATGCTTACCGATCGTAGCTAGCATGCTAGGATCCGATCGTACGTAA"
            echo "GATCGTACGTTAGCATGCCGTAACGATCGTAGCTAACGTTGCACCATGGTACGATCGTACGTTAGCTAGGCTAACCGTTAACGA"
            echo "TCGATGCGTACCATCGTTAACGTGCAACGGTACCTAGCATCGATCGTAGGCTAACGTACGATGCTTACCGATCGTAGCTAGCAT"
            echo "GCTAGGATCCGATCGTACGTAAGATCGTACGTTAGCATGCCGTAACGATCGTAGCTAACGTTGCACCATGGTACGATCGTACGT"
            echo "TAGCTAGGCTAACCGTTAACGATGCTAGCTTACGATCGATGCGTACCATCGTTAACGTGCAACGGTACCTAGCATCGATCGTAG"
            echo "GCTAACGTACGATGCTTACCGATCGTAGCTAGCATGCTAGGATCCGATCGTACGTAAGATCGTACGTTAGCATGCCGTAACGAT"
            echo "CGTAGCTAACGTTGCACCATGGTACGATCGTACGTTAGCTAGGCTAACCGTTAACGATGCTAGCTTACGATCGATGCGTACCAT"
            echo "CGTTAACGTGCAACGGTACCTAGCATCGATCGTAGGCTAACGTACGATGCTTACCGATCGTAGCTAGCATGCTAGGATCCGATC"
            echo "GTACGTAAGATCGTACGTTAGCATGCCGTAACGATCGTAGCTAACGTTGCACCATGGTACGATCGTACGTTAGCTAGGCTAAC"
        done
        echo ">chrB"
        for ((repeat_i=0; repeat_i<synthetic_repeat; repeat_i++)); do
            echo "TTGACCGTAGCTAACGATCGTACCGTTAGCATCGATGCTAACCGTAGGCTTACGATCGTAGCATGCGTTAACCGATCGTACCTA"
            echo "GCTAACGTTAGGCTACGATGCTAGCATCGTTAACGATCCGTAGCTTACGTTGCATCGATGGTACCGTAGCATGCTAACGATCGT"
            echo "ACCGTTAGCATCGATGCTAACCGTAGGCTTACGATCGTAGCATGCGTTAACCGATCGTACCTAGCTAACGTTAGGCTACGATGC"
            echo "TAGCATCGTTAACGATCCGTAGCTTACGTTGCATCGATGGTACCGTAGCATGCTAACGATCGTTTGACCGTAGCTAACGATCGT"
            echo "ACCGTTAGCATCGATGCTAACCGTAGGCTTACGATCGTAGCATGCGTTAACCGATCGTACCTAGCTAACGTTAGGCTACGATGC"
            echo "TAGCATCGTTAACGATCCGTAGCTTACGTTGCATCGATGGTACCGTAGCATGCTAACGATCGTTTGACCGTAGCTAACGATCGT"
            echo "ACCGTTAGCATCGATGCTAACCGTAGGCTTACGATCGTAGCATGCGTTAACCGATCGTACCTAGCTAACGTTAGGCTACGATGC"
            echo "TAGCATCGTTAACGATCCGTAGCTTACGTTGCATCGATGGTACCGTAGCATGCTAACGATCGTTTGACCGTAGCTAACGATCGT"
        done
    } > "${out_fasta}"
}

if [[ -z "${bench_fasta_files}" ]]; then
    make_synthetic_fasta "${fasta}"
    fasta_files=("${fasta}")
else
    read -r -a fasta_files <<< "${bench_fasta_files}"
fi

for fasta_file in "${fasta_files[@]}"; do
    if [[ ! -r "${fasta_file}" ]]; then
        echo "ERROR: FASTA is not readable: ${fasta_file}" >&2
        exit 1
    fi
done

if [[ "${#fasta_files[@]}" -eq 0 ]]; then
    echo "ERROR: no FASTA files were provided" >&2
    exit 1
fi

extract_stage_times() {
    local log_path="$1"
    local stage_path="$2"
    local log_year
    log_year="$(date +%Y)"

    awk -v year="${log_year}" '
        BEGIN {
            month["Jan"]=1; month["Feb"]=2; month["Mar"]=3; month["Apr"]=4;
            month["May"]=5; month["Jun"]=6; month["Jul"]=7; month["Aug"]=8;
            month["Sep"]=9; month["Oct"]=10; month["Nov"]=11; month["Dec"]=12;
            print "stage\tstart_time\tend_time\tseconds";
        }

        function timestamp(line, mon, day, hour, minute, second) {
            mon=substr(line,1,3);
            if (!(mon in month)) {
                return "";
            }
            day=substr(line,5,2)+0;
            hour=substr(line,8,2)+0;
            minute=substr(line,11,2)+0;
            second=substr(line,14,2)+0;
            return mktime(year " " month[mon] " " day " " hour " " minute " " second);
        }

        function mark(key,line, parsed_time) {
            if (!(key in times)) {
                parsed_time=timestamp(line);
                if (parsed_time!="") {
                    times[key]=parsed_time;
                    labels[key]=substr(line,1,15);
                }
            }
        }

        function emit(stage,start_key,end_key) {
            if ((start_key in times) && (end_key in times)) {
                printf "%s\t%s\t%s\t%d\n", stage, labels[start_key], labels[end_key], times[end_key]-times[start_key];
            } else {
                printf "%s\tNA\tNA\tNA\n", stage;
            }
        }

        /starting to sort Suffix Array/ { mark("sa_sort_start",$0); next }
        /sorting Suffix Array chunks/ { mark("sa_chunk_sort",$0); next }
        /packing SA from RAM chunks|loading chunks from disk, packing SA/ { mark("sa_pack",$0); next }
        /finished generating suffix array/ { mark("sa_done",$0); next }
        /generating Suffix Array index/ { mark("saindex_start",$0); next }
        /completed Suffix Array index/ { mark("saindex_done",$0); next }
        /Finished preparing junctions/ { mark("sj_prepare_done",$0); next }
        /inserting junctions into the genome indices/ { mark("sj_insert_start",$0); next }
        /Finished SA search:/ { mark("sj_sa_search_done",$0); next }
        /Finished sorting SA indicesL/ { mark("sj_sort_done",$0); next }
        /Finished inserting junction indices/ { mark("sj_insert_indices_done",$0); next }
        /Finished SAi$/ { mark("sj_sai_done",$0); next }
        /writing Genome to disk/ { mark("write_genome_start",$0); next }
        /writing Suffix Array to disk/ { mark("write_sa_start",$0); next }
        /writing SAindex to disk/ { mark("write_saindex_start",$0); next }
        /finished successfully/ { mark("finish_success",$0); next }

        END {
            emit("sa_prefix_plan_seconds","sa_sort_start","sa_chunk_sort");
            emit("sa_chunk_sort_seconds","sa_chunk_sort","sa_pack");
            emit("sa_pack_seconds","sa_pack","sa_done");
            emit("saindex_seconds","saindex_start","saindex_done");
            emit("sj_prepare_seconds","saindex_done","sj_prepare_done");
            emit("sj_sa_search_seconds","sj_insert_start","sj_sa_search_done");
            emit("sj_sort_seconds","sj_sa_search_done","sj_sort_done");
            emit("sj_insert_indices_seconds","sj_sort_done","sj_insert_indices_done");
            emit("sj_sai_seconds","sj_insert_indices_done","sj_sai_done");
            emit("write_genome_seconds","write_genome_start","write_sa_start");
            emit("write_sa_seconds","write_sa_start","write_saindex_start");
            emit("write_saindex_seconds","write_saindex_start","finish_success");
        }
    ' "${log_path}" > "${stage_path}"
}

stage_metric() {
    local stage_path="$1"
    local stage_name="$2"
    awk -F '\t' -v stage="${stage_name}" '$1==stage {print $4; found=1; exit} END {if (!found) print "NA"}' "${stage_path}"
}

summary="${out_root}/summary.tsv"
printf "label\tthreads\trun\treal_seconds\tuser_seconds\tsys_seconds\tmax_rss_kb\tsa_chunk_fill_strategy\tsa_chunk_sort_prefix_length\tsa_chunk_sort_granularity\tsa_chunk_storage_strategy\tsa_chunk_system_available_bytes\tsa_chunk_ram_peak_bytes\tsa_all_scatter_available_bytes\tsa_all_scatter_required_bytes\tsa_batched_fill_batches\tsaindex_traversal_strategy\tsaindex_event_count\tsa_prefix_plan_seconds\tsa_chunk_sort_seconds\tsa_pack_seconds\tsaindex_seconds\tsj_prepare_seconds\tsj_sa_search_seconds\tsj_sort_seconds\tsj_insert_indices_seconds\tsj_sai_seconds\twrite_genome_seconds\twrite_sa_seconds\twrite_saindex_seconds\tstage_times_log\tgenome_dir\tio_log\tcpu_log\n" > "${summary}"

run_dir=""

run_generate() {
    local label="$1"
    local bin="$2"
    local threads="$3"
    local run_id="$4"
    run_dir="${out_root}/${label}_t${threads}_r${run_id}"
    local log_file="${out_root}/${label}_t${threads}_r${run_id}.log"
    local time_file="${out_root}/${label}_t${threads}_r${run_id}.time"
    local io_file="${out_root}/${label}_t${threads}_r${run_id}.iostat"
    local cpu_file="${out_root}/${label}_t${threads}_r${run_id}.threads"
    local stage_file="${out_root}/${label}_t${threads}_r${run_id}.stage_times.tsv"
    local iostat_pid=""
    local cpu_sampler_pid=""
    local star_wrapper_pid=""

    if [[ -e "${run_dir}" ]]; then
        echo "ERROR: output directory already exists: ${run_dir}" >&2
        exit 1
    fi
    mkdir -p "${run_dir}"

    if [[ -n "${iostat_bin}" && -x "${iostat_bin}" ]]; then
        "${iostat_bin}" -xz "${iostat_interval}" > "${io_file}" 2>&1 &
        iostat_pid="$!"
    else
        printf "iostat unavailable\n" > "${io_file}"
    fi

    printf "timestamp\tpid\ttid\tpsr\tpcpu\tstat\tcomm\targs\n" > "${cpu_file}"

    local star_args=(
        --runMode genomeGenerate
        --runThreadN "${threads}"
        --genomeDir "${run_dir}"
        --genomeFastaFiles "${fasta_files[@]}"
        --genomeSAindexNbases "${genome_saindex_nbases}"
        --genomeChrBinNbits "${genome_chrbin_nbits}"
        --limitGenomeGenerateRAM "${limit_genome_ram}"
    )

    if [[ -n "${bench_gtf}" ]]; then
        star_args+=(--sjdbGTFfile "${bench_gtf}")
    fi

    if [[ -n "${sjdb_overhang}" ]]; then
        star_args+=(--sjdbOverhang "${sjdb_overhang}")
    fi

    set +e
    (
        cd "${out_root}"
        "${time_bin}" -f "%e\t%U\t%S\t%M" -o "${time_file}" \
            "${bin}" "${star_args[@]}" \
            > "${log_file}" 2>&1
    ) &
    star_wrapper_pid="$!"

    (
        while kill -0 "${star_wrapper_pid}" 2>/dev/null; do
            local_now="$(date -Ins)"
            ps -eLo pid,tid,psr,pcpu,stat,comm,args \
                | awk -v now="${local_now}" -v dir="${run_dir}" \
                    'NR>1 && $6=="STAR" && index($0,dir)>0 {printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t", now, $1, $2, $3, $4, $5, $6; for (i=7; i<=NF; i++) printf "%s%s", $i, (i<NF ? OFS : ORS)}' \
                >> "${cpu_file}"
            sleep "${cpu_sample_interval}"
        done
    ) &
    cpu_sampler_pid="$!"

    wait "${star_wrapper_pid}"
    local star_status="$?"
    set -e

    if [[ -n "${cpu_sampler_pid}" ]]; then
        kill "${cpu_sampler_pid}" 2>/dev/null || true
        wait "${cpu_sampler_pid}" 2>/dev/null || true
    fi

    if [[ -n "${iostat_pid}" ]]; then
        kill "${iostat_pid}" 2>/dev/null || true
        wait "${iostat_pid}" 2>/dev/null || true
    fi

    if [[ "${star_status}" -ne 0 ]]; then
        echo "ERROR: STAR genomeGenerate failed for ${label} t${threads} run ${run_id}; see ${log_file}" >&2
        exit "${star_status}"
    fi

    extract_stage_times "${run_dir}/Log.out" "${stage_file}"

    local real_s user_s sys_s max_rss sa_chunk_fill_strategy sa_chunk_sort_prefix_length sa_chunk_sort_granularity sa_chunk_storage_strategy sa_chunk_system_available_bytes sa_chunk_ram_peak_bytes sa_all_scatter_available_bytes sa_all_scatter_required_bytes sa_batched_fill_batches saindex_traversal_strategy saindex_event_count
    local sa_prefix_plan_seconds sa_chunk_sort_seconds sa_pack_seconds saindex_seconds sj_prepare_seconds sj_sa_search_seconds sj_sort_seconds sj_insert_indices_seconds sj_sai_seconds write_genome_seconds write_sa_seconds write_saindex_seconds
    read -r real_s user_s sys_s max_rss < "${time_file}"
    sa_chunk_fill_strategy="$(awk -F': ' '/SA chunk fill strategy:/ {print $2; found=1; exit} END {if (!found) print "NA"}' "${run_dir}/Log.out")"
    sa_chunk_sort_prefix_length="$(awk -F': ' '/SA chunk sort prefix length:/ {print $2; found=1; exit} END {if (!found) print "NA"}' "${run_dir}/Log.out")"
    sa_chunk_sort_granularity="$(awk -F': ' '/SA chunk sort granularity:/ {print $2; found=1; exit} END {if (!found) print "NA"}' "${run_dir}/Log.out")"
    sa_chunk_storage_strategy="$(awk -F': ' '/SA chunk storage strategy:/ {print $2; found=1; exit} END {if (!found) print "NA"}' "${run_dir}/Log.out")"
    sa_chunk_system_available_bytes="$(awk '/SA chunk system available bytes:/ {gsub(/;/,""); print $6; found=1; exit} END {if (!found) print "NA"}' "${run_dir}/Log.out")"
    sa_chunk_ram_peak_bytes="$(awk '/SA chunk retained bytes:/ {gsub(/;/,""); print $9; found=1; exit} END {if (!found) print "NA"}' "${run_dir}/Log.out")"
    read -r sa_all_scatter_available_bytes sa_all_scatter_required_bytes < <(awk '/SA chunk all-scatter available bytes:/ {gsub(/;/,""); print $6, $10; found=1; exit} END {if (!found) print "NA NA"}' "${run_dir}/Log.out")
    sa_batched_fill_batches="$(awk -F': ' '/SA chunk estimated batched-fill batches:/ {print $2; found=1; exit} END {if (!found) print "NA"}' "${run_dir}/Log.out")"
    saindex_traversal_strategy="$(awk -F': ' '/SAindex traversal strategy:/ {print $2; found=1; exit} END {if (!found) print "NA"}' "${run_dir}/Log.out")"
    saindex_event_count="$(awk -F': ' '/SAindex event count:/ {print $2; found=1; exit} END {if (!found) print "NA"}' "${run_dir}/Log.out")"
    sa_prefix_plan_seconds="$(stage_metric "${stage_file}" "sa_prefix_plan_seconds")"
    sa_chunk_sort_seconds="$(stage_metric "${stage_file}" "sa_chunk_sort_seconds")"
    sa_pack_seconds="$(stage_metric "${stage_file}" "sa_pack_seconds")"
    saindex_seconds="$(stage_metric "${stage_file}" "saindex_seconds")"
    sj_prepare_seconds="$(stage_metric "${stage_file}" "sj_prepare_seconds")"
    sj_sa_search_seconds="$(stage_metric "${stage_file}" "sj_sa_search_seconds")"
    sj_sort_seconds="$(stage_metric "${stage_file}" "sj_sort_seconds")"
    sj_insert_indices_seconds="$(stage_metric "${stage_file}" "sj_insert_indices_seconds")"
    sj_sai_seconds="$(stage_metric "${stage_file}" "sj_sai_seconds")"
    write_genome_seconds="$(stage_metric "${stage_file}" "write_genome_seconds")"
    write_sa_seconds="$(stage_metric "${stage_file}" "write_sa_seconds")"
    write_saindex_seconds="$(stage_metric "${stage_file}" "write_saindex_seconds")"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${label}" "${threads}" "${run_id}" "${real_s}" "${user_s}" "${sys_s}" "${max_rss}" \
        "${sa_chunk_fill_strategy}" "${sa_chunk_sort_prefix_length}" "${sa_chunk_sort_granularity}" \
        "${sa_chunk_storage_strategy}" "${sa_chunk_system_available_bytes}" "${sa_chunk_ram_peak_bytes}" \
        "${sa_all_scatter_available_bytes}" "${sa_all_scatter_required_bytes}" "${sa_batched_fill_batches}" "${saindex_traversal_strategy}" "${saindex_event_count}" \
        "${sa_prefix_plan_seconds}" "${sa_chunk_sort_seconds}" "${sa_pack_seconds}" "${saindex_seconds}" \
        "${sj_prepare_seconds}" "${sj_sa_search_seconds}" "${sj_sort_seconds}" "${sj_insert_indices_seconds}" "${sj_sai_seconds}" \
        "${write_genome_seconds}" "${write_sa_seconds}" "${write_saindex_seconds}" \
        "${stage_file}" "${run_dir}" "${io_file}" "${cpu_file}" >> "${summary}"
}

compare_core_outputs() {
    local dir_a="$1"
    local dir_b="$2"
    local label="$3"
    local file

    for file in Genome SA SAindex chrName.txt chrStart.txt chrLength.txt chrNameLength.txt sjdbList.out.tab; do
        if [[ -e "${dir_a}/${file}" || -e "${dir_b}/${file}" ]]; then
            if ! cmp -s "${dir_a}/${file}" "${dir_b}/${file}"; then
                echo "ERROR: ${label}: ${file} differs between ${dir_a} and ${dir_b}" >&2
                exit 1
            fi
        fi
    done
}

first_candidate_dir=""

for threads in ${threads_list}; do
    candidate_reference_dir=""

    run_id=1
    while [[ "${run_id}" -le "${runs_per_thread}" ]]; do
        run_generate candidate "${star_bin}" "${threads}" "${run_id}"
        candidate_dir="${run_dir}"

        if [[ -z "${candidate_reference_dir}" ]]; then
            candidate_reference_dir="${candidate_dir}"
        else
            compare_core_outputs "${candidate_reference_dir}" "${candidate_dir}" "candidate repeat t${threads}"
        fi

        run_id=$((run_id+1))
    done

    if [[ -z "${first_candidate_dir}" ]]; then
        first_candidate_dir="${candidate_reference_dir}"
    else
        compare_core_outputs "${first_candidate_dir}" "${candidate_reference_dir}" "candidate thread-count t${threads}"
    fi

    if [[ -n "${baseline_star_bin}" ]]; then
        baseline_reference_dir=""
        run_id=1
        while [[ "${run_id}" -le "${baseline_runs_per_thread}" ]]; do
            run_generate baseline "${baseline_star_bin}" "${threads}" "${run_id}"
            baseline_dir="${run_dir}"

            if [[ -z "${baseline_reference_dir}" ]]; then
                baseline_reference_dir="${baseline_dir}"
            else
                compare_core_outputs "${baseline_reference_dir}" "${baseline_dir}" "baseline repeat t${threads}"
            fi

            run_id=$((run_id+1))
        done

        compare_core_outputs "${baseline_reference_dir}" "${candidate_reference_dir}" "baseline vs candidate t${threads}"
    fi

    if [[ -n "${base_star_bin}" ]]; then
        base_reference_dir=""
        run_id=1
        while [[ "${run_id}" -le "${base_runs_per_thread}" ]]; do
            run_generate base "${base_star_bin}" "${threads}" "${run_id}"
            base_dir="${run_dir}"

            if [[ -z "${base_reference_dir}" ]]; then
                base_reference_dir="${base_dir}"
            else
                compare_core_outputs "${base_reference_dir}" "${base_dir}" "base repeat t${threads}"
            fi

            run_id=$((run_id+1))
        done

        compare_core_outputs "${base_reference_dir}" "${candidate_reference_dir}" "base vs candidate t${threads}"
    fi
done

echo "Benchmark complete: ${out_root}"
echo "Summary: ${summary}"
