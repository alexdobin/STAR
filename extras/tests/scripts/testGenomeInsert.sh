#!/usr/bin/env bash
set -euo pipefail

# Fast regression test for --runMode genomeInsert.
#
# The test builds a tiny annotated base genome, inserts two extra named
# sequences persistently with insert-only annotations, and compares the result
# against a full rebuild from base+insert FASTA/GTF. SA is intentionally not
# required to match the full rebuild byte-for-byte because equivalent suffix tie
# ordering can differ.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"

star_bin="${STAR_BIN:-${repo_root}/source/STAR}"
threads="${THREADS:-4}"
keep_output="${KEEP_TEST_OUTPUT:-0}"
out_root="${OUT_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/star-genome-insert-test.XXXXXX")}"

if [[ ! -x "${star_bin}" ]]; then
    echo "ERROR: STAR_BIN is not executable: ${star_bin}" >&2
    exit 1
fi

cleanup() {
    if [[ "${keep_output}" == "1" ]]; then
        echo "Keeping test output: ${out_root}" >&2
    else
        rm -rf "${out_root}"
    fi
}
trap cleanup EXIT

mkdir -p "${out_root}/inputs"

base_fasta="${out_root}/inputs/base.fa"
insert_fasta="${out_root}/inputs/insert.fa"
base_gtf="${out_root}/inputs/base.gtf"
insert_gtf="${out_root}/inputs/insert.gtf"
insert_nojunction_gtf="${out_root}/inputs/insert_nojunction.gtf"
combined_gtf="${out_root}/inputs/combined.gtf"
combined_nojunction_gtf="${out_root}/inputs/combined_nojunction.gtf"
reads_fastq="${out_root}/inputs/reads.fq"

cat > "${base_fasta}" <<'EOF_FASTA'
>chrA
ACGTTGCACCATGGTACGATCGTACGTTAGCTAGGCTAACCGTTAACGATGCTAGCTTACGATCGATGCGTACCATCGTTAACGTGCAACGGTACCTAGCATCGATCGTAGGCTAACGTACGATGCTTACCGATCGTAGCTAGCATGCTAGGATCCGATCGTACGTAAGATCGTACGTTAGCATGCCGTAACGATCGTAGCTAACGTTGCACCATGGTACGATCGTACGTTAGCTAGGCTAACCGTTAACGA
>chrB
TTGACCGTAGCTAACGATCGTACCGTTAGCATCGATGCTAACCGTAGGCTTACGATCGTAGCATGCGTTAACCGATCGTACCTAGCTAACGTTAGGCTACGATGCTAGCATCGTTAACGATCCGTAGCTTACGTTGCATCGATGGTACCGTAGCATGCTAACGATCGTTTGACCGTAGCTAACGATCGTACCGTTAGCATCGATGCTAACCGTAGGCTTACGATCGTAGCATGCGTTAACCGATCGTACCTA
>chrBaseOnly
ATGCGTACAACTGACTGGTACCGTACGTTAGGACCTGTAACGTTACCGGATCTAGTCCGATGATCGTACGATCGATAGCTAGTCCGTAGCATTCGATGCAA
EOF_FASTA

cat > "${insert_fasta}" <<'EOF_FASTA'
>addGFP
GCTAGTACCGATGACTGATCGGATCTACGATCGTACCGGATGCTAGTTCGATCGATGACCTAGGCTAATCGGATGCTACCGATGATCCGATGCTAGTACCGATG
>addGST
TACCGATGGTACCTAGGATCGTTAACCGTAGCTAGGCTTACGATCGTAGGATCCGTAACGATGCTAGCTACCGATGGTACCGATGCTTAGGCTACCGATGTAA
EOF_FASTA

{
    printf 'chrA\ttest\texon\t1\t80\t.\t+\t.\tgene_id "geneA"; transcript_id "txA";\n'
    printf 'chrA\ttest\texon\t141\t220\t.\t+\t.\tgene_id "geneA"; transcript_id "txA";\n'
    printf 'chrB\ttest\texon\t20\t120\t.\t-\t.\tgene_id "geneB"; transcript_id "txB";\n'
} > "${base_gtf}"

{
    printf 'addGFP\ttest\texon\t1\t40\t.\t+\t.\tgene_id "geneGFP"; transcript_id "txGFP";\n'
    printf 'addGFP\ttest\texon\t71\t104\t.\t+\t.\tgene_id "geneGFP"; transcript_id "txGFP";\n'
    printf 'addGST\ttest\texon\t1\t103\t.\t+\t.\tgene_id "geneGST"; transcript_id "txGST";\n'
} > "${insert_gtf}"

{
    printf 'addGFP\ttest\texon\t1\t104\t.\t+\t.\tgene_id "geneGFP"; transcript_id "txGFP";\n'
    printf 'addGST\ttest\texon\t1\t103\t.\t+\t.\tgene_id "geneGST"; transcript_id "txGST";\n'
} > "${insert_nojunction_gtf}"

cat "${base_gtf}" "${insert_gtf}" > "${combined_gtf}"
cat "${base_gtf}" "${insert_nojunction_gtf}" > "${combined_nojunction_gtf}"

emit_fastq() {
    local read_name="$1"
    local seq="$2"
    local qual
    qual="$(printf '%*s' "${#seq}" '' | tr ' ' 'I')"
    printf '@%s\n%s\n+\n%s\n' "${read_name}" "${seq}" "${qual}"
}

{
    emit_fastq chrBaseOnly_read "ATGCGTACAACTGACTGGTACCGTACGTTAGGACCTGTAACGTTACCGGATCT"
    emit_fastq addGFP_read "GCTAGTACCGATGACTGATCGGATCTACGATCGTACCGGATGCTAGTTCGATCGATG"
    emit_fastq addGFP_spliced_read "GGATCTACGATCGTACCGGAGGATGCTACCGATGATCCGATGCTAGTACC"
    emit_fastq addGST_read "TACCGATGGTACCTAGGATCGTTAACCGTAGCTAGGCTTACGATCGTAGGATCCG"
} > "${reads_fastq}"

base_index="${out_root}/base_index"
incremental_index="${out_root}/incremental_index"
incremental_repeat_index="${out_root}/incremental_repeat_index"
overlay_index="${out_root}/overlay_index"
delta_index="${out_root}/delta_index"
delta_nojunction_index="${out_root}/delta_nojunction_index"
full_rebuild_index="${out_root}/full_rebuild_index"
full_rebuild_nojunction_index="${out_root}/full_rebuild_nojunction_index"

common_generate_args=(
    --runMode genomeGenerate
    --runThreadN "${threads}"
    --limitGenomeGenerateRAM 10000000
    --genomeSAindexNbases 2
    --genomeChrBinNbits 4
    --sjdbOverhang 9
)

"${star_bin}" \
    "${common_generate_args[@]}" \
    --sjdbGTFfile "${base_gtf}" \
    --genomeDir "${base_index}" \
    --genomeFastaFiles "${base_fasta}" \
    --outFileNamePrefix "${out_root}/base_" \
    > "${out_root}/base.stdout" 2>&1

base_sjdb_n="$(wc -l < "${base_index}/sjdbList.out.tab")"
if [[ "${base_sjdb_n}" -eq 0 ]]; then
    echo "ERROR: base test genome did not build any SJDB junctions" >&2
    exit 1
fi

set +e
(
    cd "${out_root}"
    "${star_bin}" \
        --runMode genomeInsert \
        --runThreadN "${threads}" \
        --genomeDir base_index \
        --genomeFastaFiles inputs/insert.fa \
        --genomeInsertOutDir ./base_index \
        --outFileNamePrefix same_dir_
) > "${out_root}/same_dir.stdout" 2>&1
same_dir_status="$?"
set -e

if [[ "${same_dir_status}" -eq 0 ]]; then
    echo "ERROR: genomeInsert unexpectedly allowed --genomeInsertOutDir to alias --genomeDir" >&2
    exit 1
fi

if ! grep -q -- "--genomeInsertOutDir cannot be the same directory as --genomeDir" "${out_root}/same_dir.stdout"; then
    echo "ERROR: genomeInsert same-directory failure did not report the expected message" >&2
    exit 1
fi

set +e
"${star_bin}" \
    --runMode genomeInsert \
    --runThreadN "${threads}" \
    --genomeDir "${base_index}" \
    --genomeFastaFiles "${insert_fasta}" \
    --sjdbGTFfile "${base_gtf}" \
    --genomeInsertOutDir "${out_root}/bad_base_gtf_index" \
    --outFileNamePrefix "${out_root}/bad_base_gtf_" \
    > "${out_root}/bad_base_gtf.stdout" 2>&1
bad_base_gtf_status="$?"
set -e

if [[ "${bad_base_gtf_status}" -eq 0 ]]; then
    echo "ERROR: genomeInsert unexpectedly allowed a GTF containing base genome annotations" >&2
    exit 1
fi

if ! grep -q -- "only annotations for inserted sequences" "${out_root}/bad_base_gtf.stdout"; then
    echo "ERROR: genomeInsert base-GTF failure did not report the expected message" >&2
    exit 1
fi

"${star_bin}" \
    --runMode genomeInsert \
    --runThreadN "${threads}" \
    --genomeDir "${base_index}" \
    --genomeFastaFiles "${insert_fasta}" \
    --sjdbGTFfile "${insert_gtf}" \
    --genomeInsertOutDir "${incremental_index}" \
    --outFileNamePrefix "${out_root}/incremental_" \
    > "${out_root}/incremental.stdout" 2>&1

"${star_bin}" \
    --runMode genomeInsert \
    --runThreadN "${threads}" \
    --genomeDir "${base_index}" \
    --genomeFastaFiles "${insert_fasta}" \
    --sjdbGTFfile "${insert_gtf}" \
    --genomeInsertOutDir "${incremental_repeat_index}" \
    --outFileNamePrefix "${out_root}/incremental_repeat_" \
    > "${out_root}/incremental_repeat.stdout" 2>&1

"${star_bin}" \
    --runMode genomeInsert \
    --runThreadN "${threads}" \
    --genomeDir "${base_index}" \
    --genomeFastaFiles "${insert_fasta}" \
    --sjdbGTFfile "${insert_gtf}" \
    --genomeInsertOutMode Overlay \
    --genomeInsertOutDir "${overlay_index}" \
    --outFileNamePrefix "${out_root}/overlay_" \
    > "${out_root}/overlay.stdout" 2>&1

if [[ ! -s "${overlay_index}/genomeInsertOverlay.tsv" ]]; then
    echo "ERROR: genomeInsert overlay did not write genomeInsertOverlay.tsv" >&2
    exit 1
fi

for file in Genome SA SAindex; do
    if [[ -e "${overlay_index}/${file}" ]]; then
        echo "ERROR: genomeInsert overlay unexpectedly wrote ${file}" >&2
        exit 1
    fi
done

"${star_bin}" \
    --runMode genomeInsert \
    --runThreadN "${threads}" \
    --genomeDir "${base_index}" \
    --genomeFastaFiles "${insert_fasta}" \
    --sjdbGTFfile "${insert_gtf}" \
    --genomeInsertOutMode Delta \
    --genomeInsertOutDir "${delta_index}" \
    --outFileNamePrefix "${out_root}/delta_" \
    > "${out_root}/delta.stdout" 2>&1

"${star_bin}" \
    --runMode genomeInsert \
    --runThreadN "${threads}" \
    --genomeDir "${base_index}" \
    --genomeFastaFiles "${insert_fasta}" \
    --sjdbGTFfile "${insert_nojunction_gtf}" \
    --genomeInsertOutMode Delta \
    --genomeInsertOutDir "${delta_nojunction_index}" \
    --outFileNamePrefix "${out_root}/delta_nojunction_" \
    > "${out_root}/delta_nojunction.stdout" 2>&1

if [[ ! -s "${delta_index}/genomeInsertOverlay.tsv" ]]; then
    echo "ERROR: genomeInsert delta did not write genomeInsertOverlay.tsv" >&2
    exit 1
fi

if [[ ! -s "${delta_index}/genomeInsertDelta.bin" ]]; then
    echo "ERROR: genomeInsert delta did not write genomeInsertDelta.bin" >&2
    exit 1
fi

if [[ ! -s "${delta_nojunction_index}/genomeInsertOverlay.tsv" ]]; then
    echo "ERROR: genomeInsert no-junction delta did not write genomeInsertOverlay.tsv" >&2
    exit 1
fi

if [[ ! -s "${delta_nojunction_index}/genomeInsertDelta.bin" ]]; then
    echo "ERROR: genomeInsert no-junction delta did not write genomeInsertDelta.bin" >&2
    exit 1
fi

for file in Genome SA SAindex; do
    if [[ -e "${delta_index}/${file}" ]]; then
        echo "ERROR: genomeInsert delta unexpectedly wrote ${file}" >&2
        exit 1
    fi
    if [[ -e "${delta_nojunction_index}/${file}" ]]; then
        echo "ERROR: genomeInsert no-junction delta unexpectedly wrote ${file}" >&2
        exit 1
    fi
done

"${star_bin}" \
    "${common_generate_args[@]}" \
    --sjdbGTFfile "${combined_gtf}" \
    --genomeDir "${full_rebuild_index}" \
    --genomeFastaFiles "${base_fasta}" "${insert_fasta}" \
    --outFileNamePrefix "${out_root}/full_rebuild_" \
    > "${out_root}/full_rebuild.stdout" 2>&1

"${star_bin}" \
    "${common_generate_args[@]}" \
    --sjdbGTFfile "${combined_nojunction_gtf}" \
    --genomeDir "${full_rebuild_nojunction_index}" \
    --genomeFastaFiles "${base_fasta}" "${insert_fasta}" \
    --outFileNamePrefix "${out_root}/full_rebuild_nojunction_" \
    > "${out_root}/full_rebuild_nojunction.stdout" 2>&1

required_equivalent_files=(
    Genome
    SAindex
    chrName.txt
    chrStart.txt
    chrLength.txt
    chrNameLength.txt
    sjdbInfo.txt
    sjdbList.out.tab
    sjdbList.fromGTF.out.tab
    exonInfo.tab
    exonGeTrInfo.tab
    geneInfo.tab
    transcriptInfo.tab
)

idempotent_files=(
    "${required_equivalent_files[@]}"
    SA
)

compare_files() {
    local label="$1"
    local dir_a="$2"
    local dir_b="$3"
    shift 3
    local file

    for file in "$@"; do
        if [[ ! -e "${dir_a}/${file}" || ! -e "${dir_b}/${file}" ]]; then
            echo "ERROR: ${label}: missing ${file}" >&2
            exit 1
        fi
        if ! cmp -s "${dir_a}/${file}" "${dir_b}/${file}"; then
            echo "ERROR: ${label}: ${file} differs" >&2
            exit 1
        fi
    done
}

compare_files "genomeInsert idempotence" "${incremental_index}" "${incremental_repeat_index}" "${idempotent_files[@]}"
compare_files "genomeInsert versus full rebuild" "${incremental_index}" "${full_rebuild_index}" "${required_equivalent_files[@]}"

if cmp -s "${incremental_index}/SA" "${full_rebuild_index}/SA"; then
    sa_full_rebuild_status="byte-identical"
else
    sa_full_rebuild_status="different-equivalent-ordering"
fi

align_and_extract_body() {
    local genome_dir="$1"
    local prefix="$2"
    mkdir -p "${prefix%/}"
    "${star_bin}" \
        --runThreadN "${threads}" \
        --genomeDir "${genome_dir}" \
        --readFilesIn "${reads_fastq}" \
        --quantMode GeneCounts \
        --outSAMtype SAM \
        --outFileNamePrefix "${prefix}" \
        > "${prefix}stdout" 2>&1
    grep -v '^@' "${prefix}Aligned.out.sam" > "${prefix}Aligned.body.sam"
}

align_and_extract_body "${incremental_index}" "${out_root}/align_incremental/"
align_and_extract_body "${full_rebuild_index}" "${out_root}/align_full/"
align_and_extract_body "${overlay_index}" "${out_root}/align_overlay/"
align_and_extract_body "${delta_index}" "${out_root}/align_delta/"
align_and_extract_body "${full_rebuild_nojunction_index}" "${out_root}/align_full_nojunction/"
align_and_extract_body "${delta_nojunction_index}" "${out_root}/align_delta_nojunction/"

diff -u "${out_root}/align_incremental/Aligned.body.sam" "${out_root}/align_full/Aligned.body.sam" > "${out_root}/alignment_body.diff"
diff -u "${out_root}/align_incremental/SJ.out.tab" "${out_root}/align_full/SJ.out.tab" > "${out_root}/alignment_sj.diff"
diff -u "${out_root}/align_overlay/Aligned.body.sam" "${out_root}/align_full/Aligned.body.sam" > "${out_root}/alignment_overlay_body.diff"
diff -u "${out_root}/align_overlay/SJ.out.tab" "${out_root}/align_full/SJ.out.tab" > "${out_root}/alignment_overlay_sj.diff"
diff -u "${out_root}/align_delta/Aligned.body.sam" "${out_root}/align_full/Aligned.body.sam" > "${out_root}/alignment_delta_body.diff"
diff -u "${out_root}/align_delta/SJ.out.tab" "${out_root}/align_full/SJ.out.tab" > "${out_root}/alignment_delta_sj.diff"
diff -u "${out_root}/align_delta_nojunction/Aligned.body.sam" "${out_root}/align_full_nojunction/Aligned.body.sam" > "${out_root}/alignment_delta_nojunction_body.diff"
diff -u "${out_root}/align_delta_nojunction/SJ.out.tab" "${out_root}/align_full_nojunction/SJ.out.tab" > "${out_root}/alignment_delta_nojunction_sj.diff"
diff -u "${out_root}/align_incremental/ReadsPerGene.out.tab" "${out_root}/align_full/ReadsPerGene.out.tab" > "${out_root}/gene_counts.diff"
diff -u "${out_root}/align_overlay/ReadsPerGene.out.tab" "${out_root}/align_full/ReadsPerGene.out.tab" > "${out_root}/gene_counts_overlay.diff"
diff -u "${out_root}/align_delta/ReadsPerGene.out.tab" "${out_root}/align_full/ReadsPerGene.out.tab" > "${out_root}/gene_counts_delta.diff"
diff -u "${out_root}/align_delta_nojunction/ReadsPerGene.out.tab" "${out_root}/align_full_nojunction/ReadsPerGene.out.tab" > "${out_root}/gene_counts_delta_nojunction.diff"

if ! grep -q "Loaded genome insert delta" "${out_root}/align_delta/Log.out"; then
    echo "ERROR: delta overlay alignment did not load genomeInsertDelta.bin" >&2
    exit 1
fi

if ! grep -q "Loaded genome insert delta" "${out_root}/align_delta_nojunction/Log.out"; then
    echo "ERROR: no-junction delta overlay alignment did not load genomeInsertDelta.bin" >&2
    exit 1
fi

if ! grep -q "Using virtual genome insert SA overlay" "${out_root}/align_delta_nojunction/Log.out"; then
    echo "ERROR: no-junction delta overlay alignment did not use the virtual SA overlay" >&2
    exit 1
fi

{
    printf 'addGFP\n'
    printf 'addGST\n'
    printf 'chrBaseOnly\n'
} > "${out_root}/expected_alignment_references.txt"

awk '$3!="*" {print $3}' "${out_root}/align_incremental/Aligned.body.sam" \
    | sort -u > "${out_root}/observed_alignment_references.txt"
diff -u "${out_root}/expected_alignment_references.txt" "${out_root}/observed_alignment_references.txt" > "${out_root}/alignment_references.diff"

{
    printf 'check\tstatus\n'
    printf 'same_directory_guard\tpass\n'
    printf 'insert_only_gtf_guard\tpass\n'
    printf 'base_sjdb_present\tpass\n'
    printf 'genomeInsert_overlay_manifest\tpass\n'
    printf 'genomeInsert_delta_manifest\tpass\n'
    printf 'genomeInsert_delta_file\tpass\n'
    printf 'genomeInsert_idempotence\tpass\n'
    printf 'genomeInsert_vs_full_rebuild_core_files\tpass\n'
    printf 'SA_vs_full_rebuild\t%s\n' "${sa_full_rebuild_status}"
    printf 'alignment_body_vs_full_rebuild\tpass\n'
    printf 'alignment_SJ_vs_full_rebuild\tpass\n'
    printf 'lazy_shift_sjdb_alignment_equivalence\tpass\n'
    printf 'gene_counts_vs_full_rebuild\tpass\n'
    printf 'overlay_alignment_body_vs_full_rebuild\tpass\n'
    printf 'overlay_alignment_SJ_vs_full_rebuild\tpass\n'
    printf 'overlay_gene_counts_vs_full_rebuild\tpass\n'
    printf 'delta_alignment_body_vs_full_rebuild\tpass\n'
    printf 'delta_alignment_SJ_vs_full_rebuild\tpass\n'
    printf 'delta_gene_counts_vs_full_rebuild\tpass\n'
    printf 'delta_nojunction_alignment_body_vs_full_rebuild\tpass\n'
    printf 'delta_nojunction_alignment_SJ_vs_full_rebuild\tpass\n'
    printf 'delta_nojunction_gene_counts_vs_full_rebuild\tpass\n'
    printf 'delta_replay_used\tpass\n'
    printf 'delta_virtual_SA_overlay_used\tpass\n'
    printf 'alignment_references_expected\tpass\n'
} > "${out_root}/checks.tsv"

cat "${out_root}/checks.tsv"
