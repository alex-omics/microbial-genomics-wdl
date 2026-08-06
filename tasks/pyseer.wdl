version 1.0

# tasks/pyseer.wdl
#
# Standalone pyseer tasks for microbial pangenome-wide association studies
# (GWAS): LMM association with the likelihood-ratio test ("vanilla" pyseer),
# optional lineage effects, optional per-covariate-column scanning, and a
# gene/module-annotation join that never touches pyseer's own output files.
#
# `presence_absence_rtab` is deliberately generic: pyseer's own docs describe
# the Rtab format it accepts via --pres as usable "flexibly to represent
# variants from other sources" — a block_id x sample 0/1 matrix, nothing
# gene-specific about it. A Panaroo gene_presence_absence.Rtab and a
# pangenome-network module Rtab (module_id x isolate) both work here
# unchanged; only the annotation table handed to pyseer_annotate_results
# needs to match whichever block type was tested.
#
# Deliberately independent of tasks/gwas/task_pyseer.wdl (feature/multimodal-
# gwas), which is a much larger, tightly-coupled translation of the microGWAS
# Snakemake DAG (unitigs, structural variants, panfeed, whole-genome elastic
# net). That task expects its own upstream prepare_pyseer/mash/mlst inputs and
# is not meant to run on its own. This file is for calling pyseer by itself.

task pyseer_similarity_from_phylogeny {

    input {
        File    phylogeny_newick
        Int     cpu     = 1
        Int     mem_gb  = 4
        Int     disk_gb = 10
        String  docker  = "aarvani1/pyseer:1.4.2@sha256:20ba84511a4a7ebca154292b8a2af5e873556dbacd27999c5045c76f4f71d10c"
    }

    parameter_meta {
        phylogeny_newick: "Midpoint-rooted core-genome phylogeny (Newick)"
        cpu:               "Number of CPUs delegated to task (default = 1)"
        mem_gb:            "Amount of memory in GB delegated to task (default = 4)"
        disk_gb:           "Amount of disk space in GB delegated to task (default = 10)"
        docker:            "Container image"
    }

    meta {
        description: "Generate both a kinship (similarity) matrix and a plain patristic distance matrix from a phylogeny, via pyseer's own phylogeny_distance.py with and without --lmm. The distance matrix is not optional overhead: pyseer's --lineage refuses to run without one (\"Must also provide a distance matrix to report lineage effects\"), even when --lineage-clusters is supplied - it is not just the fallback for MDS-derived lineages. Both come from the same script and the same input tree, so producing both here is nearly free."
    }

    command <<<
        set -euo pipefail

        phylogeny_distance.py \
            --lmm \
            ~{phylogeny_newick} \
            > phylogeny_K.tsv

        phylogeny_distance.py \
            ~{phylogeny_newick} \
            > phylogeny_distances.tsv

        echo "Kinship matrix dimensions:"
        head -1 phylogeny_K.tsv | awk '{print NF " columns"}'
        wc -l < phylogeny_K.tsv | awk '{print $1 " rows"}'
        echo "Distance matrix dimensions:"
        head -1 phylogeny_distances.tsv | awk '{print NF " columns"}'
        wc -l < phylogeny_distances.tsv | awk '{print $1 " rows"}'
    >>>

    output {
        File kinship_matrix  = "phylogeny_K.tsv"
        File distance_matrix = "phylogeny_distances.tsv"
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk_gb} SSD"
        preemptible:    1
        maxRetries:     2
    }
}

task pyseer_lineage_effects {

    input {
        File            phenotype_tsv
        File            presence_absence_rtab
        File            distance_matrix
        File?           lineage_clusters
        File?           covariates_file
        String?         use_covariates

        Int             cpu     = 1
        Int             mem_gb  = 4
        Int             disk_gb = 10
        String          docker  = "aarvani1/pyseer:1.4.2@sha256:20ba84511a4a7ebca154292b8a2af5e873556dbacd27999c5045c76f4f71d10c"
    }

    parameter_meta {
        phenotype_tsv:         "Two columns: sample_id\\tphenotype_value"
        presence_absence_rtab: "block_id x isolate 0/1 matrix - only the first few blocks are actually used, see meta.description"
        distance_matrix:       "From pyseer_similarity_from_phylogeny. pyseer's --lineage refuses to run without one, even with lineage_clusters supplied."
        lineage_clusters:      "Two columns sample_id\\tcluster_id (e.g. BAPS), passed as --lineage-clusters. Falls back to pyseer's MDS-derived lineages if omitted."
        covariates_file:       "Tab-separated: sample_id, then one named column per covariate. Passed through so lineage effects are adjusted the same way as the main association."
        use_covariates:        "pyseer --use-covariates value, matching whatever was applied to the main association's covariates_file"
        cpu:                   "Number of CPUs delegated to task (default = 1)"
        mem_gb:                "Amount of memory in GB delegated to task (default = 4)"
        disk_gb:               "Amount of disk space in GB delegated to task (default = 10)"
        docker:                "Container image"
    }

    meta {
        description: "Report per-lineage effects as their own minimal-input pyseer call, deliberately never combined with the full variant-testing association. Mirrors the pattern already validated in the microGWAS translation (tasks/gwas/task_pyseer.wdl run_pyseer): that task's own comment says the lineage pass 'exists to produce the per-lineage effect table, not association statistics, so pyseer is given only enough input to initialise' - it feeds pyseer a 10-line slice of the real variant file, not the whole thing. This task does the same against presence_absence_rtab. Folding --lineage into the main --lmm call against the *full* Rtab is what caused an unexplained OOM against real 200-isolate BBSS data (6,427 variants, 2.5 MB Rtab) that killed a 32 GB machine in under a minute - too fast and too small to be a real memory shortage, and consistent with pyseer's lineage-effects code path and the full per-variant LMM loop compounding in one process."
    }

    command <<<
        set -euo pipefail

        head -n 10 ~{presence_absence_rtab} > small.Rtab

        pyseer \
            --phenotypes ~{phenotype_tsv} \
            --pres small.Rtab \
            --distances ~{distance_matrix} \
            --lineage --lineage-file lineage_effects.txt \
            --cpu ~{cpu} \
            ~{if defined(lineage_clusters) then "--lineage-clusters " + lineage_clusters else ""} \
            ~{if defined(covariates_file) then "--covariates " + covariates_file else ""} \
            ~{if defined(use_covariates) then "--use-covariates " + use_covariates else ""} \
            > lineage_pass.log \
            2> lineage_pass_stderr.log

        echo "Lineage effects complete:"
        wc -l lineage_effects.txt
    >>>

    output {
        File lineage_effects     = "lineage_effects.txt"
        File lineage_pass_log    = "lineage_pass.log"
        File lineage_pass_stderr = "lineage_pass_stderr.log"
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk_gb} SSD"
        preemptible:    1
        maxRetries:     2
    }
}

task pyseer_association {

    input {
        File            phenotype_tsv
        File            presence_absence_rtab
        File            kinship_matrix
        File?           variant_vcf

        Float           min_af = 0.05
        Float           max_af = 0.95

        # Joint adjustment: applied to every result this task produces.
        File?           covariates_file
        String?         use_covariates

        # Combination scan: for each entry here, run a *separate* association
        # using exactly that group of covariates_file columns (joined with
        # "+", e.g. "BAPS" or "RST+OspC+BAPS") as the sole covariates, so a
        # module's association can be checked against individual lineage
        # markers or deliberately chosen small groups of them - without
        # jointly loading every column at once, which burns degrees of
        # freedom fast (see pyseer_gwas README: an all-four joint run on an
        # 8-sample panel left zero real variants able to fit; the same panel
        # with just MLST plus one other column already lost half). Categorical
        # columns only — quantitative covariates need the use_covariates
        # escape hatch above. Not a request to enumerate every combination:
        # each entry here is one deliberately-chosen, auditable pyseer run,
        # not a black-box powerset search.
        Array[String]   covariate_combinations = []

        Int             cpu     = 4
        Int             mem_gb  = 8
        Int             disk_gb = 30
        String          docker  = "aarvani1/pyseer:1.4.2@sha256:20ba84511a4a7ebca154292b8a2af5e873556dbacd27999c5045c76f4f71d10c"
    }

    parameter_meta {
        phenotype_tsv:         "Two columns: sample_id\\tphenotype_value"
        presence_absence_rtab: "block_id x isolate 0/1 matrix - Panaroo's gene_presence_absence.Rtab or an equivalent module Rtab from a pangenome-network decomposition"
        kinship_matrix:        "From pyseer_similarity_from_phylogeny"
        variant_vcf:           "Optional VCF for an additional SNP-based association, run alongside the presence/absence test"
        min_af:                "Minimum allele/block frequency filter (default = 0.05)"
        max_af:                "Maximum allele/block frequency filter (default = 0.95)"
        covariates_file:       "Tab-separated: sample_id, then one named column per covariate (e.g. RST, OspC, MLST, BAPS)"
        use_covariates:        "pyseer --use-covariates value (column indices into covariates_file, 'q' suffix for quantitative) applied jointly to the main association"
        covariate_combinations: "Groups of covariates_file column names to test together, one group per pyseer run. Each entry is one or more column names joined with '+' (e.g. 'BAPS' or 'RST+OspC+BAPS'); a bare name is a single-covariate run. Column names must not contain '+', whitespace, or tabs."
        cpu:                   "Number of CPUs delegated to task (default = 4)"
        mem_gb:                "Amount of memory in GB delegated to task (default = 8)"
        disk_gb:               "Amount of disk space in GB delegated to task (default = 30)"
        docker:                "Container image"
    }

    meta {
        description: "Run pyseer's LMM association (likelihood-ratio test) of presence/absence blocks against a phenotype, correcting for population structure via a phylogeny-derived kinship matrix. Scans deliberately-chosen covariate combinations and runs an additional SNP-based pass if a VCF is supplied. Lineage effects are a separate task (pyseer_lineage_effects) - see its meta.description for why folding --lineage into this call is dangerous at real scale."
    }

    command <<<
        set -euo pipefail

        echo "=== pyseer association ==="
        wc -l ~{phenotype_tsv}
        head -1 ~{presence_absence_rtab} | awk '{print NF-1 " isolates"}'

        # --- Presence/absence association (primary analysis) ---
        pyseer \
            --lmm \
            --phenotypes ~{phenotype_tsv} \
            --pres ~{presence_absence_rtab} \
            --similarity ~{kinship_matrix} \
            --min-af ~{min_af} \
            --max-af ~{max_af} \
            --cpu ~{cpu} \
            --output-patterns gene_patterns.txt \
            ~{if defined(covariates_file) then "--covariates " + covariates_file else ""} \
            ~{if defined(use_covariates) then "--use-covariates " + use_covariates else ""} \
            > pyseer_gene_results.tsv \
            2> pyseer_gene_stderr.log

        echo "Gene association complete. Results:"
        wc -l pyseer_gene_results.tsv

        # Pattern-counted significance threshold (0.05 / unique presence-
        # absence patterns), via pyseer's own count_patterns.py rather than a
        # hand-rolled count - far less conservative than counting every raw
        # variant, since co-occurring blocks share a pattern.
        THRESHOLD="$(count_patterns.py --threshold gene_patterns.txt 2>> pyseer_gene_stderr.log)"
        echo "${THRESHOLD}" > significance_threshold.txt
        echo "Bonferroni threshold: ${THRESHOLD}"

        python3 -c "
import pandas as pd
threshold = float(open('significance_threshold.txt').read().strip())
df = pd.read_csv('pyseer_gene_results.tsv', sep='\t')
sig = df[pd.to_numeric(df['lrt-pvalue'], errors='coerce') < threshold].sort_values('lrt-pvalue')
sig.to_csv('pyseer_gene_significant.tsv', sep='\t', index=False)
print(f'Significant hits (p < {threshold:.3e}): {len(sig)}')
" 2>&1 | tee -a pyseer_gene_stderr.log

        # --- Optional: SNP-based association if VCF provided ---
        if [ -n "~{default="" variant_vcf}" ]; then
            echo "=== Running SNP-based association ==="
            pyseer \
                --lmm \
                --phenotypes ~{phenotype_tsv} \
                --vcf ~{default="" variant_vcf} \
                --similarity ~{kinship_matrix} \
                --min-af ~{min_af} \
                --max-af ~{max_af} \
                --cpu ~{cpu} \
                --output-patterns snp_patterns.txt \
                > pyseer_snp_results.tsv \
                2> pyseer_snp_stderr.log
            echo "SNP association complete:"
            wc -l pyseer_snp_results.tsv
        else
            echo "No VCF provided, skipping SNP-based association."
        fi

        # --- Optional: scan each deliberately-chosen covariate combination ---
        mkdir -p covariate_scan
        mapfile -t COV_COMBOS < ~{write_lines(covariate_combinations)}
        if [ "${#COV_COMBOS[@]}" -gt 0 ] && [ -z "~{default="" covariates_file}" ]; then
            echo "ERROR: covariate_combinations given but covariates_file was not supplied" >&2
            exit 1
        fi
        if [ "${#COV_COMBOS[@]}" -gt 0 ]; then
            HEADER="$(head -1 ~{default="" covariates_file})"
            for COMBO in "${COV_COMBOS[@]}"; do
                IFS='+' read -ra MEMBERS <<< "${COMBO}"
                IDXS=()
                for NAME in "${MEMBERS[@]}"; do
                    IDX="$(awk -F'\t' -v col="${NAME}" '{for (i=1;i<=NF;i++) if ($i==col) {print i; exit}}' <<< "${HEADER}")"
                    if [ -z "${IDX}" ]; then
                        echo "ERROR: covariate column '${NAME}' (from combination '${COMBO}') not found in covariates_file header: ${HEADER}" >&2
                        exit 1
                    fi
                    IDXS+=("${IDX}")
                done
                echo "Scanning covariate combination '${COMBO}' (columns ${IDXS[*]} of covariates_file)"
                pyseer \
                    --lmm \
                    --phenotypes ~{phenotype_tsv} \
                    --pres ~{presence_absence_rtab} \
                    --similarity ~{kinship_matrix} \
                    --min-af ~{min_af} \
                    --max-af ~{max_af} \
                    --cpu ~{cpu} \
                    --covariates ~{default="" covariates_file} \
                    --use-covariates "${IDXS[@]}" \
                    > "covariate_scan/${COMBO}.tsv" \
                    2> "covariate_scan/${COMBO}.log"
            done

            # Long-format comparison table: one row per (variant, combination),
            # native pyseer column names preserved (this is a provenance side
            # table, not the renamed/reformatted output - that is
            # pyseer_annotate_results's job on a single file of interest).
            # Names are passed through a file, not shell-interpolated into the
            # python source, so quoting is never a concern here.
            printf '%s\n' "${COV_COMBOS[@]}" > covariate_combinations.txt
            python3 -c "
import pandas as pd
combos = [l.strip() for l in open('covariate_combinations.txt') if l.strip()]
frames = []
for combo in combos:
    df = pd.read_csv(f'covariate_scan/{combo}.tsv', sep='\t')
    df.insert(1, 'covariates', combo)
    frames.append(df)
pd.concat(frames, ignore_index=True).to_csv('covariate_scan_combined.tsv', sep='\t', index=False)
" 2>&1 | tee -a pyseer_gene_stderr.log
        else
            touch covariate_scan_combined.tsv
        fi

        echo "=== Done ==="
    >>>

    output {
        File            gene_results           = "pyseer_gene_results.tsv"
        File            gene_significant       = "pyseer_gene_significant.tsv"
        File            gene_patterns          = "gene_patterns.txt"
        File            significance_threshold = "significance_threshold.txt"
        File            gene_log               = "pyseer_gene_stderr.log"
        File?           snp_results            = "pyseer_snp_results.tsv"
        File?           snp_patterns           = "snp_patterns.txt"
        File?           snp_log                = "pyseer_snp_stderr.log"
        Array[File]     covariate_scan_results = glob("covariate_scan/*.tsv")
        File            covariate_scan_combined = "covariate_scan_combined.tsv"
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk_gb} SSD"
        preemptible:    1
        maxRetries:     2
    }
}

task pyseer_annotate_results {

    input {
        File     pyseer_results
        File?    annotation_table
        String   id_column         = "Gene"
        String   name_column       = "Non-unique Gene name"
        String   annotation_column = "Annotation"
        String   basename

        Int      cpu     = 1
        Int      mem_gb  = 4
        Int      disk_gb = 20
        String   docker  = "aarvani1/pyseer:1.4.2@sha256:20ba84511a4a7ebca154292b8a2af5e873556dbacd27999c5045c76f4f71d10c"
    }

    parameter_meta {
        pyseer_results:    "One native pyseer LMM results table (e.g. gene_results, gene_significant, or a single covariate_scan_results entry)"
        annotation_table:  "Panaroo/Roary gene_presence_absence.csv, or an equivalent module-annotation table. Left join on id_column == pyseer_results' variant column; unmatched variants are kept and labelled 'no match', never dropped."
        id_column:         "Header name in annotation_table holding the values that match pyseer's variant column (default = 'Gene')"
        name_column:       "Header name in annotation_table for a short descriptive label. Set to '' to skip. (default = 'Non-unique Gene name')"
        annotation_column: "Header name in annotation_table for a longer free-text annotation. Set to '' to skip. (default = 'Annotation')"
        basename:          "Basename for the two output TSVs, without extension"
        cpu:               "Number of CPUs delegated to task (default = 1)"
        mem_gb:            "Amount of memory in GB delegated to task (default = 4)"
        disk_gb:           "Amount of disk space in GB delegated to task (default = 20)"
        docker:            "Container image"
    }

    meta {
        description: "Join gene/module names and annotations onto a pyseer results table, and separately emit a human-readable copy with renamed headers and fixed-decimal formatting. Never modifies pyseer's own output - both files here are additional, not replacements, so anything expecting pyseer's exact column names still gets them unannotated."
    }

    command <<<
        set -euo pipefail

        cat > join.py <<'PY'
import csv
import sys

pyseer_path = sys.argv[1]
annotation_path = sys.argv[2] or None
id_column = sys.argv[3]
name_column = sys.argv[4]
annotation_column = sys.argv[5]
basename = sys.argv[6]

with open(pyseer_path, newline="", encoding="utf-8", errors="replace") as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    rows = list(reader)

n_var = header.index("variant")

lookup = {}
if annotation_path:
    # A real CSV parser is required: Panaroo's Annotation column carries
    # free-text descriptions containing commas, and splitting on commas would
    # silently shift every column after it.
    with open(annotation_path, newline="", encoding="utf-8", errors="replace") as fh:
        areader = csv.DictReader(fh)
        afields = areader.fieldnames or []
        if id_column not in afields:
            sys.stderr.write(
                f"ERROR: id_column '{id_column}' not found in annotation_table. "
                f"Columns present: {', '.join(afields[:12])}\n"
            )
            sys.exit(1)
        for row in areader:
            key = row.get(id_column)
            if not key:
                continue
            lookup[key] = (
                row.get(name_column, "") if name_column else "",
                row.get(annotation_column, "") if annotation_column else "",
            )

n_matched = 0
n_unmatched = 0

RENAME = {
    "variant": "gene_family_id",
    "af": "allele_frequency",
    "filter-pvalue": "prefilter_pvalue",
    "lrt-pvalue": "lrt_pvalue",
    "beta": "effect_size",
    "beta-std-err": "effect_size_stderr",
    "variant_h2": "variance_explained",
    "notes": "qc_flags",
}
FIXED6 = {"af", "beta", "beta-std-err", "variant_h2"}
PVAL = {"filter-pvalue", "lrt-pvalue"}

def fmt(col, value):
    if col in FIXED6:
        try:
            return f"{float(value):.6f}"
        except ValueError:
            return value
    if col in PVAL:
        try:
            v = float(value)
        except ValueError:
            return value
        return f"{v:.6f}" if abs(v) >= 1e-4 else f"{v:.3e}"
    return value

annotated_rows = []
readable_rows = []
for row in rows:
    variant = row[n_var]
    if annotation_path:
        gene_name, annotation = lookup.get(variant, ("no match", "no match"))
        if variant in lookup:
            n_matched += 1
        else:
            n_unmatched += 1
        annotated_rows.append(row[: n_var + 1] + [gene_name, annotation] + row[n_var + 1 :])
        readable_rows.append(
            [fmt(header[i], v) for i, v in enumerate(row[: n_var + 1])]
            + [gene_name, annotation]
            + [fmt(header[i], v) for i, v in enumerate(row[n_var + 1 :], start=n_var + 1)]
        )
    else:
        readable_rows.append([fmt(header[i], v) for i, v in enumerate(row)])

if annotation_path:
    ann_header = header[: n_var + 1] + ["gene_name", "annotation"] + header[n_var + 1 :]
    with open(f"{basename}_annotated.tsv", "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(ann_header)
        w.writerows(annotated_rows)
    readable_header = (
        [RENAME.get(h, h) for h in header[: n_var + 1]]
        + ["gene_name", "annotation"]
        + [RENAME.get(h, h) for h in header[n_var + 1 :]]
    )
else:
    readable_header = [RENAME.get(h, h) for h in header]

with open(f"{basename}_readable.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(readable_header)
    w.writerows(readable_rows)

sys.stderr.write(f"Matched {n_matched}, unmatched (labelled 'no match') {n_unmatched}\n")
PY

        python3 join.py \
            ~{pyseer_results} \
            "~{default="" annotation_table}" \
            "~{id_column}" \
            "~{name_column}" \
            "~{annotation_column}" \
            "~{basename}"
    >>>

    output {
        File? annotated = "~{basename}_annotated.tsv"
        File  readable  = "~{basename}_readable.tsv"
    }

    runtime {
        docker:         docker
        memory:         "~{mem_gb} GB"
        cpu:            cpu
        disks:          "local-disk ~{disk_gb} SSD"
        preemptible:    1
        maxRetries:     2
    }
}
