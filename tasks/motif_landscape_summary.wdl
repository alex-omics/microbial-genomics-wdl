version 1.0

task motif_landscape_summary {

    input {
        Array[File]     per_isolate_landscape
        File?           ortholog_long_table
        Array[String]?  highlight_genes
        String          basename = "landscape_summary"
        String          docker   = "python:3.11-slim@sha256:db3ff2e1800a8581e2c48a27c3995339d47bdf046da21c7627accd3d51053a93"
    }

    parameter_meta {
        per_isolate_landscape: "One motif_landscape.tsv per isolate (tier 1 + tier 2, already computed within each genome)"
        ortholog_long_table:   "methylation_orthologs' _long.tsv, already gene x isolate with genic_sites_per_kb. Reused rather than recomputed -- the pangenome join already did the hard part."
        highlight_genes:       "Case-insensitive substrings to flag in gene_name or annotation, e.g. ['mex','opr','amp','nal']. Purely a sort/flag convenience over the full ranked table, not a filter -- everything is still reported."
        basename:              "Basename for the two summary tables"
        docker:                "Container image"
    }

    meta {
        description: "Tier 3: cross-isolate variability, at two resolutions. A motif or gene that is uniform across every isolate is very likely RM housekeeping -- the same enzyme doing the same thing everywhere. One that varies sharply across isolates is where isolate-specific biology, including candidate regulatory methylation, would actually show up. Neither table implies causation; both are a triage ranking for follow-up, not a phenotype association test."
    }

    command <<<
        set -euo pipefail

        cat > summarise.py <<'PY'
        import csv, statistics, sys
        from collections import defaultdict

        landscape_paths = sys.argv[1].split(",") if sys.argv[1] else []
        ortholog_path   = sys.argv[2] if sys.argv[2] != "-" else None
        highlight       = [h.strip().lower() for h in sys.argv[3].split(",") if h.strip()]
        motif_out, gene_out = sys.argv[4], sys.argv[5]

        # ---- motif-level -----------------------------------------------------
        by_motif = defaultdict(list)  # (motif, code) -> [in_motif_mean_pct, ...]
        n_isolates_seen = 0
        seen_samples = set()
        for path in landscape_paths:
            with open(path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    seen_samples.add(row["sample"])
                    pct = row.get("in_motif_mean_pct", "NA")
                    if pct == "NA":
                        continue
                    by_motif[(row["motif"], row["mod_code"])].append(float(pct))
        n_isolates_seen = len(seen_samples)

        motif_rows = []
        for (motif, code), vals in by_motif.items():
            mean = statistics.mean(vals)
            cv = (statistics.pstdev(vals) / mean) if mean > 0 and len(vals) > 1 else 0.0
            motif_rows.append([motif, code, len(vals), n_isolates_seen,
                                "%.2f" % mean, "%.2f" % min(vals), "%.2f" % max(vals),
                                "%.3f" % cv])
        motif_rows.sort(key=lambda r: -float(r[7]))

        with open(motif_out, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t", lineterminator="\n")
            w.writerow(["motif", "mod_code", "n_isolates_tested", "n_isolates_total",
                        "mean_in_motif_pct", "min_in_motif_pct", "max_in_motif_pct",
                        "cross_isolate_cv"])
            w.writerows(motif_rows)

        # ---- gene-level, from the existing ortholog long table ---------------
        # NA (gene absent) must stay excluded from variability -- an isolate
        # that lacks the gene contributes no information about how methylated
        # it is, and folding NA in as 0 would manufacture variability that is
        # really just presence/absence, the same trap the ortholog matrix
        # itself was built to avoid.
        gene_vals = defaultdict(list)
        gene_meta = {}
        if ortholog_path:
            with open(ortholog_path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    dens = row.get("genic_sites_per_kb", "NA")
                    g = row["ortholog_group"]
                    gene_meta[g] = (row.get("gene_name", ""), row.get("annotation", ""),
                                     row.get("n_isolates_with_gene", ""))
                    if dens not in ("NA", ""):
                        try:
                            gene_vals[g].append(float(dens))
                        except ValueError:
                            pass

        gene_rows = []
        for g, vals in gene_vals.items():
            gene, ann, n_with = gene_meta.get(g, ("", "", ""))
            mean = statistics.mean(vals)
            cv = (statistics.pstdev(vals) / mean) if mean > 0 and len(vals) > 1 else 0.0
            text = (gene + " " + ann).lower()
            flagged = any(h in text for h in highlight) if highlight else False
            gene_rows.append([g, gene, ann, n_with, len(vals),
                               "%.3f" % mean, "%.3f" % cv, "yes" if flagged else ""])

        # Highlighted genes surface first regardless of CV, then CV, then mean
        # density. CV alone ties constantly on sparse data: it's scale-invariant,
        # so any gene where the same handful of isolates out of the panel carry
        # a nonzero value gets an identical CV regardless of how large that value
        # actually is -- e.g. a gene spiking in 1/14 isolates ties with every
        # other gene that also spikes in exactly 1/14, no matter the magnitude.
        # Breaking those ties by mean density at least orders them by how much
        # signal is actually there, rather than by dict-iteration order.
        gene_rows.sort(key=lambda r: (r[7] != "yes", -float(r[6]), -float(r[5])))

        with open(gene_out, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t", lineterminator="\n")
            w.writerow(["ortholog_group", "gene_name", "annotation",
                        "n_isolates_with_gene", "n_isolates_with_signal",
                        "mean_genic_sites_per_kb", "cross_isolate_cv", "highlighted"])
            w.writerows(gene_rows)

        sys.stderr.write("%d motifs, %d genes summarised across %d isolates\n"
                         % (len(motif_rows), len(gene_rows), n_isolates_seen))
        PY
        # No de-indenting here: WDL strips the command block's common
        # leading whitespace before the shell ever sees it, so the heredoc
        # lands with the python's relative indentation already correct.

        python3 summarise.py \
            "~{sep=',' per_isolate_landscape}" \
            "~{select_first([ortholog_long_table, '-'])}" \
            "~{sep=',' select_first([highlight_genes, []])}" \
            ~{basename}_by_motif.tsv \
            ~{basename}_by_gene.tsv
    >>>

    output {
        File motif_summary = "~{basename}_by_motif.tsv"
        File gene_summary   = "~{basename}_by_gene.tsv"
    }

    runtime {
        docker:         docker
        memory:         "8 GB"
        cpu:            1
        disks:          "local-disk 50 SSD"
        preemptible:    1
        maxRetries:     2
    }
}
