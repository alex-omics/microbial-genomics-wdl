version 1.0

task methylation_orthologs {

    input {
        File         gene_presence_absence
        Array[File]  annotated_tables
        Array[String] sample_names
        String       basename    = "methylation_by_ortholog"
        Int          cpu         = 2
        Int          mem_gb      = 16
        Int          disk_gb     = 50
        String       docker      = "python:3.11-slim"
    }

    parameter_meta {
        gene_presence_absence: "Panaroo gene_presence_absence.csv. Its per-isolate columns hold that isolate's own locus tags for each ortholog group, which is what lets independently-annotated genomes be compared."
        annotated_tables:      "Per-isolate methylation tables from annotate_methylation, keyed on each isolate's own locus tags"
        sample_names:          "Isolate names, matching both the annotated_tables order and the column headers Panaroo derived from the input GFF filenames"
        basename:              "Basename for the output tables"
        cpu:                   "Number of CPUs delegated to task (default = 2)"
        mem_gb:                "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:               "Amount of disk space in GB delegated to task (default = 50)"
        docker:                "Container image"
    }

    meta {
        description: "Collapse per-isolate methylation calls onto pangenome ortholog groups, producing a gene-by-isolate matrix. Absent genes are emitted as NA rather than 0, because a gene an isolate does not carry has undefined methylation — conflating that with an unmethylated gene would manufacture differences that are really presence/absence."
    }

    command <<<
        set -euo pipefail

        cat > join.py <<'PY'
        import csv, sys, os
        from collections import defaultdict

        gpa_path = sys.argv[1]
        basename = sys.argv[2]
        samples  = sys.argv[3].split(",")
        tables   = sys.argv[4:]

        # gene_presence_absence.csv carries free-text annotations containing
        # commas, so it must go through a real CSV parser. Splitting on commas
        # silently shifts every column after the annotation field.
        with open(gpa_path, newline="", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh)
            header = reader.fieldnames or []
            missing = [s for s in samples if s not in header]
            if missing:
                sys.stderr.write(
                    "ERROR: these samples have no column in gene_presence_absence.csv: "
                    + ", ".join(missing) + "\n"
                    "Panaroo names columns after the input GFF filenames; they must\n"
                    "match the sample names used elsewhere in the workflow.\n"
                    "Columns present: " + ", ".join(header[:12]) + "\n"
                )
                sys.exit(1)

            # locus tag -> ortholog group, plus per-group metadata and the set of
            # isolates actually carrying the gene.
            tag2group = {}
            meta      = {}
            present   = defaultdict(set)

            for row in reader:
                group = row["Gene"]
                meta[group] = (
                    row.get("Non-unique Gene name", "") or "",
                    row.get("Annotation", "") or "",
                )
                for s in samples:
                    cell = (row.get(s) or "").strip()
                    if not cell:
                        continue
                    present[group].add(s)
                    # Paralogues and split genes appear semicolon-separated.
                    for tag in cell.replace("\t", " ").split(";"):
                        tag = tag.strip()
                        if tag:
                            tag2group[(s, tag)] = group

        sys.stderr.write("Mapped %d locus tags into %d ortholog groups\n"
                         % (len(tag2group), len(meta)))

        # Aggregate methylation per (group, sample).
        agg = defaultdict(lambda: {"genic": 0, "upstream": 0,
                                   "6mA": 0, "4mC": 0, "5mC": 0,
                                   "pct_sum": 0.0, "n": 0, "tags": set()})
        CODE = {"a": "6mA", "21839": "4mC", "m": "5mC"}
        unmatched = 0

        for path in tables:
            with open(path, newline="", encoding="utf-8", errors="replace") as fh:
                rd = csv.reader(fh, delimiter="\t")
                next(rd, None)
                for r in rd:
                    if len(r) < 13:
                        continue
                    sample, mod_code, pct, tag, region = r[0], r[5], r[7], r[8], r[12]
                    group = tag2group.get((sample, tag))
                    if group is None:
                        unmatched += 1
                        continue
                    a = agg[(group, sample)]
                    a[region if region in ("genic", "upstream") else "genic"] += 1
                    label = CODE.get(mod_code)
                    if label:
                        a[label] += 1
                    try:
                        a["pct_sum"] += float(pct); a["n"] += 1
                    except ValueError:
                        pass
                    a["tags"].add(tag)

        if unmatched:
            sys.stderr.write("WARNING: %d methylation rows had locus tags absent from the "
                             "pangenome (usually genes Panaroo filtered)\n" % unmatched)

        groups = sorted(meta)

        with open(basename + "_long.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["ortholog_group", "gene_name", "annotation",
                        "n_isolates_with_gene", "sample", "locus_tags",
                        "n_genic_sites", "n_upstream_sites",
                        "n_6mA", "n_4mC", "n_5mC", "mean_percent_modified"])
            for g in groups:
                gene, ann = meta[g]
                for s in samples:
                    if s not in present[g]:
                        continue
                    a = agg.get((g, s))
                    if a is None:
                        w.writerow([g, gene, ann, len(present[g]), s, "",
                                    0, 0, 0, 0, 0, "NA"])
                    else:
                        mean = ("%.2f" % (a["pct_sum"] / a["n"])) if a["n"] else "NA"
                        w.writerow([g, gene, ann, len(present[g]), s,
                                    ";".join(sorted(a["tags"])),
                                    a["genic"], a["upstream"],
                                    a["6mA"], a["4mC"], a["5mC"], mean])

        # Wide matrix of total methylated sites per gene per isolate.
        #
        # NA and 0 mean different things here and must not be merged: NA is
        # "this isolate does not carry the gene", 0 is "it carries it and no
        # methylation was called". Treating absence as zero would invent
        # methylation differences that are really gene presence/absence.
        with open(basename + "_matrix.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["ortholog_group", "gene_name", "annotation",
                        "n_isolates_with_gene"] + samples)
            for g in groups:
                gene, ann = meta[g]
                row = [g, gene, ann, len(present[g])]
                for s in samples:
                    if s not in present[g]:
                        row.append("NA")
                    else:
                        a = agg.get((g, s))
                        row.append(0 if a is None else a["genic"] + a["upstream"])
                w.writerow(row)

        n_core = sum(1 for g in groups if len(present[g]) == len(samples))
        with open("N_GROUPS", "w") as fh:
            fh.write(str(len(groups)) + "\n")
        with open("N_CORE", "w") as fh:
            fh.write(str(n_core) + "\n")
        with open("N_METHYLATED_GROUPS", "w") as fh:
            fh.write(str(len({g for (g, _s) in agg})) + "\n")
        PY

        # The heredoc above is indented for readability inside the WDL block;
        # strip that indentation before handing the file to python.
        sed -i 's/^        //' join.py

        python3 join.py \
            ~{gene_presence_absence} \
            ~{basename} \
            "~{sep=',' sample_names}" \
            ~{sep=' ' annotated_tables}

        gzip -c ~{basename}_long.tsv > ~{basename}_long.tsv.gz
    >>>

    output {
        File  long_table            = "~{basename}_long.tsv"
        File  long_table_gz         = "~{basename}_long.tsv.gz"
        File  matrix                = "~{basename}_matrix.tsv"
        Int   n_ortholog_groups     = read_int("N_GROUPS")
        Int   n_core_groups         = read_int("N_CORE")
        Int   n_methylated_groups   = read_int("N_METHYLATED_GROUPS")
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
