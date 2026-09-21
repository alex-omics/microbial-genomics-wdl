version 1.0

task featurecounts {

    input {
        Array[File]    input_bams
        Array[String]  sample_ids
        File           annotation
        String         reference_name      = "counts"
        String         strandness          = "2"
        String         feature_type        = "CDS"
        String         attribute_type      = "locus_tag"
        String         annotation_format   = "GFF"
        Boolean        paired_end          = true
        Boolean        count_read_pairs    = true
        Boolean        require_both_mates  = true
        Boolean        count_chimeric      = false
        Boolean        ignore_duplicates   = false
        Boolean        fraction_counting   = false
        Int            min_overlap         = 1
        Int            cpu                 = 8
        Int            mem_gb              = 16
        Int            disk_gb             = 50
        String         docker              = "quay.io/biocontainers/subread@sha256:114390a783c77f7739d86e474bedfa5a4e65309a2f71d4db430803fb04601f5d"
    }

    parameter_meta {
        input_bams:         "BAMs to count, all aligned to the reference that `annotation` describes"
        sample_ids:         "Sample names, positionally matched to input_bams. Used to rename the count matrix columns, which featureCounts otherwise labels with the BAM path."
        annotation:         "GFF3/GTF for the reference these BAMs were aligned to"
        reference_name:     "Names the outputs; for per-isolate counting, the isolate"
        strandness:         "0 = unstranded, 1 = stranded, 2 = reverse-stranded (default = 2)"
        feature_type:       "Feature type to count: CDS for bacteria, exon for eukaryotes (default = CDS)"
        attribute_type:     "Attribute to group features by: locus_tag for bacterial GFF, gene_id for GTF (default = locus_tag)"
        annotation_format:  "GFF or GTF (default = GFF)"
        paired_end:         "Input is paired-end (default = true)"
        count_read_pairs:   "Count fragments rather than reads (--countReadPairs). Since subread 2.0.2, -p alone counts each mate separately, double-counting every pair (default = true)."
        require_both_mates: "Only count pairs where both mates aligned (-B; paired-end only) (default = true)"
        count_chimeric:     "Count chimeric pairs. If false they are discarded (-C) (default = false)"
        ignore_duplicates:  "Ignore reads flagged as duplicates (--ignoreDup) (default = false). Off because, without UMIs, RNA-seq duplicates are mostly independent molecules from short, highly expressed transcripts; excluding them removes the most reads from the most expressed genes."
        fraction_counting:  "Split multi-overlapping reads fractionally (--fraction) (default = false)"
        min_overlap:        "Minimum bases overlapping a feature (default = 1)"
        cpu:                "Number of CPUs delegated to task (default = 8)"
        mem_gb:             "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:            "Amount of disk space in GB delegated to task (default = 50)"
        docker:             "Container image"
    }

    meta {
        description: "Count reads per feature across a group of BAMs that share one reference, and return a features x samples matrix with columns named by sample rather than by BAM path."
    }

    command <<<
        set -euo pipefail
        featureCounts -v 2>&1 | grep -i featurecounts | tee FEATURECOUNTS_VERSION || true

        # Read a file that may or may not be gzipped. `zcat -f` is not portable
        # across the images used here (some zcats reject plain text outright).
        plain() {
            if [ "$(head -c2 "$1" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]; then gzip -dc "$1"; else cat "$1"; fi
        }

        # A Bakta GFF3 ends with a ##FASTA block of raw sequence, which is not
        # feature lines. Drop it rather than rely on the counter skipping it.
        plain ~{annotation} | awk '/^##FASTA/{exit} {print}' > annotation.gff

        fc_flags=()
        if [ "~{paired_end}" == "true" ]; then fc_flags+=("-p"); fi
        # subread >= 2.0.2: -p alone counts reads, so every mate is counted.
        if [ "~{paired_end}" == "true" ] && [ "~{count_read_pairs}" == "true" ]; then fc_flags+=("--countReadPairs"); fi
        if [ "~{paired_end}" == "true" ] && [ "~{require_both_mates}" == "true" ]; then fc_flags+=("-B"); fi
        if [ "~{count_chimeric}" == "false" ]; then fc_flags+=("-C"); fi
        if [ "~{ignore_duplicates}" == "true" ]; then fc_flags+=("--ignoreDup"); fi
        if [ "~{fraction_counting}" == "true" ]; then fc_flags+=("--fraction"); fi

        featureCounts \
            -T ~{cpu} \
            "${fc_flags[@]+"${fc_flags[@]}"}" \
            -s ~{strandness} \
            -t ~{feature_type} \
            -g ~{attribute_type} \
            -F ~{annotation_format} \
            --minOverlap ~{min_overlap} \
            -a annotation.gff \
            -o raw_counts.txt \
            ~{sep=" " input_bams}

        # featureCounts labels sample columns with the BAM path, in command-line
        # order. Rename them by position to the sample ids, and check each path
        # really contains its sample id so a reordering cannot silently swap
        # samples. The '# Program' comment line is dropped.
        rename() {
            awk -F'\t' -v OFS='\t' '
                NR==FNR { id[FNR]=$1; n=FNR; next }
                /^#/    { next }
                !hdr    { hdr=1
                          if (NF-fixed != n) { print "ERROR: expected " n " sample columns, found " (NF-fixed) > "/dev/stderr"; exit 1 }
                          for (i=1; i<=n; i++) {
                              path=$(fixed+i); sub(/.*\//, "", path)
                              if (index(path, id[i] ".") != 1) { print "ERROR: column " i " is " path ", expected sample " id[i] > "/dev/stderr"; exit 1 }
                              $(fixed+i)=id[i]
                          }
                        }
                { print }
            ' fixed="$2" ~{write_lines(sample_ids)} "$1"
        }
        rename raw_counts.txt 6 > "~{reference_name}.counts.tsv"
        rename raw_counts.txt.summary 1 > "~{reference_name}.counts.summary.tsv"
    >>>

    output {
        File count_matrix = "~{reference_name}.counts.tsv"
        File summary      = "~{reference_name}.counts.summary.tsv"
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


task merge_count_matrices {

    input {
        Array[File]    matrices
        Array[String]  reference_ids
        File?          ortholog_table
        String         basename    = "counts"
        Int            cpu         = 2
        Int            mem_gb      = 16
        Int            disk_gb     = 20
        String         docker      = "python:3.11-slim@sha256:db3ff2e1800a8581e2c48a27c3995339d47bdf046da21c7627accd3d51053a93"
    }

    parameter_meta {
        matrices:       "Per-reference count matrices from featurecounts, in the same order as reference_ids"
        reference_ids:  "Reference (isolate) each matrix was counted against"
        ortholog_table: "Optional Panaroo gene_presence_absence.csv. Puts samples from different isolates on common rows, since each isolate's features carry its own locus tags. Its isolate column headers must equal reference_ids."
        basename:       "Basename for the outputs"
        cpu:            "Number of CPUs delegated to task (default = 2)"
        mem_gb:         "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:        "Amount of disk space in GB delegated to task (default = 20)"
        docker:         "Container image"
    }

    meta {
        description: "Combine per-isolate count matrices. Always emits a long-format table (one row per isolate x sample x feature). Emits a single wide matrix when one exists: the matrix itself for a single reference, or, for several isolates, an ortholog-group matrix when an ortholog table is supplied. A gene an isolate lacks is NA, not 0."
    }

    command <<<
        set -euo pipefail

        cat > merge.py <<'PY'
        import csv, sys
        from collections import defaultdict

        basename, gpa_path = sys.argv[1], sys.argv[2]
        refs   = sys.argv[3].split(",")
        paths  = sys.argv[4:]
        if len(refs) != len(paths):
            sys.exit("ERROR: %d reference ids but %d matrices" % (len(refs), len(paths)))

        FIXED = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
        samples_of = {}                  # ref -> [sample ids]
        counts     = {}                  # (ref, sample) -> {feature: count}
        features   = {}                  # ref -> [(feature, chr, start, end, strand, length)]
        seen       = set()

        with open(basename + "_long.tsv", "w") as out:
            out.write("reference_id\tsample_id\tfeature_id\tchr\tstart\tend\tstrand\tlength\tcount\n")
            for ref, path in zip(refs, paths):
                with open(path, newline="") as fh:
                    rd = csv.reader(fh, delimiter="\t")
                    hdr = next(rd)
                    if hdr[:6] != FIXED:
                        sys.exit("ERROR: %s does not look like a featureCounts matrix (header %r)" % (path, hdr[:6]))
                    samples_of[ref] = hdr[6:]
                    features[ref] = []
                    for s in hdr[6:]:
                        if s in seen:
                            sys.exit("ERROR: sample %r appears in more than one matrix" % s)
                        seen.add(s)
                        counts[(ref, s)] = {}
                    for row in rd:
                        f = tuple(row[:6])
                        features[ref].append(f)
                        for s, c in zip(hdr[6:], row[6:]):
                            counts[(ref, s)][f[0]] = c
                            out.write("\t".join([ref, s, f[0], f[1], f[2], f[3], f[4], f[5], c]) + "\n")

        n_samples = sum(len(v) for v in samples_of.values())
        sys.stderr.write("%d samples across %d reference(s)\n" % (n_samples, len(refs)))

        if gpa_path == "":
            if len(refs) == 1:
                ref = refs[0]
                with open(basename + "_matrix.tsv", "w") as out:
                    out.write("\t".join(FIXED + samples_of[ref]) + "\n")
                    for f in features[ref]:
                        out.write("\t".join(list(f) + [counts[(ref, s)][f[0]] for s in samples_of[ref]]) + "\n")
            else:
                sys.stderr.write(
                    "NOTE: %d references, no ortholog_table: no single wide matrix emitted. Locus tags\n"
                    "      differ between isolates, so samples from different isolates have no shared\n"
                    "      rows. Use the long table, the per-isolate matrices, or supply ortholog_table.\n" % len(refs))
            sys.exit(0)

        # ---- ortholog-level matrix
        # Panaroo's csv has free-text annotations with commas: use a real parser.
        with open(gpa_path, newline="", encoding="utf-8", errors="replace") as fh:
            rd = csv.DictReader(fh)
            header = rd.fieldnames or []
            missing = [r for r in refs if r not in header]
            if missing:
                sys.exit("ERROR: these references have no column in the ortholog table: %s\n"
                         "Panaroo names columns after the input GFF filenames, which must equal the\n"
                         "reference ids. Columns present: %s" % (", ".join(missing), ", ".join(header[:12])))
            groups = []
            for row in rd:
                groups.append(row)

        # Panaroo was run on some annotation of each isolate; if that is not the GFF the
        # reads were counted against, no locus tag matches and every cell would come out
        # NA with no error. Count matches per isolate, and refuse if an isolate has none.
        in_table = {}
        for r in refs:
            tags = {t.strip() for row in groups for t in (row.get(r) or "").split(";") if t.strip()}
            in_table[r] = tags
            counted = set(counts[(r, samples_of[r][0])])
            hit = len(tags & counted)
            sys.stderr.write("%s: %d of %d ortholog-table locus tags are counted features\n" % (r, hit, len(tags)))
            if tags and hit == 0:
                sys.exit("ERROR: none of %s's %d locus tags in the ortholog table appear in its count matrix.\n"
                         "The ortholog table and the annotation used for counting must come from the same\n"
                         "GFFs (e.g. locus tags such as %s vs %s)." % (r, len(tags), sorted(tags)[0], sorted(counted)[0]))

        all_samples = [s for r in refs for s in samples_of[r]]
        n_absent = n_unquant = n_cells = 0
        with open(basename + "_by_ortholog_matrix.tsv", "w") as out:
            out.write("\t".join(["ortholog_group", "gene_name", "annotation"] + all_samples) + "\n")
            for row in groups:
                vals = []
                for r in refs:
                    cell = (row.get(r) or "").strip()
                    tags = [t.strip() for t in cell.split(";") if t.strip()]
                    for s in samples_of[r]:
                        n_cells += 1
                        if not tags:
                            vals.append("NA"); n_absent += 1; continue
                        found = [counts[(r, s)][t] for t in tags if t in counts[(r, s)]]
                        if not found:
                            # Present in the pangenome but not a counted feature in this
                            # isolate (wrong feature_type, or a Panaroo "refound" gene
                            # with no GFF line): not quantified, so not 0.
                            vals.append("NA"); n_unquant += 1; continue
                        # Paralogues in one group are summed.
                        tot = sum(float(x) for x in found)
                        vals.append(str(int(tot)) if tot == int(tot) else "%.6g" % tot)
                out.write("\t".join([row.get("Gene", ""), row.get("Non-unique Gene name", "") or "",
                                     row.get("Annotation", "") or ""] + vals) + "\n")
            # A counted feature that is in NO ortholog group has not been shown to be absent
            # anywhere, it has just not been clustered - Panaroo's stricter clean modes drop
            # genes seen in few genomes as likely annotation error, and those are precisely
            # the accessory genes this workflow exists to keep. Dropping them here would lose
            # their expression silently, so each gets its own row (NA for other isolates).
            n_uncl = {}
            for r in refs:
                n_uncl[r] = 0
                for f in features[r]:
                    tag = f[0]
                    if tag in in_table[r]:
                        continue
                    n_uncl[r] += 1
                    vals = []
                    for rr in refs:
                        for s in samples_of[rr]:
                            vals.append(counts[(rr, s)][tag] if rr == r else "NA")
                    out.write("\t".join(["unclustered:%s:%s" % (r, tag), "", "(not in ortholog table)"] + vals) + "\n")
        sys.stderr.write("ortholog matrix: %d groups x %d samples; %d cells NA (gene absent), "
                         "%d NA (present but not a counted feature)\n"
                         % (len(groups), len(all_samples), n_absent, n_unquant))
        for r in refs:
            if n_uncl[r]:
                sys.stderr.write("NOTE: %s: %d of %d counted features are in no ortholog group and are emitted as "
                                 "'unclustered' rows. If this is many, the table was built with a strict Panaroo "
                                 "clean mode (or from different GFFs).\n" % (r, n_uncl[r], len(features[r])))
        PY

        python3 merge.py "~{basename}" "~{select_first([ortholog_table, ''])}" "~{sep=',' reference_ids}" ~{sep=' ' matrices}
    >>>

    output {
        File   counts_long   = "~{basename}_long.tsv"
        File?  counts_matrix = if defined(ortholog_table) then "~{basename}_by_ortholog_matrix.tsv" else "~{basename}_matrix.tsv"
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
