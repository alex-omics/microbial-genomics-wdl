version 1.0

task match_samples_to_references {

    input {
        Array[String]  sample_ids
        Array[String]  reference_ids
        Array[Int]     sample_array_lengths
        Array[Int]     reference_array_lengths
        String         mode                   = "match"
        File?          sample_reference_map
        Int            max_strip_depth        = 1
        String?        replicate_regex
        File?          ortholog_table
        String         docker                 = "python:3.11-slim@sha256:db3ff2e1800a8581e2c48a27c3995339d47bdf046da21c7627accd3d51053a93"
    }

    parameter_meta {
        sample_ids:              "Sequencing sample names, e.g. one per FASTQ pair. Replicates of one isolate are expected to be named as the isolate plus a replicate suffix."
        reference_ids:           "Names of the reference genomes / assemblies the samples may be aligned to, positionally matched to the reference FASTA (and GFF) arrays."
        sample_array_lengths:    "Lengths of every per-sample array (read1, read2, ...) that must equal length(sample_ids)"
        reference_array_lengths: "Lengths of every per-reference array (FASTAs, GFFs) that must equal length(reference_ids)"
        mode:                    "'match' (default) assigns each sample to the reference its name resolves to. 'all_to_one' assigns every sample to the single supplied reference (the shared-reference case, e.g. everything to PA14)."
        sample_reference_map:    "Optional two-column TSV (sample_id, reference_id; no header) that overrides name matching for the samples it lists. Anything not listed still goes through name matching."
        max_strip_depth:         "How many replicate suffixes may be peeled off a sample name while looking for its reference (default = 1). Deliberately low: each extra level makes a wrong-but-existing match more likely."
        replicate_regex:         "Optional Python regex, matched at the END of a sample name, that replaces the built-in replicate-suffix grammar. Use this for a lab-specific convention the defaults do not cover."
        ortholog_table:          "Optional Panaroo gene_presence_absence.csv that was supplied rather than built. Checked here, up front, so a bad table fails the run before any alignment rather than at the very end. Not needed for the workflow to build its own."
        docker:                  "Container image"
    }

    meta {
        description: "Resolve each RNA-seq sample to the assembly it should be aligned against, and fail fast if any cannot be resolved. Replicates are usually named as the isolate plus a suffix (1/2/3, a/b/c, rep2, ...), so a sample is matched by exact name first, then by peeling replicate suffixes off its name until it lands on a name that is actually in reference_ids. Matching is only ever against real reference names and never falls through to a guess: an unmatched sample stops the run before any alignment is billed. Emits the mapping as a table so every assignment can be audited, plus the indices the workflow uses to scatter over only the references that are actually used."
    }

    command <<<
        set -euo pipefail

        cat > match.py <<'PY'
        import csv, difflib, re, sys

        (samples_f, refs_f, mode, depth, override_f, user_regex,
         sample_lens, ref_lens, gpa_f) = sys.argv[1:10]

        samples = [l.rstrip("\n") for l in open(samples_f)]
        refs    = [l.rstrip("\n") for l in open(refs_f)]
        depth   = int(depth)
        errors  = []

        # ---- structural checks (cheap, and the usual cause of an index error deep in a scatter)
        for what, n, lens in (("sample", len(samples), sample_lens), ("reference", len(refs), ref_lens)):
            for c in [x for x in lens.split(",") if x != ""]:
                if int(c) != n:
                    errors.append("a per-%s input array has %s entries but there are %d %s names; "
                                  "these arrays are matched positionally" % (what, c, n, what))
        if not samples:
            errors.append("no sample_ids supplied")
        if not refs:
            errors.append("no reference_ids supplied")
        if mode not in ("match", "all_to_one"):
            errors.append("mode must be 'match' or 'all_to_one', got %r" % mode)
        if mode == "all_to_one" and len(refs) != 1:
            errors.append("mode 'all_to_one' needs exactly one reference, got %d" % len(refs))

        # Names end up in filenames, read groups and matrix headers.
        bad = re.compile(r"[^A-Za-z0-9._-]")
        for label, names in (("sample", samples), ("reference", refs)):
            for n in names:
                if bad.search(n) or n == "":
                    errors.append("%s name %r contains characters outside [A-Za-z0-9._-]" % (label, n))
            dup = sorted({n for n in names if names.count(n) > 1})
            if dup:
                errors.append("duplicate %s names: %s" % (label, ", ".join(dup)))
        if errors:
            sys.stderr.write("ERROR:\n  " + "\n  ".join(errors) + "\n")
            sys.exit(1)

        # ---- replicate-suffix grammar
        # A digit-ending suffix must carry a separator or a keyword. Without that
        # rule, isolate numbers would be read as replicates ("PSA_11" -> "PSA_1"),
        # which produces a wrong-but-existing match with no error.
        # A letter-ending suffix is a single letter. A bare letter (no separator,
        # no keyword) only counts when it directly follows a digit ("PSA_1b").
        KW = r"(?:replicate|rep|bio|tech)"
        if user_regex:
            rules = [("custom", re.compile("(?:" + user_regex + ")$"))]
        else:
            rules = [
                # "_R2" is accepted with a separator only; unseparated, "R" + digits is too
                # likely to be part of an isolate name.
                ("numeric", re.compile(r"(?:[._-]+(?:" + KW + r"|r)?[._-]*|" + KW + r"[._-]*)[0-9]{1,2}$", re.I)),
                ("letter",  re.compile(r"(?:[._-]+" + KW + r"?[._-]*|" + KW + r"[._-]*|(?<=[0-9]))[A-Za-z]$", re.I)),
            ]

        def strip_once(name):
            """The (rule, stripped name) candidates one suffix-strip away. Digit-ending
            and letter-ending rules cannot both match, so this is a single chain."""
            out = []
            for rname, rx in rules:
                m = rx.search(name)
                if m and m.start() > 0:
                    out.append((rname, name[:m.start()]))
            return out

        # ---- explicit overrides
        override = {}
        if override_f:
            for ln, line in enumerate(open(override_f), 1):
                line = line.rstrip("\r\n")
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.split("\t")
                if len(f) != 2:
                    errors.append("sample_reference_map line %d: expected 2 tab-separated columns" % ln)
                    continue
                if f[0] in override and override[f[0]] != f[1]:
                    errors.append("sample_reference_map: %r is listed twice with different references" % f[0])
                override[f[0]] = f[1]
            for s, r in override.items():
                if s not in samples:
                    errors.append("sample_reference_map: %r is not one of the sample_ids" % s)
                if r not in refs:
                    errors.append("sample_reference_map: reference %r (for %r) is not one of the reference_ids" % (r, s))

        # ---- match
        refset = set(refs)
        rows, unmatched = [], []
        for s in samples:
            note = ""
            if s in override:
                ref, rule = override[s], "override"
            elif mode == "all_to_one":
                ref, rule = refs[0], "all_to_one"
            elif s in refset:
                ref, rule = s, "exact"
                # Both "PSA_1b" and "PSA_1" being assemblies is legitimate, and exact
                # wins, but it is the one case where a replicate could be misfiled.
                stripped = [c for _, c in strip_once(s) if c in refset]
                if stripped:
                    note = "also_strips_to:" + stripped[0]
            else:
                ref = rule = None
                cur = s
                for level in range(1, depth + 1):
                    cands = strip_once(cur)
                    if not cands:
                        break
                    rname, cur = cands[0]
                    if cur in refset:
                        ref, rule = cur, "strip%d:%s" % (level, rname)
                        break
            if ref is None:
                near = difflib.get_close_matches(s, refs, n=3, cutoff=0.5)
                unmatched.append((s, near))
                continue
            rows.append((s, ref, refs.index(ref), rule, note))

        if unmatched or errors:
            msg = list(errors)
            if unmatched:
                msg.append("%d of %d samples did not resolve to any reference:" % (len(unmatched), len(samples)))
                for s, near in unmatched:
                    msg.append("    %s   (closest reference names: %s)" %
                               (s, ", ".join(near) if near else "none"))
                msg.append("Fix the names, raise max_strip_depth, supply replicate_regex, or list the "
                           "sample in sample_reference_map.")
            sys.stderr.write("ERROR:\n  " + "\n  ".join(msg) + "\n")
            sys.exit(1)

        used = sorted({r[2] for r in rows})
        pos  = {idx: p for p, idx in enumerate(used)}

        # ---- ortholog table: fail here, not after the alignments have been paid for
        used_names = [refs[i] for i in used]
        if gpa_f:
            with open(gpa_f, newline="", encoding="utf-8", errors="replace") as fh:
                header = next(csv.reader(fh), [])
            problems = []
            if len(header) < 2 and header and "\t" in header[0]:
                problems.append("this looks like the tab-separated gene_presence_absence.Rtab. The workflow "
                                "needs gene_presence_absence.csv, which holds each isolate's locus tags.")
            else:
                if "Gene" not in header:
                    problems.append("no 'Gene' column: this does not look like Panaroo's gene_presence_absence.csv")
                missing = [r for r in used_names if r not in header]
                if missing:
                    problems.append("these isolates have no column in the ortholog table: " + ", ".join(missing) +
                                    "\n    Panaroo names columns after the input GFF filenames, which must equal reference_ids."
                                    "\n    Columns present: " + ", ".join(header[:14]) + (" ..." if len(header) > 14 else ""))
            if problems:
                sys.stderr.write("ERROR: ortholog_table:\n  " + "\n  ".join(problems) + "\n")
                sys.exit(1)

        with open("mapping.tsv", "w") as fh:
            fh.write("sample_id\treference_id\treference_index\trule\tnote\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        with open("used_reference_idx.txt", "w") as fh:
            fh.write("".join("%d\n" % i for i in used))
        with open("sample_to_used_pos.txt", "w") as fh:
            fh.write("".join("%d\n" % pos[r[2]] for r in rows))

        unused = [refs[i] for i in range(len(refs)) if i not in pos]
        if unused:
            sys.stderr.write("NOTE: %d reference(s) have no samples and will be skipped: %s\n"
                             % (len(unused), ", ".join(unused)))
        for r in rows:
            if r[4]:
                sys.stderr.write("WARNING: %s matched exactly, but %s\n" % (r[0], r[4]))
        PY

        python3 match.py \
            ~{write_lines(sample_ids)} \
            ~{write_lines(reference_ids)} \
            "~{mode}" \
            ~{max_strip_depth} \
            "~{select_first([sample_reference_map, ''])}" \
            "~{select_first([replicate_regex, ''])}" \
            "~{sep=',' sample_array_lengths}" \
            "~{sep=',' reference_array_lengths}" \
            "~{select_first([ortholog_table, ''])}"
    >>>

    output {
        File        mapping             = "mapping.tsv"
        Array[Int]  used_reference_idx  = read_lines("used_reference_idx.txt")
        Array[Int]  sample_to_used_pos  = read_lines("sample_to_used_pos.txt")
    }

    runtime {
        docker:         docker
        memory:         "2 GB"
        cpu:            1
        disks:          "local-disk 10 SSD"
        preemptible:    1
        maxRetries:     0
    }
}
