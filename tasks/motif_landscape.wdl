version 1.0

task build_motif_list {

    input {
        File?  find_motifs_tsv
        File?  rebase_mtases_tsv
        String basename = "motif_list"
        String docker    = "python:3.11-slim@sha256:db3ff2e1800a8581e2c48a27c3995339d47bdf046da21c7627accd3d51053a93"
    }

    parameter_meta {
        find_motifs_tsv:   "modkit find-motifs output (mod_code, motif, offset, frac_mod, ...). Either input may be omitted, not both."
        rebase_mtases_tsv: "rebase_mtase_search output (sample, locus_tag, rebase_enzyme, ..., predicted_motif, modification_type)"
        basename:          "Basename for the merged motif list"
        docker:            "Container image"
    }

    meta {
        description: "Merge de novo (find-motifs) and homology-predicted (REBASE) motifs into one deduplicated (motif, mod_code) list for motif_landscape to test. A motif discovered both ways appears once. Carries find-motifs' offset (the index of the modified base within the motif) so motif_landscape can score that base alone; a motif with no usable offset, such as one known only from REBASE, is written with offset NA and scored by window instead."
    }

    command <<<
        set -euo pipefail

        cat > merge.py <<'PY'
        import csv, sys

        find_motifs_path   = sys.argv[1] if sys.argv[1] != "-" else None
        rebase_path        = sys.argv[2] if sys.argv[2] != "-" else None
        out_path           = sys.argv[3]

        # bedMethyl / find-motifs use short codes; REBASE's TSV spells the
        # human name. Accept either so this does not silently drop rows if a
        # future modkit version changes its own convention.
        TO_CODE = {"a": "a", "6ma": "a", "6mA": "a",
                   "m": "m", "5mc": "m", "5mC": "m",
                   "21839": "21839", "4mc": "21839", "4mC": "21839"}

        # IUPAC codes that can stand for the base each modification sits on:
        # A for 6mA, C for 4mC and 5mC.
        CAN_BE = {"a": set("ARWMDHVN"), "m": set("CYSMBHVN"), "21839": set("CYSMBHVN")}

        # (motif, code) -> offsets known for it. An empty set means none.
        offsets = {}

        def note(motif, code, offset=None):
            offsets.setdefault((motif, code), set())
            if offset is None:
                return
            if offset < 0 or offset >= len(motif) or motif[offset] not in CAN_BE[code]:
                # modkit's offset is the index of the modified base within the
                # motif as written. If it does not land on a base that
                # modification can sit on, the convention is not what this was
                # written against, and scoring on it would silently test the
                # wrong base. Fall back to window scoring for this motif.
                sys.stderr.write(
                    "WARNING: offset %d is not a plausible '%s' base in motif %s; "
                    "scoring this motif by window instead\n" % (offset, code, motif))
                return
            offsets[(motif, code)].add(offset)

        if find_motifs_path:
            with open(find_motifs_path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    motif = (row.get("motif") or "").strip().upper()
                    code  = TO_CODE.get((row.get("mod_code") or "").strip())
                    if not (motif and code):
                        continue
                    try:
                        off = int((row.get("offset") or "").strip())
                    except ValueError:
                        off = None
                    note(motif, code, off)

        if rebase_path:
            with open(rebase_path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    motif = (row.get("predicted_motif") or "").strip().upper()
                    code  = TO_CODE.get((row.get("modification_type") or "").strip())
                    if motif and motif != "NA" and code:
                        note(motif, code)

        with open(out_path, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t", lineterminator="\n")
            w.writerow(["motif", "mod_code", "offset"])
            for (motif, code), offs in sorted(offsets.items()):
                w.writerow([motif, code, ",".join(str(o) for o in sorted(offs)) if offs else "NA"])

        n_off = sum(1 for o in offsets.values() if o)
        sys.stderr.write("%d distinct (motif, mod_code) pairs, %d with a modified-base offset\n"
                         % (len(offsets), n_off))
        PY
        # No de-indenting here: WDL strips the command block's common
        # leading whitespace before the shell ever sees it, so the heredoc
        # lands with the python's relative indentation already correct.

        python3 merge.py \
            "~{select_first([find_motifs_tsv, '-'])}" \
            "~{select_first([rebase_mtases_tsv, '-'])}" \
            ~{basename}.tsv

        awk 'NR>1' ~{basename}.tsv | wc -l > N_MOTIFS
    >>>

    output {
        File motif_list = "~{basename}.tsv"
        Int  n_motifs    = read_int("N_MOTIFS")
    }

    runtime {
        docker:         docker
        memory:         "2 GB"
        cpu:            1
        disks:          "local-disk 10 SSD"
        preemptible:    1
        maxRetries:     2
    }
}


task motif_landscape {

    input {
        File     bedmethyl
        File     motif_list
        String   sample_name
        File     reference_fasta
        Int      min_coverage       = 10
        Float    heterogeneous_low_cutoff = 50.0
        Int      cpu                = 4
        Int      mem_gb             = 16
        Int      disk_gb            = 50
        String   docker             = "python:3.11-slim@sha256:db3ff2e1800a8581e2c48a27c3995339d47bdf046da21c7627accd3d51053a93"
    }

    parameter_meta {
        bedmethyl:                "Unfiltered bedMethyl from modkit_pileup. The enrichment test needs the true genome-wide background rate, which a coverage/percent-filtered table would have discarded."
        motif_list:               "Output of build_motif_list: one (motif, mod_code, offset) row per motif. offset is optional; without it the motif is scored by window."
        sample_name:              "Some identifier for naming outputs"
        reference_fasta:          "This isolate's own assembly, in the same coordinates as bedmethyl"
        min_coverage:             "Minimum Nvalid_cov for a site to be included in either the in-motif or background set (default = 10)"
        heterogeneous_low_cutoff: "A motif occurrence counts as 'low' in the heterogeneity summary if its representative percent-modified falls below this. Interpret relative to the panel's typical housekeeping level, usually near 100 (default = 50.0)"
        cpu:                      "Number of CPUs delegated to task (default = 4)"
        mem_gb:                   "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:                  "Amount of disk space in GB delegated to task (default = 50)"
        docker:                   "Container image"
    }

    meta {
        description: "Test each candidate motif against this isolate's own methylation data. Tier 1: whether methylation at the motif's occurrences is elevated over the genome-wide background for that base (log2 enrichment, chi-square). Tier 2: whether those occurrences are uniformly methylated or split into high and low populations, which suggests phase variation or a regulator competing for specific copies. Occurrences are located on both strands and IUPAC codes are supported. Scoring has two modes, reported per motif in the scoring column. 'offset': the motif has a known modified-base offset, so only that base (on the strand the motif reads 5'->3' on) is scored, which is what the motif's methylation actually is. 'window': no offset is known, so every base of the matching canonical type (A for 6mA, C for 4mC/5mC) within an occurrence counts. Window scoring dilutes the in-motif rate with the unmethylated bases the window also covers, spacer N's of a split motif above all, and so reads lower in a genome where that base is common (A in a low-GC genome) than the same methylation would elsewhere. The background excludes every candidate motif's scored sites, so methylated motifs do not inflate it."
    }

    command <<<
        set -euo pipefail

        cat > landscape.py <<'PY'
        import bisect, csv, math, re, statistics, sys
        from collections import defaultdict

        bedmethyl_path, motif_path, ref_path, sample, mincov, low_cut, out_path = sys.argv[1:8]
        mincov = int(mincov)
        low_cut = float(low_cut)

        IUPAC = {'A':'A','C':'C','G':'G','T':'T',
                 'R':'[AG]','Y':'[CT]','S':'[GC]','W':'[AT]','K':'[GT]','M':'[AC]',
                 'B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
        COMP = {'A':'T','T':'A','C':'G','G':'C','R':'Y','Y':'R','S':'S','W':'W',
                'K':'M','M':'K','B':'V','V':'B','D':'H','H':'D','N':'N'}

        def revcomp(m):
            return ''.join(COMP[b] for b in reversed(m))

        def to_regex(motif):
            return ''.join(IUPAC[b] for b in motif)

        def find_occurrences(seq, motif):
            # Start positions of the motif read on the plus strand, and of its
            # reverse complement, i.e. of the motif read on the minus strand.
            # A palindromic motif is its own reverse complement, so each of its
            # sites is both: GATC is methylated on both strands.
            fwd = sorted(m.start() for m in
                         re.finditer('(?=(' + to_regex(motif) + '))', seq))
            rc = revcomp(motif)
            if rc == motif:
                rev = list(fwd)
            else:
                rev = sorted(m.start() for m in
                             re.finditer('(?=(' + to_regex(rc) + '))', seq))
            return fwd, rev

        # Reference may hold multiple contigs; keep per-contig sequence and scan
        # each independently so motif windows never span a contig boundary.
        seqs = {}
        name, buf = None, []
        with open(ref_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(buf).upper()
                    name = line[1:].split()[0]
                    buf = []
                else:
                    buf.append(line)
            if name is not None:
                seqs[name] = "".join(buf).upper()

        def parse_offsets(field):
            field = (field or "").strip()
            if not field or field == "NA":
                return None
            return [int(x) for x in field.split(",")]

        motifs = []
        with open(motif_path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                motifs.append((row["motif"].strip().upper(), row["mod_code"].strip(),
                               parse_offsets(row.get("offset"))))

        # Per motif, how its in-motif rows are recognised.
        #   window mode: starts[contig] = sorted start positions; a row is in
        #     the motif if it falls inside any [start, start+len) window.
        #   offset mode: scored[(contig, strand, pos)] = the occurrence's start;
        #     a row is in the motif only if it is the modified base of one.
        # Every window for one motif shares the same length, so storing just
        # the start is enough to reconstruct membership.
        entries = []
        for motif, code, offs in motifs:
            motif_len = len(motif)
            starts = {}
            scored = {}
            n_occ = 0
            for contig, seq in seqs.items():
                fwd, rev = find_occurrences(seq, motif)
                allstarts = sorted(set(fwd) | set(rev))
                if allstarts:
                    starts[contig] = allstarts
                n_occ += len(allstarts)
                if offs is not None:
                    # Offset k is the modified base's index within the motif as
                    # written, on the strand the motif reads 5'->3' on. Read on
                    # the minus strand it lies at the mirrored index, and its
                    # bedMethyl row is the minus-strand one.
                    for s in fwd:
                        for k in offs:
                            scored[(contig, '+', s + k)] = s
                    for s in rev:
                        for k in offs:
                            scored[(contig, '-', s + motif_len - 1 - k)] = s
            entries.append({"motif": motif, "code": code, "offs": offs,
                            "len": motif_len, "starts": starts, "scored": scored,
                            "n_occ": n_occ})

        def in_any_window(contig, pos, starts_by_contig, motif_len):
            starts = starts_by_contig.get(contig)
            if not starts:
                return None
            # A window covers pos iff its start falls in
            # [pos-motif_len+1, pos]; binary-search that range instead of
            # scanning every occurrence. With a motif occurring thousands
            # of times and millions of bedMethyl rows to test, the linear
            # scan was O(sites x occurrences) -- this is what made
            # motif_landscape take up to two hours on some isolates.
            lo = pos - motif_len + 1
            i = bisect.bisect_left(starts, lo)
            if i < len(starts) and starts[i] <= pos:
                return starts[i]
            return None

        # bedMethyl: 0 chrom,1 start,3 mod_code,5 strand,9 Nvalid_cov,10 pct,11 Nmod
        by_code = defaultdict(list)
        for_code_entries = defaultdict(list)
        for e in entries:
            for_code_entries[e["code"]].append(e)
        with open(bedmethyl_path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                cov = int(f[9])
                if cov < mincov:
                    continue
                by_code[f[3]].append((f[0], int(f[1]), float(f[10]), int(f[11]), cov, f[5]))

        # A row that is a scored site of ANY candidate motif is not background
        # for any of them. Without this, a genome with several methylated
        # motifs tests each against a background carrying the others' sites.
        # One pass per modification code, so it costs about one extra motif.
        excluded = {}
        for code, sites in by_code.items():
            offset_sites = set()
            window_entries = []
            for e in for_code_entries.get(code, []):
                if e["offs"] is not None:
                    offset_sites.update(e["scored"])
                else:
                    window_entries.append(e)
            flags = []
            for contig, pos, pct, nmod, cov, strand in sites:
                hit = (contig, strand, pos) in offset_sites
                if not hit:
                    for e in window_entries:
                        if in_any_window(contig, pos, e["starts"], e["len"]) is not None:
                            hit = True
                            break
                flags.append(hit)
            excluded[code] = flags

        results = []
        for e in entries:
            motif, code, offs = e["motif"], e["code"], e["offs"]
            scoring = "offset" if offs is not None else "window"
            offset_str = ",".join(str(o) for o in offs) if offs is not None else "NA"
            sites = by_code.get(code, [])
            flags = excluded.get(code, [])
            in_pct, bg_pct = [], []
            in_mod = in_cov = bg_mod = bg_cov = 0
            occ_rep = {}  # (contig, occurrence start) -> max pct seen in that occurrence
            for (contig, pos, pct, nmod, cov, strand), is_excluded in zip(sites, flags):
                if offs is not None:
                    w = e["scored"].get((contig, strand, pos))
                else:
                    w = in_any_window(contig, pos, e["starts"], e["len"])
                if w is not None:
                    in_pct.append(pct); in_mod += nmod; in_cov += cov
                    key = (contig, w)
                    occ_rep[key] = max(occ_rep.get(key, 0.0), pct)
                elif not is_excluded:
                    bg_pct.append(pct); bg_mod += nmod; bg_cov += cov

            if not in_pct or not bg_pct:
                results.append([sample, motif, code, e["n_occ"],
                                 len(in_pct), "NA", "NA", "NA", "NA",
                                 "NA", "NA", "NA", offset_str, scoring])
                continue

            in_mean = statistics.mean(in_pct)
            bg_mean = statistics.mean(bg_pct)
            log2fc = math.log2((in_mean + 0.1) / (bg_mean + 0.1))

            a, b = in_mod, in_cov - in_mod
            c, d = bg_mod, bg_cov - bg_mod
            n = a + b + c + d
            denom = (a + b) * (c + d) * (a + c) * (b + d)
            if denom > 0:
                chi2 = n * (a * d - b * c) ** 2 / denom
                pval = math.erfc(math.sqrt(chi2 / 2))
            else:
                chi2 = 0.0
                pval = 1.0

            reps = list(occ_rep.values())
            rep_mean = statistics.mean(reps)
            rep_cv = (statistics.pstdev(reps) / rep_mean) if rep_mean > 0 and len(reps) > 1 else 0.0
            n_low = sum(1 for r in reps if r < low_cut)

            results.append([sample, motif, code, e["n_occ"],
                             len(in_pct), "%.2f" % in_mean, "%.2f" % bg_mean,
                             "%.3f" % log2fc, "%.3e" % pval,
                             "%.3f" % rep_cv, len(reps), n_low, offset_str, scoring])

        with open(out_path, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t", lineterminator="\n")
            w.writerow(["sample", "motif", "mod_code", "n_genomic_occurrences",
                        "n_sites_tested", "in_motif_mean_pct", "background_mean_pct",
                        "log2_enrichment", "chi2_pvalue",
                        "occurrence_cv", "n_occurrences_scored", "n_occurrences_low",
                        "offset", "scoring"])
            w.writerows(results)
        PY
        # No de-indenting here: WDL strips the command block's common
        # leading whitespace before the shell ever sees it, so the heredoc
        # lands with the python's relative indentation already correct.

        python3 landscape.py \
            ~{bedmethyl} ~{motif_list} ~{reference_fasta} ~{sample_name} \
            ~{min_coverage} ~{heterogeneous_low_cutoff} \
            ~{sample_name}_motif_landscape.tsv

        awk 'NR>1' ~{sample_name}_motif_landscape.tsv | wc -l > N_TESTED
    >>>

    output {
        File landscape_tsv = "~{sample_name}_motif_landscape.tsv"
        Int  n_motifs_tested = read_int("N_TESTED")
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
