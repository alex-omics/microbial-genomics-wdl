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
        description: "Merge de novo (find-motifs) and homology-predicted (REBASE) motifs into one deduplicated (motif, mod_code) list for motif_landscape to test. A motif discovered both ways appears once."
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

        motifs = set()

        if find_motifs_path:
            with open(find_motifs_path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    motif = (row.get("motif") or "").strip().upper()
                    code  = TO_CODE.get((row.get("mod_code") or "").strip())
                    if motif and code:
                        motifs.add((motif, code))

        if rebase_path:
            with open(rebase_path, newline="") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    motif = (row.get("predicted_motif") or "").strip().upper()
                    code  = TO_CODE.get((row.get("modification_type") or "").strip())
                    if motif and motif != "NA" and code:
                        motifs.add((motif, code))

        with open(out_path, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t", lineterminator="\n")
            w.writerow(["motif", "mod_code"])
            for motif, code in sorted(motifs):
                w.writerow([motif, code])

        sys.stderr.write("%d distinct (motif, mod_code) pairs\n" % len(motifs))
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
        bedmethyl:                "UNFILTERED bedMethyl from modkit_pileup. Deliberately not the coverage/percent-filtered table -- the enrichment test needs the true genome-wide background rate, which a pre-filtered table would already have thrown away."
        motif_list:               "Output of build_motif_list: one (motif, mod_code) pair per row"
        sample_name:              "Some identifier for naming outputs"
        reference_fasta:          "This isolate's own assembly, in the same coordinates as bedmethyl"
        min_coverage:             "Minimum Nvalid_cov for a site to be included in either the in-motif or background set (default = 10)"
        heterogeneous_low_cutoff: "A motif occurrence is counted as 'low' in the heterogeneity summary if its representative percent-modified falls below this. Meaningful only relative to the panel's typical housekeeping level, which is usually near 100 (default = 50.0)"
        cpu:                      "Number of CPUs delegated to task (default = 4)"
        mem_gb:                   "Amount of memory in GB delegated to task (default = 16)"
        disk_gb:                  "Amount of disk space in GB delegated to task (default = 50)"
        docker:                   "Container image"
    }

    meta {
        description: "Test each candidate motif against this isolate's own methylation data. Tier 1: is methylation at the motif's genomic occurrences elevated over the genome-wide background for that base (log2 enrichment, chi-square). Tier 2: are those occurrences uniformly methylated (the RM-housekeeping null) or split into high/low populations, which is the signature of phase variation or a regulator competing for specific copies. Occurrences are located on both strands and IUPAC ambiguity codes are supported; a motif is tested once per position even where it is palindromic, and its reverse complement is not tested separately if identical to the motif itself. This is intentionally coarse at the base level: every base of the matching canonical type (A for 6mA, C for 4mC/5mC) within an occurrence window counts as part of that occurrence, rather than pinning one exact intra-motif position, because REBASE's own position convention could not be verified against source documentation."
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
        CANONICAL_BASE = {'a': 'A', 'm': 'C', '21839': 'C'}

        def revcomp(m):
            return ''.join(COMP[b] for b in reversed(m))

        def to_regex(motif):
            return ''.join(IUPAC[b] for b in motif)

        def find_occurrences(seq, motif):
            fwd = re.compile('(?=(' + to_regex(motif) + '))')
            rc = revcomp(motif)
            hits = {m.start() for m in fwd.finditer(seq)}
            if rc != motif:
                rev = re.compile('(?=(' + to_regex(rc) + '))')
                hits |= {m.start() for m in rev.finditer(seq)}
            return sorted(hits)

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

        motifs = []
        with open(motif_path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                motifs.append((row["motif"].strip().upper(), row["mod_code"].strip()))

        # windows[(motif, mod_code)] = {contig: [start, start, ...]}, each
        # contig's starts sorted ascending (find_occurrences already returns
        # them sorted). Every window for one motif shares the same length,
        # so storing just the start is enough to reconstruct membership.
        windows = {}
        for motif, code in motifs:
            per_contig = {}
            for contig, seq in seqs.items():
                occ = find_occurrences(seq, motif)
                if occ:
                    per_contig[contig] = occ
            windows[(motif, code)] = per_contig

        def in_any_window(contig, pos, per_contig, motif_len):
            starts = per_contig.get(contig)
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

        # bedMethyl: 0 chrom,1 start,3 mod_code,9 Nvalid_cov,10 pct,11 Nmod
        by_code = defaultdict(list)
        with open(bedmethyl_path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 12:
                    continue
                cov = int(f[9])
                if cov < mincov:
                    continue
                by_code[f[3]].append((f[0], int(f[1]), float(f[10]), int(f[11]), cov))

        results = []
        for motif, code in motifs:
            per_contig = windows[(motif, code)]
            motif_len = len(motif)
            sites = by_code.get(code, [])
            in_pct, bg_pct = [], []
            in_mod = in_cov = bg_mod = bg_cov = 0
            occ_rep = {}  # (contig, window_start) -> max pct seen in that window
            for contig, pos, pct, nmod, cov in sites:
                w = in_any_window(contig, pos, per_contig, motif_len)
                if w is not None:
                    in_pct.append(pct); in_mod += nmod; in_cov += cov
                    key = (contig, w)
                    occ_rep[key] = max(occ_rep.get(key, 0.0), pct)
                else:
                    bg_pct.append(pct); bg_mod += nmod; bg_cov += cov

            n_genomic_occurrences = sum(len(v) for v in per_contig.values())

            if not in_pct or not bg_pct:
                results.append([sample, motif, code, n_genomic_occurrences,
                                 len(in_pct), "NA", "NA", "NA", "NA",
                                 "NA", "NA", "NA"])
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

            results.append([sample, motif, code, n_genomic_occurrences,
                             len(in_pct), "%.2f" % in_mean, "%.2f" % bg_mean,
                             "%.3f" % log2fc, "%.3e" % pval,
                             "%.3f" % rep_cv, len(reps), n_low])

        with open(out_path, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t", lineterminator="\n")
            w.writerow(["sample", "motif", "mod_code", "n_genomic_occurrences",
                        "n_sites_tested", "in_motif_mean_pct", "background_mean_pct",
                        "log2_enrichment", "chi2_pvalue",
                        "occurrence_cv", "n_occurrences_scored", "n_occurrences_low"])
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
