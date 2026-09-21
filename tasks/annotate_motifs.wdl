version 1.0

task annotate_motif_membership {

    input {
        File     annotated_tsv
        File     motif_list
        String   sample_name
        File     reference_fasta
        Int      cpu     = 2
        Int      mem_gb  = 8
        Int      disk_gb = 20
        String   docker  = "python:3.11-slim@sha256:db3ff2e1800a8581e2c48a27c3995339d47bdf046da21c7627accd3d51053a93"
    }

    parameter_meta {
        annotated_tsv:    "Per-site methylation table from annotate_methylation (locus_tag/gene/product/region already joined)"
        motif_list:       "Output of build_motif_list: one (motif, mod_code) pair per row, this isolate's de novo + REBASE candidates"
        sample_name:      "Some identifier for naming outputs"
        reference_fasta:  "This isolate's own assembly, in the same coordinates as annotated_tsv"
        cpu:              "Number of CPUs delegated to task (default = 2)"
        mem_gb:           "Amount of memory in GB delegated to task (default = 8)"
        disk_gb:          "Amount of disk space in GB delegated to task (default = 20)"
        docker:           "Container image"
    }

    meta {
        description: "Adds a motif column to annotate_methylation's per-site table: which of this isolate's candidate motifs (de novo + REBASE), if any, each methylated site falls inside. This is the single per-isolate deliverable combining gene/product/region annotation, promoter windows, and motif context alongside the actual methylation call, in one row per site."
    }

    command <<<
        set -euo pipefail

        cat > annotate_motifs.py <<'PY'
        import bisect, csv, re, sys
        from collections import defaultdict

        annotated_path, motif_path, ref_path, out_path = sys.argv[1:5]

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
            fwd = re.compile('(?=(' + to_regex(motif) + '))')
            rc = revcomp(motif)
            hits = {m.start() for m in fwd.finditer(seq)}
            if rc != motif:
                rev = re.compile('(?=(' + to_regex(rc) + '))')
                hits |= {m.start() for m in rev.finditer(seq)}
            return sorted(hits)

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

        # Grouped by mod_code, since a site can only ever be inside a motif
        # tested for the SAME modification it carries -- an 'a' (6mA) site
        # is never a hit for a motif keyed to code 'm' (5mC).
        by_code = defaultdict(list)  # code -> [(motif, {contig: sorted_starts})]
        with open(motif_path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                motif = row["motif"].strip().upper()
                code  = row["mod_code"].strip()
                per_contig = {}
                for contig, seq in seqs.items():
                    occ = find_occurrences(seq, motif)
                    if occ:
                        per_contig[contig] = occ
                by_code[code].append((motif, per_contig))

        def motifs_at(code, contig, pos):
            hits = []
            for motif, per_contig in by_code.get(code, []):
                starts = per_contig.get(contig)
                if not starts:
                    continue
                # Same bisect trick as motif_landscape: a window covers pos
                # iff its start falls in [pos-len(motif)+1, pos].
                lo = pos - len(motif) + 1
                i = bisect.bisect_left(starts, lo)
                if i < len(starts) and starts[i] <= pos:
                    hits.append(motif)
            return hits

        with open(annotated_path, newline="") as fh, open(out_path, "w", newline="") as out:
            rd = csv.reader(fh, delimiter="\t")
            wr = csv.writer(out, delimiter="\t", lineterminator="\n")
            header = next(rd)
            wr.writerow(header + ["motif"])
            # annotate_methylation's own column order: sample,chrom,start,end,
            # site_strand,mod_code,n_valid_cov,percent_modified,locus_tag,gene,
            # feature_strand,region_length,product,region
            for row in rd:
                if len(row) < 14:
                    wr.writerow(row + ["NA"])
                    continue
                hits = motifs_at(row[5], row[1], int(row[2]))
                wr.writerow(row + [",".join(hits) if hits else "NA"])
        PY
        # No de-indenting here: WDL strips the command block's common
        # leading whitespace before the shell ever sees it, so the heredoc
        # lands with the python's relative indentation already correct.

        python3 annotate_motifs.py \
            ~{annotated_tsv} ~{motif_list} ~{reference_fasta} \
            ~{sample_name}_methylation_annotated_with_motifs.tsv
    >>>

    output {
        File combined_tsv = "~{sample_name}_methylation_annotated_with_motifs.tsv"
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
