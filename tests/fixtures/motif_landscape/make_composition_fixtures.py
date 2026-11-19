#!/usr/bin/env python3
"""Generate the composition-independence fixtures for motif_landscape. Pure python.

Two 20 kb genomes of different GC (38% and 66%) carry the SAME methylation truth:
every occurrence of the bipartite Type I motif GAAGNNNNNNTCC (and of GATC) is
methylated at ~95% on its modified base(s); every other A is ~2% noise. The
spacer N's are drawn from the genome's own composition, so an A-rich genome puts
many more unmethylated A's inside a motif window than a GC-rich one does.

A scorer that counts every A in the window therefore reports a different
in-motif rate for the two genomes; one that scores only the modified base does
not. The in-motif rate for the planted motif should be ~95 in both.

Not committed pre-built: the bedMethyl rows are megabytes. The test builds them
into a scratch directory.

Writes, into OUTDIR:  lowgc.fa, lowgc.bedmethyl.bed, highgc.fa, highgc.bedmethyl.bed,
                      motifs_offset.tsv (motif, mod_code, offset), motifs_window.tsv

Usage: make_composition_fixtures.py OUTDIR
"""
import os
import random
import re
import sys

L = 20_000
N_PLANTED = 60
# (motif as written 5'->3' on the modified strand, offset of the modified A).
# GGANNNNNNCTTC is the minus-strand half of GAAGNNNNNNTCC: the same sites,
# methylated on the other strand.
MOTIFS = [("GATC", 1), ("GAAGNNNNNNTCC", 2), ("GGANNNNNNCTTC", 2)]
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def occurrences(seq, motif):
    pat = lambda m: "(?=(" + m.replace("N", "[ACGT]") + "))"
    fwd = [x.start() for x in re.finditer(pat(motif), seq)]
    rev = [x.start() for x in re.finditer(pat(rc(motif)), seq)] if rc(motif) != motif else list(fwd)
    return fwd, rev


def build(gc, seed, outdir, tag):
    rng = random.Random(seed)
    w = [(1 - gc) / 2, (1 - gc) / 2, gc / 2, gc / 2]  # A T C G
    seq = rng.choices("ATCG", weights=w, k=L)
    spacer = MOTIFS[1][0]
    for _ in range(N_PLANTED):
        p = rng.randrange(0, L - len(spacer))
        site = "".join(rng.choices("ATCG", weights=w)[0] if c == "N" else c for c in spacer)
        seq[p:p + len(spacer)] = rc(site) if rng.random() < 0.5 else site
    seq = "".join(seq)

    methylated = set()
    for motif, k in MOTIFS:
        fwd, rev = occurrences(seq, motif)
        for s in fwd:
            methylated.add(("+", s + k))
        for s in rev:
            methylated.add(("-", s + len(motif) - 1 - k))

    with open(os.path.join(outdir, f"{tag}.fa"), "w") as fh:
        fh.write(">chr1\n" + seq + "\n")
    with open(os.path.join(outdir, f"{tag}.bedmethyl.bed"), "w") as fh:
        for p, base in enumerate(seq):
            if base not in "AT":
                continue
            strand = "+" if base == "A" else "-"
            if (strand, p) in methylated:
                pct = min(100.0, max(60.0, rng.gauss(95, 3)))
            else:
                pct = max(0.0, rng.gauss(2, 1.5))
            cov = 30
            row = ["chr1", p, p + 1, "a", cov, strand, p, p + 1, "0,0,0", cov,
                   f"{pct:.2f}", round(pct * cov / 100)] + [0] * 6
            fh.write("\t".join(map(str, row)) + "\n")


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    build(0.38, 21, outdir, "lowgc")
    build(0.66, 22, outdir, "highgc")
    with open(os.path.join(outdir, "motifs_offset.tsv"), "w") as fh:
        fh.write("motif\tmod_code\toffset\n")
        for m, k in MOTIFS:
            fh.write(f"{m}\ta\t{k}\n")
    with open(os.path.join(outdir, "motifs_window.tsv"), "w") as fh:
        fh.write("motif\tmod_code\n")
        for m, _ in MOTIFS:
            fh.write(f"{m}\ta\n")


if __name__ == "__main__":
    main(sys.argv[1])
