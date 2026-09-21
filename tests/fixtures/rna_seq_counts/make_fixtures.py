#!/usr/bin/env python3
"""Regenerate the rna_seq_counts fixtures. Deterministic (fixed seed).

Two synthetic isolates with ground-truth expression, so the workflow's output
can be asserted exactly rather than eyeballed:

  IsoA  6 genes.  Samples IsoA_a, IsoA_b        (letter replicates)
  IsoB  5 genes.  Samples IsoB_1                (numeric replicate)

IsoA and IsoB share three genes (identical sequence, different locus tags), IsoA
carries two the other lacks, and IsoB carries one IsoA lacks - so the ortholog
matrix has genuine absences. Two IsoA genes are paralogues in one ortholog group
(counts must be summed). Every gene is fully covered by fragments and none is
repeated, so every read maps uniquely and expected counts are exact.

Reads are dUTP-stranded (read 1 antisense to the transcript), matching the
workflow's default strandness of 2.
"""
import gzip, json, os, random

HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(20260918)
COMP = str.maketrans("ACGT", "TGCA")

def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))

# Codons other than the three stops, so a gene is a real ORF. Panaroo (and Bakta)
# discard anything that is not one, so random sequence would vanish from the
# pangenome and the built-ortholog-table test would have nothing to join.
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
          if a + b + c not in ("TAA", "TAG", "TGA")]

def rand_orf(n_codons):
    return "ATG" + "".join(rng.choice(CODONS) for _ in range(n_codons)) + "TAA"

def revcomp(s):
    return s.translate(COMP)[::-1]

# ---- shared gene sequences
shared = {name: rand_orf(rng.randint(200, 260)) for name in ("core1", "core2", "core3")}

# gene list per isolate: (group, sequence, strand)
isolates = {
    "IsoA": [
        ("core1",  shared["core1"], "+"),
        ("core2",  shared["core2"], "-"),
        ("core3",  shared["core3"], "+"),
        ("onlyA",  rand_orf(230), "-"),
        ("para1",  rand_orf(215), "+"),    # paralogue pair: one ortholog group
        ("para2",  rand_orf(215), "+"),
    ],
    "IsoB": [
        ("core1",  shared["core1"], "-"),
        ("core2",  shared["core2"], "-"),
        ("core3",  shared["core3"], "+"),
        ("onlyB",  rand_orf(238), "+"),
        ("silent", rand_orf(198), "+"),    # present, never expressed -> real 0
    ],
}

# ground-truth fragments per gene, per sample
truth = {
    "IsoA_a": {"core1": 200, "core2": 100, "core3": 50, "onlyA": 40, "para1": 30, "para2": 10},
    "IsoA_b": {"core1": 150, "core2": 120, "core3": 60, "onlyA": 35, "para1": 20, "para2": 25},
    "IsoB_1": {"core1": 180, "core2": 90,  "core3": 70, "onlyB": 60, "silent": 0},
}

def tag(iso, i):
    return "%s_%05d" % (iso, i + 1)

fasta, gff_rows, genes = {}, {}, {}
for iso, glist in isolates.items():
    seq, feats, pos = [], [], 0
    spacer = lambda: rand_seq(rng.randint(300, 500))
    seq.append(spacer()); pos = len(seq[0])
    for i, (grp, gseq, strand) in enumerate(glist):
        start = pos + 1
        placed = gseq if strand == "+" else revcomp(gseq)
        seq.append(placed); pos += len(placed)
        feats.append((tag(iso, i), grp, start, pos, strand))
        sp = spacer(); seq.append(sp); pos += len(sp)
    full = "".join(seq)
    fasta[iso] = full
    genes[iso] = feats

    with open(os.path.join(HERE, iso + ".fasta"), "w") as fh:
        fh.write(">%s_contig1 test contig\n" % iso)
        for i in range(0, len(full), 80):
            fh.write(full[i:i+80] + "\n")
    # A Bakta-style GFF3: features, then a ##FASTA block that featureCounts must not choke on.
    with open(os.path.join(HERE, iso + ".gff3"), "w") as fh:
        fh.write("##gff-version 3\n##sequence-region %s_contig1 1 %d\n" % (iso, len(full)))
        fh.write("%s_contig1\tBakta\tregion\t1\t%d\t.\t+\t.\tID=%s_contig1\n" % (iso, len(full), iso))
        for t, grp, s, e, strand in feats:
            fh.write("%s_contig1\tPyrodigal\tCDS\t%d\t%d\t.\t%s\t0\tID=%s;locus_tag=%s;product=hypothetical protein %s\n"
                     % (iso, s, e, strand, t, t, grp))
        fh.write("##FASTA\n>%s_contig1\n" % iso)
        for i in range(0, len(full), 80):
            fh.write(full[i:i+80] + "\n")

# ---- reads
READ, FRAG = 75, 200
def write_reads(sample, iso, gene_counts):
    lookup = {grp: (s, e, strand) for _, grp, s, e, strand in genes[iso]}
    g = fasta[iso]
    r1, r2 = [], []
    n = 0
    for grp, cnt in gene_counts.items():
        s, e, strand = lookup[grp]
        for _ in range(cnt):
            a = rng.randint(s - 1, e - FRAG)          # fragment lies wholly inside the gene
            frag = g[a:a + FRAG]                       # forward-strand genomic sequence
            if strand == "-":
                frag = revcomp(frag)                   # transcript sense sequence
            # dUTP: read 1 is antisense to the transcript, read 2 is sense.
            a1 = revcomp(frag)[:READ]
            a2 = frag[:READ]
            n += 1
            name = "%s_read%d" % (sample, n)
            r1.append((name, a1)); r2.append((name, a2))
    order = list(range(len(r1))); rng.shuffle(order)
    for suffix, recs in (("R1", r1), ("R2", r2)):
        with gzip.open(os.path.join(HERE, "%s_%s.fastq.gz" % (sample, suffix)), "wt") as fh:
            for j in order:
                nm, sq = recs[j]
                fh.write("@%s/%s\n%s\n+\n%s\n" % (nm, suffix[1], sq, "I" * len(sq)))

for sample, gc in truth.items():
    write_reads(sample, sample.split("_")[0], gc)

# ---- ortholog table (Panaroo gene_presence_absence.csv layout; annotation has a comma)
groups = ["core1", "core2", "core3", "onlyA", "para", "onlyB", "silent"]
def cell(iso, grp):
    if grp == "para":
        return ";".join(tag(iso, i) for i, (g, *_r) in enumerate(isolates[iso]) if g in ("para1", "para2")) if iso == "IsoA" else ""
    for i, (g, *_r) in enumerate(isolates[iso]):
        if g == grp:
            return tag(iso, i)
    return ""
import csv
with open(os.path.join(HERE, "gene_presence_absence.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["Gene", "Non-unique Gene name", "Annotation", "IsoA", "IsoB"])
    for grp in groups:
        w.writerow(["group_" + grp, "", "hypothetical protein, %s family" % grp, cell("IsoA", grp), cell("IsoB", grp)])

# ---- expected counts, keyed by locus tag, for the test to assert against
expected = {}
for sample, gc in truth.items():
    iso = sample.split("_")[0]
    expected[sample] = {tag(iso, i): gc.get(grp, 0) for i, (grp, *_r) in enumerate(isolates[iso])}
with open(os.path.join(HERE, "expected_counts.json"), "w") as fh:
    json.dump(expected, fh, indent=1, sort_keys=True)
print("fixtures written to", HERE)
