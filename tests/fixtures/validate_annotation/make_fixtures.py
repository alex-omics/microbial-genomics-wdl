#!/usr/bin/env python3
"""Regenerate the validate_annotation fixtures. Pure python; needs no docker.

Two isolates, each an assembly plus the annotator's files for it (fna, gff3, faa),
and the ways a mismatched set can look:

  isoA.*       chrA 1200 bp + plA 400 bp. Clean set; the fna is wrapped differently
               from the assembly, and the GFF carries ##sequence-region lines, a
               tRNA, a short ORF (sORF, in the faa but not a CDS) and a ##FASTA block
  isoB.*       chrB 1000 bp + plB 500 bp, a different isolate
  renamed.*    isoA with the annotator's contig_1/contig_2 names (same lengths)
  partial.gff3 isoA with only the plasmid renamed
  badregion.gff3  isoA with chrA's ##sequence-region off by one
  genericA.*, genericB.*  two isolates whose annotator named the contigs
               contig_1/contig_2 and numbered the genes GEN_00001... the same way, so
               only the contig LENGTHS tell them apart

expected.json holds the statistics the clean set must report.

Usage: make_fixtures.py
"""
import json
import os
import random
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(5)


def seq(n, gc):
    return "".join(rng.choices("ATCG", weights=[(1 - gc) / 2, (1 - gc) / 2, gc / 2, gc / 2], k=n))


def fasta(path, records, width):
    with open(os.path.join(HERE, path), "w") as fh:
        for name, s in records:
            fh.write(f">{name} some description\n")
            fh.write("\n".join(textwrap.wrap(s, width)) + "\n")


def gff(path, records, features, region_override=None, ids_prefix=""):
    region_override = region_override or {}
    with open(os.path.join(HERE, path), "w") as fh:
        fh.write("##gff-version 3\n")
        for name, s in records:
            fh.write(f"##sequence-region {name} 1 {region_override.get(name, len(s))}\n")
        for contig, ftype, start, end, strand, fid in features:
            fh.write(f"{contig}\tBakta\t{ftype}\t{start}\t{end}\t.\t{strand}\t.\t"
                     f"ID={fid};locus_tag={fid};product=hypothetical protein\n")
        fh.write("##FASTA\n")
        for name, s in records:
            fh.write(f">{name}\n" + "\n".join(textwrap.wrap(s, 80)) + "\n")


def faa(path, ids):
    with open(os.path.join(HERE, path), "w") as fh:
        for i in ids:
            fh.write(f">{i} hypothetical protein\n{seq(40, 0.5).replace('A', 'M').replace('T', 'K')}\n")


def isolate(prefix, contigs, names, gc):
    recs = [(n, seq(L, gc)) for n, L in zip(names, contigs)]
    first, second = names
    feats = [
        (first, "CDS", 10, 300, "+", f"{prefix}_00001"),
        (first, "CDS", 400, 700, "-", f"{prefix}_00002"),
        (first, "CDS", 800, 1000, "+", f"{prefix}_00003"),
        (second, "CDS", 20, 200, "+", f"{prefix}_00004"),
        (second, "sORF", 250, 300, "+", f"{prefix}_00005"),
        (first, "tRNA", 1050, 1120, "+", f"{prefix}_t0001"),
    ]
    return recs, feats, [f"{prefix}_00001", f"{prefix}_00002", f"{prefix}_00003",
                         f"{prefix}_00004", f"{prefix}_00005"]


# ---- isoA: the clean set
recsA, featsA, idsA = isolate("ISOA", [1200, 400], ["chrA", "plA"], 0.38)
fasta("isoA.fa", recsA, 60)
fasta("isoA.fna", recsA, 80)
gff("isoA.gff3", recsA, featsA)
faa("isoA.faa", idsA)

# ---- isoB: a different isolate
recsB, featsB, idsB = isolate("ISOB", [1000, 500], ["chrB", "plB"], 0.66)
fasta("isoB.fa", recsB, 60)
fasta("isoB.fna", recsB, 80)
gff("isoB.gff3", recsB, featsB)
faa("isoB.faa", idsB)

# ---- isoA with the annotator's generic names
ren = {"chrA": "contig_1", "plA": "contig_2"}
recsR = [(ren[n], s) for n, s in recsA]
featsR = [(ren[c], t, a, b, st, i) for c, t, a, b, st, i in featsA]
fasta("renamed.fna", recsR, 80)
gff("renamed.gff3", recsR, featsR)

# ---- isoA with only the plasmid renamed in the GFF
recsP = [("chrA", recsA[0][1]), ("plasmid_x", recsA[1][1])]
featsP = [("plasmid_x" if c == "plA" else c, t, a, b, st, i) for c, t, a, b, st, i in featsA]
gff("partial.gff3", recsP, featsP)

# ---- isoA's GFF declaring chrA one base short
gff("badregion.gff3", recsA, featsA, region_override={"chrA": len(recsA[0][1]) - 1})

# ---- two isolates that look alike to everything but contig length
for tag, lens, gc in (("A", [1200, 400], 0.4), ("B", [1000, 500], 0.6)):
    recs = [(f"contig_{k + 1}", seq(L, gc)) for k, L in enumerate(lens)]
    feats = [
        ("contig_1", "CDS", 10, 300, "+", "GEN_00001"),
        ("contig_1", "CDS", 400, 700, "-", "GEN_00002"),
        ("contig_2", "CDS", 20, 200, "+", "GEN_00003"),
    ]
    fasta(f"generic{tag}.fa", recs, 60)
    fasta(f"generic{tag}.fna", recs, 80)
    gff(f"generic{tag}.gff3", recs, feats)
    faa(f"generic{tag}.faa", ["GEN_00001", "GEN_00002", "GEN_00003"])

gcA = sum(c in "GC" for _, s in recsA for c in s) / sum(len(s) for _, s in recsA) * 100
with open(os.path.join(HERE, "expected.json"), "w") as fh:
    json.dump({"isoA": {"n_contigs": 2, "genome_size": 1600, "n_cds": 4,
                        "gc_percent": f"{gcA:.2f}"}}, fh, indent=2)
    fh.write("\n")
print("isoA GC %.2f%%" % gcA)
