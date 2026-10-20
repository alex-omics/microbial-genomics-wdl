#!/usr/bin/env python3
"""Regenerate the call_ospA_type fixtures.

Everything is derived from a LISTT checkout (Pfizer, Apache-2.0), so the
fixtures can be rebuilt and their provenance checked:

    python3 make_fixtures.py --listt-dir /path/to/LISTT --out .

Two fixtures are real: a ~4 kb window around ospA cut from each of LISTT's two
PubMLST test genomes (1_B31, IST1; 155_PBi, IST4 - the second on the reverse
strand). The others are synthetic ospA genes built by editing the codons of a
LISTT reference allele until its protein matches a chosen target, so BLAST
still sees one ~97% identical hit, set in random flanking sequence:

  below_threshold  variant 1 with ~11% of residues changed: closest to IST1 but
                   well under its 95.8 identity threshold
  provisional      the best-buildable variant of an IST above 17 (not in the paper)
  variant68        LISTT's variant 68, which its serotypes.csv does not list, so
                   LISTT itself crashes on it
  no_ospa          a random contig with no ospA at all

Reads are 2x150 bp pairs at 40x with a 0.2% error rate. Seeded, so output is
reproducible; gzip headers carry no timestamp.
"""
import argparse
import csv
import gzip
import json
import os
import random

CODONS = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, c in enumerate(f"{a}{b}{d}" for a in _bases for b in _bases for d in _bases):
    CODONS[c] = _aas[i]
BY_AA = {}
for c, a in CODONS.items():
    BY_AA.setdefault(a, []).append(c)
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def translate(nt):
    return "".join(CODONS[nt[i:i + 3]] for i in range(0, len(nt) - 2, 3))


def read_fasta(path):
    out, name, seq = {}, None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out[name] = "".join(seq).upper()
            name, seq = line[1:].split()[0], []
        elif line:
            seq.append(line)
    if name:
        out[name] = "".join(seq).upper()
    return out


def write_fasta(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


def random_dna(rng, n):
    return "".join(rng.choices("ACGT", weights=[36, 14, 14, 36], k=n))  # AT-rich, like Borrelia


def reads(rng, seq, path_stub, depth=40, rl=150, err=0.002):
    n = int(depth * len(seq) / (2 * rl))
    files = [gzip.GzipFile(f"{path_stub}_R{i}.fastq.gz", "wb", mtime=0) for i in (1, 2)]
    name = os.path.basename(path_stub)
    for k in range(n):
        flen = rng.randint(250, 400)
        if len(seq) < flen:
            continue
        p = rng.randint(0, len(seq) - flen)
        frag = seq[p:p + flen]
        if rng.random() < 0.5:
            frag = rc(frag)
        mates = (frag[:rl], rc(frag)[:rl])
        for fh, mate, i in zip(files, mates, (1, 2)):
            m = "".join(rng.choice([b for b in "ACGT" if b != c]) if rng.random() < err else c for c in mate)
            fh.write(f"@{name}.{k}/{i}\n{m}\n+\n{'I' * len(m)}\n".encode())
    for fh in files:
        fh.close()


def build_nt(target_aa, allele_nt):
    """Edit allele_nt codon by codon to encode target_aa, changing as few bases as possible."""
    assert len(allele_nt) == 3 * (len(target_aa) + 1), "needs a same-length allele"
    out = []
    for i, aa in enumerate(target_aa):
        old = allele_nt[3 * i:3 * i + 3]
        if CODONS[old] == aa:
            out.append(old)
        else:
            out.append(min(BY_AA[aa], key=lambda c: sum(x != y for x, y in zip(c, old))))
    return "".join(out) + allele_nt[-3:]  # keep the allele's stop codon


def best_allele(target_aa, alleles):
    best = None
    for name, nt in alleles.items():
        if len(nt) != 3 * (len(target_aa) + 1):
            continue
        ref = translate(nt)[:-1]
        ident = sum(a == b for a, b in zip(ref, target_aa)) / len(target_aa)
        if best is None or ident > best[0]:
            best = (ident, name)
    return best


def locate_ospa(contigs, alleles):
    """Find ospA in a genome by 28-mer seeds from every reference allele, either strand."""
    hits = {}
    for nt in alleles.values():
        for off in (60, 400, 760):
            seed = nt[off:off + 28]
            for strand, s in (("+", seed), ("-", rc(seed))):
                for cname, cseq in contigs.items():
                    p = cseq.find(s)
                    if p >= 0:
                        start = p - off if strand == "+" else p + len(seed) + off - len(nt)
                        hits.setdefault((cname, strand), []).append(start)
    (cname, strand), starts = max(hits.items(), key=lambda kv: len(kv[1]))
    starts.sort()
    return cname, strand, starts[len(starts) // 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--listt-dir", required=True)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    rng = random.Random(20251020)
    L = a.listt_dir
    alleles = {k.split("_")[-1]: v for k, v in read_fasta(f"{L}/alleles/ref_alleles.fasta").items()}
    protein = {int(k.split("_")[-1]): v.rstrip("*") for k, v in read_fasta(f"{L}/variants/ospA_protein.fas").items()}
    ist = {}
    with open(f"{L}/variants/serotypes.csv", newline="", encoding="utf-8-sig") as fh:
        for r in csv.reader(fh):
            if r and r[0].isdigit():
                ist.setdefault(int(r[0]), r[1])
    expected = {}

    # ---- real ospA windows from LISTT's PubMLST test genomes
    for stem, label in (("1_B31", "ist1_b31"), ("155_PBi", "ist4_pbi")):
        contigs = read_fasta(f"{L}/test/pubmlst/{stem}.fas")
        cname, strand, start = locate_ospa(contigs, alleles)
        seq = contigs[cname][max(0, start - 1700):start + 822 + 1700]
        write_fasta(f"{a.out}/{label}.fasta", label, seq)
        reads(rng, seq, f"{a.out}/{label}")
        expected[label] = {"source": stem, "contig": cname, "strand": strand, "length": len(seq)}

    def flank():
        return random_dna(rng, 1500)

    # ---- synthetic: below the IST1 threshold
    p = list(protein[1])
    for i in rng.sample(range(30, len(p)), 30):  # keep the signal peptide
        p[i] = rng.choice([x for x in "ACDEFGHIKLMNPQRSTVWY" if x != p[i]])
    below = "".join(p)
    ident, allele = best_allele(below, alleles)
    write_fasta(f"{a.out}/below_threshold.fasta", "below_threshold", flank() + build_nt(below, alleles["1"]) + flank())
    expected["below_threshold"] = {"derived_from_variant": 1, "approx_identity_to_variant_1": round(
        sum(x == y for x, y in zip(below, protein[1])) / len(below), 3)}

    # ---- synthetic: closest variant belongs to an IST above 17
    cands = []
    for v, aa in protein.items():
        if ist.get(v, "").isdigit() and int(ist[v]) > 17:
            r = best_allele(aa, alleles)
            if r:
                cands.append((r[0], v, r[1]))
    ident, v, allele = max(cands)
    write_fasta(f"{a.out}/provisional.fasta", "provisional", flank() + build_nt(protein[v], alleles[allele]) + flank())
    expected["provisional"] = {"variant": v, "ist": int(ist[v]), "built_from_allele": allele, "protein_identity_to_allele": round(ident, 3)}

    # ---- synthetic: variant 68, absent from serotypes.csv
    ident, allele = best_allele(protein[68], alleles)
    write_fasta(f"{a.out}/variant68.fasta", "variant68", flank() + build_nt(protein[68], alleles[allele]) + flank())
    expected["variant68"] = {"variant": 68, "built_from_allele": allele, "protein_identity_to_allele": round(ident, 3)}

    # ---- no ospA at all
    seq = random_dna(rng, 4200)
    write_fasta(f"{a.out}/no_ospa.fasta", "no_ospa", seq)
    reads(rng, seq, f"{a.out}/no_ospa")

    json.dump(expected, open(f"{a.out}/expected.json", "w"), indent=2, sort_keys=True)
    print(json.dumps(expected, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
