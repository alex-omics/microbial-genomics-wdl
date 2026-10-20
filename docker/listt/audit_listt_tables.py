#!/usr/bin/env python3
"""Audit LISTT's reference tables as shipped in the image.

    docker run --rm aarvani1/listt:1.2.0 python /opt/audit/audit_listt_tables.py

Checks the three files LISTT classifies from (variants/ospA_protein.fas,
serotypes.csv, thresholds.csv) against each other, then re-derives each IST's
identity threshold by the method in Lee et al. 2024 (midpoint between the lowest
within-IST identity and the highest out-of-IST identity, taking each variant's
best match) on the current membership, and prints it beside the shipped value.

Identity is LISTT's own: matching columns over the length of a MAFFT alignment of
every reference variant. A threshold that was derived on today's membership
should reproduce; one that does not was set on an older membership or by hand.
"""
import collections
import csv
import io
import subprocess
import sys
import tempfile

V = "/opt/LISTT/variants"


def read_fasta(text):
    out, name, seq = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name:
        out[name] = "".join(seq)
    return out


protein = {int(k.split("_")[-1]): v.rstrip("*") for k, v in read_fasta(open(f"{V}/ospA_protein.fas").read()).items()}

rows = []
with open(f"{V}/serotypes.csv", newline="", encoding="utf-8-sig") as fh:
    for r in csv.reader(fh):
        if r and r[0].isdigit():
            rows.append((int(r[0]), r[1], r[2]))
label = {}
for v, ist, sp in rows:
    label.setdefault(v, (ist, sp))
thr = {}
with open(f"{V}/thresholds.csv", newline="", encoding="utf-8-sig") as fh:
    for r in csv.reader(fh):
        if r and r[0].isdigit():
            thr[int(r[0])] = float(r[1])

print("== table consistency")
print(f"protein variants: {len(protein)}; serotypes.csv rows: {len(rows)}; thresholds rows: {len(thr)}")
print("in the protein FASTA but not serotypes.csv (LISTT crashes if it best-matches one):",
      sorted(set(protein) - set(label)) or "none")
print("in serotypes.csv but not the protein FASTA:", sorted(set(label) - set(protein)) or "none")
dups = [v for v, n in collections.Counter(v for v, _, _ in rows).items() if n > 1]
print("variants repeated in serotypes.csv:", {v: sum(1 for x, _, _ in rows if x == v) for v in dups} or "none")
numeric = sorted({int(i) for i, _ in label.values() if i.isdigit()})
print("ISTs with a threshold row:", sorted(thr) == numeric and "all" or f"missing {sorted(set(numeric) - set(thr))}")
print("variants with no IST ('unk'):", sorted(v for v, (i, _) in label.items() if not i.isdigit()))
print("an IST whose members carry more than one genospecies name:")
by_ist = collections.defaultdict(set)
for v, (i, sp) in label.items():
    by_ist[i].add(sp)
for i in sorted((k for k in by_ist if k.isdigit()), key=int):
    if len(by_ist[i]) > 1:
        print(f"  IST{i}: {sorted(by_ist[i])}")

aln = subprocess.run(["mafft", "--quiet", "/dev/stdin"], input="".join(f">v{v}\n{s}\n" for v, s in protein.items()),
                     capture_output=True, text=True).stdout
seqs = {int(k[1:]): v for k, v in read_fasta(aln).items()}
L = len(next(iter(seqs.values())))


def ident(a, b):
    return 100.0 * sum(x == y for x, y in zip(seqs[a], seqs[b])) / L


ids = sorted(seqs)
M = {(a, b): ident(a, b) for a in ids for b in ids if a != b}
members = collections.defaultdict(list)
for v in ids:
    if v in label:
        members[label[v][0]].append(v)

print(f"\n== thresholds re-derived on current membership (alignment length {L})")
print(f"{'IST':>4} {'variants':>8} {'shipped':>8} {'rederived':>9} {'within_min':>10} {'out_max':>8}  separable")
for ist in sorted((k for k in members if k.isdigit()), key=int):
    vs = members[ist]
    inside = [max([M[(v, u)] for u in vs if u != v] or [100.0]) for v in vs]
    outside = [max(M[(v, u)] for u in ids if u not in vs and u in label) for v in vs]
    lo, hi = min(inside), max(outside)
    print(f"{ist:>4} {len(vs):>8} {thr.get(int(ist), float('nan')):>8.2f} {(lo + hi) / 2:>9.2f} {lo:>10.2f} {hi:>8.2f}  "
          f"{'yes' if lo > hi else 'NO: within-IST identity falls below an out-of-IST match'}")
