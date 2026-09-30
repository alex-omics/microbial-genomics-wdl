#!/usr/bin/env python3
"""Regenerate the classify_kraken2_batched fixtures.

A hand-made taxonomy of two genera (A with species A1 and A2, B with species B1)
plus a host species, a random genome for each, and four unaligned BAMs whose
reads are exact substrings of those genomes or random sequence absent from the
database, so every read's classification is known in advance. The expected
counts are written to expected.json.

Needs docker (kraken2 and samtools run in the image the workflow uses).
Usage: make_fixtures.py
"""
import json
import os
import random
import subprocess
import tempfile

IMAGE = "quay.io/broadinstitute/viral-classify:2.5.21.0"
HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(7)

# taxid: (parent, rank, name)
TAXA = {
    1: (1, "no rank", "root"),
    2: (1, "superkingdom", "Bacteria"),
    10: (2, "genus", "GenusA"),
    11: (10, "species", "SpeciesA1"),
    12: (10, "species", "SpeciesA2"),
    20: (2, "genus", "GenusB"),
    21: (20, "species", "SpeciesB1"),
    9606: (1, "species", "Host"),
}
GENOMES = {t: "".join(rng.choice("ACGT") for _ in range(6000)) for t in (11, 12, 21, 9606)}
COMP = str.maketrans("ACGT", "TGCA")
READ = 100

# sample -> (paired, planted read or pair counts per source; "U" is random sequence)
SAMPLES = {
    "sampleA": (True, {11: 40, 12: 20, 21: 30, 9606: 10, "U": 25}),
    "sampleB": (True, {21: 50, 9606: 15, "U": 5}),
    "sampleC": (False, {12: 35, "U": 5}),
    "sampleD": (True, {}),
}


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_from(seq):
    i = rng.randrange(0, len(seq) - 400)
    frag = seq[i:i + 300]
    return frag[:READ], revcomp(frag)[:READ]


def sam_records(sample, paired, plan):
    n = 0
    for source, count in plan.items():
        for _ in range(count):
            n += 1
            seq = "".join(rng.choice("ACGT") for _ in range(1000)) if source == "U" else GENOMES[source]
            r1, r2 = read_from(seq)
            name = f"{sample}_{n}"
            q = "I" * READ
            if paired:
                yield f"{name}\t77\t*\t0\t0\t*\t*\t0\t0\t{r1}\t{q}\tRG:Z:{sample}"
                yield f"{name}\t141\t*\t0\t0\t*\t*\t0\t0\t{r2}\t{q}\tRG:Z:{sample}"
            else:
                yield f"{name}\t4\t*\t0\t0\t*\t*\t0\t0\t{r1}\t{q}\tRG:Z:{sample}"


def docker(work, *cmd):
    subprocess.run(["docker", "run", "--rm", "--platform", "linux/amd64", "-v", f"{work}:/w", "-w", "/w",
                    "--entrypoint", cmd[0], IMAGE, *cmd[1:]], check=True)


with tempfile.TemporaryDirectory() as work:
    os.makedirs(f"{work}/lib/taxonomy")
    with open(f"{work}/lib/taxonomy/nodes.dmp", "w") as nodes, open(f"{work}/lib/taxonomy/names.dmp", "w") as names:
        for tid, (parent, rank, name) in TAXA.items():
            nodes.write(f"{tid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n")
            names.write(f"{tid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")
    with open(f"{work}/genomes.fa", "w") as fa:
        for tid, seq in GENOMES.items():
            fa.write(f">seq{tid}|kraken:taxid|{tid}\n{seq}\n")

    docker(work, "bash", "-c",
           "kraken2-build --add-to-library genomes.fa --db lib --no-masking "
           "&& kraken2-build --build --db lib --kmer-len 31 --minimizer-len 15 --minimizer-spaces 0 "
           "&& tar -czf tiny_db.tar.gz -C lib hash.k2d opts.k2d taxo.k2d")
    os.replace(f"{work}/tiny_db.tar.gz", f"{HERE}/tiny_db.tar.gz")

    for sample, (paired, plan) in SAMPLES.items():
        with open(f"{work}/{sample}.sam", "w") as sam:
            sam.write("@HD\tVN:1.6\tSO:unsorted\n")
            sam.write(f"@RG\tID:{sample}\tSM:{sample}\n")
            sam.write("\n".join(sam_records(sample, paired, plan)) + "\n" if plan else "")
        docker(work, "samtools", "view", "-b", "-o", f"{sample}.bam", f"{sample}.sam")
        os.replace(f"{work}/{sample}.bam", f"{HERE}/{sample}.bam")

expected = {s: {str(k): v for k, v in plan.items()} for s, (_, plan) in SAMPLES.items()}
json.dump(expected, open(f"{HERE}/expected.json", "w"), indent=1)
