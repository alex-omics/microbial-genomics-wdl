#!/usr/bin/env python3
"""Regenerate the filter_reads fixtures.

Eight synthetic ONT-style reads whose length and quality are planted, so what any
filter should keep is known in advance (see READS). Each is written three ways:

  reads.bam         unaligned modBAM: MM/ML/MN tags and a dorado-style qs tag
  reads_noqs.bam    the same records without the qs tag
  reads.fastq.gz    the same reads as FASTQ

plus aligned.bam (one mapped record; the BAM filter must refuse it) and
expected.json (which reads each parameter set keeps).

Needs docker (samtools runs in the image the workflow uses).
Usage: make_fixtures.py
"""
import gzip
import json
import math
import os
import random
import subprocess
import tempfile

IMAGE = "aarvani1/read-filter:1.0.0"
HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(11)

# name: (length, quality). A number is the same Phred score at every base; a pair
# (a, b) is the first half at a and the second half at b.
READS = {
    "short_hq":  (300,  20),         # fails min_length 500
    "ok_a":      (1000, 15),         # kept by everything below
    "ok_b":      (2000, 12),         # kept by everything below
    # No read sits exactly on a quality threshold: a mean Q is a float, and ties
    # are decided by rounding noise (chopper -q 15 kept a 1000 bp Q15 read and
    # dropped a 6000 bp one).
    "lowq":      (1500, 6),          # fails min_mean_q 10
    "long_hq":   (6000, 15),         # fails max_length 5000
    "boundary":  (500,  11),         # exactly on min_length 500 (both tools are inclusive)
    "mixed":     (2000, (30, 2)),    # arithmetic mean Q is 16, but its mean error rate is Q5
    "short_lq":  (200,  5),          # fails both
}


def mean_q(quals):
    """Mean Phred score taken over error probabilities, which is what read-level
    quality means for ONT; averaging the scores themselves overstates a read with
    a good stretch and a bad one."""
    p = sum(10 ** (-q / 10) for q in quals) / len(quals)
    return -10 * math.log10(p)


def phred(spec, n):
    if isinstance(spec, tuple):
        return [spec[0]] * (n // 2) + [spec[1]] * (n - n // 2)
    return [spec] * n


# Drawn once, so the BAMs and the FASTQ hold identical sequences
RECORDS = [("".join(rng.choice("ACGT") for _ in range(n)), name, phred(spec, n))
           for name, (n, spec) in READS.items()]
RECORDS = [(name, seq, q) for seq, name, q in RECORDS]


def sam_text(with_qs):
    lines = ["@HD\tVN:1.6\tSO:unsorted",
             "@RG\tID:fix\tPL:ONT\tSM:fixture",
             "@PG\tID:basecaller\tPN:dorado"]
    for name, seq, q in RECORDS:
        qual = "".join(chr(33 + x) for x in q)
        tags = ["MM:Z:C+m?,0,3,1;", f"ML:B:C,200,180,150", f"MN:i:{len(seq)}", "RG:Z:fix"]
        if with_qs:
            tags.append(f"qs:f:{mean_q(q):.2f}")
        lines.append("\t".join([name, "4", "*", "0", "0", "*", "*", "0", "0", seq, qual] + tags))
    return "\n".join(lines) + "\n"


def docker(work, *cmd):
    subprocess.run(["docker", "run", "--rm", "--platform", "linux/amd64", "-v", f"{work}:/w", "-w", "/w",
                    IMAGE, *cmd], check=True)


def main():
    with tempfile.TemporaryDirectory() as work:
        for fname, with_qs in (("reads", True), ("reads_noqs", False)):
            with open(os.path.join(work, fname + ".sam"), "w") as fh:
                fh.write(sam_text(with_qs))
            docker(work, "samtools", "view", "-b", "-o", fname + ".bam", fname + ".sam")
            os.replace(os.path.join(work, fname + ".bam"), os.path.join(HERE, fname + ".bam"))

        # One mapped record: flag 0, on a contig named ref. Enough for the "this BAM is aligned" guard.
        with open(os.path.join(work, "aligned.sam"), "w") as fh:
            seq = "".join(rng.choice("ACGT") for _ in range(1000))
            fh.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:5000\n")
            fh.write("\t".join(["mapped", "0", "ref", "100", "60", "1000M", "*", "0", "0", seq, "I" * 1000]) + "\n")
        docker(work, "samtools", "view", "-b", "-o", "aligned.bam", "aligned.sam")
        os.replace(os.path.join(work, "aligned.bam"), os.path.join(HERE, "aligned.bam"))

    with gzip.open(os.path.join(HERE, "reads.fastq.gz"), "wt") as fh:
        for name, seq, q in RECORDS:
            fh.write(f"@{name} runid=fixture\n{seq}\n+\n{''.join(chr(33 + x) for x in q)}\n")

    summary = {name: {"length": n, "mean_q": round(mean_q(phred(spec, n)), 2)} for name, (n, spec) in READS.items()}
    with open(os.path.join(HERE, "expected.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
