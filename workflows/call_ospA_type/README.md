# call_ospA_type

OspA in silico type (IST) for one *Borrelia burgdorferi* s.l. isolate, from its assembly and,
optionally, its Illumina paired-end reads, using Pfizer's
[LISTT](https://github.com/pfizer-opensource/LISTT) (Lee et al. 2024). Built for Terra: one
sample per execution, no `Array` inputs, so a thousand isolates are a thousand independent
executions.

The assembly is always typed. If reads are supplied they are typed independently, and the two
results are reconciled into one `ospA_call` with a confidence `ospA_flag`.

## Terra usage

Run at the data-table level on rows that hold an assembly and, if you have them, read paths. Map
`assembly_fasta` to `this.assembly_fasta`, `read1`/`read2` to `this.read1`/`this.read2`, select all
rows and launch. Rows without reads simply run assembly-only. Outputs write back as `ospA_*`
columns, prefixed so they do not collide with other workflows' columns in a shared table.

## What it does

| Step | Tool | Notes |
| ---- | ---- | ----- |
| Validate | | `read1` and `read2` together or neither. A lone `read1` is an error, not silently dropped |
| Type the assembly | LISTT `-m a` | BLAST against 10 reference ospA alleles; needs an 813+ bp hit |
| Type the reads *(optional)* | LISTT `-m r` | `bwa mem` to each reference, best by mean-minus-SD depth, QC (≥100 aligned reads, median depth ≥10, ≥99% of ospA at ≥5×), `bcftools` consensus |
| Classify each | this workflow | Closest of 118 reference OspA variants, its identity, and **the identity threshold of that variant's IST** |
| Reconcile | this workflow | One call, one flag |

LISTT's code runs unmodified. What this workflow adds is described under *Differences from LISTT*.

## Inputs

| Input | Type | Notes |
| ----- | ---- | ----- |
| `sample_name` | `String` | **Required.** Label for the outputs |
| `assembly_fasta` | `File` | **Required.** Short- or long-read assembly, plain or gzip. Required because it is the one input every isolate has |
| `read1`, `read2` | `File?` | Illumina paired-end FASTQ. Optional, but give both or neither. The read pileup is designed for Illumina PE data |
| `max_trusted_ist` | `Int` | Default `17`. See *Which ISTs to trust* |
| `cpu`, `mem_gb`, `disk_gb`, `preemptible` | | Defaults `2`, `4`, 10 + 2× input size, `1` |

## Outputs

| Output | Type | Notes |
| ------ | ---- | ----- |
| `ospA_call` | `String` | `IST<n>`, `Unknown` (flag `UK` or `QF`) or `Ambig` (flag `DS`) |
| `ospA_flag` | `String` | `OK`, `SM`, `UK`, `DS` or `QF` (below) |
| `ospA_species` | `String` | Genospecies of the closest reference variant; `a\|b` if the modes differ; `NA` for `QF` |
| `ospA_identity` | `String` | Amino-acid identity to the closest variant per mode, e.g. `assembly=100.00;reads=99.64` |
| `ospA_details` | `String` | Per-mode status, closest IST and variant, threshold, runner-up IST, reasons, LISTT version |
| `ospA_results_tsv` | `File` | One row per mode plus a `final` row, including ospA locus QC |
| `ospA_protein_fasta` | `File` | The consensus OspA protein from each mode |
| `ospA_raw_tar` | `File` | LISTT's own outputs (BAMs, VCFs, logs), for debugging. Leave unmapped unless needed |
| `ospA_listt_version` | `String` | LISTT release and commit |

### Flags

| Flag | `ospA_call` | Meaning |
| ---- | ----------- | ------- |
| `OK` | `IST<n>` | Assembly and reads both gave a confident call, and the same one |
| `SM` | `IST<n>` | One mode gave a confident call and nothing contradicts it: assembly-only runs, or the other mode failed QC |
| `UK` | `Unknown` | No confident IST. The closest is in `ospA_details` and `ospA_results_tsv`. Either the best match is **under its IST's identity threshold** (a candidate novel allele if the locus QC is clean), or its IST is not trusted (below) |
| `DS` | `Ambig` | The modes disagree, or one is under threshold for a different IST than the other's call |
| `QF` | `Unknown` | No usable ospA in any supplied mode. `ospA_details` says why: `no_ospa_hit`, `ospa_incomplete_<n>nt`, `reads_qc_failed:<metric>`, `no_ospa_alignment`, `premature_stop_codon` |

Assembly and reads share the same sequencing run, so agreement is a consistency check and not
independent confirmation. Disagreement is the useful signal: a collapsed or fragmented assembly,
a mixed infection, or an ospA cut by a contig break.

`UK` is not `QF`. A `UK` whose `nt_identity` is near 100 and whose `at_contig_edge` is `no`
(`locus_qc` in `ospA_results_tsv`) is a different situation from a `QF`. When both modes recover the
identical protein, `ospA_details` says so.

## Differences from LISTT

1. **Identity thresholds are enforced.** The paper assigns an IST to a non-exact match only if its
   identity meets that IST's threshold, otherwise "unknown IST". LISTT's released code never reads
   `variants/thresholds.csv`, so it assigns the nearest variant whatever the identity. Here an exact
   match is a call; a non-exact match is a call only at or above the threshold in
   `variants/thresholds.csv`, and is otherwise `UK`, with the closest IST still reported.
2. **A cut-off on which ISTs are called** (`max_trusted_ist`).
3. **LISTT's crash on variant 68 is survived.** See *Known problems in LISTT's tables*.
4. **Inputs are staged under a safe name.** LISTT takes the sample name from the filename up to the
   first `_` (reads) or `.f` (genome), so `B31_S1.fasta` would become `B31`.
5. **LISTT exits 0 on failure**, so the result is judged from the files it leaves behind.

## Which ISTs to trust

LISTT 1.2.0 knows 23 ISTs. The paper (Lee et al. 2024) describes IST1-17, and `max_trusted_ist`
defaults to 17. Nothing public that I could find documents ISTs 18 and up beyond the LISTT
commit history: they were added in LISTT 1.1.0 (Dec 2024, as ISTs 18-25), and renumbered to 18-23
in 1.2.0 (May 2025), when several genospecies moved (*B. finlandensis* 21→14, *B. valaisiana* 25→16,
*B. turdi* 24→22, *B. japonica* 22→20, *B. carolinensis* 20→19). A closest match to IST18+ is
reported as `UK`, with the IST shown. Raise `max_trusted_ist` to call them.

`docker/listt/audit_listt_tables.py` re-derives each IST's threshold by the paper's method on the
current 118 variants. For ISTs 1-9 and 12 the result is within 0.4 of the shipped value, and for 11
and 13 within 1.3. The rest differ materially:

| IST | Variants | Shipped | Re-derived | Note |
| --- | -------- | ------- | ---------- | ---- |
| 10 | 9 | 87.40 | 84.15 | Within-IST and out-of-IST identity are 0.7 points apart; includes a *B. garinii* variant |
| 14 | 3 | 96.70 | 93.84 | Set when IST14 was *B. mayonii* only; *B. finlandensis* was added later. A *B. finlandensis*-like allele may fall under the shipped threshold |
| 15 | 4 | 90.29 | n/a | Not separable: its members match other ISTs at up to 91.55% but each other as low as 85.21% |
| 16 | 5 | 90.29 | 92.78 | *B. valaisiana* variants were added after the paper |
| 17 | 3 | 91.90 | 93.84 | |
| 18-22 | 1-5 | 93.00 (98.17 for 20) | 88.6-94.4 | 18, 19, 21 and 22 share one value, which looks like a placeholder; IST21 has one variant |
| 23 | 3 | 80.14 | 76.06 | *B. lusitaniae* is distant from everything else; the low value is plausible |

The shipped thresholds are the ones used. The re-derived values are for reference.

## Known problems in LISTT's tables

Present in 1.1.0 and 1.2.0; `docker/listt/audit_listt_tables.py` lists them.

- **Variant 68** is in `ospA_protein.fas` but has no row in `serotypes.csv`. An isolate whose
  closest variant is 68 makes LISTT stop with a `KeyError`. Variant 68 is closest to variant 27
  (IST10, 98.6%). Reported here as `UK` with `variant_not_in_listt_table`.
- **Variant 107** (*B. tanukii*, IST21) is listed three times.
- **IST3 includes a variant named *B. afzelii*** (73) and **IST10 a *B. garinii*** (85). The
  reported species follows the closest variant, not the IST.
- **Eight variants have no IST** (`unk`), including two named *B. garinii* (90, 117).
- The README says 25 ISTs; the data has 23.

## Validation

The paper tested reads from 22 Valneva isolates covering serotypes 1-6: all concordant with
serology, plus one *B. garinii* isolate with a novel ospA allele at 99% identity to IST3. It does
not report an assembly-mode validation. This repository's tests are regression tests, not
validation: two PubMLST genomes (IST1, IST4) and synthetic genes with planted properties.
Comparing assembly and reads calls on a subset of real isolates is the intended check of
assembly mode.

Limits stated by the authors: a co-infection is typed by its dominant IST only, and genospecies
inferred from IST could be wrong after horizontal transfer.

## Running locally

```bash
miniwdl run workflows/call_ospA_type/call_ospA_type.wdl -i tests/call_ospA_type_inputs.json
# the fixture is a 4 kb window around ospA from PubMLST's B31 genome with simulated reads: expect IST1, flag OK
tests/run_unit_tests.sh      # includes this workflow's regression section
```

Reads mode makes ten `bwa mem` passes over all the reads (one per reference allele). On the
fixtures' simulated 25x genome-wide reads (130-190k pairs) it took 1-1.5 minutes under amd64
emulation on an Apple-silicon Mac; real MiSeq data is deeper, so expect it to be the slower mode.

## Credit

LISTT is by Jonathan T. Lee and colleagues at Pfizer, released under Apache-2.0 at
<https://github.com/pfizer-opensource/LISTT>; the image pins release 1.2.0 (commit `b9c6569`) and
carries its licence. If you use this workflow, cite:

> Lee JT, Li Z, Nunez LD, et al. Development of a sequence-based *in silico* OspA typing method
> for *Borrelia burgdorferi* sensu lato. *Microbial Genomics* 2024;10:001252.
> doi:[10.1099/mgen.0.001252](https://doi.org/10.1099/mgen.0.001252)

The IST scheme is also part of the [PubMLST](https://pubmlst.org/organisms/borrelia-spp) *Borrelia*
database.
