version 1.0

task validate_annotation {

    input {
        File    assembly
        File    gff3
        File    fna
        File    faa
        String  sample_name
        String  docker = "ubuntu:22.04@sha256:0e0a0fc6d18feda9db1590da249ac93e8d5abfea8f4c3c0c849ce512b5ef8982"
    }

    parameter_meta {
        assembly:    "The assembly the reads will be aligned to (plain or gzipped FASTA)"
        gff3:        "GFF3 annotation of that same assembly, from whatever annotator produced it (Bakta, Prokka, ...)"
        fna:         "The annotator's nucleotide FASTA of that assembly. Must carry the same contig names and lengths as `assembly`"
        faa:         "The annotator's protein FASTA. Each protein must be a feature in `gff3`"
        sample_name: "Isolate name. Returned unchanged as checked_name, so that anything consuming checked_name cannot start before these checks have passed"
        docker:      "Container image"
    }

    meta {
        description: "Fail fast if an isolate's assembly and annotation files do not belong together. They are matched positionally across five arrays, which is easy to get wrong, and the later stages cannot tell: a GFF from another isolate with generic contig names (contig_1, contig_2, ...) still intersects the pileup and silently labels sites with the wrong isolate's genes. Compares contig names and lengths between the assembly and the annotator's FASTA, that the GFF only names those contigs and agrees with its own ##sequence-region lengths, and that the proteins are features of the GFF. Also reports the assembly statistics the summary table needs, since those are otherwise only available if this workflow runs the annotator itself."
    }

    command <<<
        set -euo pipefail

        fail() { echo "ERROR [~{sample_name}]: $*" >&2; exit 1; }

        # name<TAB>length per record. The name is the first whitespace-delimited
        # token of the header, which is how minimap2, samtools faidx and the
        # annotators' GFFs all name a contig.
        contig_lengths() {
            zcat -f "$1" | awk '
                /^>/ { if (n != "") print n "\t" l; n = substr($1, 2); l = 0; next }
                { gsub(/[ \t\r]/, ""); l += length($0) }
                END  { if (n != "") print n "\t" l }'
        }

        contig_lengths ~{assembly} > assembly.lengths
        contig_lengths ~{fna}      > fna.lengths

        [ -s assembly.lengths ] || fail "the assembly contains no sequences."
        [ -s fna.lengths ]      || fail "the annotation FASTA (fna) contains no sequences."

        for f in assembly fna; do
            DUP="$(cut -f1 ${f}.lengths | sort | uniq -d | head -n 3 | tr '\n' ' ')"
            [ -z "${DUP}" ] || fail "duplicate contig names in the ${f} (${DUP}). Contig names must be unique."
        done

        # Assembly vs the annotator's FASTA: same contigs, same lengths.
        if ! cmp -s <(LC_ALL=C sort assembly.lengths) <(LC_ALL=C sort fna.lengths); then
            echo "ERROR [~{sample_name}]: the assembly and the annotation FASTA (fna) are not the same sequences." >&2
            echo "  Only in the assembly (name, length):" >&2
            comm -23 <(LC_ALL=C sort assembly.lengths) <(LC_ALL=C sort fna.lengths) | head -n 5 | sed 's/^/    /' >&2
            echo "  Only in the fna (name, length):" >&2
            comm -13 <(LC_ALL=C sort assembly.lengths) <(LC_ALL=C sort fna.lengths) | head -n 5 | sed 's/^/    /' >&2
            if cmp -s <(cut -f2 assembly.lengths) <(cut -f2 fna.lengths); then
                echo "  The contig lengths agree in file order but the names differ, so the annotator" >&2
                echo "  probably renamed the contigs. Supply its fna as the assembly, or re-annotate" >&2
                echo "  keeping the original headers." >&2
            else
                echo "  If the lengths differ too, the assembly and annotation arrays are probably in" >&2
                echo "  different orders, or the annotation is of a different assembly." >&2
            fi
            exit 1
        fi

        # GFF: it must only name contigs of this assembly. Stop at ##FASTA, after
        # which the lines are sequence, not features.
        awk -F'\t' '/^##FASTA/ {exit}  /^#/ {next}  NF {print $1}' ~{gff3} | LC_ALL=C sort -u > gff.seqids
        [ -s gff.seqids ] || fail "the GFF3 contains no features."
        UNKNOWN="$(comm -23 gff.seqids <(cut -f1 assembly.lengths | LC_ALL=C sort) | head -n 5 | tr '\n' ' ')"
        [ -z "${UNKNOWN}" ] || fail "the GFF3 names contigs that are not in the assembly (${UNKNOWN}). It is probably from a different isolate, or the annotator renamed the contigs."

        # ...and agree with its own sequence-region lengths where it declares them.
        awk '/^##FASTA/ {exit}  /^##sequence-region/ {print $2 "\t" ($4 - $3 + 1)}' ~{gff3} > gff.regions
        awk 'BEGIN {FS = OFS = "\t"}
             NR == FNR {len[$1] = $2; next}
             ($1 in len) && len[$1] != $2 {print $1, "assembly=" len[$1], "gff=" $2}' \
            assembly.lengths gff.regions > region_mismatch.txt
        if [ -s region_mismatch.txt ]; then
            echo "ERROR [~{sample_name}]: the GFF3 declares contig lengths that differ from the assembly:" >&2
            head -n 5 region_mismatch.txt | sed 's/^/    /' >&2
            echo "  The annotation is of a different assembly of this isolate, or another isolate's." >&2
            exit 1
        fi

        # Proteins: each header must be a feature in the GFF, by ID or locus_tag.
        # Any feature type counts, since annotators put short ORFs in the protein
        # FASTA under a type other than CDS. Tolerates a few strays rather than
        # failing on an annotator convention that was not anticipated; a protein
        # file from another isolate matches almost nothing.
        awk -F'\t' '/^##FASTA/ {exit}  /^#/ {next}
             {n = split($9, kv, ";")
              for (i = 1; i <= n; i++) {
                  if (kv[i] ~ /^ID=/)        print substr(kv[i], 4)
                  if (kv[i] ~ /^locus_tag=/) print substr(kv[i], 11)
              }}' ~{gff3} | LC_ALL=C sort -u > gff.ids
        awk '/^>/ {print substr($1, 2)}' ~{faa} | LC_ALL=C sort -u > faa.ids
        N_FAA="$(wc -l < faa.ids | tr -d ' ')"
        [ "${N_FAA}" -gt 0 ] || fail "the protein FASTA (faa) contains no sequences."
        N_STRAY="$(comm -23 faa.ids gff.ids | wc -l | tr -d ' ')"
        if [ "${N_STRAY}" -gt 0 ]; then
            echo "WARNING [~{sample_name}]: ${N_STRAY} of ${N_FAA} proteins in the faa are not features in the GFF3 (e.g. $(comm -23 faa.ids gff.ids | head -n 3 | tr '\n' ' '))" >&2
        fi
        if [ $((N_STRAY * 100)) -gt $((N_FAA * 5)) ]; then
            fail "more than 5% of the proteins in the faa are not features in the GFF3; the protein FASTA is probably from a different isolate."
        fi

        # Statistics, from the assembly itself.
        wc -l < assembly.lengths | tr -d ' ' > N_CONTIGS
        awk -F'\t' '{s += $2} END {print s + 0}' assembly.lengths > GENOME_SIZE
        zcat -f ~{assembly} | awk '
            /^>/ {next}
            {gc += gsub(/[GCgc]/, "&"); at += gsub(/[ATat]/, "&")}
            END  {if (gc + at > 0) printf "%.2f\n", 100 * gc / (gc + at); else print "NA"}' > GC_PERCENT
        awk -F'\t' '/^##FASTA/ {exit}  !/^#/ && $3 == "CDS" {n++}  END {print n + 0}' ~{gff3} > N_CDS
        if [ "$(cat N_CDS)" -eq 0 ]; then
            echo "WARNING [~{sample_name}]: the GFF3 has no CDS features." >&2
        fi

        echo "~{sample_name}: $(cat N_CONTIGS) contigs, $(cat GENOME_SIZE) bp, GC $(cat GC_PERCENT)%, $(cat N_CDS) CDS; assembly and annotation agree"
    >>>

    output {
        String  checked_name = sample_name
        Int     n_contigs    = read_int("N_CONTIGS")
        Int     genome_size  = read_int("GENOME_SIZE")
        String  gc_percent   = read_string("GC_PERCENT")
        Int     n_cds        = read_int("N_CDS")
    }

    runtime {
        docker:         docker
        memory:         "2 GB"
        cpu:            1
        disks:          "local-disk 20 SSD"
        preemptible:    1
        # A failed check is deterministic, so retrying it only delays the error.
        maxRetries:     0
    }
}
