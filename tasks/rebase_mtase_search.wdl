version 1.0

task rebase_blastp {

    input {
        File     faa
        String   sample_name
        File     rebase_goldset_fasta
        String   evalue          = "1e-25"
        Int      cpu             = 4
        Int      mem_gb          = 8
        Int      disk_gb         = 20
        String   docker          = "staphb/blast:2.16.0@sha256:48cb5071ad646bd062fbb6f170bc9950912db0428d4175a30ac4706cd36c3240"
    }

    parameter_meta {
        faa:                   "This isolate's predicted proteins, from bakta.faa"
        sample_name:           "Some identifier for naming outputs"
        rebase_goldset_fasta:  "REBASE Gold Standard protein set: experimentally characterised MTases only, not the full REBASE database of putative homologs. Stage as a File input rather than baking into the image or committing to the repo -- both MPore (GPL-3.0) and REBASE itself (NEB's own terms) carry redistribution terms that a runtime input sidesteps."
        evalue:                "BLASTP e-value cutoff for a homology call, passed straight through to blastp -evalue. Deliberately a String, not a Float: WDL renders a Float this small as a fixed 6-decimal string, so 1e-25 becomes the literal text \"0.000000\" and blastp rejects it as non-positive. 1e-25 matches the threshold MPore itself uses for candidate MTase identification (default = \"1e-25\")"
        cpu:                   "Number of CPUs delegated to task (default = 4)"
        mem_gb:                "Amount of memory in GB delegated to task (default = 8)"
        disk_gb:                "Amount of disk space in GB delegated to task (default = 20)"
        docker:                "Container image"
    }

    meta {
        description: "BLASTP this isolate's predicted proteins against the REBASE Gold Standard set. Split from the motif join into its own task deliberately: the BLAST image carries no Python, and reaching for a heavier combined image would be worse than one more small, single-purpose task in the same pattern as the rest of this pipeline."
    }

    command <<<
        set -euo pipefail

        makeblastdb -in ~{rebase_goldset_fasta} -dbtype prot -out goldset_db

        # -max_target_seqs 1 keeps one row per query gene: the single best
        # homology call, not every REBASE paralogue it resembles.
        blastp \
            -query ~{faa} \
            -db goldset_db \
            -evalue ~{evalue} \
            -max_target_seqs 1 \
            -outfmt '6 qseqid sseqid pident length evalue bitscore' \
            -num_threads ~{cpu} \
            > ~{sample_name}_blast_hits.tsv

        echo "$(wc -l < ~{sample_name}_blast_hits.tsv) candidate MTase hits at e-value <= ~{evalue}"
    >>>

    output {
        File blast_hits = "~{sample_name}_blast_hits.tsv"
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


task rebase_join_motifs {

    input {
        File     blast_hits
        String   sample_name
        File     rebase_motif_tsv
        String   docker = "python:3.11-slim@sha256:db3ff2e1800a8581e2c48a27c3995339d47bdf046da21c7627accd3d51053a93"
    }

    parameter_meta {
        blast_hits:        "Output of rebase_blastp: qseqid, sseqid, pident, length, evalue, bitscore (no header)"
        sample_name:       "Some identifier for naming outputs"
        rebase_motif_tsv:  "Enzyme -> recognition motif -> modification type -> strand table, keyed on the same enzyme names as the BLAST subject database"
        docker:            "Container image"
    }

    meta {
        description: "Join BLASTP hits against REBASE's motif table to report each candidate MTase's predicted recognition motif and modification type. This is inventory, not a statistical test: it says which restriction-modification systems the isolate is likely to carry and which motifs they SHOULD be protecting with methylation, which is what a genome's own motif occurrences can then be checked against."
    }

    command <<<
        set -euo pipefail

        cat > join.py <<'PY'
        import csv, sys
        from collections import defaultdict

        hits_path, motif_path, sample, out_path = sys.argv[1:5]

        # TSV_REBASE_data.tsv carries one row per methylated position within a motif
        # (e.g. three rows for one enzyme recognising CCGCGG, positions 1/2/4), so the
        # same (enzyme, motif, mod type) triple repeats. Collapse to distinct triples;
        # the exact intra-motif position is REBASE's own bookkeeping, not something
        # needed to test whether THIS genome's copies of the motif are methylated.
        enzyme_motifs = defaultdict(set)
        with open(motif_path, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                enzyme = row.get("Enzyme", "").strip()
                motif  = row.get("Motif", "").strip()
                modt   = row.get("MethylationType", "").strip()
                if enzyme and motif and motif != "?":
                    enzyme_motifs[enzyme].add((motif, modt))

        with open(hits_path, newline="") as fh, open(out_path, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t", lineterminator="\n")
            w.writerow(["sample", "locus_tag", "rebase_enzyme", "pident", "evalue",
                        "bitscore", "predicted_motif", "modification_type"])
            n_with_motif = 0
            for line in fh:
                qseqid, sseqid, pident, length, ev, bitscore = line.rstrip("\n").split("\t")
                motifs = enzyme_motifs.get(sseqid, set())
                if motifs:
                    n_with_motif += 1
                    for motif, modt in sorted(motifs):
                        w.writerow([sample, qseqid, sseqid, pident, ev, bitscore, motif, modt])
                else:
                    # Homology confirmed but REBASE has no characterised recognition
                    # sequence for this exact enzyme entry. Still real signal -- an
                    # MTase is present -- so keep the row with the motif columns blank
                    # rather than dropping it.
                    w.writerow([sample, qseqid, sseqid, pident, ev, bitscore, "NA", "NA"])

        sys.stderr.write("%d hits carried a characterised motif\n" % n_with_motif)
        PY
        # No de-indenting here: WDL strips the command block's common
        # leading whitespace before the shell ever sees it, so the heredoc
        # lands with the python's relative indentation already correct.

        python3 join.py ~{blast_hits} ~{rebase_motif_tsv} ~{sample_name} ~{sample_name}_rebase_mtases.tsv

        awk 'NR>1' ~{sample_name}_rebase_mtases.tsv | cut -f2 | sort -u | wc -l > N_CANDIDATE_GENES
        awk -F'\t' 'NR>1 && $7!="NA" {print $7}' ~{sample_name}_rebase_mtases.tsv | sort -u | wc -l > N_MOTIFS
    >>>

    output {
        File   rebase_mtases      = "~{sample_name}_rebase_mtases.tsv"
        Int    n_candidate_genes  = read_int("N_CANDIDATE_GENES")
        Int    n_predicted_motifs = read_int("N_MOTIFS")
    }

    runtime {
        docker:         docker
        memory:         "2 GB"
        cpu:            1
        disks:          "local-disk 10 SSD"
        preemptible:    1
        maxRetries:     2
    }
}
