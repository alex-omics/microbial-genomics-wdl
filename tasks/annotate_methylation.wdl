version 1.0

task annotate_methylation {

    input {
        File     bedmethyl
        String   sample_name
        File     reference_gff
        File     reference_fai
        Int      flank_upstream    = 300
        Boolean  trim_to_intergenic = true
        Int      min_coverage      = 10
        String   feature_type      = "CDS"
        Int      cpu               = 2
        Int      mem_gb            = 8
        Int      disk_gb           = 50
        String   docker            = "staphb/bedtools:2.31.1"
    }

    parameter_meta {
        bedmethyl:          "bedMethyl from modkit_pileup, in the reference's coordinate space"
        sample_name:        "Some identifier for naming outputs"
        reference_gff:      "Annotation for the reference. Prefer the reference's own curated GFF (for PAO1, the RefSeq/Pseudomonas Genome DB annotation carrying PA numbers) over a fresh Bakta run — the locus tags become the join key to transcriptomics and proteomics, and Bakta mints new ones that match nothing."
        reference_fai:      "samtools faidx index of the reference, used to build the genome file bedtools flank requires"
        flank_upstream:     "Bases upstream of each feature to treat as putative promoter/regulatory region (default = 300)"
        trim_to_intergenic: "Subtract annotated feature bodies from the upstream windows. Bacterial genomes are operonic and densely packed, so a fixed upstream window routinely lands inside the neighbouring gene (default = true)"
        min_coverage:       "Minimum Nvalid_cov for a site to be carried into the annotated table (default = 10)"
        feature_type:       "GFF feature type to annotate against (default = CDS)"
        cpu:                "Number of CPUs delegated to task (default = 2)"
        mem_gb:             "Amount of memory in GB delegated to task (default = 8)"
        disk_gb:            "Amount of disk space in GB delegated to task (default = 50)"
        docker:             "Container image"
    }

    command <<<
        set -euo pipefail

        # bedtools flank needs chrom sizes; the .fai already has them.
        cut -f1,2 ~{reference_fai} > genome.txt

        # GFF is 1-based inclusive, BED is 0-based half-open. Converting here
        # explicitly (start-1, end unchanged) rather than leaning on bedtools'
        # extension sniffing means the coordinate shift is visible and auditable
        # instead of being a silent one-base error in every downstream motif.
        #
        # locus_tag is pulled preferentially because it is the stable join key
        # to expression and abundance data; gene and product ride along for
        # readability.
        awk -F'\t' -v ftype="~{feature_type}" 'BEGIN {OFS="\t"}
            /^#/ {next}
            $3==ftype {
                lt="NA"; gn="NA"; pr="NA"
                n=split($9, kv, ";")
                for (i=1; i<=n; i++) {
                    if (kv[i] ~ /^locus_tag=/) { lt=substr(kv[i], 11) }
                    else if (kv[i] ~ /^gene=/)  { gn=substr(kv[i], 6)  }
                    else if (kv[i] ~ /^product=/) { pr=substr(kv[i], 9) }
                }
                gsub(/\t/, " ", pr)
                print $1, $4-1, $5, lt, gn, $7, pr
            }' ~{reference_gff} | sort -k1,1 -k2,2n > features.bed

        if [ ! -s features.bed ]; then
            echo "ERROR: no ~{feature_type} features parsed from ~{reference_gff}." >&2
            echo "       Check the feature type and that this is a GFF3." >&2
            exit 1
        fi
        echo "Parsed $(wc -l < features.bed) ~{feature_type} features"

        # Coverage floor before annotating, so the table is not padded with
        # positions that were never confidently observed.
        awk -F'\t' -v mincov=~{min_coverage} 'BEGIN {OFS="\t"} $10>=mincov' ~{bedmethyl} \
            | sort -k1,1 -k2,2n > sites.bed

        # Contig-name mismatch between the pileup and the annotation is the
        # classic silent failure: intersect returns zero rows and looks exactly
        # like "no methylation here". Check for shared names before trusting an
        # empty result.
        SHARED=$(comm -12 \
            <(cut -f1 sites.bed    | sort -u) \
            <(cut -f1 features.bed | sort -u) | wc -l)
        if [ "${SHARED}" -eq 0 ]; then
            echo "ERROR: no contig names shared between the bedMethyl and the annotation." >&2
            echo "  bedMethyl: $(cut -f1 sites.bed    | sort -u | head -3 | tr '\n' ' ')" >&2
            echo "  annotation: $(cut -f1 features.bed | sort -u | head -3 | tr '\n' ' ')" >&2
            exit 1
        fi

        # --- Sites inside feature bodies -------------------------------------
        #
        # Deliberately NOT using -s. Methylation is strand-specific and many
        # bacterial RM motifs are palindromic, methylated on both strands; a
        # site on the minus strand inside a plus-strand gene is still relevant
        # to that gene. Both strands are kept and the distinction is preserved
        # in the output columns for filtering later.
        # Column map after -wa -wb: bedMethyl is 18 columns and features.bed is
        # 7, so the feature fields land at 19-25 (22=locus_tag, 23=gene,
        # 24=strand, 25=product). Verified against bedtools 2.31.1 rather than
        # counted by eye.
        bedtools intersect -a sites.bed -b features.bed -wa -wb \
            | awk -F'\t' 'BEGIN {OFS="\t"} {
                print $1, $2, $3, $6, $4, $10, $11, $22, $23, $24, $25, "genic"
              }' > genic.tsv || true

        # --- Sites in upstream / putative promoter regions --------------------
        bedtools flank -i features.bed -g genome.txt -l ~{flank_upstream} -r 0 -s \
            | sort -k1,1 -k2,2n > upstream_raw.bed

        ~{if trim_to_intergenic
            then "bedtools subtract -a upstream_raw.bed -b features.bed | sort -k1,1 -k2,2n > upstream.bed"
            else "cp upstream_raw.bed upstream.bed"}

        if [ -s upstream.bed ]; then
            bedtools intersect -a sites.bed -b upstream.bed -wa -wb \
                | awk -F'\t' 'BEGIN {OFS="\t"} {
                    print $1, $2, $3, $6, $4, $10, $11, $22, $23, $24, $25, "upstream"
                  }' > upstream_sites.tsv || true
        else
            : > upstream_sites.tsv
        fi

        {
            printf "sample\tchrom\tstart\tend\tsite_strand\tmod_code\tn_valid_cov\tpercent_modified\tlocus_tag\tgene\tfeature_strand\tproduct\tregion\n"
            cat genic.tsv upstream_sites.tsv \
                | sort -k1,1 -k2,2n \
                | awk -F'\t' -v s="~{sample_name}" 'BEGIN {OFS="\t"} {print s, $0}'
        } > ~{sample_name}_methylation_annotated.tsv

        awk 'NR>1' ~{sample_name}_methylation_annotated.tsv | wc -l > N_ANNOTATED
        awk -F'\t' 'NR>1 && $13=="genic"    {n++} END {print n+0}' ~{sample_name}_methylation_annotated.tsv > N_GENIC
        awk -F'\t' 'NR>1 && $13=="upstream" {n++} END {print n+0}' ~{sample_name}_methylation_annotated.tsv > N_UPSTREAM

        # Distinct loci carrying at least one methylation call, which is the
        # count that actually maps onto a gene list for multi-omics joining.
        awk -F'\t' 'NR>1 && $9!="NA" {seen[$9]=1} END {print length(seen)}' \
            ~{sample_name}_methylation_annotated.tsv > N_LOCI
    >>>

    output {
        File    annotated_tsv   = "~{sample_name}_methylation_annotated.tsv"
        File    features_bed    = "features.bed"
        File    upstream_bed    = "upstream.bed"
        Int     n_annotated     = read_int("N_ANNOTATED")
        Int     n_genic         = read_int("N_GENIC")
        Int     n_upstream      = read_int("N_UPSTREAM")
        Int     n_loci_with_methylation = read_int("N_LOCI")
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
