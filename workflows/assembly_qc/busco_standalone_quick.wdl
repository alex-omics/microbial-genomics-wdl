version 1.0

import "../../tasks/busco.wdl" as busco_task

workflow busco_standalone {
    input {
        Array[File]   assemblies
        Array[String] sample_names
        String        busco_lineage = "spirochaetia_odb10"
    }

    scatter (pair in zip(assemblies, sample_names)) {
        call busco_task.busco {
            input:
                assembly      = pair.left,
                sample_name   = pair.right,
                busco_lineage = busco_lineage
        }
    }

    output {
        Array[File]   summaries     = busco.summary_txt
        Array[String] busco_results = busco.busco_result
        Array[String] lineages_used = busco.lineage_used
    }
}