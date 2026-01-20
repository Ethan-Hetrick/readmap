/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_readmap_pipeline'

include { FASTQ_FASTQC_UMITOOLS_FASTP                 } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp/main'
include { BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS } from '../subworkflows/nf-core/bam_variant_calling_sort_freebayes_bcftools/main'
include { MINIMAP2_ALIGN                              } from '../modules/nf-core/minimap2/align/main.nf'
include { SAMTOOLS_COVERAGE                           } from '../modules/nf-core/samtools/coverage/main'
include { CUSTOMSTATS                                 } from '../modules/local/customstats/main'
include { SAMTOOLS_SORMADUP                           } from '../modules/nf-core/samtools/sormadup/main'
include { BAM_STATS_SAMTOOLS                          } from '../subworkflows/nf-core/bam_stats_samtools/main'
include { SAMTOOLS_INDEX                              } from '../modules/nf-core/samtools/index/main'
include { BCFTOOLS_STATS                              } from '../modules/nf-core/bcftools/stats/main'

include { KRAKEN2_KRAKEN2                             } from '../modules/nf-core/kraken2/kraken2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow READMAP {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_reference = ch_samplesheet.map { meta, it -> [ meta, params.reference ] }
    ch_fai = ch_samplesheet.map { meta, it -> [ meta, params.fai ] }
    ch_reads = ch_samplesheet.map { meta, it -> [ meta, it, [] ] }

    //
    // SUBWORKFLOW: Run FastQC and FastP
    //
    FASTQ_FASTQC_UMITOOLS_FASTP (
        ch_reads,
        false,
        false,
        true,
        [],
        false,
        false,
        false,
        1
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions.first())

    //
    // MODULE: Kraken2
    //

    KRAKEN2_KRAKEN2 (
        FASTQ_FASTQC_UMITOOLS_FASTP.out.reads,
        params.kraken2_db,
        false,
        true
    )

    MINIMAP2_ALIGN (
        FASTQ_FASTQC_UMITOOLS_FASTP.out.reads,
        ch_reference,
        true,
        'bai',
        false,
        false
    )

    SAMTOOLS_SORMADUP (
        MINIMAP2_ALIGN.out.bam,
        ch_reference
    )

    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_SORMADUP.out.metrics.collect{it[1]}.ifEmpty([]))

    ch_sorted_bam = SAMTOOLS_SORMADUP.out.bam

    SAMTOOLS_INDEX (
        ch_sorted_bam
    )

    ch_bam_index = SAMTOOLS_INDEX.out.bai

    ch_bam_combined = ch_sorted_bam.join(ch_bam_index)

    BAM_STATS_SAMTOOLS (
        ch_bam_combined,
        ch_reference
    )

    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.flagstat.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]))

    //
    // MODULE: Samtools coverage
    //

    SAMTOOLS_COVERAGE (
        ch_bam_combined,
        ch_reference,
        ch_fai
    )
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_COVERAGE.out.coverage.collect{it[1]}.ifEmpty([]))

    if (params.run_variant_calling) {

        //
        // SUBWORKFLOW: Run variant calling with FreeBayes and BCFtools
        //

        ch_input = ch_bam_combined.map { meta, bam, index -> [ meta, bam, index, [], [], [] ] }
        ch_fasta_fai = ch_reference.join(ch_fai)

        BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS (
            ch_input,       // channel: [mandatory] [ val(meta), path(input1), path(index1), path(input2), path(index2), path(bed) ]
            ch_fasta_fai,   // channel: [mandatory] [ val(meta2), path(fasta), path(fai) ]
            [ [], [] ],     // channel: [optional]  [ val(meta3), path(samples) ]
            [ [], [] ],     // channel: [optional]  [ val(meta4), path(populations) ]
            [ [], [] ]      // channel: [optional]  [ val(meta5), path(cnv) ]
        )

        ch_vcf = BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.vcf
        ch_tbi = BAM_VARIANT_CALLING_SORT_FREEBAYES_BCFTOOLS.out.tbi

        //
        // MODULE: BCFtools stats
        //

        ch_bcftools_input = ch_vcf.join(ch_tbi)

        BCFTOOLS_STATS (
            ch_bcftools_input,
            [ [], [] ],
            [ [], [] ],
            [ [], [] ],
            [ [], [] ],
            [ [], [] ]
        )

        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.collect{it[1]}.ifEmpty([]))

        //
        // MODULE: Custom stats
        //

        ch_bcf_stats = BCFTOOLS_STATS.out.stats.ifEmpty([ [], [] ])

        CUSTOMSTATS (
            BAM_STATS_SAMTOOLS.out.stats,
            SAMTOOLS_COVERAGE.out.coverage,
            ch_bcf_stats
        )

    } else {

        CUSTOMSTATS (
            BAM_STATS_SAMTOOLS.out.stats,
            SAMTOOLS_COVERAGE.out.coverage,
            [ [], [] ]
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'readmap_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
