import "https://gitlab.com/intelliseq/workflows/raw/resources-kit@1.1.0/src/main/wdl/tasks/resources-kit/resources-kit.wdl" as resources_kit_task
import "https://gitlab.com/intelliseq/workflows/raw/fq-organize@2.1.1/src/main/wdl/tasks/fq-organize/fq-organize.wdl" as fq_organize_task
import "https://gitlab.com/intelliseq/workflows/raw/fq-bwa-align@1.5.3/src/main/wdl/modules/fq-bwa-align/fq-bwa-align.wdl" as fq_bwa_align_module
import "https://gitlab.com/intelliseq/workflows/raw/bam-varcalling@1.3.5/src/main/wdl/modules/bam-varcalling/bam-varcalling.wdl" as bam_varcalling_module
import "https://gitlab.com/intelliseq/workflows/raw/gvcf-genotype-by-vcf@2.0.2/src/main/wdl/tasks/gvcf-genotype-by-vcf/gvcf-genotype-by-vcf.wdl" as gvcf_genotype_by_vcf_task
import "https://gitlab.com/intelliseq/workflows/raw/pgx-genotyping-report@2.1.1/src/main/wdl/modules/pgx-genotyping-report/pgx-genotyping-report.wdl" as pgx_genotyping_report_module
import "https://gitlab.com/intelliseq/workflows/raw/bco@1.0.0/src/main/wdl/modules/bco/bco.wdl" as bco_module

workflow pgx {

  String sample_id = "no_id_provided"
  String reference_genome = "grch38-no-alt"
  # Java options for bam_metrics and bam_gatk_hc
  String java_mem_options =  "-Xms4g -Xmx63g" # use "-Xms4g -Xmx8g" on anakin

  # Start from fastq files
  Array[File]? fastqs
  Boolean is_fastqs_defined = defined(fastqs)
  Array[File]? fastqs_left
  Boolean is_fastqs_left_defined = defined(fastqs_left)
  Array[File]? fastqs_right

  # Start from bam file (with or without gvcf)
  File? bam
  File? bai

  # Start from gvcf file (only with bam)
  File? gvcf_gz
  File? gvcf_gz_tbi

  File? interval_list

  String kit = "custom_pgx" # It is not actually kit but only intervals to genotyping, it can be also "genome"
  String pipeline_name = "pgx"
  String pipeline_version = "1.1.4"

  # 1. Prepare interval_list if not defined
if (!defined(interval_list)) { 
    call resources_kit_task.resources_kit {
         input:
                kit = kit,
                reference_genome = reference_genome
        }
    File interval_list_resources = resources_kit.interval_list[0]
    }

  # 2a. Start analysis from fastq file
    if(is_fastqs_defined || is_fastqs_left_defined) {
        if(is_fastqs_defined) {
            # Organise fastq files
            call fq_organize_task.fq_organize {
                input:
                    fastqs = fastqs,
                    paired = true
            }
       }

    Array[File] fastqs_1 = select_first([fq_organize.fastqs_1, fastqs_left])
    Array[File] fastqs_2 = select_first([fq_organize.fastqs_2, fastqs_right])

  # Align reads
    call fq_bwa_align_module.fq_bwa_align {
                input:
                    fastqs_left = fastqs_1,
                    fastqs_right = fastqs_2,
                    sample_id = sample_id,
                    reference_genome = reference_genome
           }
    }

# 2b. Start from bam file
    File bam_to_genotype = select_first([fq_bwa_align.recalibrated_markdup_bam, bam])
    File bai_to_genotype = select_first([fq_bwa_align.recalibrated_markdup_bai, bai])


    File interval = select_first([interval_list_resources, interval_list])

# 2c. Call variants if no gvcf provided
    if(!defined(gvcf_gz)) {
        call bam_varcalling_module.bam_varcalling {
            input:
                input_bam = bam_to_genotype,
                input_bai = bai_to_genotype,
                sample_id = sample_id,
                interval_list = interval,
                max_no_pieces_to_scatter_an_interval_file = 6,
                run_bam_concat_task = false,
                run_gvcf_call_variants_task = false,
                haplotype_caller_java_mem = java_mem_options         }
     }
  
  File gvcf_to_genotype = select_first([bam_varcalling.gvcf_gz, gvcf_gz])
  File gvcf_tbi_to_genotype = select_first([bam_varcalling.gvcf_gz_tbi, gvcf_gz_tbi])

# 3. Pharmacogenomics module
    call pgx_genotyping_report_module.pgx_genotyping_report {
        input:
          gvcf_gz = gvcf_to_genotype, 
          gvcf_gz_tbi = gvcf_tbi_to_genotype,
          sample_id = sample_id,
          bam = bam_to_genotype,
          bai = bai_to_genotype
          }
    
    Array[File] bcos_module = select_all([resources_kit.bco, fq_organize.bco, fq_bwa_align.bco, bam_varcalling.bco, pgx_genotyping_report.bco])
    Array[File] stdout_module = select_all([resources_kit.stdout_log, fq_organize.stdout_log, fq_bwa_align.stdout_log, pgx_genotyping_report.stdout_log])
    Array[File] stderr_module = select_all([resources_kit.stderr_log, fq_organize.stderr_log, fq_bwa_align.stderr_log, bam_varcalling.stderr_log, pgx_genotyping_report.stderr_log])

     call bco_module.bco as gather_bco_info {
      input:
        bco_array = bcos_module,
        stdout_array = stdout_module,
        stderr_array = stderr_module,
        pipeline_name = pipeline_name,
        pipeline_version = pipeline_version,
        sample_id = sample_id
     }

  output {

    File? gvcf_out = bam_varcalling.gvcf_gz
    File? tsv = pgx_genotyping_report.tsv
    File? json = pgx_genotyping_report.json
    File bco = gather_bco_info.bco_merged
    File stdout_log = gather_bco_info.stdout_log
    File stderr_log = gather_bco_info.stderr_log

  }
}

