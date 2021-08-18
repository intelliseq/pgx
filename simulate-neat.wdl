workflow simulate_neat_workflow {

  meta {
    keywords: '{"keywords": ["some","keywords"]}'
    name: 'simulate_neat'
    author: 'https://gitlab.com/lltw'
    copyright: 'Copyright 2019 Intelliseq'
    description: ''
    changes: '{"1.0.0": "new docker"}'

    input_simulated_sample_name: '{"name": "Simulated sample name", "type": "String", "description": "Simulated sample name."}'

    input_reference_genome: '{"name": "Reference genome","type": "String","constraints": {"values": ["hg38","grch38-no-alt"]},"description": "Version of the reference genome to simulate reads from."}'
    input_chromosome: '{"name": "Chromosome","type": "String","required": "false","constraints": {"values": ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY-and-the-rest"]},"description": "Restrict reference genome to specified chromosoms (read will be simulated only from this chromosome)."}'

    input_ploidy: '{"name": "Ploidy","type": "Int","default": "2","constraints": {"min": "1"},"description": "Ploidy"}'

    input_read_length: '{"name": "Read length","type": "Int","default": "150","constraints": {"min": "20"},"description": "Read length."}'
    input_average_coverage: '{"name": "Average coverage","type": "Float","default": "150","constraints": {"min": "0"},"description": "Average coverage."}'
    input_simulate_paired_end_reads: '{"name": "Simulate paired end reads","type": "Boolean","default": "true","description": "Simulate paired end reads."}'
    input_paired_end_fragment_lenghth_mean: '{"name": "Paired end fragment lenghth mean","type": "Int","required": "false","default": "false","description": "Paired end fragment lenghth mean. Applicable only if pair-end option is set true."}'
    input_paired_end_fragment_lenghth_std: '{"name": "Paired end fragment lenghth standard deviation","type": "Float","required": "false","default": "false","description": "Paired end fragment lenghth standard deviation. Applicable only if pair-end option is set true."}'

    input_empirical_fragment_length_distribution: '{"name": "Empirical fragment length distribution (Advanced)", "type": "String", "constraints": {"values": ["wes-pe150-20m-150gb", "wgs-50x-pcr-free-150gb"]}, "required": "false", "description": "Select model of empirical fragment length distribution. Applicable only if pair-end option is set true. If specified, paired end fragment length mean and deviation will be ignored."}'
    input_empirical_GC_coverage_bias_distribution: '{"name": "Empirical GC coverage bias distribution (Advanced)", "type": "String", "constraints": {"values": ["wes-pe150-20m-150gb", "wgs-50x-pcr-free-150gb"]}, "required": "false", "description": "Select model of empirical GC coverage bias distribution."}'

    input_predefined_targeted_regions: '{"name": "Predefined targeted regions","type": "String","constraints": {"values": ["whole-reference-sequence","sureselect-human-all-exon-v6-r2","sureselect-human-all-exon-v6-utr-r2","sureselect-human-all-exon-v7"]},"default": "sureselect-human-all-exon-v6-r2","description": "Predefined sets of regions to to generate reads from. If provided, Targeted regions BED file (user defined) will be used instead of predefined set. Default coverage for targeted regions is 98% of avarage coverage, default coverage outside targeted regions is 2% of avarage coverage."}'
    input_user_defined_targeted_regions_bed_file: '{"name": "Targeted regions BED file (user defined)","type": "File","required": "false","description": "User defined BED file containing regions to generate reads from. If not specified, predefined targeted regions will be used. Default coverage for targeted regions is 98% of avarage coverage, default coverage outside targeted regions is 2% of avarage coverage."}'
    input_off_target_coverage_scalar: '{"name": "Off-target coverage","type": "Float","default": "0.02","constraints": {"min": "0","max": "1"},"description": "Average off-target coverage as a percent of average coverage"}'

    input_vcf_gz_with_variants_to_insert: '{"name": "Bgzippped VCF file with mandatory variants","type": "File","required": "false", "extension": [".vcf.gz"],"description": "Bgzippped VCF file with mandatory variants. variants specified in the VCF filee will be inserted into simulated sample"}'

    # input_sequencing_error_model: '{"name": "Sequencing error model (Advanced)","type": "String","required": "false","description": "Sequencing error model - can be generated with genSeqErrorModel.py from NEAT-genReads"}'
    # input_rescale_avg_sequencing_error_rate_to_this: '{"name": "Value to rescale avg sequencing error rate (Advanced)","type": "Float","required": "false","description": "Rescale average sequencing error rate to this value."}'
    # input_mutation_model_pickle_file
    # input_rescale_avg_mutation_rate_to_this
    # input_positional_mut_rates_bed_file: '{"name": "Positional mutation rates BED file (Advanced)","type": "File","required": "false","description": "Positional mutation rates BED file."}'
    # input_qual_threshold_n:

    input_output_golden_bam_file: '{"name": "Output golden BAM file","type": "Boolean","default": "false","description": "Output golden BAM file"}'
    input_output_golden_vcf_file: '{"name": "Output golden VCF file","type": "Boolean","default": "true","description": "Output golden VCF file"}'
    input_output_fasta_instead_of_fastq: '{"name": "Output FASTA instead of FASTQ","type": "Boolean","default": "false","description": "Output FASTA instead of FASTQ"}'
    input_bypass_fastq_generation: '{"name": "Bypass FASTQ generation","type": "Boolean","default": "false","description": "Bypass FASTQ generation"}'
    input_rng_seed_value: '{"name": "RNG seed value","type": "Int","default": "1","description": "RNG seed value (Identical RNG value should produce identical runs of the program)."}'

    input_no_simulate_neat_jobs: '{"name": "(Advanced) Number of NEAT-genRead jobs","type": "Int","default": "1","description": "(Advanced) Number of NEAT-genRead jobs."}'
    input_average_mutation_rate: '{"name": "average_mutation_rate", "type": "String", "required": "false", "constraints": {"min": "0", "max": "0.3"}, "description": "Mutation rate. Can be between 0 - 0.3. If 0 no random mutations will be simulated in reads. Integer or float" }'

    output_stdout_log: '{"name": "Standard out","type": "File","copy": "True","description": "Standard out"}'
    output_stderr_err: '{"name": "Standard err","type": "File","copy": "True","description": "Standard error"}'
    output_bco: '{"name": "Biocompute object","type": "File","copy": "True","description": "Biocompute object"}'
  }

  call simulate_neat

}

task simulate_neat {

  Int? index
  String simulated_sample_name = 'simulated_sample_name'

  String? chromosome
  String reference_genome = "hg38"
  String reference_genome_scope =  if defined(chromosome) then reference_genome + "-" + chromosome else reference_genome

  Int ploidy = 2

  Int read_length = 150
  Float average_coverage = 30
  Boolean simulate_paired_end_reads = true
  Int paired_end_fragment_lenghth_mean = 300
  Int paired_end_fragment_lenghth_std = 30

  File? user_defined_targeted_regions_bed_file
  String predefined_targeted_regions = "agilent-sureselect-human-all-exon-v6-r2"
  String predefined_targeted_regions_command = if defined(user_defined_targeted_regions_bed_file) then " " else "-t /resources/intervals/" + predefined_targeted_regions + "/covered.bed"
  Float off_target_coverage_scalar = 0.00

  File? vcf_gz_with_variants_to_insert
  Boolean is_vcf_gz_with_variants_to_insert_present = defined(vcf_gz_with_variants_to_insert)

  String? empirical_fragment_length_distribution
  String paired_end_mean_and_std_command = if simulate_paired_end_reads && ! defined(empirical_fragment_length_distribution) then "--pe " + paired_end_fragment_lenghth_mean + " " + paired_end_fragment_lenghth_mean else ""
  String empirical_fragment_length_distribution_command = if defined(empirical_fragment_length_distribution) then "--pe-model /resources/neat-genreads/" + empirical_fragment_length_distribution + "/fraglen.p" else ""

  String? empirical_GC_coverage_bias_distribution
  String empirical_GC_coverage_bias_distribution_command =  if defined(empirical_GC_coverage_bias_distribution) then "--gc-model /resources/neat-genreads/" + empirical_GC_coverage_bias_distribution + "/model.p" else ""

  ########################################
  # Advanced options - to be added later #
  ########################################
  #
  #   String? sequencing_error_model
  #   Float? rescale_avg_sequencing_error_rate_to_this
  #   String? mutation_model_pickle_file
  #   Float? rescale_avg_mutation_rate_to_this
  #   File? positional_mut_rates_bed_file
  #   Int? qual_threshold_n
  #
  #### To be added to ${genreads_py} command in the commad section of wdl
  #
  #   ${default='' '-e ' + sequencing_error_model} \
  #   ${default='' '-E ' + rescale_avg_sequencing_error_rate_to_this}
  #   ${default='' '-m ' + mutation_model_pickle_file}
  #   ${default='' '-M ' + rescale_avg_mutation_rate_to_this}
  #   ${default='' '-Mb ' + positional_mut_rates_bed_file} \
  #   ${default='' '-N ' + qual_threshold_n} \

  String average_mutation_rate = 0 
  Int rng_seed_value = 1
  Boolean output_golden_bam_file = false
  Boolean output_golden_vcf_file = true
  Boolean output_fasta_instead_of_fastq = false
  Boolean bypass_fastq_generation = false

  # Tools runtime settings, paths etc.
  Int no_jobs = 20
  Boolean is_single_job = if no_jobs == 1 then true else false

  String genreads_py = "genReads.py"
  String mergejobs_py = "mergeJobs.py"
  String reference_genome_fasta_gz = "/resources/reference-genomes/" + reference_genome_scope + "/" + reference_genome_scope + ".fa.gz"
  String reference_genome_fasta = "/resources/reference-genomes/" + reference_genome_scope + "/" + reference_genome_scope + ".fa"
  String reference_genome_fasta_fai = "/resources/reference-genomes/" + reference_genome_scope + "/" + reference_genome_scope + ".fa.fai"

  String task_name = "simulate_neat"
  String task_name_with_index = if defined(index) then task_name + "_" + index else task_name
  String task_version = "1.0.0"
  String docker_image = if defined(chromosome) then "intelliseqngs/task_simulate-neat-" + reference_genome + "-chr-wise:1.0.0-" + chromosome else "intelliseqngs/task_simulate-neat-" + reference_genome + ":1.0.0"

  command <<<
    set -e -o pipefail
    bash /intelliseqtools/bco-after-start.sh --task-name-with-index ${task_name_with_index}

    # unbgzip vcf file with variants to insert...
    vcf_file_with_madatory_variants_command=""
    if ${is_vcf_gz_with_variants_to_insert_present}; then
      zcat ${vcf_gz_with_variants_to_insert} > variants.vcf
      vcf_file_with_madatory_variants_command="-v variants.vcf"
    fi

    # unbgzip reference genome...
    gunzip ${reference_genome_fasta_gz}


    if ${is_single_job}; then

        ###########################################
       ## run genReads.py in SINGLE THREAD MODE ##
      ###########################################

      INFIX=""

      ${genreads_py} \
        -o $PWD/${simulated_sample_name} \
        -r ${reference_genome_fasta} \
        -p ${ploidy} \
        -R ${read_length} \
        -c ${average_coverage} \
        -M ${average_mutation_rate} \
        ${default='' '-t ' + user_defined_targeted_regions_bed_file} ${predefined_targeted_regions_command} \
        -to ${off_target_coverage_scalar} \
        ${empirical_fragment_length_distribution_command} \
        ${paired_end_mean_and_std_command} \
        ${empirical_GC_coverage_bias_distribution_command} \
        $vcf_file_with_madatory_variants_command \
        ${true='--bam' false='' output_golden_bam_file} \
        ${true='--vcf' false='' output_golden_vcf_file} \
        ${true='--fa' false='' output_fasta_instead_of_fastq} \
        ${true='--no-fastq' false='' bypass_fastq_generation} \
        --rng ${rng_seed_value}

    else

        ######################################
       ## run genReads.py in PARRALEL MODE ##
      ######################################

      INFIX="_merged"

      parallel -j ${no_jobs} \
        ${genreads_py} \
          -o $PWD/${simulated_sample_name} \
          -r ${reference_genome_fasta} \
          -p ${ploidy} \
          -R ${read_length} \
          -c ${average_coverage} \
          -M ${average_mutation_rate} \
          ${default='' '-t ' + user_defined_targeted_regions_bed_file} ${predefined_targeted_regions_command} \
          -to ${off_target_coverage_scalar} \
          ${empirical_fragment_length_distribution_command} \
          ${paired_end_mean_and_std_command} \
          ${empirical_GC_coverage_bias_distribution_command} \
          $vcf_file_with_madatory_variants_command \
          ${true='--bam' false='' output_golden_bam_file} \
          ${true='--vcf' false='' output_golden_vcf_file} \
          ${true='--fa' false='' output_fasta_instead_of_fastq} \
          ${true='--no-fastq' false='' bypass_fastq_generation} \
          --rng ${rng_seed_value} \
          --job ::: {1..${no_jobs}} ::: ${no_jobs}

      # merge outputs with the use of mergeJobs.py
      ${mergejobs_py} -i $PWD/${simulated_sample_name} -o $PWD/${simulated_sample_name}"$INFIX" -s samtools

    fi

    # prepare final outputs
    if ${output_golden_bam_file}; then
      mv $PWD/${simulated_sample_name}"$INFIX"_golden.bam $PWD/${simulated_sample_name}-golden.bam
      samtools index $PWD/${simulated_sample_name}-golden.bam
    fi

    if ${output_golden_vcf_file}; then
      cat $PWD/${simulated_sample_name}"$INFIX"_golden.vcf | bgzip > $PWD/${simulated_sample_name}-golden.vcf.gz
      tabix -p vcf $PWD/${simulated_sample_name}-golden.vcf.gz
      rm $PWD/${simulated_sample_name}"$INFIX"_golden.vcf
    fi

    if ! ${bypass_fastq_generation}; then
      if ${output_fasta_instead_of_fastq}; then
        cat $PWD/${simulated_sample_name}"$INFIX"_read1.fa | bgzip > $PWD/${simulated_sample_name}_1.fa.gz
        cat $PWD/${simulated_sample_name}"$INFIX"_read2.fa | bgzip > $PWD/${simulated_sample_name}_2.fa.gz
        rm $PWD/${simulated_sample_name}"$INFIX"_read1.fa $PWD/${simulated_sample_name}"$INFIX"_read2.fa
      else
        cat $PWD/${simulated_sample_name}"$INFIX"_read1.fq | bgzip > $PWD/${simulated_sample_name}_1.fq.gz
        cat $PWD/${simulated_sample_name}"$INFIX"_read2.fq | bgzip > $PWD/${simulated_sample_name}_2.fq.gz
        rm $PWD/${simulated_sample_name}"$INFIX"_read1.fq $PWD/${simulated_sample_name}"$INFIX"_read2.fq
      fi
    fi


    bash /intelliseqtools/bco-before-finish.sh --task-name ${task_name} \
                                               --task-name-with-index ${task_name_with_index} \
                                               --task-version ${task_version} \
                                               --task-docker ${docker_image}
  >>>

  runtime {

    docker: docker_image
    memory: "80G"
    cpu: "20"
    maxRetries: 2

  }

  output {

    File? fastq_1_gz = "${simulated_sample_name}_1.fq.gz"
    File? fastq_2_gz = "${simulated_sample_name}_2.fq.gz"

    File? golden_bam = "${simulated_sample_name}-golden.bam"
    File? golden_bam_bai = "${simulated_sample_name}-golden.bam.bai"

    File? golden_vcf_gz = "${simulated_sample_name}-golden.vcf.gz"
    File? golden_vcf_gz_tbi = "${simulated_sample_name}-golden.vcf.gz.tbi"

    File? fasta_1_gz = "${simulated_sample_name}_1.fa.gz"
    File? fasta_2_gz = "${simulated_sample_name}_2.fa.gz"

    #File stdout_log = stdout()
    #File stderr_log = stderr()
    #File bco = "bco.json"

  }

}
