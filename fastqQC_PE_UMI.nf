/* --------------------------
   Parameters
--------------------------- */

// File locations
params.input      = "${baseDir}/fastq/with/UMIs/*_{R1,R2}.fastq.gz"  	// paired-end pattern
params.salmonIdx = "${baseDir}/path/to/salmonIdx" 						//pre-built salmon index
params.outdir = "${baseDir}/out" 										// main output directory
params.smrIndex = '/path/to/idx'
params.riboFasta = '/path/to/rrna_sequences.fa'

// Fastq trimming
params.paired_end = true                       		// set true for PE
params.threads    = 8
params.min_len    = 30
params.q_right    = 20
params.max_ns     = 5
params.trim_poly  = true     						// poly-G/X trimming for NovaSeq/NextSeq
params.correction = true     						// PE overlap correction

// UMI parameters
params.umi_length = 8								// UMI length. Default is 8
params.umi_separator = ':'							// UMI separator. Default is ":"
params.subsampler = 'seqtk'							// Default is 'seqtk'


// Index for sortmerna *params.smrIndex must contain folder "idx" with the index files
IDX_CACHE_CH = Channel.value( file(params.smrIndex, type: 'dir') )

/*
   *   params.sample_prob  : e.g., 0.10  (mutually exclusive with sample_n)
   *   params.sample_n     : e.g., 100000 (mutually exclusive with sample_prob)
   *   params.seed         : RNG seed for deterministic subsampling (default 13)
   *   params.validateMaps : true/false for Salmon --validateMappings (default true)
   */
params.seed = 13
params.sample_n = 100000
//params.sample_prob = 0.1 // Not an option in this workflow, but can replace sample_n
params.validateMaps = "--validateMappings" // "--validateMappings" or ""


/* --------------------------
   Inputs
--------------------------- */
Channel
    .fromFilePairs( params.input, size: 2, flat: true )
    .ifEmpty { error "No paired FASTQs found at: ${params.input}" }
    .set { READS_PE }

/* 1) Baseline fastq metrics */
process INTEGRITY_STATS {
    tag "${base}"

    publishDir "${params.outdir}/qc/basic_qc", mode: 'copy'

    input:
    tuple val(base), path(r1), path(r2)

    output:
    path("${base}.integrity.log")
    path("${base}.seqkit.tsv")
	path("*.log")

    script:
    """
    set -euo pipefail
    echo "== Integrity check for ${base} ==" > ${base}.integrity.log
    (gunzip -t ${r1} && echo "R1 OK") >> ${base}.integrity.log 2>&1
    (gunzip -t ${r2} && echo "R2 OK") >> ${base}.integrity.log 2>&1
    echo -e "\\n== seqkit stats ==" >> ${base}.integrity.log
    seqkit stats -a ${r1} ${r2} | tee ${base}.seqkit.tsv >> ${base}.integrity.log
    """
}

/* 
Infer library stranded from a subset of fastq file pairs
- Subsampled reads are not saved
- Outputs are salmon analysis files
*/
process SUBSAMPLE_AND_INFER_STRAND {
  tag "$sample"

  publishDir "${params.outdir}/strandedness", mode: 'copy'

  input:
  tuple val(sample), path(r1), path(r2)
  val paired
  path salmon_index

  output:
  path "${sample}.lib_format_counts.json", emit: LIB_COUNTS
  path "${sample}.meta_info.json",        emit: META
  path "${sample}.salmon_infer.tsv",      emit: SUMMARY, optional: true

  /*
    Unnecessarily complicated tunable parameters:
      params.subsampler   : 'seqtk' | 'reformat'  (default 'seqtk')
      params.sample_prob  : e.g., 0.10  (ignored if sample_n is set)
      params.sample_n     : e.g., 100000
      params.seed         : e.g., 13
      params.validateMaps : true/false for Salmon --validateMappings
  */
  script:
  def seed     = params.seed ?: 13
  def frac     = params.sample_prob ?: 0.10
  def nReads   = params.sample_n ?: null
  def method   = (params.subsampler ?: 'seqtk').toString().toLowerCase()
  def validate = (params.validateMaps == false) ? "" : "--validateMappings"

  // Build subsample command
  def subsampleCmd
  if (method == 'reformat') {
    // BBTools reformat.sh (pair-aware; supports samplerate= or reads=)
    def amount = nReads ? "reads=${nReads}" : "samplerate=${frac}"
    subsampleCmd = paired ?
      """
      reformat.sh in1=${r1} in2=${r2} out1=subsampled/${sample}_R1.sub.fq.gz out2=subsampled/${sample}_R2.sub.fq.gz \
        ${amount} sampleseed=${seed} overwrite=t gz=t
      """
      :
      """
      reformat.sh in=${r1} out=subsampled/${sample}.sub.fq.gz \
        ${amount} sampleseed=${seed} overwrite=t gz=t
      """
  } else {
    // seqtk sample (use same seed for R1/R2; assumes synchronized mates)
    if (nReads) {
      subsampleCmd = paired ?
        """
        seqtk sample -s${seed} ${r1} ${nReads} | gzip -c > subsampled/${sample}_R1.sub.fq.gz
        seqtk sample -s${seed} ${r2} ${nReads} | gzip -c > subsampled/${sample}_R2.sub.fq.gz
        """
        :
        """
        seqtk sample -s${seed} ${r1} ${nReads} | gzip -c > subsampled/${sample}.sub.fq.gz
        """
    } else {
      subsampleCmd = paired ?
        """
        seqtk sample -s${seed} ${r1} ${frac} | gzip -c > subsampled/${sample}_R1.sub.fq.gz
        seqtk sample -s${seed} ${r2} ${frac} | gzip -c > subsampled/${sample}_R2.sub.fq.gz
        """
        :
        """
        seqtk sample -s${seed} ${r1} ${frac} | gzip -c > subsampled/${sample}.sub.fq.gz
        """
    }
  }

  def salmonCmd = paired ?
    """
    salmon quant -i ${salmon_index} -l A \
      -1 subsampled/${sample}_R1.sub.fq.gz -2 subsampled/${sample}_R2.sub.fq.gz \
      -o quant --threads ${task.cpus} ${validate}
    """
    :
    """
    salmon quant -i ${salmon_index} -l A \
      -r subsampled/${sample}.sub.fq.gz \
      -o quant --threads ${task.cpus} ${validate}
    """

  """
  set -euo pipefail
  mkdir -p subsampled quant

  ${subsampleCmd}
  ${salmonCmd}

  cp -f quant/lib_format_counts.json      ${sample}.lib_format_counts.json || true
  cp -f quant/aux_info/meta_info.json     ${sample}.meta_info.json

  if command -v jq >/dev/null 2>&1; then
    jq -r '
      ["sample","inferred_library_type","num_mapped","percent_mapped"],
      ["${sample}", (.library_type // "NA"), (.num_mapped // 0), (.percent_mapped // 0)]
      | @tsv
    ' quant/aux_info/meta_info.json > ${sample}.salmon_infer.tsv || true
  fi
  """
}

/* 2) Initial FastQC quality scores */
process FASTQC_RAW {
    tag "${base}"
    publishDir "${params.outdir}/qc/raw_fastqc", mode: 'copy'

    input:
    tuple val(base), path(r1), path(r2)

	output:
	path "*"

    script:
    """
    set -euo pipefail
    fastqc -t ${params.threads} -o . ${r1} ${r2}
    """
}

/* Extract UMI barcodes - PE only, no cell barcodes, only UMI barcodes */
process UMI_TOOLS_EXTRACT_PE {
  tag "$sample"
  errorStrategy 'terminate'
  publishDir "${params.outdir}/umi_ext", pattern: "*.log", mode: 'copy'

  /*
   * INPUTS
   *  - sample : sample id (string)
   *  - r1     : FASTQ R1 with UMI prefix
   *  - r2     : FASTQ R2 (insert only)
   *
   * PARAMS (set in nextflow.config or via CLI)
   *  - params.umi_length : length of UMI (default 8)
   *  - params.umi_separator : separator added to read name (default ":")
   */

  input:
  tuple val(sample), path(r1), path(r2)

  output:
  tuple val(sample), 
        path("${sample}_R1.extracted.fq.gz"), 
		path("${sample}_R2.extracted.fq.gz"), 
		emit: UMI_EXT_FQ
  path "${sample}.umi_extract.log",  emit: LOG

  script:
  def umiLen = params.umi_length ?: 8
  def bcPattern = "N" * umiLen                  // e.g. "NNNNNNNN"
  def sep = params.umi_separator ?: ":"

  """
  set -euo pipefail

  umi_tools extract \
    --bc-pattern='${bcPattern}' \
    --stdin=${r1} \
    --read2-in=${r2} \
    --stdout=${sample}_R1.extracted.fq.gz \
    --read2-out=${sample}_R2.extracted.fq.gz \
    --umi-separator='${sep}' \
    --log=${sample}.umi_extract.log
  """
}

/* 3) Cleaning with fastp */
process FASTP_CLEAN_PE {
  tag "${base}"

  // Publish cleaned FASTQs to .../clean, QC to .../qc/fastp
  publishDir "${params.outdir}/clean",    mode: 'copy', pattern: '*_R{1,2}.clean.fastq.gz'
  publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: '*.fastp.*'

  input:
  tuple val(base), path(r1), path(r2)

  output:
  //cleaned read pair (keyed by base)
  tuple val(base),
		path("${base}_R1.clean.fastq.gz"),
        path("${base}_R2.clean.fastq.gz"), emit: files

  // QC artifacts (each keyed by base so you can join later if needed)
  	tuple val(base), path("${base}.fastp.html"), emit: html
  	tuple val(base), path("${base}.fastp.json"), emit: json

  script:
  // Optional flags assembled via Groovy ternary
  def trimPoly = params.trim_poly   ? '--trim_poly_g --trim_poly_x' : ''
  def corr     = params.correction  ? '--correction'                : ''
  """
  set -euo pipefail

  fastp \
    -i ${r1} -I ${r2} \
    -o ${base}_R1.clean.fastq.gz -O ${base}_R2.clean.fastq.gz \
    --detect_adapter_for_pe \
    --cut_right --cut_right_mean_quality ${params.q_right} \
    ${trimPoly} \
    -n ${params.max_ns} -l ${params.min_len} \
    ${corr} \
    --thread ${params.threads} \
    --html ${base}.fastp.html \
    --json ${base}.fastp.json
  """
}

/* Remove contaiminating RNAs like rRNA (PE only), SortMeRNA ≥4.x */
process SORTMERNA_FILTER_PE {
  tag "${sample}"

  publishDir "${params.outdir}/sortmerna", mode: 'copy'

  input:
  tuple val(sample), path(read1), path(read2)
  path workdir_with_idx
  path (rrna_fasta)

  output:
  tuple val(sample),
        path("${sample}.rClean_R1.fastq.gz"), 
        path("${sample}.rClean_R2.fastq.gz"), 
		emit: non_rrna_reads
  tuple val(sample),
        path("${sample}.rRNA_R1.fastq.gz"), 
        path("${sample}.rRNA_R2.fastq.gz"), 
		emit: rRNA_reads
  path("*.log"), emit: log

  script:
  """
  set -euo pipefail

  # prep workdir + optional cached index
  mkdir -p ./work/idx
  if [ -d ./idx ] && find ./idx -mindepth 1 -maxdepth 1 -print -quit >/dev/null; then
    cp -a ./idx/. ./work/idx/
  fi

  # run sortmerna
  sortmerna \
  --ref "${rrna_fasta}" \
  --reads "${read1}" --reads "${read2}" \
  --workdir ./work \
  --aligned "${sample}.rRNA" \
  --other   "${sample}.non_rRNA" \
  --fastx --paired_in --out2 --num_alignments 1 \
  2> "${sample}.sortmerna.log"

  # Fix outfile names
  mv ${sample}.non_rRNA_fwd.fq.gz ./${sample}.rClean_R1.fastq.gz
  mv ${sample}.non_rRNA_rev.fq.gz ./${sample}.rClean_R2.fastq.gz
  mv ${sample}.rRNA_fwd.fq.gz ./${sample}.rRNA_R1.fastq.gz
  mv ${sample}.rRNA_rev.fq.gz ./${sample}.rRNA_R2.fastq.gz
  """
}

/* 4) FastQC quality scores for cleaned fastq files */
process FASTQC_CLEAN_PE {
    tag "${base}"
    publishDir "${params.outdir}/qc/clean_fastqc", mode: 'copy'

    input:
    tuple val(base), path(r1c), path(r2c)

	output:
	path "*.html"
	path "*.zip"

    script:
    """
    set -euo pipefail
    fastqc -t ${params.threads} -o . ${r1c} ${r2c}
    """
}

workflow {
	INTEGRITY_STATS( READS_PE )
	SUBSAMPLE_AND_INFER_STRAND( READS_PE, params.paired_end, params.salmonIdx )
	FASTQC_RAW( READS_PE )
	UMI_EXT_OUT = UMI_TOOLS_EXTRACT_PE( READS_PE )
	CLEANED = FASTP_CLEAN_PE( UMI_EXT_OUT.UMI_EXT_FQ ) // Clean the UMI-extracted fastqs
	SMR = SORTMERNA_FILTER_PE( CLEANED.files, IDX_CACHE_CH, params.riboFasta )
	FASTQC_CLEAN_PE( SMR.non_rrna_reads )
    }