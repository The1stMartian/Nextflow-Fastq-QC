/* --------------------------
   Parameters
--------------------------- */
params.input      = '/your/input/path/*_{R1,R2}.fastq.gz'  	// paired-end pattern
params.single_end = false                                  	// set true for SE
params.threads    = 8
params.min_len    = 30
params.q_right    = 20
params.max_ns     = 5
params.trim_poly  = true     								// poly-G/X trimming for NovaSeq/NextSeq
params.correction = true     								// PE overlap correction
params.outdir = "/your/output/path" 						// folder must already exist


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

    publishDir "${params.outdir}/qc", mode: 'copy'

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

/* 2) Initial FastQC quality scores */
process FASTQC_RAW {
    tag "${base}"
    publishDir "${params.outdir}/qc/raw_fastqc", mode: 'copy'

    input:
    tuple val(base), path(r1), path(r2)

	output:
	path "*"

    when:
    !params.single_end

    script:
    """
    set -euo pipefail
    fastqc -t ${params.threads} -o . ${r1} ${r2}
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

  when:
  !params.single_end

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

/* 4) FastQC quality scores for cleaned fastq files */
process FASTQC_CLEAN_PE {
    tag "${base}"
    publishDir "${params.outdir}/qc/clean_fastqc", mode: 'copy'

    input:
    tuple val(base), path(r1c), path(r2c)

	output:
	path "*.html"
	path "*.zip"

    when:
    !params.single_end

    script:
    """
    set -euo pipefail
    fastqc -t ${params.threads} -o . ${r1c} ${r2c}
    """
}

workflow {
	INTEGRITY_STATS(READS_PE)
	FASTQC_RAW(READS_PE)
	CLEANED = FASTP_CLEAN_PE(READS_PE)
	CLEANED.files.view()
	FASTQC_CLEAN_PE( CLEANED.files )
    }