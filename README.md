![Banner](./pix/banner.png)
# Nextflow Fastq Quality Control Pipelines
<i> Fastq metrics and trimming for paired-end reads</i>

## Available Pipelines: (2)
<i>Separate pipelines are available for different circumstances:</i><br>
#### <span style="color:skyblue;">General purpose: </span> fastqQC<br> 
- Trims on quality and adaptors 
- Collects pre/post cleaning metrics

#### <span style="color:skyblue;">RNA-Seq reads w/ UMI barcodes: </span>fastqQC_PE_UMI
- Trims on quality and adaptors
- Determines library strandedness
- Moves UMIs to read headers
- Removes reads that map to rRNA
- Collects pre and post-cleaning QC metrics

## Running the pipelines:
- Download the scripts and nextflow.config files
- Configure your input and output folders as parameters
- Script can be run in <b>docker</b> by setting nextflow.config file's docker setting to "true". That will execute the pipeline in container cbreuer/fastqc_preprocess:latest which has the required programs for all nextflow scripts.
- Running locally without Docker - required command line tools:
	- python 3
	- Java 17
	- seqkit
	- fastp
	- FastQC
	- seqtk
	- AWS CLI (optional)

## Processing Steps - <i>FastQC_PE</i>:
1. <b>Fastq statistics</b> - seqkit stats is used to collect basic metrics about the fastq files including the number of sequences, minimum and maximum sequence lengths, average length, sum of lengths, and the format and type of sequences in the file. 
Process: "INTEGRITY_STATS".<br>
2. <b>Fastq quality control (pre-cleaning)</b> - FastQC is used to collect metrics on per-base quality scores. 
Process "FASTQC_RAW"
3. <b>Fastq cleaning</b> - fastp is used to trim on adaptors and quality. 
Process: "FASTP_CLEAN_PE".
4. <b>Fastq quality control (post-cleaning)</b> - FastQC is used again to collect metrics on cleaned fastq files.
Process: "FASTQC_CLEAN_PE". 

## Processing Steps - <i>FastQC_PE_UMI</i>:
1. <b>Fastq statistics</b> - seqkit stats is used to collect basic metrics about the fastq files including the number of sequences, minimum and maximum sequence lengths, average length, sum of lengths, and the format and type of sequences in the file. 
Process: "INTEGRITY_STATS".<br>
2. <b>Determine Strandedness</b> - salmon determines the library strandedness 
Process: "SUBSAMPLE_AND_INFER_STRAND"
3. <b>Fastq quality control (pre-cleaning)</b> - FastQC is used to collect metrics on per-base quality scores. 
Process "FASTQC_RAW"
4. <b>UMI Extraction</b> - UMI-tools moves UMIs to read headers 
Process "UMI_TOOLS_EXTRACT_PE"
5. <b>Fastq cleaning</b> - fastp is used to trim on adaptors and quality. 
Process: "FASTP_CLEAN_PE".
6. <b>rRNA Read Removal</b> - SortMeRNA removes rRNA reads (user adjusts the target sequence)
Process: "SORTMERNA_FILTER_PE".
7. <b>Fastq quality control (post-cleaning)</b> - FastQC is used again
Process "FASTQC_CLEAN_PE"

## Quality Metrics
- seqkit outputs general statistics
 ![seqkit](./pix/stats.jpg)

 - fastp trimming results
 ![fastqp](./pix/fastp.jpg)

 - fastp fragment length analysis
 ![fastp fragment analysis](./pix/fastpFrag.jpg)

- FastQC per-base PHRED scores (pre-cleaning)
![fastqc phred pre](./pix/fastqcRaw.jpg)

- FastQC per-base PHRED scores (post-cleaning)
![fastqc phred post](./pix/fastqcClean.jpg)

# Notes:
- My dockerfile for container cbreuer/fastq_quality_control is included if needed.