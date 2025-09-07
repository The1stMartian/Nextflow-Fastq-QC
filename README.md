![Banner](./pix/banner.jpg)
# Nextflow Fastq Pre-Processing Pipelines
<i> Fastq metrics and trimming for paired-end reads</i>

## Available Pipelines: (2)

#### 1) FastqQC (PE)
- Trims on quality and adaptors
- Determines library strandedness
- Removes reads that map to rRNA
- Collects pre/post-cleaning QC metrics

#### 2) FastqQC (PE + UMI)
- Trims on quality and adaptors
- Determines library strandedness
- <b>Moves UMIs</b> to read headers
- Removes reads that map to rRNA
- Collects pre/post-cleaning QC metrics

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
- <b>Downstream:</b> Proceed to my [RNA-Seq analysis](https://github.com/The1stMartian/Nextflow-RNA-Seq) pipelines

## Processing Steps - <i>FastQC-PE</i>:
1. <b>Fastq statistics</b> - seqkit stats is used to collect basic metrics about the fastq files including the number of sequences, minimum and maximum sequence lengths, average length, sum of lengths, and the format and type of sequences in the file. 
2. <b>Determine Strandedness</b> - salmon determines the library strandedness 
3. <b>QC (pre-cleaning)</b> - FastQC is used to collect metrics on per-base quality scores. 
4. <b>Fastq cleaning</b> - fastp is used to trim on adaptors and quality. 
5. <b>rRNA Read Removal</b> - SortMeRNA removes rRNA reads (user adjusts the target sequence)
6. <b>QC (post-cleaning)</b> - FastQC is used again to collect metrics on cleaned fastq files. 

## Processing Steps - <i>FastQC-PE-UMI</i>:
1. <b>Fastq statistics</b> - seqkit stats is used to collect basic metrics about the fastq files including the number of sequences, minimum and maximum sequence lengths, average length, sum of lengths, and the format and type of sequences in the file. 
Process: "INTEGRITY_STATS".<br>
2. <b>Determine Strandedness</b> - salmon determines the library strandedness 
3. <b>QC (pre-cleaning)</b> - FastQC is used to collect metrics on per-base quality scores. 
4. <b>UMI Extraction</b> - UMI-tools moves UMIs to read headers 
5. <b>Fastq cleaning</b> - fastp is used to trim on adaptors and quality. 
6. <b>rRNA Read Removal</b> - SortMeRNA removes rRNA reads (users can adjusts the rRNA sequence to eukaryotic or prokaryotic sequences)
7. <b>QC (post-cleaning)</b> - FastQC is used again

# Understanding Strandedness
- Salmon detects the library type and outputs a analysis that can be hard to understand - read the <b>/strandedness/meta_info.json</b> for the detected library type information. <br>
- The libary types are described using the following terms. e.g "IU" would be paired-end unstranded.
- All combinations are listed in the output file with the number of reads fitting the description beside it. Typically, one code will predominate, allowing users to make an informed decision about the library type. This information can be fed into downstream mapping pipelines.

|Layout Term|Meaning|
|:-----|:----|
|S|single-end| 
|I|paired-end inward (typical)|
|O|paired-end outward (atypical)|
|M|paired-end matching|
|U|unstranded|
|S |stranded|
|F|Forward|
|R|Reverse|

## Quality Metrics
- seqkit outputs general statistics
 ![seqkit](./pix/stats.jpg)

 - fastp trimming results
 ![fastqp](./pix/fastp.jpg)

 - fastp fragment length analysis
 ![fastp fragment analysis](./pix/fastpFrag.jpg)

- FastQC per-base PHRED scores (pre-cleaning)
![fastqc phred pre](./pix/quality.png)


# Background work:
- My dockerfile for container cbreuer/fastq_quality_control is included if needed.- My dockerfile for container cbreuer/fastq_quality_control is included if needed.<br>

Create a sortmerna index:
- Open rrna-db-defaults.txt and [download](https://www.arb-silva.de/download/arb-files/) the silva database containing the rRNA sequences appropriate for your library. Keeping the number of sequences to a minimum is advised because the processing time of sortmerna can be high. Combine the desired sequences into a single .fa (your 'ribo_fasta.fa' sequence) file and run a one sortmerna example analysis with dummy fastq files. This will cause sortmerna to create the index in the work folder/idx. It's the /idx folder you want to keep and reuse. Point the nextflow paramter.
- Save the .fa file with the rRNA sequences used to make the sortmerna library. Point the nextflow paramter at the file.<br>
<br>
Create a salmon index:
- A tutorial is [here](https://combine-lab.github.io/salmon/getting_started/#indexing-txome)). 
- Briefly: use "salmon index -t genome.fa.gz -i salmon_gnm_name" and point the nextflow parameter to the index folder.