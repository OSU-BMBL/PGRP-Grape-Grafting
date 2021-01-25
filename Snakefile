####### Libraries #######
from utils import findLibraries, loadGenome, verifyGenome, which

####### Global variables #######
EXTENSION = config["reads"]["extension"]
PREFIX = config["reads"]["prefix"]
READS_PATH = config["reads"]["path"]
FORWARD_READ_ID = config["reads"]["forward_read_id"]
SUFFIX = "_" + FORWARD_READ_ID + "." + EXTENSION
LIBS = findLibraries(READS_PATH,PREFIX,SUFFIX)
RAW_ENDS = "_R1"
ADAPTER_PATH = which("trimmomatic")
if (os.environ.get('ADAPTERS', '') == ''):
    os.environ['ADAPTERS'] = os.path.join(ADAPTER_PATH,"../share/trimmomatic/adapters")
try:
    TRIMMOMATIC_OPTIONS = config["trimmomatic"]["options"]
except:
    raise ValueError("trimmomatic > options not found in the configuration file")
# BBDUK_VERSION = config["bbduk"]["version"]

###### Multithread configuration #####
CPUS_FASTQC = 4
CPUS_TRIMMING = 5
CPUS_BWA = 20
CPUS_READCOUNTS = 5

####### Output directories #######
REF_GENOME = "GENOME/"
LOGS = "0.LOGS/"
RAW_FASTQC = "1.QC.RAW/"
TRIMMED_READS = "2.TRIMMED/"
TRIMMED_READS_FASTQC = "3.QC.TRIMMED/"
ALIGNMENT = "4.ALIGNMENT/"
ALIGNMENT_QC = "5.QC.ALIGNMENT/"
REPORTS = "999.REPORTS/"

####### Reference datasets #######
FA,GTF = loadGenome(config["genome"])
GENOME_FILENAMES = {"FA":FA,"GTF":GTF}
verifyGenome(config["genome"],REF_GENOME + FA, REF_GENOME + GTF)

####### Rules #######
rule all:
    input:
        expand(RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"]),
        expand(TRIMMED_READS_FASTQC + "{raw_reads}{raw_ends}_fastqc.{format}",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"]),
        expand(ALIGNMENT + "{raw_reads}{raw_ends}_sorted.bam",
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"]),
        expand(ALIGNMENT_QC + "{raw_reads}{raw_ends}.flagstat_BWA.txt",
            raw_reads = LIBS, raw_ends = RAW_ENDS),
        "readCounts.txt"
    output:
        expand(REPORTS + "Report_{step}.html", step = ["FastQC_Raw", "FastQC_Trimmed"])
    params:
        logs 	= directory(LOGS),
        reports	= directory(REPORTS)
    run:
        shell("multiqc -f -o {params.reports} -n Report_FastQC_Raw.html -d " + RAW_FASTQC)
        shell("multiqc -f -o {params.reports} -n Report_FastQC_Trimmed.html -d " + TRIMMED_READS_FASTQC)
        shell("multiqc -f -o {params.reports} -n Report_Alignment_BAM.html -d " + ALIGNMENT)
        shell("multiqc -f -o {params.reports} -n Report_Alignment_flagstat.html -d " + ALIGNMENT_QC)

rule fastqc_raw:
    input:
        reads = READS_PATH + "{raw_reads}{raw_ends}." + EXTENSION
    output:
        html = RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.html",
        zip  = RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.zip"
    message:
        "FastQC on raw data"
    log:
        RAW_FASTQC + "{raw_reads}{raw_ends}.log"
    threads:
        CPUS_FASTQC
    shell:
        "fastqc -o " + RAW_FASTQC + " -t {threads} {input.reads} 2> {log}"

rule trim_reads:
    input:
        # adapter = os.path.join(ADAPTER_PATH,"../share/trimmomatic/adapters"),
        # adapter = os.path.join(ADAPTER_PATH, "../opt/bbmap-" + TRIMMOMATIC_VERSION + "/resources/"),
        reads = READS_PATH + "{raw_reads}{raw_ends}." + EXTENSION
    output:
        TRIMMED_READS + "{raw_reads}{raw_ends}." + EXTENSION
    params:
        options = TRIMMOMATIC_OPTIONS
    log:
        err = TRIMMED_READS + "{raw_reads}{raw_ends}.err",
        txt = TRIMMED_READS + "{raw_reads}{raw_ends}.txt"
    message:
        "Using Single End Trimming"
    threads:
        CPUS_TRIMMING
    shell:
        "trimmomatic SE -threads {threads} {input.reads} {output} {params.options} -trimlog {log.txt} 2> {log.err}"
        # "bbduk.sh threads={threads} in={input.reads} out={output} {params.options} 2> {log}"

rule fastqc_trimmed:
    input:
        reads = rules.trim_reads.output
    output:
        html = TRIMMED_READS_FASTQC + "{raw_reads}{raw_ends}_fastqc.html",
        zip  = TRIMMED_READS_FASTQC + "{raw_reads}{raw_ends}_fastqc.zip"
    message:
        "FastQC on trimmed data"
    log:
        TRIMMED_READS_FASTQC + "{raw_reads}{raw_ends}.log"
    threads:
        CPUS_FASTQC
    shell:
        "fastqc -o " + TRIMMED_READS_FASTQC + " -t {threads} {input.reads} 2> {log}"

rule genome_index:
	input:
		genome_files = expand(REF_GENOME + "{genome_file}", genome_file = GENOME_FILENAMES.values())
	output:
		dir = directory(REF_GENOME + "GENOME_INDEX")
	params:
		genome_files = expand("{genome_file}", genome_file = GENOME_FILENAMES.values())
	message:
		"Generate genome index for BWA"
	log:
		REF_GENOME + "genome_index.log"
	# threads:
	# 	CPUS_BWA
	shell:
        	"mkdir -p {output.dir} && ln -sf ../{params.genome_files[0]} {output.dir} && bwa index {output.dir}/{params.genome_files[0]} 2> {log}"
		# "mkdir -p {output.dir} && BWA --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome_files[0]}  --sjdbGTFfile {input.genome_files[1]} --sjdbOverhang 50 2> {log}"

rule alignment:
	input:
		genome = rules.genome_index.output.dir,
		reads = rules.trim_reads.output
	output:
		# unmapped_m81 = ALIGNMENT + "{raw_reads}{raw_ends}_Unmapped.out.mate1",
		# unmapped_m82 = ALIGNMENT + "{raw_reads}{raw_ends}_Unmapped.out.mate2",
		# aligned_bam  = ALIGNMENT + "{raw_reads}{raw_ends}_Aligned.sortedByCoord.out.bam"
		ALIGNMENT + "{raw_reads}{raw_ends}_sorted.sam"
	params:
		ref = rules.genome_index.params.genome_files[0]
	message:
		"BWA alignment"
	log:
		ALIGNMENT + "{raw_reads}{raw_ends}_sam.log"
	params:
		prefix = ALIGNMENT + "{raw_reads}{raw_ends}_"
	threads:
		CPUS_BWA
	shell:
        	"bwa mem -t {threads} {input.genome}/{params.ref} {input.reads} -o {output} 2> {log}"
		# "STAR --runThreadN {threads} --genomeDir {input.genome} --readFilesIn {input.reads} --readFilesCommand gunzip -c --outFilterIntronMotifs RemoveNoncanonical --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped  Fastx 2> {log}"

rule sam2bam:
    input:
        rules.alignment.output
    output:
        ALIGNMENT + "{raw_reads}{raw_ends}_sorted.bam"
    log:
        ALIGNMENT + "{raw_reads}{raw_ends}_bam.log"
    message:
        "Converting SAM to BAM"
    threads:
        CPUS_BWA
    shell:
        "samtools view -@ {threads} -bS {input} | samtools sort -@ {threads} -o {output} 2> {log}"

rule alignment_quality:
    input:
        rules.alignment.output
    output:
        ALIGNMENT_QC + "{raw_reads}{raw_ends}.flagstat_BWA.txt"
    message:
        "Assessing alignment quality"
    shell:
        "samtools flagstat {input} > {output}"
        # "SAMstats --sorted_sam_file {input} --outf {output} 2> {log}"

rule read_counts:
	input:
		aligned = expand(rules.sam2bam.output, raw_reads = LIBS, raw_ends = RAW_ENDS),
		genome = rules.genome_index.input.genome_files[1]
	output:
		readCounts = "readCounts.txt"
	log:
		"read_counts.log"
	threads:
		CPUS_READCOUNTS
	shell:
		"featureCounts -a {input.genome} -o {output} -T {threads} {input.aligned} 2> {log}"
