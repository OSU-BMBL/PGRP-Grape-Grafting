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
ADAPTER_PATH = which("bbduk.sh")
try:
    BBDUK_OPTIONS = config["bbduk"]["options"]
except:
    raise ValueError("bbduk > options not found in the configuration file")
BBDUK_VERSION = config["bbduk"]["version"]

###### Multithread configuration #####
CPUS_FASTQC = 4
CPUS_TRIMMING = 5
CPUS_STAR = 20

####### Output directories #######
REF_GENOME = "GENOME/"
LOGS = "0.LOGS/"
RAW_FASTQC = "1.QC.RAW/"
TRIMMED_READS = "2.TRIMMED/"
TRIMMED_READS_FASTQC = "3.QC.TRIMMED/"
ALIGNMENT = "4.ALIGNMENT/"
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
            raw_reads = LIBS, raw_ends = RAW_ENDS, format = ["html","zip"])
        # expand(RAW_FASTQC + "{raw_reads}{raw_ends}_fastqc.{format}",
        #     raw_reads = LIBS, raw_ends = [1, 2], format = ["html","zip"])
    output:
        expand(REPORTS + "Report_{step}.html", step = ["FastQC_Raw", "FastQC_Trimmed"])
    params:
        logs 	= directory(LOGS),
        reports	= directory(REPORTS)
    run:
        shell("multiqc -f -o {params.reports} -n Report_FastQC_Raw.html -d " + RAW_FASTQC)
        shell("multiqc -f -o {params.reports} -n Report_FastQC_Trimmed.html -d " + TRIMMED_READS_FASTQC)

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
        adapter = os.path.join(ADAPTER_PATH, "../opt/bbmap-" + BBDUK_VERSION + "/resources/"),
        reads = READS_PATH + "{raw_reads}{raw_ends}." + EXTENSION
    output:
        TRIMMED_READS + "{raw_reads}{raw_ends}." + EXTENSION
    params:
        options = BBDUK_OPTIONS
    log:
        TRIMMED_READS + "{raw_reads}{raw_ends}.log"
    message:
        "Using Single End Trimming"
    threads:
        CPUS_TRIMMING
    shell:
        "bbduk.sh threads={threads} in={input.reads} out={output} {params.options} 2> {log}"

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
		genome_files = expand(REF_GENOME + "{genome_file}", genome_file = GENOME_FILENAMES)
	output:
		dir = directory(REF_GENOME + "GENOME_INDEX")
	message:
		"Generate genome index for STAR"
	log:
		REF_GENOME + "genome_index.log"
	threads:
		CPUS_STAR
	shell:
		"mkdir -p {output.dir} && STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.genome_files[0]}  --sjdbGTFfile {input.genome_files[1]} --sjdbOverhang 50 2> {log}"
