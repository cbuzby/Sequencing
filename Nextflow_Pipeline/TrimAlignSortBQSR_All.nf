/*  GATK4 Variant Calling Pipeline
 *  Usage: nextflow run /path/to/main.nf
 *
 *  Author: Mohammed Khalfan < mkhalfan@nyu.edu > adapted by Cassandra Buzby
 *  NYU Center for Genetics and System Biology 2020
 */

// Setting some defaults here,
// can be overridden in config or via command line

nextflow.enable.dsl=2
params.out = "${params.outdir}/out"
params.tmpdir = "${params.outdir}/gatk_temp"

// Define modules here
BWA = 'bwa/intel/0.7.17'
PICARD = 'picard/2.17.11'
GATK = 'gatk/4.1.9.0'
R = 'r/intel/4.0.3'
SAMTOOLS = 'samtools/intel/1.14/'
BAMTOOLS = 'bamtools/intel/2.5.1'
BCFTOOLS = 'bcftools/intel/1.14'
MULTIQC = 'multiqc/1.9'
TRIMMOMATIC = 'trimmomatic/0.36'
PICARD_JAR = '/share/apps/picard/2.17.11/picard.jar'

// Define Reference Here
REF = '/scratch/cb4097/Sequencing/Reference/GCF_000146045.2_R64_genomic.fna'
KNOWN_SITES = '/scratch/cb4097/Sequencing/ParentSequences/OakWine_variants.vcf.gz'

// Alternative to setup the reference file
ref = file(params.ref)

// Print some stuff here
println "reads: $params.reads"
println "ref: $params.ref"
println "output: $params.out"
println "gatk temp dir: $params.tmpdir"


process trim {
// comments
        publishDir "${params.out}/trimmed", mode:'copy'

        input:
        tuple val(id),
                file(reads)

        output:
        tuple val(id),
                file(trimmed_file)

        script:
        trimmed_file = "${id}.trimmed.fq.gz"
        """
        module load $TRIMMOMATIC
        module load $GATK
        module load $BWA
        module load $PICARD

        java -jar /share/apps/trimmomatic/0.36/trimmomatic-0.36.jar SE -phred33 \
               $reads $trimmed_file ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

        """
}


process align {
    publishDir "${params.out}/aligned_reads", mode:'copy'

    input:
    tuple val(id),
        file(read_1)

    output:
    tuple val(id), file("${id}.aln-se.sam")

    script:
    readGroup = \
        "@RG\\tID:${id}\\tLB:${id}\\tPL:ILLUMINA\\tPM:NOVASEQ\\tSM:${id}"
    """
    module load $BWA
    bwa mem -Y -K 100000000 -R \"${readGroup}\" \
               $REF $read_1 > ${id}.aln-se.sam
    """
}

process sort {
        // comments
        publishDir "${params.out}/sorted", mode:'copy'

        input:
        tuple val(id),
            file(read_1)

        output:
        tuple val(id), file("${id}.aln-se.bam") 

        script:

        """
        module load $GATK
        module load $BWA
        module load $PICARD

        java -jar $PICARD_JAR SortSam \
                  INPUT=$read_1 \
                  OUTPUT=${id}.aln-se.bam \
                  SORT_ORDER=coordinate
        """
}

process MakeBQSRTable {
	// comments 	
	publishDir "${params.out}/bqsr_table", mode:'copy'
	
	input:
	tuple val(id),
	    file(read) 	

	output:
	tuple val(id), file(read), file("${id}_recal_data.table")

	script:
	"""
	
	module load $GATK
	module load $BWA
	module load $PICARD
	
	gatk BaseRecalibrator -R $REF -I $read --known-sites $KNOWN_SITES -O ${id}_recal_data.table
	
	"""
}

process ApplyBQSR {
	// comments 

	publishDir "${params.out}/bqsr_bams", mode:'copy'

	input:
	tuple val(id),
	    file(bamfile),
		file(table)

	output:
	tuple file("${id}.bqsr.output.bam")

	script:
	"""
	
	module load $GATK
	module load $BWA
	module load $PICARD
	module load $SAMTOOLS 

	gatk ApplyBQSR \
	   -R $REF \
	   -I ${bamfile} \
	   --bqsr-recal-file ${table} \
	   -O ${id}.bqsr.output.bam
	"""
}

process MergeSplit {

	publishDir "${params.out}/mergedsplit", mode:'copy'

	input:
	path(bamfiles) // NEEDS EDITING?

	output:
	path("${params.fcid}*.bam")

	script:
	"""
	module $GATK
	module load $BWA
	module load $PICARD
	module load $SAMTOOLS 
	module load $BAMTOOLS

	bamstring = ${bamfiles}.join(" ")

	samtools merge ${params.fcid}.bam ${bamstring}
	samtools index ${params.fcid}.bam
	bamtools split -in ${params.fcid}.bam -reference

	"""
}


workflow {
    main:
        files_ch = Channel.fromPath( "${params.reads}" ).map{ [it.getBaseName(), it ] }
        files_ch.view()
	
        trim(files_ch) | set{ trimmed_ch }
	align(trimmed_ch) | set{ aligned_ch }

	// aligned_ch = Channel.fromPath( "${params.out}/aligned_reads/*" ).map{ [it.getBaseName(), it ]}
	// aligned_ch.view()
	
	sort(aligned_ch) | set{ sorted_ch }

	// sorted_ch = Channel.fromPath( "${params.out}/sorted/*" ).map{ [it.getBaseName(), it ]}
	// sorted_ch.view()

	MakeBQSRTable(sorted_ch) | set{ bqsrtables_ch }
	ApplyBQSR(bqsrtables_ch) | set{ bqsr_bam_ch }

	// collect( bqsr_bam_ch ) | MergeSplit | set{ splitchr_ch }
	// splitchr_ch.view()
}

