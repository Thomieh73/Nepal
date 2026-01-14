// this module runs the latest version of DORADO
// Process DORADO is used for basecalling the POD5 raw data files.


process DORADO_SIMPLEX {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/dorado_gpu_1.3.0"
	
	//publishDir "${params.out_dir}/dorado_simplex/", pattern: "*", mode: "copy"
	//publishDir "${params.out_dir}/01_guppy/", pattern: "fastq", mode: "copy"
	//publishDir "${params.out_dir}/01_dorado_simplex/", pattern: "sequencing_logs/sequencing_*.*", mode: "copy"

	label 'gpu_A100'

	input:
	file("*")

	output:
	path "basecalls.bam", emit: simplex_ch
	path "sequencing_summary.txt", emit: summary_ch

	script:
	"""
	dorado basecaller -x "cuda:all" --min-qscore 7 --no-trim $params.dorado.moddir/$params.dorado.model pod5 > basecalls.bam

	# creating a summary file for calculating stats
	dorado summary basecalls.bam > sequencing_summary.txt

	"""

}

// This process does demultiplexing with dorado demux.
// It takes a samplesheet and the information for the kit that was used.
// it also does trimming of the barcode. 

process DORADO_DEMUX {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/dorado_gpu_1.3.0"
	
	//publishDir "${params.out_dir}/dorado_demux/", pattern: "demultiplexed/*", mode: "copy"
	//publishDir "${params.out_dir}/01_guppy/", pattern: "fastq", mode: "copy"
	//publishDir "${params.out_dir}/01_dorado_simplex/", pattern: "sequencing_logs/sequencing_*.*", mode: "copy"
	publishDir "${params.out_dir}/data_overview/", pattern: "demultiplexed/barcoding_summary.txt", mode: "copy"

	label 'heavy'

	input:
	file(x)

	output:
	// Captures every BAM in the nested structure
	path "demultiplexed/**/*.bam", emit: demux_ch

	script:
	"""
    dorado demux --kit-name $params.dorado.barcode --sample-sheet $params.samplesheet.location --emit-summary --output-dir demultiplexed $x

	"""

}


// This process takes each demultiplex bam file and then determines which reads can be combined
// into duplex reads.
// It only takes the bam file as input

process DORADO_DUPLEX {
	conda "/cluster/projects/nn9305k/src/miniconda/envs/dorado_gpu_1.3.0"
	
	//publishDir "${params.out_dir}/dorado_duplex/", pattern: "*", mode: "copy"
	//publishDir "${params.out_dir}/01_guppy/", pattern: "fastq", mode: "copy"
	//publishDir "${params.out_dir}/01_dorado_simplex/", pattern: "sequencing_logs/sequencing_*.*", mode: "copy"

	label 'gpu_A100'

	input:
	file("*")
	

	output:
	//tuple val(samplename), path('*.bam'), emit: duplex_ch
	path('*.bam'), emit: duplex_ch

	script:
	//samplename = x.toString() - ~/.read_ids.txt$/
	"""
	## creating a variable from the name of the text file
	SAMPLENAME=(\$(ls *.txt|rev|cut -c 14- | rev))

	echo \$SAMPLENAME

	dorado duplex -x "cuda:all" --min-qscore 7 -t 32 $params.dorado.moddir/$params.dorado.model --read-ids \$SAMPLENAME.read_ids.txt pod5 > \$SAMPLENAME.bam

	"""

}


// the bellow process is to rename the bam files after the demuxing step.
// by doing this we get bam files with the same sample name as the alias in the samplesheet.

process MERGE_BAMS {
    tag "$alias"

    input:
    tuple val(alias), path(chunks)

    output:
    path "${alias}.bam", emit: merged_bam

    script:
    if (chunks instanceof List && chunks.size() > 1)
        """
        samtools merge ${alias}.bam $chunks
        """
    else
        """
        mv $chunks ${alias}.bam
        """
}
