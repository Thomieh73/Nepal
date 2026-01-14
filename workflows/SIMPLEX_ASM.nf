// modules for DUPLEX_ASM workflow
	
    include { DORADO_DEMUX } from "../modules/DORADO.nf"
    include { DORADO_SIMPLEX } from "../modules/DORADO.nf"
    include { NANOFILT_SIMPLEX } from "../modules/NANOFILT.nf"
  	include { NANOPLOT_SIMPLEX } from "../modules/NANOPLOT.nf"
    include { NANOPLOT_FASTQ } from "../modules/NANOPLOT.nf"
    include { PYCOQC_SIMPLEX } from "../modules/PYCOQC.nf"
    include { SAMTOOLS_BAM2FQ } from "../modules/SAMTOOLS.nf"
    include { SEQKIT_SIMPLEX } from "../modules/SEQKIT.nf"
    include { SEQKIT_NFILT } from "../modules/SEQKIT.nf"
    include { SEQKIT_FLYE } from "../modules/SEQKIT.nf"
    include { FLYE_ASM } from "../modules/FLYE.nf"
	include { MERGE_BAMS } from "../modules/DORADO.nf"


// workflows

workflow SIMPLEX_ASM {
	pod5_ch=channel.fromPath(params.reads, checkIfExists: true)
                        .collect()

    // process reads to get demultiplexed reads IDs
	DORADO_SIMPLEX(pod5_ch)   
    DORADO_DEMUX(DORADO_SIMPLEX.out.simplex_ch.flatten())

	// 3. CHANNEL TRANSFORMATION: Grouping by the parent folder (the alias)
    ch_grouped_bams = DORADO_DEMUX.out.demux_ch
        .flatten()
        .map { file -> 
            def alias = file.parent.name 
            return [ alias, file ] 
        }
        .filter { alias, file -> alias != 'unclassified' } // Optional: ignore unclassified
        .groupTuple()

    // 4. Merging chunks into sample-named BAMs (e.g., my_first_sample.bam)
    MERGE_BAMS(ch_grouped_bams)

    // 5. Convert to FastQ (now using the clean merged BAMs)
    // We no longer need .flatten() here because MERGE_BAMS emits one file at a time
    SAMTOOLS_BAM2FQ(MERGE_BAMS.out.merged_bam)

    // filtering the reads to remove poor reads
    NANOFILT_SIMPLEX(SAMTOOLS_BAM2FQ.out.filter_ch.flatten())

    // Doing an assembly with FLYE on all samples
    FLYE_ASM(NANOFILT_SIMPLEX.out.nfilt_ch.flatten())

    // Generating the stats of the sequence data
	NANOPLOT_SIMPLEX(DORADO_SIMPLEX.out.summary_ch.collect())
    PYCOQC_SIMPLEX(DORADO_SIMPLEX.out.summary_ch.collect())
    NANOPLOT_FASTQ(SAMTOOLS_BAM2FQ.out.filter_ch.flatten())
    SEQKIT_SIMPLEX(SAMTOOLS_BAM2FQ.out.filter_ch.flatten())
    SEQKIT_NFILT(NANOFILT_SIMPLEX.out.nfilt_ch.collect())
    SEQKIT_FLYE(FLYE_ASM.out.assembly_ch.collect())
    
}
