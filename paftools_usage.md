# `paftools` usage

There is only one shell command, `paftools`, that provides access to all the tools. The first free parameter to the `paftools` command
specifies the tool to be used.

Throughout paftools, the `-h` option provides help. Specifically,
* `paftools -h` shows help about using the `paftools` command, including a list of all tools (i.e. subcommands),
* `paftools <cmd> -h` shows help about the tool `<cmd>`.

```shell
paftools -h
paftools recoverSeqs -h
```

## Recover Target Sequences from Fastq Files

The `recoverSeqs` tool is designed to recover target sequences, provided via a 'targets file'. It provides a complex process comprised
of the following stages:

* **Trimming** of reads, currently only trimmomatic is supported  (`--trimmer`)
* **Mapping** of reads, using either BWA or tblastn (`--mapper`)
* **Assembly** of reads into contigs, using spades or the built-in overlap based assembler,
* **Recovery** of coding sequences of the target genes.


### Run Basic Recovery

```shell
targetsFile = Angiosperms353_targetSequences_organism-gene_format_corrected.fasta
sampleID = "11943"
unzippedR1FastqFile = 11943_R1.fastq
unzippedR1FastqFile = 11943_R2.fastq
adapterFasta = illumina_adapters.fasta

paftools recoverSeqs $targetsFile ${sampleId}.fasta \
	-f $unzippedR1FastqFile -r $unzippedR2FastqFile \
	--trimmer trimmomatic --trimmomaticLeadingQuality 10 --trimmomaticTrailingQuality 10 \
	--trimmomaticMinLength 40 --trimmomaticSlidingWindowSize 4 --trimmomaticSlidingWindowQuality 20 \
	--trimmomaticAdapterFname $adapterFasta  \
	--mapper tblastn --assembler overlapSerial --blastNumThreads 4 --allowInvalidBases \
	--windowSizeReference 50 --relIdentityThresholdReference 0.7 --windowSizeReadOverlap 30 \
	--relIdentityThresholdReadOverlap 0.9 --summaryCsv ${sampleId}_summary.csv
```

The output consists of :

* The generated fasta file containing the recovered sequences
* a comma separated values (csv) file containing statistics of the recovery process,
* a `tar.gz` archive of the temporary directory used during the recovery process.


## Extracting Target Genes from Genomes or Transcriptomes

The tools for processing genomes and transcriptomes internally use a concept of a reference genome, which should be searchable by BLAST tools, specifically `blastn`. Therefore, FASTA formatted files are needed and must be indexed using `makeblastdb` before running `retrievetargets`. 

In addition to these two sequence files, a name for the reference genome, specified using the `--genomeName` option, is also required.

```shell
paftools retrievetargets targets353.fasta PRJDB1747.fasta --fastaFname Oryza_sativa.IRGSP-1.0.cds.all.fa --genomeName O_sativa
```


## Database Related Tools

* `addTargetsFile`
* `addPaftolFastq`



## Additional Notes

### Outline of the `recoverySeqs` process

The process is comprised of four major stages:

* **Trimming**, which is optional and is requested by the `--trimmer`
  option. The only trimmer currently supported is `trimmomatic`.
* **Read mapping**, which associates ("maps") reads to target genes.
  Mapping is based on sequence similarity and uses an external aligner
  to find similar sequences. Currently, `tblastn` and `bwa` are supported.
* **Assembly** of reads into contigs. Currently, support for `spades`
  and a built-in overlap based assembler are provided.
* **Coding sequence recovery**. This final stage attempts to select a
  subset of contigs that minimises redundancy (i.e. contigs pertaining
  to the same part of a target gene), and pieces a coding sequence together
  from these. This stage heavily relies on the external `exonerate` program.


### Outline of the `overlapSerial` Assembly Process

This process is comprised of the following stages:

* The target to which the largest number of reads were mapped is
  selected as the **reference target sequence**.

* Reads are **ordered** by aligning them to the reference sequence.
  Reads are pairwisely aligned to the reference sequence, and if there
  is a window in the alignment in which relative sequence identity
  meets a given threshold, the alignment start position is used as
  the ordering criterion. This is parameterised by the
  `--windowSizeReference` and `--relIdentityThresholdReference`
  options.

* Consecutive reads within this order are aligned to find **read
  overlaps**. If there is a window in this alignment in which
  relative identity meets a given threshold, the aligned reads are
  added to the same contig. This is parameterised by the
  `--windowSizeReadOverlap` and `--relIdentityThresholdReadOverlap`
  parameters.
