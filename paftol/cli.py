import sys
import argparse
import logging

import Bio
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord

import paftol
    
    
def addBwaRunnerToParser(p):
    p.add_argument('--bwaNumThreads', type=int, help='set number of threads for BWA (see bwa mem -t)')
    p.add_argument('--bwaMinSeedLength', type=int, help='set minimum seed length for BWA (see bwa mem -k)')
    p.add_argument('--bwaScoreThreshold', type=int, help='set minimum score for BWA (see bwa mem -T)')
    p.add_argument('--bwaReseedTrigger', type=float, help='set re-seed trigger BWA (see bwa mem -r)')


def argToBwaRunner(argNamespace):
    bwaRunner = paftol.tools.BwaRunner()
    bwaRunner.numThreads = argNamespace.bwaNumThreads
    bwaRunner.minSeedLength = argNamespace.bwaMinSeedLength
    bwaRunner.scoreThreshold = argNamespace.bwaScoreThreshold
    bwaRunner.reseedTrigger = argNamespace.bwaReseedTrigger
    return bwaRunner


def runHybseq(argNamespace):
    """Experimental.
"""
    hybseqAnalyser = paftol.HybseqAnalyser(argNamespace.targetsfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.tgz)
    sys.stderr.write('%s\n' % str(hybseqAnalyser))
    
    
def runHybseqstats(argNamespace):
    bwaRunner = argToBwaRunner(argNamespace)
    fastqPairList = []
    sys.stderr.write('fastqList: %s\n' % str(argNamespace.fastqList))
    for i in range(0, len(argNamespace.fastqList), 2):
        fastqPairList.append((argNamespace.fastqList[i], argNamespace.fastqList[i + 1], ))
    sdf = paftol.paftolSummary(argNamespace.targetseqs, fastqPairList, bwaRunner)
    if argNamespace.outfile is None:
        sdf.writeCsv(sys.stdout)
    else:
        with open(argNamespace.outfile, 'w') as f:
            sdf.writeCsv(f)
    

def runHybpiper(argNamespace):
    """Run an analysis (currently CDS reconstruction) using a HybPiper like approach.
"""
    bwaRunner = argToBwaRunner(argNamespace)
    hybpiperAnalyser = paftol.HybpiperAnalyser(argNamespace.tgz, bwaRunner=bwaRunner)
    if argNamespace.csv is not None:
        hybpiperAnalyser.statsCsvFilename = argNamespace.csv
    hybpiperAnalyser.keepTmpDir = True
    hybpiperResult = hybpiperAnalyser.analyse(argNamespace.targetsfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.allowInvalidBases)
    if argNamespace.outfile is not None:
        Bio.SeqIO.write([sr for sr in hybpiperResult.reconstructedCdsDict.values() if sr is not None], argNamespace.outfile, 'fasta')
    else:
        Bio.SeqIO.write([sr for sr in hybpiperResult.reconstructedCdsDict.values() if sr is not None], sys.stdout, 'fasta')
    if argNamespace.summaryCsv is not None:
        summaryStats = hybpiperResult.summaryStats()
        with open(argNamespace.summaryCsv, 'w') as f:
            summaryStats.writeCsv(f)

        
def runTargetGeneScan(argNamespace):
    paftolTargetSet = paftol.PaftolTargetSet()
    if argNamespace.targetsfile is None:
        paftolTargetSet.readFasta(sys.stdin)
    else:
        paftolTargetSet.readFasta(argNamespace.targetsfile)
    sys.stderr.write('read target set with %d genes and %d organisms\n' % (len(paftolTargetSet.paftolGeneDict), len(paftolTargetSet.organismDict)))
    # FIXME: hack to use scanMethod as the genome name as well
    referenceGenome = paftol.ReferenceGenome(argNamespace.scanMethod, argNamespace.refFasta, argNamespace.refGenbank)
    referenceGenome.scanGenes(argNamespace.scanMethod)
    sys.stderr.write('read reference genome and scanned %d genes\n' % len(referenceGenome.geneList))
    targetGeneTable = referenceGenome.blastTargetSet(paftolTargetSet)
    if argNamespace.outfile is None:
        targetGeneTable.writeCsv(sys.stdout)
    else:
        with open(argNamespace.outfile, 'w') as csvFile:
            targetGeneTable.writeCsv(csvFile)
            
            
def runGenomeReadScan(argNamespace):
    bwaRunner = argToBwaRunner(argNamespace)
    referenceGenome = paftol.ReferenceGenome(argNamespace.scanMethod, argNamespace.refFasta, argNamespace.refGenbank)
    referenceGenome.scanGenes(argNamespace.scanMethod)
    statsTable, rawmapTable = referenceGenome.mapReadsStatsBwaMem(argNamespace.forwardreads, argNamespace.reversereads, bwaRunner=bwaRunner)
    if argNamespace.rawmapTable is not None:
        with open(argNamespace.rawmapTable, 'w') as csvFile:
            rawmapTable.writeCsv(csvFile)
    if argNamespace.outfile is None:
        statsTable.writeCsv(sys.stdout)
    else:
        with open(argNamespace.outfile, 'w') as csvFile:
            statsTable.writeCsv(csvFile)


def addDevParser(subparsers):
    p = subparsers.add_parser('dev')
    p.add_argument('-s', '--str', help='set a parameter')
    p.add_argument('-i', '--int', type=int, help='set an integer')
    p.add_argument('-f', '--flag', action='store_true', help='set a flag')
    p.set_defaults(func=showArgs)

    
def addHybseqstatsParser(subparsers):
    p = subparsers.add_parser('hybseqstats')
    p.add_argument('-t', '--targetseqs', help='target sequences FASTA file')
    p.add_argument('-o', '--outfile', help='output file (default: stdout)')
    p.add_argument('fastqList', nargs='*', help='fastq file list')
    addBwaRunnerToParser(p)
    p.set_defaults(func=runHybseqstats)

    
def addHybseqParser(subparsers):
    p = subparsers.add_parser('hybseq')
    # p.add_argument('-t', '--targetseqs', help='target sequences (FASTA)')
    p.add_argument('-f', '--forwardreads', help='forward reads (FASTQ)', required=True)
    p.add_argument('-r', '--reversereads', help='reverse reads (FASTQ), omit to use single end mode')
    p.add_argument('--tgz', help='put temporary working directory into tgz')
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), default stdin')
    p.add_argument('outfile', nargs='?', help='output file (FASTA), default stdout')
    p.set_defaults(func=runHybseq)

    
def addHybpiperParser(subparsers):
    p = subparsers.add_parser('hybpiper')
    # p.add_argument('-t', '--targetseqs', help='target sequences (FASTA)')
    p.add_argument('-f', '--forwardreads', help='forward reads (FASTQ)', required=True)
    p.add_argument('-r', '--reversereads', help='reverse reads (FASTQ), omit to use single end mode')
    p.add_argument('--allowInvalidBases', action='store_true', help='allow any symbol in reference sequence (e.g. IUPAC ambiguity but also entirely invalid ones)')
    # FIXME: obsolete, to be superseded by summaryCsv
    p.add_argument('--csv', help='write basic analysis stats in CSV format')
    p.add_argument('--summaryCsv', help='write analysis stats in CSV format')
    p.add_argument('--tgz', help='put temporary working directory into tgz')
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), default stdin')
    addBwaRunnerToParser(p)
    p.add_argument('outfile', nargs='?', help='output file (FASTA), default stdout')
    p.set_defaults(func=runHybpiper)
    
    
def addTargetGeneScanParser(subparsers):
    p = subparsers.add_parser('genescan')
    p.add_argument('-r', '--refFasta', help='FASTA file of reference genome, must be BLAST indexed', required=True)
    p.add_argument('-g', '--refGenbank', help='GenBank file of reference genome, used to find genes from gene features', required=True)
    p.add_argument('-m', '--scanMethod', help='method for scanning for genes in reference genome', required=True)
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), default stdin')
    p.add_argument('outfile', nargs='?', help='output file (CSV), default stdout')
    p.set_defaults(func=runTargetGeneScan)
    
    
def addGenomeReadScanParser(subparsers):
    p = subparsers.add_parser('readscan')
    p.add_argument('--refFasta', help='FASTA file of reference genome, must be BLAST indexed', required=True)
    p.add_argument('--refGenbank', help='GenBank file of reference genome, used to find genes from gene features', required=True)
    p.add_argument('-m', '--scanMethod', help='method for scanning for genes in reference genome', required=True)
    p.add_argument('-f', '--forwardreads', help='forward reads (FASTQ)', required=True)
    p.add_argument('-r', '--reversereads', help='reverse reads (FASTQ), omit to use single end mode')
    p.add_argument('--rawmapTable', help='genomic coordinates of mapped reads')
    addBwaRunnerToParser(p)
    p.add_argument('outfile', nargs='?', help='output file (CSV), default stdout')
    p.set_defaults(func=runGenomeReadScan)
    
    
def showArgs(args):
    sys.stderr.write('%s\n' % str(args))


def paftoolsMain():
    """Entry point for the C{paftools} script.
"""
    logging.basicConfig(format='%(levelname)s: %(module)s:%(lineno)d, %(funcName)s, %(asctime)s: %(message)s')
    # logger = logging.getLogger(__name__)
    p = argparse.ArgumentParser(description='paftools -- tools for the Plant and Fungal Trees of Life (PAFTOL) project')
    p.add_argument('--loglevel', help='set logging level [DEBUG, INFO, WARNING, ERROR, CRITICAL]')
    subparsers = p.add_subparsers(title='paftools subcommands')
    addDevParser(subparsers)
    addHybseqParser(subparsers)
    addHybpiperParser(subparsers)
    addTargetGeneScanParser(subparsers)
    addGenomeReadScanParser(subparsers)
    addHybseqstatsParser(subparsers)
    args = p.parse_args()
    if args.loglevel is not None:
        loglevel = getattr(logging, args.loglevel.upper(), None)
        if loglevel is None:
            raise ValueError('invalid log level: %s' % args.loglevel)
        logging.getLogger().setLevel(loglevel)
    args.func(args)
