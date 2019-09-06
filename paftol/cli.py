import sys
import argparse
import logging

import Bio
import Bio.SeqIO
import Bio.Seq
import Bio.SeqRecord
import Bio.AlignIO
import Bio.Phylo

import paftol
import paftol.database


logger = logging.getLogger(__name__)


def addTrimmomaticRunnerToParser(p):
    p.add_argument('--trimmomaticNumThreads', type=int, help='set number of threads for trimmomatic (see trimmomatic -threads)')
    p.add_argument('--trimmomaticLeadingQuality', type=int, help='set leading quality threshold (see trimmomatic LEADING step)')
    p.add_argument('--trimmomaticTrailingQuality', type=int, help='set trailing quality threshold (see trimmomatic TRAILING step)')
    p.add_argument('--trimmomaticMinLength', type=int, help='set minimum length after trimming (see trimmomatic MINLEN step)')
    p.add_argument('--trimmomaticSlidingWindowSize', type=int, help='set sliding window size (see trimmomatic SLIDINGWINDOW step, first parameter)')
    p.add_argument('--trimmomaticSlidingWindowQuality', type=int, help='set sliding window average quality (see trimmomatic SLIDINGWINDOW step, second parameter)')
    p.add_argument('--trimmomaticAdapterFname', help='set adapter file name (see trimmomatic ILLUMINACLIP step)')
    
    
def addBwaRunnerToParser(p):
    p.add_argument('--bwaNumThreads', type=int, help='set number of threads for BWA (see bwa mem -t)')
    p.add_argument('--bwaMinSeedLength', type=int, help='set minimum seed length for BWA (see bwa mem -k)')
    p.add_argument('--bwaScoreThreshold', type=int, help='set minimum score for BWA (see bwa mem -T)')
    p.add_argument('--bwaReseedTrigger', type=float, help='set re-seed trigger BWA (see bwa mem -r)')

    
def addBlastRunnerToParser(p):
    p.add_argument('--blastNumThreads', type=int, help='set number of threads for blast program (see -num_threads)')
    p.add_argument('--blastGapOpen', type=int, help='cost of opening a gap (see -gapopen)')
    p.add_argument('--blastGapExtend', type=int, help='cost of extending a gap (see -gapextend)')
    p.add_argument('--blastEvalue', type=float, help='E value threshold (see -evalue)')
    p.add_argument('--blastWindowSize', type=int, help='multiple hits window size (see -window_size)')

    
def addOverlapAssemblerToParser(p):
    p.add_argument('--windowSizeReference', type=int, help='window size for reference to read alignment')
    p.add_argument('--relIdentityThresholdReference', type=float, help='percent identity threshold for reference to read alignment')
    p.add_argument('--windowSizeReadOverlap', type=int, help='window size for read overlap alignment')
    p.add_argument('--relIdentityThresholdReadOverlap', type=float, help='percent identity threshold for read overlap alignment')
    
    
def addHybseqToParser(p):
    # p.add_argument('-t', '--targetseqs', help='target sequences (FASTA)')
    p.add_argument('--allowInvalidBases', action='store_true', help='allow any symbol in reference sequence (e.g. IUPAC ambiguity but also entirely invalid ones)')
    p.add_argument('--strictOverlapFiltering', action='store_true', help='filter contigs so that no region of the reference target is covered by multiple (overlapping) contigs')
    p.add_argument('--maxNumReadsPerGene', type=int, help='maximal number of reads used for assembly (in development)')
    p.add_argument('-f', '--forwardreads', help='forward reads (FASTQ)', required=True)
    p.add_argument('-r', '--reversereads', help='reverse reads (FASTQ), omit to use single end mode')
    p.add_argument('--exoneratePercentIdentityThreshold', type=float, help='percent identity threshold for reference to contig alignment')
    p.add_argument('--tgz', help='put temporary working directory into tgz')
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), default stdin')
    p.add_argument('outfile', nargs='?', help='output file (FASTA), default stdout')
    p.add_argument('--summaryCsv', help='write analysis stats in CSV format')


def addHybpiperToParser(p):
    addHybseqToParser(p)
    addSpadesRunnerToParser(p)


def addTblastnRunnerToParser(p):
    addBlastRunnerToParser(p)


def addBlastnRunnerToParser(p):
    addBlastRunnerToParser(p)

    
def addSpadesRunnerToParser(p):
    p.add_argument('--spadesNumThreads', type=int, help='set number of threads for SPAdes (see spades -t)')
    p.add_argument('--spadesCovCutoff', help='set coverage cutoff for SPAdes (see spades --cov-cutoff, string to allow "auto" and "off")')


def checkAbsenceOfOptions(illegalOptNameList, argNamespace, msg):
    foundIllegalOptNameList = []
    for optName in illegalOptNameList:
        if argNamespace.__dict__[optName] is not None :
            foundIllegalOptNameList.append(optName)
    if len(foundIllegalOptNameList) > 0:
        raise StandardError, '%s: %s' % (msg, ', '.join(foundIllegalOptNameList))


def requiredArg(arg, msg):
    if arg is None:
        raise StandardError, msg
    return arg

def argToOverlapAssemblerSerial(argNamespace):
    targetAssemblerOverlapSerial = paftol.TargetAssemblerOverlapSerial()
    targetAssemblerOverlapSerial.windowSizeReference = requiredArg(argNamespace.windowSizeReference, 'windowSizeReference is required')
    targetAssemblerOverlapSerial.relIdentityThresholdReference = requiredArg(argNamespace.relIdentityThresholdReference, 'relIdentityThresholdReference is required')
    targetAssemblerOverlapSerial.windowSizeReadOverlap = requiredArg(argNamespace.windowSizeReadOverlap, 'windowSizeReadOverlap is required')
    targetAssemblerOverlapSerial.relIdentityThresholdReadOverlap = requiredArg(argNamespace.relIdentityThresholdReadOverlap, 'relIdentityThresholdReadOverlap is required')
    targetAssemblerOverlapSerial.semiglobalAlignmentRunner = paftol.tools.SemiglobalAlignmentRunner()
    return targetAssemblerOverlapSerial


def argToTrimmomaticRunner(argNamespace):
    trimmomaticRunner = paftol.tools.TrimmomaticRunner()
    trimmomaticRunner.numThreads = argNamespace.trimmomaticNumThreads
    trimmomaticRunner.leadingQuality = argNamespace.trimmomaticLeadingQuality
    trimmomaticRunner.trailingQuality = argNamespace.trimmomaticTrailingQuality
    trimmomaticRunner.minLength = argNamespace.trimmomaticMinLength
    trimmomaticRunner.slidingWindowSize = argNamespace.trimmomaticSlidingWindowSize
    trimmomaticRunner.slidingWindowQuality = argNamespace.trimmomaticSlidingWindowQuality
    trimmomaticRunner.adapterFname = argNamespace.trimmomaticAdapterFname
    return trimmomaticRunner


def checkAbsenceOfTrimmomaticOptions(argNamespace, msg):
    checkAbsenceOfOptions(['trimmomaticNumThreads', 'trimmomaticLeadingQuality', 'trimmomaticTrailingQuality', 'trimmomaticMinLength', 'trimmomaticSlidingWindowSize', 'trimmomaticSlidingWindowQuality', 'trimmomaticAdapterFname'], argNamespace, msg)


def argToBwaRunner(argNamespace):
    bwaRunner = paftol.tools.BwaRunner()
    bwaRunner.numThreads = argNamespace.bwaNumThreads
    bwaRunner.minSeedLength = argNamespace.bwaMinSeedLength
    bwaRunner.scoreThreshold = argNamespace.bwaScoreThreshold
    bwaRunner.reseedTrigger = argNamespace.bwaReseedTrigger
    return bwaRunner


def checkAbsenceOfBwaOptions(argNamespace, msg):
    checkAbsenceOfOptions(['bwaNumThreads', 'bwaMinSeedLength', 'bwaScoreThreshold', 'bwaReseedTrigger'], argNamespace, msg)


def argToBlastRunnerParams(argNamespace, blastRunner):
    blastRunner.numThreads = argNamespace.blastNumThreads
    blastRunner.gapOpen = argNamespace.blastGapOpen
    blastRunner.gapExtend = argNamespace.blastGapExtend
    blastRunner.evalue = argNamespace.blastEvalue
    blastRunner.windowSize = argNamespace.blastWindowSize
    return blastRunner


def argToTblastnRunner(argNamespace):
    tblastnRunner = paftol.tools.TblastnRunner()
    argToBlastRunnerParams(argNamespace, tblastnRunner)
    return tblastnRunner


def checkAbsenceOfTblastnOptions(argNamespace, msg):
    checkAbsenceOfOptions(['bastNumThreads', 'blastGapOpen', 'blastGapExtend', 'blastEvalue', 'blastWindowSize'], argNamespace, msg)


def argToBlastnRunner(argNamespace):
    blastnRunner = paftol.tools.BlastnRunner()
    argToBlastRunnerParams(argNamespace, blastnRunner)
    return blastnRunner


def checkAbsenceOfTblastnOptions(argNamespace, msg):
    checkAbsenceOfOptions(['bastNumThreads', 'blastGapOpen', 'blastGapExtend', 'blastEvalue', 'blastWindowSize'], argNamespace, msg)


def argToSpadesRunner(argNamespace):
    spadesRunner = paftol.tools.SpadesRunner()
    spadesRunner.numThreads = argNamespace.spadesNumThreads
    spadesRunner.covCutoff = argNamespace.spadesCovCutoff
    return spadesRunner


def checkAbsenceOfSpadesnOptions(argNamespace, msg):
    checkAbsenceOfOptions(['spadesNumThreads', 'spadesCovCutoff'], argNamespace, msg)


# FIXME: obsolete -- summary stats are now provided by result (so not adding spadesRunner here)
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


def runHybpiperBwa(argNamespace):
    """Run an analysis (currently CDS reconstruction) using a HybPiper
like approach, unsing BWA for mapping reads to targets.

    """
    bwaRunner = argToBwaRunner(argNamespace)
    spadesRunner = argToSpadesRunner(argNamespace)
    # FIXME: backwards compatibility to previous implementation and to hybpiper (??)
    if spadesRunner.covCutoff is None:
        spadesRunner.covCutoff = 8
        logger.warning('SPAdes coverage cutoff not specified, set to %d for backwards compatibility', spadesRunner.covCutoff)
    hybpiperBwaAnalyser = paftol.HybpiperBwaAnalyser(argNamespace.tgz, bwaRunner=bwaRunner, spadesRunner=spadesRunner)
    # hybpiperBwaAnalyser.keepTmpDir = True
    targetsfile = sys.stdin if argNamespace.targetsfile is None else argNamespace.targetsfile
    hybpiperResult = hybpiperBwaAnalyser.analyse(targetsFile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.allowInvalidBases, argNamespace.strictOverlapFiltering, argNamespace.maxNumReadsPerGene)
    if argNamespace.outfile is not None:
        Bio.SeqIO.write([sr for sr in hybpiperResult.reconstructedCdsDict.values() if sr is not None], argNamespace.outfile, 'fasta')
    else:
        Bio.SeqIO.write([sr for sr in hybpiperResult.reconstructedCdsDict.values() if sr is not None], sys.stdout, 'fasta')
    if argNamespace.summaryCsv is not None:
        summaryStats = hybpiperResult.summaryStats()
        with open(argNamespace.summaryCsv, 'w') as f:
            summaryStats.writeCsv(f)


def runHybpiperTblastn(argNamespace):
    """Run an analysis (currently CDS reconstruction) using a HybPiper
like approach, unsing tblastn for mapping reads to targets.

    """
    tblastnRunner = argToTblastnRunner(argNamespace)
    spadesRunner = argToSpadesRunner(argNamespace)
    # FIXME: backwards compatibility to previous implementation and to hybpiper (??)
    if spadesRunner.covCutoff is None:
        spadesRunner.covCutoff = 8
        logger.warning('SPAdes coverage cutoff not specified, set to %d for backwards compatibility', spadesRunner.covCutoff)
    hybpiperTblastnAnalyser = paftol.HybpiperTblastnAnalyser(argNamespace.tgz, tblastnRunner=tblastnRunner, spadesRunner=spadesRunner)
    targetsfile = sys.stdin if argNamespace.targetsfile is None else argNamespace.targetsfile
    hybpiperResult = hybpiperTblastnAnalyser.analyse(targetsfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.allowInvalidBases, argNamespace.strictOverlapFiltering, argNamespace.maxNumReadsPerGene)
    if argNamespace.outfile is not None:
        Bio.SeqIO.write([sr for sr in hybpiperResult.reconstructedCdsDict.values() if sr is not None], argNamespace.outfile, 'fasta')
    else:
        Bio.SeqIO.write([sr for sr in hybpiperResult.reconstructedCdsDict.values() if sr is not None], sys.stdout, 'fasta')
    if argNamespace.summaryCsv is not None:
        summaryStats = hybpiperResult.summaryStats()
        with open(argNamespace.summaryCsv, 'w') as f:
            summaryStats.writeCsv(f)


def runTargetRecovery(argNamespace):
    if argNamespace.usePaftolDb:
        paftol.database.preRecoveryCheck(argNamespace.forwardreads, argNamespace.reversereads)
    trimmomaticRunner = None
    if argNamespace.trimmer == 'trimmomatic':
        trimmomaticRunner = argToTrimmomaticRunner(argNamespace)
        logger.debug('got trimmomaticRunner')
    else:
        checkAbsenceOfTrimmomaticOptions(argNamespace, 'options illegal without "--trim trimmomatic"')
    if argNamespace.mapper == 'tblastn':
        checkAbsenceOfBwaOptions(argNamespace, 'options illegal without "--mapper bwa"')
        tblastnRunner = argToTblastnRunner(argNamespace)
        targetMapper = paftol.TargetMapperTblastn(tblastnRunner)
    elif argNamespace.mapper == 'bwa':
        checkAbsenceOfTblastnOptions(argNamespace, 'options illegal without "--mapper tblastn"')
        bwaRunner = argToBwaRunner(argNamespace)
        targetMapper = paftol.TargetMapperBwa(bwaRunner)
    else:
        raise StandardError, 'unsupported / unknown mapper %s' % argNamespace.mapper
    if argNamespace.assembler == 'spades':
        raise StandardError, 'spades assembly not yet refactored'
    elif argNamespace.assembler == 'overlapSerial':
        targetAssembler = argToOverlapAssemblerSerial(argNamespace)
    targetRecoverer = paftol.TargetRecoverer(argNamespace.tgz, 'targetrecover', trimmomaticRunner=trimmomaticRunner, targetMapper=targetMapper, targetAssembler=targetAssembler)
    targetsfile = sys.stdin if argNamespace.targetsfile is None else argNamespace.targetsfile
    result = targetRecoverer.recoverTargets(targetsfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.allowInvalidBases, argNamespace.strictOverlapFiltering, argNamespace.maxNumReadsPerGene)
    result.cmdLine = argNamespace.rawCmdLine
    if argNamespace.outfile is not None:
        Bio.SeqIO.write([sr for sr in result.reconstructedCdsDict.values() if sr is not None], argNamespace.outfile, 'fasta')
    else:
        Bio.SeqIO.write([sr for sr in result.reconstructedCdsDict.values() if sr is not None], sys.stdout, 'fasta')
    if argNamespace.contigFname is not None:
        result.writeContigFastaFile(argNamespace.contigFname)
    if argNamespace.summaryCsv is not None:
        summaryStats = result.summaryStats()
        with open(argNamespace.summaryCsv, 'w') as f:
            summaryStats.writeCsv(f)
    if argNamespace.usePaftolDb:
        paftol.database.addRecoveryResult(result)
    

def runOverlapAnalysis(argNamespace):
    """Run an analysis using tblastn for mapping to targets and overlap based assembly for gene recovery.

    """
    tblastnRunner = argToTblastnRunner(argNamespace)
    overlapAnalyser = paftol.OverlapAnalyser(argNamespace.tgz, tblastnRunner=tblastnRunner)
    overlapAnalyser.windowSizeReference = argNamespace.windowSizeReference
    overlapAnalyser.relIdentityThresholdReference = argNamespace.relIdentityThresholdReference
    overlapAnalyser.windowSizeReadOverlap = argNamespace.windowSizeReadOverlap
    overlapAnalyser.relIdentityThresholdReadOverlap = argNamespace.relIdentityThresholdReadOverlap
    targetsfile = sys.stdin if argNamespace.targetsfile is None else argNamespace.targetsfile
    result = overlapAnalyser.analyse(targetsfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.allowInvalidBases, argNamespace.strictOverlapFiltering, argNamespace.maxNumReadsPerGene)
    if argNamespace.outfile is not None:
        Bio.SeqIO.write([sr for sr in result.reconstructedCdsDict.values() if sr is not None], argNamespace.outfile, 'fasta')
    else:
        Bio.SeqIO.write([sr for sr in result.reconstructedCdsDict.values() if sr is not None], sys.stdout, 'fasta')
    if argNamespace.summaryCsv is not None:
        summaryStats = result.summaryStats()
        with open(argNamespace.summaryCsv, 'w') as f:
            summaryStats.writeCsv(f)


def runRetrieveTargets(argNamespace):
    paftolTargetSet = paftol.PaftolTargetSet()
    if argNamespace.targetsfile is None:
        paftolTargetSet.readFasta(sys.stdin)
    else:
        paftolTargetSet.readFasta(argNamespace.targetsfile)
    blastnRunner = argToBlastnRunner(argNamespace)
    paftolTargetSeqRetriever = paftol.PaftolTargetSeqRetriever()
    targetList = paftolTargetSeqRetriever.retrievePaftolTargetList(argNamespace.genomeName, argNamespace.fastaFname, paftolTargetSet, blastnRunner)
    if argNamespace.outfile is None:
        Bio.SeqIO.write(targetList, sys.stdout, 'fasta')
    else :
        with open(argNamespace.outfile, 'w') as f:
            Bio.SeqIO.write(targetList, f, 'fasta')


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
    targetGeneTable, cdsList = referenceGenome.blastTargetSet(paftolTargetSet)
    if argNamespace.outfile is None:
        targetGeneTable.writeCsv(sys.stdout)
    else:
        with open(argNamespace.outfile, 'w') as csvFile:
            targetGeneTable.writeCsv(csvFile)
    if argNamespace.cdsFasta is not None:
        Bio.SeqIO.write(cdsList, argNamespace.cdsFasta, 'fasta')
            
            
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

            
def runExtractCdsList(argNamespace):
    referenceGenome = paftol.ReferenceGenome(argNamespace.genomeName, argNamespace.refFasta, argNamespace.refGenbank)
    cdsList = referenceGenome.makeCdsList()
    if argNamespace.outfile is None:
        Bio.SeqIO.write(cdsList, sys.stdout, 'fasta')
    else :
        with open(argNamespace.outfile, 'w') as f:
            Bio.SeqIO.write(cdsList, f, 'fasta')


def runSelectgenes(argNamespace):
    organismNameSet = None
    if argNamespace.organism is not None:
        organismNameSet = set(argNamespace.organism)
    geneNameSet = []
    if argNamespace.gene is not None:
        geneNameSet = set(argNamespace.gene)
    if argNamespace.genefile is not None:
        with open(argNamespace.genefile) as f:
            for line in f:
                geneNameSet.append(line.strip())
    if len(geneNameSet) == 0:
        geneNameSet = None
    paftolTargetSet = paftol.PaftolTargetSet()
    if argNamespace.targetsfile is None:
        paftolTargetSet.readFasta(sys.stdin)
    else:
        paftolTargetSet.readFasta(argNamespace.targetsfile)
    srList = paftolTargetSet.getSeqRecordSelection(organismNameSet, geneNameSet)
    if argNamespace.outfile is None:
        Bio.SeqIO.write(srList, sys.stdout, 'fasta')
    else:
        Bio.SeqIO.write(srList, argNamespace.outfile, 'fasta')


def runNeedleComparison(argNamespace):
    sr1Dict = Bio.SeqIO.to_dict(Bio.SeqIO.parse(argNamespace.seqfile1, 'fasta'))
    sr2Dict = Bio.SeqIO.to_dict(Bio.SeqIO.parse(argNamespace.seqfile2, 'fasta'))
    needleRunner = paftol.tools.NeedleRunner()
    cmpStats = paftol.tools.pairwiseAlignmentStats(sr1Dict, sr2Dict, needleRunner)
    if argNamespace.outfile is None:
        cmpStats.writeCsv(sys.stdout)
    else:
        with open(argNamespace.outfile, 'w') as csvFile:
            cmpStats.writeCsv(csvFile)


def runGeneSetStats(argNamespace):
    paftolTargetSet = paftol.PaftolTargetSet()
    sampleId = 'unknown'
    if argNamespace.sampleId is not None:
        sampleId = argNamespace.sampleId
    if argNamespace.targetsfile is None:
        raise StandardError, 'no targets file specified'
    else:
        paftolTargetSet.readFasta(argNamespace.targetsfile)
    if argNamespace.seqfile is None:
        geneSetStatsDataFrame = paftol.makeGeneSetStatsDataFrame(sys.stdin, sampleId, paftolTargetSet)
    else:
        with open(argNamespace.seqfile, 'r') as f:
            geneSetStatsDataFrame = paftol.makeGeneSetStatsDataFrame(f, sampleId, paftolTargetSet)
    if argNamespace.outfile is None:
        geneSetStatsDataFrame.writeCsv(sys.stdout)
    else:
        with open(argNamespace.outfile, 'w') as f:
            geneSetStatsDataFrame.writeCsv(f)


def runExonerateStarAlignment(argNamespace):
    reference = Bio.SeqIO.read(argNamespace.reference, 'fasta')
    fastaFname = argNamespace.seqfile
    exonerateStarAlignment = paftol.tools.ExonerateStarAlignment(reference, fastaFname)
    if argNamespace.epsfname is not None:
        with open(argNamespace.epsfname,  'w') as epsfile:
            exonerateStarAlignment.epsSketch(epsfile)
    if argNamespace.outfile is None:
        Bio.AlignIO.write(exonerateStarAlignment.xstarAlignment, sys.stdout, 'fasta')
    else:
        with open(argNamespace.outfile, 'w') as f:
            Bio.AlignIO.write(exonerateStarAlignment.xstarAlignment, f, 'fasta')

            
def runAddOrganism(argNamespace):
    logger.info('starting, organismName = %s', argNamespace.organismName)
    if '-' in argNamespace.organismName:
        raise StandardError, 'invalid organism name: %s (contains at least one dash)' % argNamespace.organismName
    if argNamespace.inFasta is None:
        infile = sys.stdin
    else:
        infile = open(argNamespace.inFasta, 'r')
    if argNamespace.outFasta is None:
        outfile = sys.stdout
    else:
        outfile = open(argNamespace.outFasta, 'w')
    for seqRecord in Bio.SeqIO.parse(infile, 'fasta'):
        seqRecord.id = '%s-%s' % (argNamespace.organismName, seqRecord.id)
        Bio.SeqIO.write([seqRecord], outfile, 'fasta')

            
def runAddTargetsFile(argNamespace):
    targetsfile = argNamespace.targetsfile
    description = argNamespace.description
    # sys.stderr.write('insertGenes: %s, geneType: %s\n' % (str(argNamespace.insertGenes), str(argNamespace.geneType)))
    paftol.database.addTargetsFile(targetsfile, description, argNamespace.insertGenes, argNamespace.geneType)

    
def runAddPaftolFastqFiles(argNamespace):
    paftol.database.addPaftolFastqFiles(argNamespace.fastq)


def runAlignmentOutline(argNamespace):
    alignment = Bio.AlignIO.read(argNamespace.infile, 'fasta')
    paftol.tools.plotAlignmentPostscript(alignment, argNamespace.outfile, alignment.get_alignment_length() * 0.01, len(alignment) * 0.01)


def runGeneFasta(argNamespace):
    paftolTargetSet = paftol.PaftolTargetSet()
    paftolTargetSet.readFasta(argNamespace.infile)
    fastaFnameList = []
    for geneName in paftolTargetSet.paftolGeneDict:
        srList = paftolTargetSet.getSeqRecordSelection(organismNameList=None, geneNameList=[geneName])
        fastaFname = argNamespace.outFastaFormat % geneName
        Bio.SeqIO.write(srList, fastaFname, 'fasta')
        fastaFnameList.append(fastaFname)
    sys.stdout.write('%s\n' % ' '.join(fastaFnameList))


def runDelgeneNewick(argNamespace):
    treeList = []
    for tree in Bio.Phylo.NewickIO.parse(argNamespace.infile):
        terminalCladeList = tree.clade.get_terminals()
        # sys.stderr.write('%s\n' % str(terminalCladeList)))
        for clade in terminalCladeList:
            if clade.name is None:
                raise StandardError, 'clade has no name'
            organismName, geneName = paftol.extractOrganismAndGeneNames(clade.name)
            clade.name = organismName
        treeList.append(tree)
    Bio.Phylo.NewickIO.write(treeList, argNamespace.outfile)


def addDevParser(subparsers):
    p = subparsers.add_parser('dev', help='for testing argparse')
    p.add_argument('-s', '--str', help='set a parameter')
    p.add_argument('-i', '--int', type=int, help='set an integer')
    p.add_argument('-f', '--flag', action='store_true', help='set a flag')
    p.set_defaults(func=showArgs)

    
def addHybseqstatsParser(subparsers):
    p = subparsers.add_parser('hybseqstats', help='(obsolete)')
    p.add_argument('-t', '--targetseqs', help='target sequences FASTA file')
    p.add_argument('-o', '--outfile', help='output file (default: stdout)')
    p.add_argument('fastqList', nargs='*', help='fastq file list')
    addBwaRunnerToParser(p)
    addSpadesRunnerToParser(p)
    p.set_defaults(func=runHybseqstats)
    
    
def addHybpiperBwaParser(subparsers):
    p = subparsers.add_parser('hybpiperBwa', help='recover PAFTOL gene sequences using BWA and SPAdes assembly')
    addBwaRunnerToParser(p)
    addHybpiperToParser(p)
    p.set_defaults(func=runHybpiperBwa)

    
def addHybpiperTblastnParser(subparsers):
    p = subparsers.add_parser('hybpiperTblastn', help='recover PAFTOL gene sequences using tblastn and SPAdes assembly')
    addTblastnRunnerToParser(p)
    addHybpiperToParser(p)
    p.set_defaults(func=runHybpiperTblastn)

    
def addOverlapAnalyserParser(subparsers):
    p = subparsers.add_parser('overlapRecover', help='recover PAFTOL gene sequences using tblastn and overlap based assembly')
    addTblastnRunnerToParser(p)
    addHybseqToParser(p)
    addOverlapAssemblerToParser(p)
    p.set_defaults(func=runOverlapAnalysis)
    
    
def addRecoverParser(subparsers):
    p = subparsers.add_parser('recoverSeqs', help='recover target sequences from HybSeq NGS files')
    addHybseqToParser(p)
    p.add_argument('--trimmer', choices=['trimmomatic'], help='method for trimming reads')
    p.add_argument('--mapper', choices=['tblastn', 'bwa'], help='method to be used for mapping reads to target genes', required=True)
    p.add_argument('--assembler', choices=['spades', 'overlapSerial'], help='method to be used to assemble reads mapped to a gene into contigs', required=True)
    p.add_argument('--contigFname', help='filename for contigs')
    p.add_argument('--usePaftolDb', action='store_true', help='store results in PAFTOL database')
    addTrimmomaticRunnerToParser(p)
    addTblastnRunnerToParser(p)
    addBwaRunnerToParser(p)    
    addOverlapAssemblerToParser(p)
    p.set_defaults(func=runTargetRecovery)

    
def addTargetGeneScanParser(subparsers):
    p = subparsers.add_parser('genescan', help='identify PAFTOL genes within a reference genome')
    p.add_argument('-r', '--refFasta', help='FASTA file of reference genome, must be BLAST indexed', required=True)
    p.add_argument('-g', '--refGenbank', help='GenBank file of reference genome, used to find genes from gene features', required=True)
    p.add_argument('-m', '--scanMethod', help='method for scanning for genes in reference genome', required=True)
    p.add_argument('-c', '--cdsFasta', help='output FASTA of coding sequences of matched genes')
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), default stdin')
    p.add_argument('outfile', nargs='?', help='output file (CSV), default stdout')
    p.set_defaults(func=runTargetGeneScan)


def addGenomeReadScanParser(subparsers):
    p = subparsers.add_parser('readscan', help='map reads to reference genome and compute stats of unmapped / mapping to gene / mapping to intergenic region')
    p.add_argument('--refFasta', help='FASTA file of reference genome, must be BLAST indexed', required=True)
    p.add_argument('--refGenbank', help='GenBank file of reference genome, used to find genes from gene features', required=True)
    p.add_argument('-m', '--scanMethod', help='method for scanning for genes in reference genome', required=True)
    p.add_argument('-f', '--forwardreads', help='forward reads (FASTQ)', required=True)
    p.add_argument('-r', '--reversereads', help='reverse reads (FASTQ), omit to use single end mode')
    p.add_argument('--rawmapTable', help='genomic coordinates of mapped reads')
    addBwaRunnerToParser(p)
    p.add_argument('outfile', nargs='?', help='output file (CSV), default stdout')
    p.set_defaults(func=runGenomeReadScan)


def addExtractCdsListParser(subparsers):
    p = subparsers.add_parser('xcds', help='extract CDS list from a genome')
    p.add_argument('--genomeName', help='genome name', required=True)
    p.add_argument('--refFasta', help='FASTA file of reference genome', required=True)
    p.add_argument('--refGenbank', help='GenBank file of reference genome, used to find genes from gene features', required=True)
    p.add_argument('outfile', nargs='?', help='output file (fasta), default stdout')
    p.set_defaults(func=runExtractCdsList)


def addRetrieveTargetListParser(subparsers):
    p = subparsers.add_parser('retrievetargets', help='retrieve targets list from a fasta file containing a transcriptome')
    p.add_argument('--genomeName', help='genome name', required=True)
    p.add_argument('--fastaFname', help='FASTA file ', required=True)
    addBlastnRunnerToParser(p)
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), default stdin')
    p.add_argument('outfile', nargs='?', help='output file (fasta), default stdout')
    p.set_defaults(func=runRetrieveTargets)


def addSelectgenesParser(subparsers):
    p = subparsers.add_parser('selectgenes', help='select genes by name from a PAFTOL target set')
    p.add_argument('-o', '--organism', action='append', help='specify organism name on command line')
    p.add_argument('-g', '--gene', action='append', help='specify gene name on command line')
    p.add_argument('-f', '--genefile', help='specify file containing gene names')
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), default stdin')
    p.add_argument('outfile', nargs='?', help='output file (CSV), default stdout')
    p.set_defaults(func=runSelectgenes)


def addNeedleComparisonParser(subparsers):
    p = subparsers.add_parser('needlecmp', help='compare two sets of sequences by pairwise alignment using EMBOSS needle')
    p.add_argument('seqfile1', help='sequence file 1 (required)')
    p.add_argument('seqfile2', help='sequence file 2 (required)')
    p.add_argument('outfile', nargs='?', help='output file (CSV), default stdout')
    p.set_defaults(func=runNeedleComparison)
    
    
def addGeneSetStatsParser(subparsers):
    p = subparsers.add_parser('genesetstats', help='extract stats from a FASTA file containing a set of PAFTOL genes from one organism')
    p.add_argument('--targetsfile', nargs='?', help='target sequences (FASTA), required', required=True)
    p.add_argument('-s', '--sampleId', help='sample ID, default "unknown"')
    p.add_argument('seqfile', nargs='?', help='sequence file (FASTA), default stdin')
    p.add_argument('outfile', nargs='?', help='output file (CSV), default stdout')
    p.set_defaults(func=runGeneSetStats)


def addExonerateStarAlignmentParser(subparsers):
    p = subparsers.add_parser('xstar', help='compute a star multiple alignment around a reference sequence using exonerate')
    p.add_argument('reference', help='reference sequence file')
    p.add_argument('seqfile', help='further sequences for constructing the star alignment (required)')
    p.add_argument('outfile', nargs='?', help='output file, default stdout')
    p.add_argument('--epsfname', help='sketch (encapsulated postscript)')
    p.set_defaults(func=runExonerateStarAlignment)
    
    
def addAddTargetsFileParser(subparsers):
    p = subparsers.add_parser('addTargetsFile', help='add reference targets from a targets file')
    p.add_argument('--description', help='description of the targets file')
    p.add_argument('--geneType', help='specify gene type')
    p.add_argument('--insertGenes', action='store_true', help='insert new genes that are not already in database')
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), required')
    p.set_defaults(func=runAddTargetsFile)

    
def addAddPaftolFastqFilesParser(subparsers):
    p = subparsers.add_parser('addPaftolFastq', help='add PAFTOL fastq files')
    p.add_argument('fastq', nargs='+', help='fastq files (any number, at least one)')
    p.set_defaults(func=runAddPaftolFastqFiles)


def addAddOrganismParser(subparsers):
    p = subparsers.add_parser('addorganism', help='add organism name to FASTA sequence identifier (assuming existing FASTA identifier is a gene name)')
    # FIXME: make organismName mandatory
    p.add_argument('organismName', help='organism name')
    p.add_argument('inFasta', help='input file (FASTA format)')
    p.add_argument('outFasta', help='output file (FASTA format)')
    p.set_defaults(func=runAddOrganism)

    
def addAlignmentOutlineParser(subparsers):
    p = subparsers.add_parser('alignmentOutline', help='write an outline graphics file (EPS format) of an alignment')
    p.add_argument('infile', nargs='?', default=sys.stdin, type=argparse.FileType('r'), help='input file (multiple sequence alignment in FASTA format)')
    p.add_argument('outfile', nargs='?', default=sys.stdout, type=argparse.FileType('w'), help='output file (encapsulated PostScript)')
    p.set_defaults(func=runAlignmentOutline)


def addGeneFastaParser(subparsers):
    p = subparsers.add_parser('geneFasta', help='split target sequence FASTA input into multiple FASTA files by gene name')
    p.add_argument('infile', nargs='?', default=sys.stdin, type=argparse.FileType('r'), help='input file (target sequences in FASTA format with canonical organism-gene IDs)')
    p.add_argument('outFastaFormat', nargs='?', default='%s.fasta', help='format string for output fasta files, should contain one %%s conversion for gene name')
    p.set_defaults(func=runGeneFasta)


def addDelgeneNewickParser(subparsers):
    p = subparsers.add_parser('delgeneNewick', help='delete gene name part from leaf labels in a Newick file')
    p.add_argument('infile', nargs='?', default=sys.stdin, type=argparse.FileType('r'), help='input file (Newick format, leaf labels must be canonical organism-gene IDs)')
    p.add_argument('outfile', nargs='?', default=sys.stdout, type=argparse.FileType('w'), help='output file (Newick format)')
    p.set_defaults(func=runDelgeneNewick)

    
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
    addHybpiperBwaParser(subparsers)
    addHybpiperTblastnParser(subparsers)
    addOverlapAnalyserParser(subparsers)
    addRecoverParser(subparsers)
    addTargetGeneScanParser(subparsers)
    addGenomeReadScanParser(subparsers)
    addExtractCdsListParser(subparsers)
    addRetrieveTargetListParser(subparsers)
    addHybseqstatsParser(subparsers)
    addSelectgenesParser(subparsers)
    addNeedleComparisonParser(subparsers)
    addGeneSetStatsParser(subparsers)
    addExonerateStarAlignmentParser(subparsers)
    addAddTargetsFileParser(subparsers)
    addAddPaftolFastqFilesParser(subparsers)
    addAddOrganismParser(subparsers)
    addAlignmentOutlineParser(subparsers)
    addGeneFastaParser(subparsers)
    addDelgeneNewickParser(subparsers)
    args = p.parse_args()
    args.rawCmdLine = ' '.join(['%s' % arg for arg in sys.argv])
    if args.loglevel is not None:
        loglevel = getattr(logging, args.loglevel.upper(), None)
        if loglevel is None:
            raise ValueError('invalid log level: %s' % args.loglevel)
        logging.getLogger().setLevel(loglevel)
    args.func(args)
    logger.info('paftools command completed')
