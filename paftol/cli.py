import sys
import argparse
import logging

import Bio
import Bio.SeqIO

import paftol


def runHybseq(argNamespace):
    """Experimental.
"""
    hybseqAnalyser = paftol.HybseqAnalyser(argNamespace.targetsfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.tgz)
    sys.stderr.write('%s\n' % str(hybseqAnalyser))


def runHybpiper(argNamespace):
    """Run an analysis (currently CDS reconstruction) using a HybPiper like approach.
"""
    hybpiperAnalyser = paftol.HybpiperAnalyser(argNamespace.targetsfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.tgz)
    if argNamespace.bwaMinSeedLength is not None:
        hybpiperAnalyser.bwaMinSeedLength = argNamespace.bwaMinSeedLength
    if argNamespace.bwaScoreThreshold is not None:
        hybpiperAnalyser.bwaScoreThreshold = argNamespace.bwaScoreThreshold
    if argNamespace.bwaReseedTrigger is not None:
        hybpiperAnalyser.bwaReseedTrigger = argNamespace.bwaReseedTrigger
    hybpiperAnalyser.allowInvalidBases = argNamespace.allowInvalidBases
    hybpiperAnalyser.keepTmpDir = True
    reconstructedCdsDict = hybpiperAnalyser.analyse()
    if argNamespace.outfile is not None:
        Bio.SeqIO.write([sr for sr in reconstructedCdsDict.values() if sr is not None], argNamespace.outfile, 'fasta')
    else:
        Bio.SeqIO.write([sr for sr in reconstructedCdsDict.values() if sr is not None], sys.stdout, 'fasta')
        

def addDevParser(subparsers):
    p = subparsers.add_parser('dev')
    p.add_argument('-s', '--str', help='set a parameter')
    p.add_argument('-i', '--int', type=int, help='set an integer')
    p.add_argument('-f', '--flag', action='store_true', help='set a flag')
    p.set_defaults(func=showArgs)

    
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
    p.add_argument('--bwaMinSeedLength', type=int, help='set minimum seed length for BWA (see bwa mem -k)')
    p.add_argument('--bwaScoreThreshold', type=int, help='set minimum score for BWA (see bwa mem -T)')
    p.add_argument('--bwaReseedTrigger', type=float, help='set re-seed trigger BWA (see bwa mem -r)')
    p.add_argument('--allowInvalidBases', action='store_true', help='allow any symbol in reference sequence (e.g. IUPAC ambiguity but also entirely invalid ones)')
    p.add_argument('--tgz', help='put temporary working directory into tgz')
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), default stdin')
    p.add_argument('outfile', nargs='?', help='output file (FASTA), default stdout')
    p.set_defaults(func=runHybpiper)
    

def showArgs(args):
    print args


def paftoolsMain():
    """Entry point for the C{paftools} script.
"""
    logging.basicConfig(format='%(levelname)s: %(module)s:%(lineno)d, %(funcName)s: %(message)s')
    # logger = logging.getLogger(__name__)
    p = argparse.ArgumentParser(description='paftools -- tools for the Plant and Fungal Trees of Life (PAFTOL) project')
    p.add_argument('--loglevel', help='set logging level [DEBUG, INFO, WARNING, ERROR, CRITICAL]')
    subparsers = p.add_subparsers(title='paftools subcommands')
    addDevParser(subparsers)
    addHybseqParser(subparsers)
    addHybpiperParser(subparsers)
    args = p.parse_args()
    if args.loglevel is not None:
        loglevel = getattr(logging, args.loglevel.upper(), None)
        if loglevel is None:
            raise ValueError, 'invalid log level: %s' % args.loglevel
        logging.getLogger().setLevel(loglevel)
    args.func(args)
