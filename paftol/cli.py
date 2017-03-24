import sys
import argparse
import logging

import paftol


def runHybseq(argNamespace):
    hybseqAnalyser = paftol.HybseqAnalyser(argNamespace.targetsfile, argNamespace.outfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.tgz)
    sys.stderr.write('%s\n' % str(hybseqAnalyser))


def runHybpiper(argNamespace):
    hybpiperAnalyser = paftol.HybpiperAnalyser(argNamespace.targetsfile, argNamespace.outfile, argNamespace.forwardreads, argNamespace.reversereads, argNamespace.tgz)
    if argNamespace.bwaMinSeedLength is not None:
        hybpiperAnalyser.bwaMinSeedLength = argNamespace.bwaMinSeedLength
    if argNamespace.bwaScoreThreshold is not None:
        hybpiperAnalyser.bwaScoreThreshold = argNamespace.bwaScoreThreshold
    if argNamespace.bwaReseedTrigger is not None:
        hybpiperAnalyser.bwaReseedTrigger = argNamespace.bwaReseedTrigger
    hybpiperAnalyser.keepTmpDir = True
    hybpiperAnalyser.analyse()


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
    p.add_argument('--tgz', help='put temporary working directory into tgz')
    p.add_argument('targetsfile', nargs='?', help='target sequences (FASTA), default stdin')
    p.add_argument('outfile', nargs='?', help='output file (FASTA), default stdout')
    p.set_defaults(func=runHybpiper)
    

def showArgs(args):
    print args


def paftoolsMain():
    logging.basicConfig(format='%(levelname)s: %(funcName)s: %(message)s')
    # logger = logging.getLogger(__name__)
    p = argparse.ArgumentParser(description='paftools -- tools for the Plant and Fungal Trees of Life (PAFTOL) project')
    p.add_argument('--loglevel', help='set logging level [DEBUG, INFO, WARNING, EROR, CRITICAL]')
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
