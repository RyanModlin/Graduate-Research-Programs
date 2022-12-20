#!/usr/bin/env python3

import sys

class FastAreader:
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''

        self.fname = fname
        self.fileH = None

    def doOpen(self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as self.fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = self.fileH.readline()
            while not line.startswith('>'):
                line = self.fileH.readline()
            header = line[1:].rstrip()

            for line in self.fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence


########################################################################
# CommandLine
########################################################################
class Range():
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )

        self.parser.add_argument('-m', '--max', nargs='?', default=8, action='store',
                                 help='max kMer size ')
        self.parser.add_argument('-c', '--cutoffProb', nargs='?', type=float, default=.01, action='store',
     #                            choices=[Range(0.0, 1.0)],
                                 help='max kMer size ')

        self.parser.add_argument('-p', '--pseudoCount', type=int, default=1, action='store',
                                 help='pseudocount to be added to each tetramer bin')

        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


def rc(s):
    return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))


import itertools, math, random
import scipy.stats as stats


class Genome():
    bases = 'ACGT'

    def __init__(self, min=1, max=8, pseudo=1):
        self.n = 0
        self.counts = {}
        for n in range(min, max + 1):
            for k in itertools.product(Genome.bases, repeat=n):
                kmer = ''.join(k)
                if rc(kmer) in self.counts:
                    self.counts[kmer] = self.counts[rc(kmer)]
                else:
                    self.counts[kmer] = [pseudo]
                    self.n += pseudo
        self.min, self.max, self.pseudo = min, max, pseudo

        self.statsN = [0 for i in range(max+1)]
        self.statsRsum = [0 for i in range (max+1)]
        self.statsRsum2 = [0 for i in range(max + 1)]

    def addSequence(self,seq):
        for l in range(self.min,self.max+1):
            for p in range(len(seq) - self.max):
                try:
                    self.counts[seq[p:p+l]][0] += 1
                except KeyError:
                    pass
        self.n += len(seq) - self.max

    def addSequenceN(self, seq):
        end = len(seq) - self.max # end position
        probes = end // 1
        for s in range(self.min, self.max+1):
            for p in range(probes):
                rp = random.randint(0,end)
                try:
                    self.counts[seq[rp:rp+s]][0] += 1
                except KeyError: # some other base was found - N typically
                    pass
        self.n += probes


    def E(self,s):
        ''' return expected count of s. For even len s, P(ab) = P(a|b) * P(b). P(abc) = P(a|b) * P(b|c) * P(c)'''
        # P(abcd) = P(a) * P(b|a) * P(c|ab) * P(d|bc)
        # E(abcd) = c(a)/len(s) * c(ab)/c(a) * c(abc)/c(ab) * c(bcd)/c(bc) *len(s) = c(abc) * c(bcd) / c(bc)
        # note: sequences 3 or less will return the counts collected

        try:
            return self.counts[s[:-1]][0] * self.counts[s[1:]][0] / self.counts[s[1:-1]][0]
        except KeyError:
            return self.counts[s][0]
        except ZeroDivisionError:
            return 0.


    def E1(self,s):
        ''' return expected count of s. For even len s, P(ab) = P(a|b) * P(b). P(abc) = P(a|b) * P(b|c) * P(c)'''

        try:
            exp = self.counts[s[0:2]][0]
            for p in range (1,len(s)-1):
                exp *= self.counts[s[p:p+2]][0] / self.counts[s[p]][0]
            return exp
        except KeyError:
            return self.counts[s][0]

        except ZeroDivisionError:
            return 0.

    def E1r(self,s):
        ''' return expected count of s. For even len s, P(ab) = P(a|b) * P(b). P(abc) = P(a|b) * P(b|c) * P(c)'''
        rs = s[::-1]
        try:
            exp = self.counts[rs[0:2]][0]
            for p in range (1,len(s)-1):
                exp *= self.counts[rs[p:p+2]][0] / self.counts[rs[p]][0]
            return exp
        except KeyError:
            return self.counts[s][0]

        except ZeroDivisionError:
            return 0.

    def E0(self,s):
        ''' return expected count of s. For even len s, P(ab) = P(a|b) * P(b). P(abc) = P(a|b) * P(b|c) * P(c)'''

        try:
            exp = self.counts[s[0]][0]
            for p in range (1,len(s)):
                exp *= self.counts[s[p]][0] / self.n
            return exp
        except KeyError:
            return self.counts[s][0]
        except ZeroDivisionError:
            return 0.



    def pValue(self,s):

        return stats.binom.cdf (self.counts[s][0], self.n, self.E(s)/self.n)
        return stats.hypergeom.cdf(k=self.counts[s], M=self.n*Ninflate, N = self.n, n = self.E(s) * Ninflate)

    def Zscore(self,s, r):

        k = len(s)
        mu = self.statsRsum[k] / self.statsN[k]
        sd = math.sqrt (
            self.statsRsum2[k] / self.statsN[k] -
            (self.statsRsum[k] / self.statsN[k]) ** 2 )
        if r == mu: return 0
        if sd == 0:
            if r < mu:
                return -10000000
            else:
                return 10000000
        return (r - mu) / sd

    def actRatio (self, s):
        exp = self.E(s)
        if exp:
            v = self.counts[s][0] / exp
        else:
            v = 0
        k = len(s)
        self.statsN[k] += 1
        self.statsRsum[k] += v
        self.statsRsum2[k] += v*v
        return v

########################################################################
# Main
# Here is the main program
#
#
########################################################################


def main(myCommandLine=None):
    cl = CommandLine()

    # get the genome
    sourceReader = FastAreader()

    # setup the Tetramer object using centerSequences as the reference
    thisGenome = Genome(1, cl.args.max, cl.args.pseudoCount)
#    thisGenome.E = thisGenome.Ecomposition
    for head, seq in sourceReader.readFasta():
        thisGenome.addSequence(seq)
    print ('N = {}'.format(thisGenome.n))
    # print (thisGenome.E('CGCG'), thisGenome.Zscore('CGCG'))
    # print (thisGenome.n)
    # for s in ('C', 'G','CG', 'GC', 'CGC', 'GCG', 'CGCG'):
    #     print (s, thisGenome.counts[s])
    # exit()

    outArray = []
    for seq,count in thisGenome.counts.items():
 #       outArray.append( (seq, count, thisGenome.E(seq), thisGenome.pValue(seq)) )
        outArray.append((seq, count[0], thisGenome.E(seq), thisGenome.actRatio(seq)))
    outArray.sort(key=lambda e: (len(e[0]), -e[3]), reverse = True )
    for seq, count, E, ratio in outArray:
        rSeq = rc(seq)
        z = thisGenome.Zscore(seq, ratio)
        if (seq in thisGenome.counts) and (z < cl.args.cutoffProb) and (seq <= rSeq) :
            #print('{0:8}\t{1:0d}\t{2:0.2f}\t{3:0.2f}'.format(seq, count, E, pVal))
            print('{0:8}:{1:8}\t{2:0d}\t{3:0.2f}\t{4:0.2f}'.format(seq, rSeq, count,E, z))
if __name__ == "__main__":
    main()
