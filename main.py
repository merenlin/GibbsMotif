#!/usr/bin/env python
"""
GibbsMotif main module

Copyright (C) {2014}  {Oxana Sachenkova}

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 2 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
"""

ITERATIONS_DEFAULT = 50

import sys
import argparse 
import numpy
from random import choice

from gibbs import GibbsSampler
from generator import RandomProtein
  
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 
from Bio.SeqUtils import ProtParam


def main():
    """
    Main program.
    """
    data = initialize()

    g = GibbsSampler(sequences=data['sequences'],
              motif_width=data['width'])

    g.find_motif(iterations=data['iterations'])

    print_sequences(data['sequences'], data['width'])
    return

def generate_synthetic_data(numseqs, length, m_position, m_width):
    f = open("data/long_disorder.fasta","rU")
    seqs = SeqIO.parse(f,"fasta")
    scores = list()
    for record in seqs: 
        score = ProtParam.ProteinAnalysis(str(record.seq)).get_amino_acids_percent().values()
        scores.append(score)

    betas_background = numpy.mean(numpy.array(scores), axis=0)
    betas_motif = [0.05]*20
   
    f_result = open("data/synthetic_orthologs.fasta", "w+") 
    
    query = RandomProtein().create_sequence(length,betas_background, betas_motif, m_position, m_width)
    query_record = SeqRecord(Seq(query),id="pr0",description="")
    SeqIO.write(query_record, f_result, "fasta")
          
    for i in range(1,numseqs): 
        seq = RandomProtein().generate_homolog(query)
        seqrecord = SeqRecord(Seq(seq),id="pr"+str(i),description="")
        SeqIO.write(seqrecord, f_result, "fasta")

    f_result.close()

def initialize():
    """
    Parse command line options, and read input Fasta file.
    """
    parser = argparse.ArgumentParser(usage = "usage: %prog -i FILE -w WIDTH [-h] "
                          "[options]",
                          version = "Gibbs Motif Finder 0.1",
                          description = "Gibbs Motif Finder is an implementation of the "
                          "Gibbs sampling algorithm for finding "
                          "shared motifs in protein sequences. "
                          )

    parser.add_argument("-i", "--input", dest="input",
                      help="Full pathname of the input file in Fasta format")
    parser.add_argument("-w", "--width", dest="width",
                      type=int, help="Motif width to search for")
    parser.add_argument("-t", "--iterations", dest="iterations",
                       default=ITERATIONS_DEFAULT, type=int,
                       help="Number of non-improving iterations "
                      "(default " + str(ITERATIONS_DEFAULT) + ")")
    parser.add_argument("-s", "--synthetic", dest="synthetic", action='store_true',
                      help="Generate and use synthetic data instead of the input file")
   
    
    options = parser.parse_args()

    if not options.input and not options.synthetic:
        parser.error("You forgot the input file")

    if not options.width:
        parser.error("You forgot the width argument")

    if options.width < 3:
        parser.error("Motif width can't be too short")

    if options.synthetic:
       generate_synthetic_data(100, 100, 5, 5)
       inputfile = "data/synthetic_orthologs.fasta"
    else:
       inputfile = options.input

    try:
        file = open(inputfile, "rU")
    except IOError:
        parser.error("Could not read file %s" % inputfile)
    
    fasta_iterator = SeqIO.parse(file, "fasta") 
    sequences = [{'name':          record.name,
                  'sequence':       record.seq,
                  'motif_position': 0}
                for record in fasta_iterator]

    if len(sequences) < 2:
        parser.error("I found %i sequences in input file %s" % (len(sequences),
                                                              options.input))

    return {'sequences':        sequences,
            'width':            options.width,
            'iterations':       options.iterations
            }

def print_sequences(sequences, motif_width):
    """
    Print the occurrence of the motif in each sequence.
    """
    print "Here is what I found:"
    for i in range(len(sequences)):
        start, end = (sequences[i]['motif_position'],
                      sequences[i]['motif_position'] + motif_width)
        print "Sequence #%2i  %s  (at position %i)" % (
            i + 1,
            sequences[i]['sequence'][start:end],
            sequences[i]['motif_position'] + 1)
    return


if __name__ == "__main__":
    main()
