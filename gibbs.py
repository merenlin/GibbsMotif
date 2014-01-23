"""
This module contains an implementation of the Gibbs sampling algorithm used
for finding shared motifs in protein sequences.

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
import random
import numpy as np
import math
import csv

class GibbsSampler:
    """
    An implementation of the Gibbs sampling algorithm.

    Methods:

        find_motif  start algorithm on provided sequences
    """
    __sequences = []
    __motif_width = 0
    __pseudocounts = {}
    __alphabet = "ACEDGFIHKMLNQPSRTWVY"

    def __init__(self, sequences, motif_width, pseudocounts_weight=0.05):
        self.__sequences = sequences
        self.__motif_width = motif_width
        self.__compute_pseudocounts(pseudocounts_weight)
        return


    def find_motif(self, iterations=50, initial_num_occurrences=0, initial_pattern_width=0):
        """
        Find a motif in sequences using Gibbs sampling.
        """
        if initial_num_occurrences > 0:
            if initial_pattern_width == 0:
                initial_pattern_width = self.__motif_width
            self.__initial_motif_positions(initial_num_occurrences,
                                           initial_pattern_width)
        else:
            self.__random_motif_positions()

        best_entropy = 0
        best_motif = [0] * len(self.__sequences)

        entropies = []
        i = j = 0
        numimprov = 0
        while i < iterations:
            i += 1
            j += 1

            # Pick a random sequence
            current_sequence = random.choice(self.__sequences)
            sequences_left = filter(
                lambda s: s != current_sequence,
                self.__sequences)

            motif = self.__compute_motif(sequences_left)
            entropies.append(self.__compute_entropy(motif))
            self.__compute_position(motif, current_sequence)
            
            # Save current iterations positions for analysis
            self.__save_positions()
            
            entropy = self.__compute_entropy(motif)
            if entropy > best_entropy:
                best_entropy = entropy
                best_motif = [s['motif_position'] for s in self.__sequences]
                i = 0

        # Restore to best entropy
        for i in range(len(self.__sequences)):
            self.__sequences[i]['motif_position'] = best_motif[i]
        motif = self.__compute_motif(self.__sequences)
        entropy = self.__compute_entropy(motif)

        # Save/Print the final result
        self.__save_entropies(entropies)
        self.__print_motif_weights(motif)
        return

    def __compute_pseudocounts(self, weight):
        """
        compute for each base a weighted pseudocount. Store results in
        __pseudocounts dictionary.

        Parameters:

            weight  weight to use for pseudocounts
        """
        # Total number of bases in sequences
        total = float(sum([len(s['sequence']) for s in self.__sequences]))

        for base in self.__alphabet:
            # Number of occurrences of this base in sequences
            n = sum([s['sequence'].count(base) for s in self.__sequences])
            self.__pseudocounts[base] = weight * (n / total)

        return

    def __random_motif_positions(self):
        """
        Populate the list of sequences with a random position of the motif for
        each sequence.
        """
        for s in self.__sequences:
            s['motif_position'] = random.randint(
                0, len(s['sequence']) - self.__motif_width)
        return

    def __initial_motif_positions(self, number_of_occurrences, pattern_width):
        """
        Populate the list of sequences with a random position of the motif for
        each sequence based on a heuristic.
        """
        # Get base with lowest frequency
        lowest_freq = self.__pseudocounts['A']
        lowest_base = 'A'
        for base in self.__alphabet:
            if self.__pseudocounts[base] < lowest_freq:
                lowest_freq = self.__pseudocounts[base]
                lowest_base = base

        for s in self.__sequences:
            positions = []
            occurrences = number_of_occurrences
            while occurrences > 0:
                # For every word of length pattern_width in sequence
                for r in range(len(s['sequence']) - pattern_width + 1):
                    # Count occurrences of lowest_base in word
                    n = s['sequence'][r : r+pattern_width].count(lowest_base)
                    if n >= occurrences:
                        positions.append(r)
                if len(positions) > 0:
                    break
                # Try the same with one occurrence of lowest_base less
                occurrences -= 1
            if len(positions) > 0:
                position = random.choice(positions)
                position += pattern_width / 2
                position -= self.__motif_width / 2
                if position < 0:
                    position = 0
                if position > len(s['sequence']) - self.__motif_width:
                    position = len(s['sequence']) - self.__motif_width
            else:
                # Choose a random position if the sequence doesn't 
                # have at least one occurrence of lowest_base
                position = random.randint(
                    0, len(s['sequence']) - self.__motif_width)
            s['motif_position'] = position
        return

    def __compute_motif(self, sequences):
        """
        compute the position weight matrix of the motif for the given
        sequences and their alignments.

        Return the matrix as a dictionary with for each key in __alphabet a list of
        length motif_width with position weights.
        """
        # Calculate pseudocounts for all positions in the motif
        q = {}
        for base in self.__alphabet:
            q[base] = [self.__pseudocounts[base]] * self.__motif_width
        pseudocounts_total = sum(self.__pseudocounts.values())
        
        # Calculate c_ij
        for i in self.__alphabet:
            for j in range(self.__motif_width):
                # For each sequence, add 1 if it has base i at position j
                for s in sequences:
                    position = s['motif_position'] + j
                    if s['sequence'][position] == i:
                        q[i][j] += 1
                # Divide by the number of sequences and pseudocounts sum
                q[i][j] /= float(len(sequences)) + pseudocounts_total

        return q


    def __compute_position(self, motif, sequence):
        """
        compute new position of motif in sequence based on motif matrix.
        """
        scores = []

        # For every word of length motif_width in sequence
        for r in range(len(sequence['sequence']) - self.__motif_width + 1):
            p_motif = p_background = 1
            for x in range(self.__motif_width):
                p_motif *= motif[ sequence['sequence'][r+x] ][x]
                p_background *= self.__pseudocounts[sequence['sequence'][r+x]]    
            scores.append(p_motif / p_background)

        # Normalize distribution
        total = sum(scores)
        distribution = map(lambda n: float(n) / total, scores)
        sequence['motif_position'] = np.random.choice(len(scores),p=distribution)

        return

    def __compute_entropy(self, motif):
        """
        compute the relative entropy of the given motif matrix.
        """
        entropy = 0
        for base in self.__alphabet:
            for i in range(self.__motif_width):
                # This is called F in the paper
                entropy += motif[base][i] * math.log(
                    motif[base][i] / self.__pseudocounts[base], 2) if self.__pseudocounts[base] else 0
        return entropy

    def __save_positions(self):
        positions = list()
        for s in self.__sequences:
            positions.append(s['motif_position'])
        with open('results/positions.csv', 'a+') as f:
            writer = csv.writer(f)
            writer.writerow(positions)

    def __save_entropies(self, entropies):
        with open('results/entropies.csv', 'w+') as f:
            writer = csv.writer(f)
            writer.writerow(entropies)
        return

    def __print_motif_weights(self, motif):
        # Print position weight matrix for motif
        for base in self.__alphabet:
            print base,
            for weight in motif[base]:
                print "%1.2f" % weight,
            print

        # Print significant motif positions
        print '',
        for i in range(self.__motif_width):
            if max([motif[base][i] for base in self.__alphabet]) > .6:
                print "   *",
            else:
                print "    ",
        print
        return

    def __print_sequences(self):
        # Print motif occurrence in each sequence and motif position
        for s in self.__sequences:
            start, end = s['motif_position'], s['motif_position'] + self.__motif_width
            print "%s motif at position %s" % (s['sequence'][start:end], s['motif_position'])
