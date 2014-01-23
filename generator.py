import numpy as np
import random 

"""
This module contains a generator of random protein sequences.

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
class RandomProtein:
    """
    Random proteins generator

    Methods:

        simulate_protein
        create_sequence
        generate_homolog
    """    
    _alphabet = 'ACEDGFIHKMLNQPSRTWVY'
    _mutation_freq = 0.05
    _indel_freq = 0.01

    def __simulate_protein(self, length, betas):
       AAs = list(self._alphabet)
       alphas = np.random.multinomial(80, betas)
       sequence = ''
       weights = np.random.dirichlet(alphas)
       for i in range(length):
           sequence += np.random.choice(AAs,p=weights)
       return sequence

    def create_sequence(self, length, b1, b2, m_position, m_width):
       sequence = ''
       sequence += self.__simulate_protein(m_position-1, b1)
       sequence += self.__simulate_protein(m_width, b2)
       sequence += self.__simulate_protein(length-(m_position+m_width-1), b1)
       print len(sequence)
       return sequence

    def generate_homolog(self, seq):
        # substitution
        seq = list(seq)
        for i, s in enumerate(seq):
           val = random.random()
           if val < self._mutation_freq:
               # choose a random AA
               seq[i] = random.choice([x for x in self._alphabet if x != s.upper()])
        
        val = random.random()
        if val < self._indel_freq:
            # simulate insertion
            newseq = list()
            ins_ind = random.randint(1, len(seq))
            for i, s in enumerate(seq):
                if i != ins_ind: 
                    newseq.append(s)
                else:
                    ins_len =  random.randint(1, len(seq)*0.2)
                    inseq = self.__simulate_protein(ins_len, [0.05]*20)
                    newseq.append(inseq)
            seq = newseq

        val = random.random()
        if val < self._indel_freq:
            # simulate deletion
            newseq = list()
            del_ind = random.randint(1, len(seq))
            del_len =  random.randint(1, len(seq)*0.1)
            for i, s in enumerate(seq):
                if i not in range(del_ind,del_len): 
                    newseq.append(s)
            seq = newseq
        return ''.join(seq)

    

