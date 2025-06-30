'''
Created on Aug 18, 2016

@author: mluksza
'''

import numpy as np
from Bio import SeqIO

from cfit.CoreObject import CoreObject


class EpitopeMask(CoreObject):
    '''
    Class implements masking of epitopes so that they are omitted from analysis

    Attributes:
        epitopemask: dict
            Dictionary mapping epitope identifiers to their mask value

        epitopes: dict
            Dictionary mapping epitope identifiers to cfit.fitness.Epitope objects.
    '''

    def __init__(self, params):
        '''
        Constructor
        :param params: list
            [typ, val]: str, str, eg. ["file", <path_to_the_file>]
            else initializes empty mask


        '''
        self.epitopemask = {}
        self.epitopes = {}

        [typ, val] = params
        if typ == "file":
            values = set()
            maskfile = val
            f = open(maskfile)
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                [epi, val] = line.strip().split("\t")
                val = float(val)
                values.add(val)
                self.epitopemask[epi] = val
                if len(values) > 2:
                    vals = filter(lambda v: v > 0, self.epitopemask.values())
                    w = 1
                    if len(vals) > 0:
                        w = len(vals) / sum(vals)
                    for epi in self.epitopemask:
                        self.epitopemask[epi] *= w
        else:
            self.epitopemask = {}

    def set_epitopes(self, epitopes):
        '''
        Sets the epitopes
        :param epitopes: dict
        :return:
        '''
        self.epitopes = epitopes

    def get_epitope_ids(self):
        '''
        Return the list of epitope identifiers.
        :return: list
        '''
        return list(self.epitopemask.keys())

    def mask_value(self, epi):
        '''
        Return the value of masking for a given epitope identifier.

        :param epi: int

        :return: int
        '''
        try:
            return self.epitopemask[epi]
        except:
            return 0

    def get_N(self):
        '''
        Returns the number of epitopes
        :return: int
        '''
        return sum(self.epitopemask.values())

    def random_masking(self, remainingFraction):
        '''
        Masks random epitopes.

        :param remainingFraction: float
            the fraction of epitopes not to be masked.
        :return:
        '''

        for epi in self.epitopemask:
            self.epitopemask[epi] = 1
            u = np.random.uniform()
            if u > remainingFraction:
                self.epitopemask[epi] = 0

    def reset(self):
        '''
        Unmasks all epitopes.

        :return:
        '''
        for epi in self.epitopemask:
            self.epitopemask[epi] = 1

    def reset_0(self):
        '''
        Masks all epitopes.
        :return:
        '''
        for epi in self.epitopemask:
            self.epitopemask[epi] = 0

    def set_mask_value(self, epi, val):
        '''
        Set the value of the mask for a given epitope.

        :param epi: int

        :param val: int

        :return:
        '''
        self.epitopemask[epi] = val

    def multiply(self, emask2):
        '''
        Multiple with another EpitopeMask object. Will mask the union of epitopes, mask in either of the two.

        :param emask2: cfit.fitness.EpitopeMask
        :return:
        '''

        epis = emask2.get_epitope_ids()
        for epi in epis:
            if epi in self.epitopemask:
                self.epitopemask[epi] * emask2.mask_value(epi)

    def mask_sequences(self, fastafile, i=-1):
        '''

        :param fastafile: str
        :param i: int
        :return:
        '''
        self.reset()
        f = open(fastafile)
        seqs = SeqIO.parse(f, "fasta")
        sseqs = set()
        j = 1
        for seq in seqs:
            if i == -1 or i == j:
                sseqs.add(str(seq.seq))
            j += 1
        f.close()
        for epi in self.epitopemask:
            epitope = self.epitopes[epi]
            if epitope.seq in sseqs:
                self.set_mask_value(epi, 0)

    def mask_epitope(self, eid):
        '''
        Mask epitope.

        :param eid: int
            identifier of the epitope to be masked

        :return:
        '''
        self.epitopemask[eid] = 0.0

    def mask_epitopes(self, eids):
        '''
        Mask epitopes.

        :param eids: list
            list of identifiers of epitopes to be masked
        :return:
        '''

        for eid in eids:
            self.epitopemask[eid] = 0.0
