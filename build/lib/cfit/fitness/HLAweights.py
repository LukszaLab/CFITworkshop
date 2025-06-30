'''
Created on Apr 13, 2017

@author: mluksza
'''

import os
from collections import defaultdict

from cfit.CoreObject import CoreObject


class HLAweights(CoreObject):
    '''
    classdocs
    '''

    def __init__(self, params):
        '''
        Constructor
        '''

        psetmap = {"1-9": set([1, 2, 3, 4, 5, 6, 7, 8, 9]),
                   "1-9H": set([1, 2, 3, 4, 5, 6, 7, 8, 9]),
                   "3-8": set([3, 4, 5, 6, 7, 8])}

        # 1-9: all neoantigens
        # 1-9H: all neoantigens, excluding those with
        # a nonhydrophobic wildtype residues on positions 2&9 
        # 3-8: neoantigens mutated positions between anchor sites (2&9)

        self.w = {"all": defaultdict(lambda: 1)}

        if len(params) == 9:
            self.w["all"] = defaultdict(lambda: 1)  # [0 for _ in range(0,10)]
            for pos, p in enumerate(params):
                self.w["all"][pos] = p
        else:
            if os.path.exists(params):
                f = open(params)
                lines = f.readlines()
                lines = [line for line in lines if len(line) > 0]
                for line in lines:
                    tab = line.strip().split()
                    hla = tab[0]
                    w = [0] + list(map(float, tab[1:]))
                    self.w[hla] = w
                f.close()
            else:
                if params in psetmap:
                    self.w["all"] = defaultdict(lambda: 0)
                    for pos in psetmap[params]:
                        self.w["all"][pos] = 1

    def get_weight(self, neo):
        hla = neo.allele
        pos = neo.position
        if hla in self.w:
            return self.w[hla][pos]
        else:
            return self.w["all"][pos]
