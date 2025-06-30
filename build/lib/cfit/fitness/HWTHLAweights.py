'''
Created on Jun 20, 2017

@author: mluksza
'''
from cfit.fitness.HLAweights import HLAweights


class HWTHLAweights(HLAweights):
    '''
    Neoantigen masking for the IO cohorts
    - positions 2&9 from a hydrophobic amino-acid
    '''

    def get_weight(self, neo):
        '''

        :param neo: cfit.tree.mutation.Neoantigen
        :return:
        '''
        w = HLAweights.get_weight(self, neo)
        if neo.residueChange[0] != "H" and (neo.position == 2 or neo.position == 9):
            w = 0
        return w
