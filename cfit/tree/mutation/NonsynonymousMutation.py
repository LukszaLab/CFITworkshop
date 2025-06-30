'''
Created on Dec 14, 2015

@author: mluksza
'''
from cfit.tree.mutation.CodingMutation import CodingMutation


class NonsynonymousMutation(CodingMutation):
    '''
    Nonsynonymous mutation class

    Attributes:

    refAA: str

    altAA: str


    '''

    def is_nonsynonymous(self):
        '''

        :return: bool
        '''
        return True

    def is_missense(self):
        '''

        :return: bool
        '''
        return True
