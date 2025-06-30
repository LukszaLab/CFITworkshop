'''
Created on Dec 14, 2015

@author: mluksza
'''
from cfit.tree.mutation.CodingMutation import CodingMutation


class SynonymousMutation(CodingMutation):
    '''
    Synonymous mutation class
    '''

    def is_synonymous(self):
        '''

        :return: bool
        '''
        return True
