'''
Created on July 12, 2021

@author: mluksza
'''

from cfit.tree.mutation.Mutation import Mutation


class CodingMutation(Mutation):
    '''
    Nonsynonymous mutation class

    Attributes:

    geneCodon: int

    __expression: float
    '''

    def __init__(self):
        Mutation.__init__(self)
        self.geneCodon = -1

    def is_coding(self):
        '''

        :return: bool
        '''
        return True

    def initialize_VCF(self, hdict):
        '''
        Initializes from a VCF file entry

        :param hdict: dict

        :return:
        '''

        Mutation.initialize_VCF(self, hdict)
        elems = hdict["INFO"].split("|")
        elems = [el for el in elems if el[:2] == "p."]

        self.refAA = ""
        self.altAA = ""
        self.geneCodon = -1

        if len(elems) > 0:
            el = elems[0]
            el = el.replace("p.", "")
            # self.logger(el)
            geneCodon = "".join([x for x in el if x.isnumeric()])
            try:
                [refAA, altAA] = el.split(geneCodon)
                geneCodon = int(geneCodon)
                self.refAA = refAA
                self.altAA = altAA
                self.geneCodon = geneCodon
            except:
                print(self.__class__.__name__)
                print("error!")
                print(hdict["INFO"])
                print(elems)

    @property
    def substitution(self):
        return self.refAA + "->" + self.altAA

    @property
    def substitution_with_codon(self):
        return self.refAA + str(self.geneCodon) + self.altAA
