'''
Created on July 30, 2021

@author: mluksza
'''

from cfit.CoreObject import CoreObject


class SampleMutation(CoreObject):
    '''
    SampleMutation class.

    Attributes:
        __sample: str
            cfit.patient.Sample

        __mutation: cfit.tree.mutation.Mutation

        __DP: int
            obsolete

        __AP: int
            obsolete

        __normalDP: int
            obsolete

        __normalAP: int
            obsolete

        __caller: str:
            the algorithm, obsolete

        __numOfCallers: int
            number of callers, obsolete

        __Ncov: int
            number of reads in the normal sample

        __Naf: float
            frequency in the normal sample

        __Taf: float:
            frequency in the tumor sample

        adjustedTaf: float:
            adjusted frequency in the tumor sample

        __X: float
            not sure

        __expression: float
    '''

    def __init__(self, mutation):
        '''
        Constructor, class representing a mutation in a sample (frequency, expression)

        :param mutation: cfit.tree.mutation.Mutation


        '''
        self.__mutation = mutation

        self.__DP = 0
        self.__AP = 0

        self.__normalDP = 100
        self.__normalAP = 0
        self.__Ncov = 0
        self.__Naf = 0
        self.__Taf = 1
        self.adjustedTaf = 0.0
        self.__X = 0.
        self.__expression = 1.0

    def __str__(self):
        return self.__mutation.id

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__mutation.id == other.__mutation.id
        else:
            return False

    @property
    def chrom(self):
        return self.__mutation.chrom

    @property
    def id(self):
        return self.__mutation.id

    @property
    def expression(self):
        return self.__expression

    @expression.setter
    def expression(self, expr):
        self.__expression = expr

    @property
    def Tcov(self):
        tcov = self.__DP  # - self.__AP
        return tcov

    @property
    def Taf(self):
        taf = 0
        if self.__DP != 0:
            taf = 1. - float(self.__AP) / float(self.__DP)
        return taf

    @property
    def Ncov(self):
        ncov = self.__normalDP  # - self.__normalAP
        return ncov

    @property
    def Naf(self):
        naf = 0.0
        if self.__normalDP != 0:
            naf = float(self.Ncov) / float(self.__normalDP)
        return naf

    @property
    def gene(self):
        return self.__mutation.gene

    @property
    def ref(self):
        return self.__mutation.ref

    @property
    def alt(self):
        return self.__mutation.alt

    @property
    def refAA(self):
        return self.__mutation.refAA

    @property
    def altAA(self):
        return self.__mutation.altAA

    @property
    def ENSG(self):
        return self.__mutation.ENSG

    @property
    def ENST(self):
        return self.__mutation.ENST

    @property
    def DP(self):
        return self.__DP

    @DP.setter
    def DP(self, DP):
        self.__DP = DP

    @property
    def substitution(self):
        return self.__mutation.substitution

    @property
    def substitution_with_codon(self):
        return self.__mutation.substitution

    def is_coding(self):
        '''

        :return: bool
        '''
        return self.__mutation.is_coding()

    def is_synonymous(self):
        '''
        :return: bool
        '''
        return self.__mutation.is_synonymous()

    def is_nonsynonymous(self):
        '''

        :return: bool
        '''
        return self.__mutation.is_nonsynonymous()

    def initialize_VCF(self, hdict):
        '''
        Initializes from a VCF file entry

        :param hdict: dict

        :return:
        '''

        self.__sample = hdict["Sample"]
        try:
            [self.__normalDP, self.__normalAP] = [int(i) for i in hdict["NORMAL"].split(":")]
            [self.__DP, self.__AP] = [int(i) for i in hdict["TUMOR"].split(":")]
        except:
            self.logger("VCF file in non-DP:AP format.", 0)

    def toVCF(self):
        '''
        Return a line of a vcf-like DataFrame for the mutation

        :return:
        '''

        js = {}
        js["#CHROM"] = "chr" + str(self.chrom)
        js["INFO"] = self.gene
        js["POS"] = self.pos
        js["ID"] = self.id
        js["REF"] = self.ref
        js["ALT"] = self.alt
        js["Sample"] = self.__sample
        js["NORMAL"] = ":".join(map(str, [self.__normalDP, self.__normalAP]))
        js["TUMOR"] = ":".join(map(str, [self.__DP, self.__AP]))
        js["QUAL"] = "."
        js["FILTER"] = "PASS"
        js["FORMAT"] = "DP:AP"
        cols = ["#CHROM", "POS", "ID", "REF", "ALT",
                "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR"]
        line = [js[col] for col in cols]
        return line
