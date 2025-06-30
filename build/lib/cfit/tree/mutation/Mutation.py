'''
Created on Dec 14, 2015

@author: mluksza
'''

from cfit.CoreObject import CoreObject


class Mutation(CoreObject):
    '''
    Mutation class.

    Attributes:

        __chrom: str
            chromosome

        __pos: int
            chromosome position

        __gene: str
            gene name

        __ENSG: str
            ensembl gene identifier

        __ENST: str
            ensembl transcript identifier

        __ref: str
            reference nucleotide

        refAA: str

        altAA: str

        __alt: str
            alternative nucleotide

        __caller: str:
            the algorithm, obsolete

        __numOfCallers: int
            number of callers, obsolete

        __impact: str:
            snpeff annotation

        __effect: str
            snpeff annotation

        __inFinalSet: bool
            obsolete

        __functionalClass: str
            to be checked

    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.__chrom = ""
        self.__pos = 0
        self.__gene = ""
        self.__ENSG = ""
        self.__ENST = ""
        self.__ref = ""
        self.__alt = ""
        self.refAA = ""
        self.altAA = ""
        self.vaf = 0

        self.__caller = None
        self.__numOfCallers = 0
        self.__impact = None
        self.__effect = None
        self.__inFinalSet = None
        self.__functionalClass = None

    def __str__(self):
        return self.id

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.id == other.id
        else:
            return False

    @property
    def chrom(self):
        return self.__chrom

    @property
    def pos(self):
        return self.__pos

    @property
    def id(self):
        __id = "_".join(map(str, [self.__chrom, self.__pos, self.__ref, self.__alt]))
        __id = __id.upper()
        return __id

    @property
    def gene(self):
        return self.__gene

    @property
    def ref(self):
        return self.__ref

    @property
    def alt(self):
        return self.__alt

    @property
    def ENSG(self):
        return self.__ENSG

    @property
    def ENST(self):
        return self.__ENST

    @gene.setter
    def gene(self, gene):
        self.__gene = gene

    @property
    def substitution(self):
        return ""

    @property
    def substitution_with_codon(self):
        return ""

    def initialize_VCF(self, hdict):
        '''
        Initializes from a VCF file entry

        :param hdict: dict

        :return:
        '''

        self.__chrom = hdict["#CHROM"]  # [2:]
        self.__chrom = self.__chrom.lower().replace("chr", "")
        self.__gene = hdict["INFO"]
        self.__ENSG = ""
        self.__ENST = ""
        self.INFO = hdict["INFO"]
        if len(self.__gene.split("|")) > 4:
            self.__gene = self.gene.split("|")[3]
            self.__ENSG = hdict["INFO"].split("|")[4]
            self.__ENST = hdict["INFO"].split("|")[5]
        self.__pos = int(hdict["POS"])
        self.__ref = hdict["REF"]
        self.__alt = hdict["ALT"]
        self.__sample = hdict["Sample"]
        try:
            [self.__normalDP, self.__normalAP] = [int(i) for i in hdict["NORMAL"].split(":")]
            [self.__DP, self.__AP] = [int(i) for i in hdict["TUMOR"].split(":")]
            self.vaf = (self.__DP-self.__AP)/self.__DP
        except:
            self.logger("VCF file in non-DP:AP format.", 0)

    def is_synonymous(self):
        '''
        :return: bool
        '''
        return False

    def is_nonsynonymous(self):
        '''

        :return: bool
        '''
        return False

    def is_missense(self):
        '''

        :return: bool
        '''
        return False

    def is_coding(self):
        '''

        :return: bool
        '''
        return False

    def toVCF(self):
        '''
        Return a line of a vcf-like DataFrame for the mutation

        :return:
        '''

        js = {}
        js["#CHROM"] = "chr" + str(self.__chrom)
        js["INFO"] = self.__gene
        js["POS"] = self.__pos
        js["ID"] = self.id
        js["REF"] = self.__ref
        js["ALT"] = self.__alt
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

    def toJSON(self):
        '''
        Creates json for Sibyl app.

        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''
        js = {'id': self.id,
#              'chr': self.__chrom,
#              'position': self.__pos,
              'gene': self.__gene,
              'missense': int(self.is_nonsynonymous())}
#              'ref_nl': self.__ref,
#              'alt_nl': self.__alt}
        return js
