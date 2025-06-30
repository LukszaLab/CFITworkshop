'''
Created on Mar 30, 2015

@author: mluksza
'''

from cfit.CoreObject import CoreObject
from cfit.util.Utils import Utils


class Neoantigen(CoreObject):
    '''
    Class implements neantigens

    Attributes:
        __id: int
            unique neontigen identifier

        __peptideid: str
            neoantigen identifier of the form
            <chrom>_<pos>_<ref>_<alt>_<mutated_position_in_the_peptide>

        __mid: str:
            mutation identifier of the form <chrom>_<pos>_<ref>_<alt>

        sample: str
            sample name

        wtPeptide: str
            wildtype peptide

        mtPeptide: str
            mutant peptide

        __peptide_length: int
            length of the neoantigen peptide

        position: int
            neoantigen position at which the peptide is mutated

        allele: str
            MHC allele

        HLA: str
            obsolete

        chopscore: int
            1 by default, if 0 neoantigen will not be used

        gene: str
            name of the mutated gene

        mutation: cfit.tree.mutation.Mutation
            the underlying mutation

        qattributres: attributes related to the quality model

        kD: float
            the dissociation constant

        wtkD: float
            the dissociation constant of the wildtype peptide

        quality: float
            neoantigen quality

    '''

    L = 1.  # concentration
    M = 1.  # mutant peptide concentration
    W = 1.  # wildtype peptide concentration
    WEPS = 0.0003
    WTCAP = Utils.INF

    def __init__(self, params):
        '''
        Constructor
        '''
        pparams = params
        if len(params) == 9:
            pparams.append("1")
        [nid, peptideid, mid, sample, wtPeptide, mtPeptide, allele, wtScore, mtScore, HLA, chopscore] = params
        self.__id = nid  # int(nid)  # unique neontigen identifier -
        self.__peptideid = peptideid
        self.__mid = None
        self.__set_mid(mid)
        self.sample = sample  # name of the sample
        self.wtPeptide = wtPeptide  # wildtype peptiden aa sequence
        self.mtPeptide = mtPeptide  # mutant peptide aa sequence
        self.__peptide_length = len(self.mtPeptide)
        try:
            [res1, res2] = [el for el in zip(self.wtPeptide, self.mtPeptide) if el[0] != el[1]][0]
            self.residueChange = Utils.residueChangeClass(res1, res2)  # HH,NN,HN,NH,...
            self.position = [i for i in range(0, self.__peptide_length) if self.mtPeptide[i] != self.wtPeptide[i]]
            self.position = self.position[0] + 1  # neoantigen position at which the peptide is mutated
        except:
            self.residueChange = -1
            self.position = -1

        self.allele = allele
        self.HLA = HLA
        self.chopscore = int(chopscore)  # 1 by default, if 0 neoantigen will not be used
        self.gene = ""
        self.quality = 0.0

        '''
        setting the wildtype and mutant dissociation constants and the MHC amplitude A
        '''
        try:
            self.kD = float(mtScore)
            self.wtkD = min(Neoantigen.WTCAP, float(wtScore))
        except:
            self.kD = Utils.INF
            self.wtkD = Utils.INF

    def __str__(self):
        return self.id

    def __hash__(self):
        return hash(str(self.id))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.id == other.id
        else:
            return False

    @property
    def mid(self):
        return self.__mid

    @mid.setter
    def mid(self, mid):
        mid = mid.split("_")
        mid = "_".join(map(str, mid[:4]))
        self.__mid = mid

    @property
    def peptide_id(self):
        return self.__peptideid

    @property
    def substitution(self):
        ref = self.wtPeptide[self.position - 1]
        alt = self.mtPeptide[self.position - 1]
        return ref + "->" + alt

    def __set_mid(self, mid):
        mid = mid.split("_")
        mid = "_".join(map(str, mid[:4]))
        self.__mid = mid

    @property
    def id(self):
        return self.__id

    @id.setter
    def id(self, id):
        self.__id = id

    def set_gene(self, gene):
        self.gene = gene

    def presentation_score(self, kd0=500., strict=False):
        '''
        Presentation score for the neoantigen

        :param kd0: float
            threshold on neoantigen dissociation constant

        :param strict: bool
            whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :return: float
        '''

        if strict:
            pr_score = int(self.kD <= kd0)
        else:
            pr_score = kd0 / (self.kD + kd0)

        return pr_score

    def set_sample_name(self, sampleName):
        '''

        :param sampleName: str
        :return:
        '''
        self.sample = sampleName

    def compute_quality(self, QModel):
        q = QModel.quality(self)
        self.quality = q

    def get_load(self):
        '''

        :return: float
        '''
        return self.load

    def toJSON(self):
        '''
        Creates json for Sibyl app.
        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''
        js = {}
        js['id'] = self.__id
        js['mutation_id'] = self.mid
        js['HLA_gene_id'] = self.allele
        js['sequence'] = self.mtPeptide
        js['WT_sequence'] = self.wtPeptide
        js['mutated_position'] = self.position
        js['Kd'] = self.kD
        js['KdWT'] = self.wtkD
#        js['score'] = self.R * self.A
        return js
