'''
Created on Apr 10, 2015

@author: mluksza
'''
# from Bio.SubsMat import MatrixInfo as matlist
from Bio.pairwise2 import align

from cfit.CoreObject import CoreObject


class AlignedEpitope(CoreObject):
    '''
    Class for an epitope aligned to a neoantigen.

    Attributes:

        wtScore: float
        mutScore: float
        wtSeq: str
        mutSeq: str

        mtfrom: int
        mtto: int
        mtaln: str

        wtfrom: int
        wtto: int

        epitopeSeq: str
        nid: str
        wtPeptide: str
        mutPeptide: str
        epitope: cfit.fitness.Epitope
        epitopeName: str
        wtScoreSW: float
        mutScoreSW: float
        wtSW: str
        mutSW: str
        normalizedScoreSW: float
        mutEpitopeSW: str

        matrix: dict
            default matlist.blosum62,
            {('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ...}

        gap_open: float
            alignment parameters
        gap_extend: float
            alignment parameters


    '''

    @staticmethod
    def align(seq1, seq2, matrix):
        #        matrix = matlist.blosum62
        gap_open = -11
        gap_extend = -1
        aln = align.localds(seq1.upper(), seq2.upper(), matrix, gap_open, gap_extend)
        return aln[0]

    def __init__(self, epitope, neoantigen, substmatrix=None, mt_score=0, wt_score=0, recompute_scores=True):
        '''

        :param epitope: cfit.fitness.Epitope

        :param neoantigen: cfit.tree.mutation.Neoantigen

        :param substmatrix: dict

        :param mt_score: float

        :param wt_score: float

        :param recompute_scores: bool
        '''

        self.wtScore = wt_score
        self.mutScore = mt_score
        self.wtSeq = ""
        self.mutSeq = ""

        self.mtfrom = -1
        self.mtto = -1
        self.mtaln = ""

        self.wtfrom = -1
        self.wtto = -1

        self.epitopeSeq = epitope.seq.upper()
        self.nid = neoantigen.id
        self.wtPeptide = neoantigen.wtPeptide.upper()
        self.mutPeptide = neoantigen.mtPeptide.upper()
        self.epitope = epitope
        self.epitopeName = epitope.get_eheader()
        self.wtScoreSW = wt_score
        self.mutScoreSW = mt_score
        self.wtSW = "-"
        self.mutSW = "-"
        self.normalizedScoreSW = 0
        self.mutEpitopeSW = "-"

        if recompute_scores:
            gap_open = -11
            gap_extend = -1
            if len(self.wtPeptide) >= 8:
                alnWt = align.localds(self.wtPeptide.upper(), self.epitopeSeq.upper(), substmatrix, gap_open,
                                      gap_extend)
                if len(alnWt) > 0:
                    alnWt = alnWt[0]
                    self.wtScoreSW = alnWt[2]
                    self.wtfrom = alnWt[3]
                    self.wtto = alnWt[4]
                    self.wtaln = alnWt[0][self.wtfrom:self.wtto]
                    self.wtSW = alnWt[0].replace("-", "")

            alnMutSelf = align.localds(self.mutPeptide.upper(), self.mutPeptide.upper(), substmatrix, gap_open,
                                       gap_extend)
            self.selfScoreSW = alnMutSelf[0][2]

            mtseq = self.mutPeptide.upper()
            alnMut = align.localds(mtseq, self.epitopeSeq.upper(), substmatrix, gap_open, gap_extend)
            if len(alnMut) > 0:
                alnMut = alnMut[0]
                mtfrom = alnMut[3]
                mtto = alnMut[4]
                self.mtaln = alnMut[0][mtfrom:mtto]
                self.mutSW = alnMut[0].replace("-", "")
                pos = self.mutSW.find(self.mtaln)
                self.mtfrom = pos + 1
                self.mtto = pos + len(self.mtaln)
                self.mutEpitopeSW = alnMut[1].replace("-", "")
                self.mutScoreSW = alnMut[2]
                self.normalizedScoreSW = self.mutScore / self.selfScoreSW

    def get_epitope_id(self):
        '''

        :return: int
        '''
        return self.epitope.id

    def set_mutant_score(self, score, seq):
        '''

        :param score: float
        :param seq: str
        :return:
        '''
        self.mutScore = score
        self.mutSeq = seq

    def set_wildtype_score(self, score, seq):
        '''

        :param score: float
        :param seq: str
        :return:
        '''
        self.wtScore = score
        self.wtSeq = seq

    def get_score(self, wt=False):
        '''

        :param wt: bool

        :return: float
        '''
        if wt:
            return self.wtScoreSW
        else:
            return self.mutScoreSW
