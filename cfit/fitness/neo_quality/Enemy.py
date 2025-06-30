from cfit.fitness.neo_quality.EpitopeDistance import EpitopeDistance
import vptree

from Bio import SeqIO

class Enemy(CoreObject):

    '''
    Class for an Foreign proteome.

    Attributes:

        edist: EpitopeDistance

        peptides: list of 9-mers

        vptrees: dict: int->vptree

    '''

    def __init__(self, enemydb, model_name="all_tcr_all_combos_model", ns=[9]):
        '''
        :param model_name: str
            name of the metric

        :param enemydb: str

        :param ns: list of integers
        '''

        self.edist = EpitopeDistance(model_name=model_name)
        nseqs = {}
        for n in ns:
            nseqs[n] = set()
        with open(enemydb) as f:
            seqs = SeqIO.parse(f, "fasta")
            for seq in seqs:
                for n in ns:
                    for i in range(len(seq)-n):
                        nseq = seq[i:(i+n)]
                        nseqs[n].add(nseq)

        self.vptrees = {}
        for n in ns:
            points = list(nseqs[n])
            self.vptrees[n] = vptree.VPTree(points, edist.epitope_dist)



