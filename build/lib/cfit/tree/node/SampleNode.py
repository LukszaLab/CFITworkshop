'''
Created on July 31, 2021

@author: mluksza
'''
import numpy as np
import pandas as pd
from collections import  defaultdict
from cfit.CoreObject import CoreObject


class SampleNode(CoreObject):
    '''
    Class for representing a clone in a sample

    Attributes:
        __node: cfit.tree.Node
            The node on the tree

        __X: float
            inclusive frequency of the clone (including subclones)

        __Y: float
            exclusive frequency of the clone (excluding subclones)

        __sharedY: float
            exclusive frequency of the clone (excluding subclones) of the part of the clone shared with another tumor (set externally)

        __privateY: float
            exclusive frequency of the clone (excluding subclones) of the part of the clone shared with another tumor (set externally)

        __shared_clone: int
            1 or 0

        __preservedY: float
            exclusive frequency of the clone (excluding subclones) of the part of the clone shared with another tumor (set externally)

        __lostY: float
            exclusive frequency of the clone (excluding subclones) of the part of the clone shared with another tumor (set externally)

        __predictedY: float
            predicted exclusive frequency

        __predictedX: float
            predicted inclusive frequency

        nload: float
            neoantigen load TODO: check if it belongs here or to Node

        fattributes: dict TODO: check if it is used
            fitness model specific attributes

        fitness_components: dict
            maps fitness component names to float fitness values

        fitness: float
            absolute fitness - sum over compoenents

        rfitness: float
            normalized relative fitness to predict frequency dynamics

        Y2: float
            predicted exclusive frequency

        X2: float
            predicted inclusive frequency

        neoantigen_immunogenicities: dict (str -> float)
            maps neoantigen identifiers to their externally estimated immunogenicities.
            The computation can for example be done at the tumor level, taking into account competition
            with other neoantigens


    '''

    def __init__(self, node):
        '''
        Constructor

        :param node: cfit.tree.node.Node

        '''

        self.__node = node

        self.__X = 0.0  # inclusive frequency of the clone (including subclones)
        self.__Y = 0.0  # exclusive frequency of the clone (excluding subclones)
        self.__sharedY = 0.0  # exclusive frequency of the clone (excluding subclones) of the part of the clone shared with another tumor (set externally)
        self.__privateY = 0.0  # exclusive frequency of the clone (excluding subclones) of the part of the clone shared with another tumor (set externally)
        self.__preservedY = 0.0  # exclusive frequency of the clone (excluding subclones) of the part of the clone shared with another tumor (set externally)
        self.__lostY = 0.0  # exclusive frequency of the clone (excluding subclones) of the part of the clone shared with another tumor (set externally)
        self.__predictedY = 0.0  # predicted exclusive frequency of the clone (excluding subclones)
        self.__predictedX = 0.0  # predicted inclusive frequency of the clone (including subclones)
        self.__shared_clone = 0
        self.fattributes = {}  # fitness model specific attrubutes
        self.fitness_components = {}
        self.neoantigen_immunogenicities = {} #externally set neoantigen immunogenicities
        self.fitness = 0.0  # absolute fintess
        self.rfitness = 0.0  # relative fitness

        self.nload = 0.0  # neoantigen load

        self.Y2 = 0.0
        self.X2 = 0.0

    @property
    def shared_clone(self):
        return self.__shared_clone

    @shared_clone.setter
    def shared_clone(self, val):
        self.__shared_clone = val

    @property
    def node(self):
        return self.__node

    @property
    def id(self):
        return self.__node.id

    @property
    def children(self):
        '''

        :return: list of cfit.tree.node.Node
        '''
        return self.__node.children

    @property
    def neoantigens(self):
        '''

        :return: list of Neoantigen objects
        '''
        return self.__node.neoantigens

    @property
    def exclusiveNeoantigens(self):
        '''

        :return: list of Neoantigen objects
        '''
        return self.__node.exclusiveNeoantigens

    @property
    def mutations(self):
        '''

        :return: list of Mutation objects
        '''
        return self.__node.mutations

    @property
    def exclusiveMutations(self):
        '''

        :return: list of Mutation objects
        '''
        return self.__node.exclusiveMutations


    @property
    def cnvs(self):
        '''

        :return: list of str objects
        '''
        return self.__node.cnvs

    @property
    def exclusiveCnvs(self):
        '''

        :return: list of Mutation objects
        '''
        return self.__node.exclustiveCnvs
    @property
    def X(self):
        return self.__X

    @X.setter
    def X(self, x):
        self.__X = x

    @property
    def Y(self):
        return self.__Y

    @Y.setter
    def Y(self, y):
        self.__Y = max(y, 0.0)

    @property
    def sharedY(self):
        return self.__sharedY

    @sharedY.setter
    def sharedY(self, y):
        self.__sharedY = max(y, 0.0)

    @property
    def privateY(self):
        return self.__privateY

    @privateY.setter
    def privateY(self, y):
        self.__privateY = max(y, 0.0)

    @property
    def preservedY(self):
        return self.__preservedY

    @preservedY.setter
    def preservedY(self, y):
        self.__preservedY = max(y, 0.0)

    @property
    def lostY(self):
        return self.__lostY

    @lostY.setter
    def lostY(self, y):
        self.__lostY = max(y, 0.0)

    @property
    def predictedY(self):
        return self.__predictedY

    @predictedY.setter
    def predictedY(self, y):
        self.__predictedY = y  # = max(y, 0.0)

    @property
    def predictedX(self):
        return self.__predictedX

    @predictedX.setter
    def predictedX(self, y):
        self.__predictedX = y  # max(y, 0.0)

    def get_info(self):
        '''
        Return a dictionary with information about the node: node id, X, Y, fitness and the number of mutations
        :return: dict
        '''
        info = {
            'ID': self.id,
            'X': self.X,
            'Y': self.Y,
            'F': self.fitness,
            'nmut': len(self.node.mutations)
        }
        return info

    def toJSON(self):  # TODO  correct with SamleTree
        '''
        Creates json representation of the node.
        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''

        js = self.node.toJSON()
#        js['id'] = self.id
#        js['mutations'] = [str(mut.id) for mut in self.node.mutations]
#        js['neoantigens'] = [str(neo.id) for neo in self.node.neoantigens]
#        js['fitness'] = self.fitness
#        if len(self.node.children) > 0:
#            children = [cnode.toJSON() for cnode in self.children]
#            js['children'] = children
        js['X'] = self.X
        js['Y'] = self.Y
        return js

    def set_X2(self):
        '''
        Sets the predicted inclusive frequency of the clone, based on
        the precomputed exclusive frequencies node.Y2
        :return:
        '''
        self.X2 = self.Y2
        for cnode in self.children:
            cnode.set_X2()
            self.X2 += cnode.X2


    def set_neoantigen_immunogenicities(self, neoantigen_immunogenicities):
        '''
        :param neoantigen_immunogenicities: dict (str -> float)
        '''

        self.neoantigen_immunogenicities = neoantigen_immunogenicities

    def get_fitness(self):
        '''
        Return precomputed node fitness

        :return: float
        '''
        return self.fitness

    def mutation_presentation_score(self, mid, kd0=500., strict=True, dominant=True):
        '''

        :param mid: str
            mutation identifier, <chrom>_<pos>_<ref>_<alt>

        :param kd0: float
            threshold on neoantigen dissociation constant

        :param strict: bool
            whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :param dominant: bool
            report the max presentation score over the neoantigens for mutation mid, or the sum of the scores.

        :return: float
            the presentation score for the mutation
        '''
        return self.node.mutation_presentation_score(mid, kd0=kd0, strict=strict, dominant=dominant)

    def get_relative_fitness(self):
        '''
        Return the precomputed relative fitness (node.nfitness, the fitness relative to the average
        tumor fitness)

        :return: float
        '''

        return self.rfitness

    def get_mutation_data(self, sampleName="", exclusive=False):
        '''

        Returns a DataFrame with data about mutations in the clone.
        The columns are:
        "Sample", "CloneNumber", "Mutation", "Gene",
                  "X", "Y", "Taf", "Tcov", "Naf", "Ncov"

        :param sampleName: str
            name the sample

        :param exclusive: bool
            whether to include only exclusive neoantigens.

        :return: pd.DataFrame
        '''

        header = ["Sample", "CloneNumber", "Mutation", "Gene",
                  "X", "Y"]
        muts = self.mutations
        if exclusive:
            muts = self.exclusiveMutations
        if len(muts) == 0:
            return None
        lines = [[sampleName, self.node.id, mut.id, mut.gene, self.X, self.Y] for mut
                 in muts]
        lines = pd.DataFrame(lines)
        lines.columns = header
        return lines

    def write_neoantigens(self, sampleName="", exclusive=False, kdthr=500., nq={}):
        '''
        Writes neoantigen data to the output file.

        :param sampleName: str
            name the sample

        :param exclusive: bool
            whether to include only exclusive neoantigens.

        :param kdthr: float
            threshold on Kd of neoantigens to include (<=kdthr)

        :return pandas.DataFrame
        '''
        neos = self.node.neoantigens
#        neos_indel = [neo for neo in self.node.fsneoantigens]
        if exclusive:
            neos = self.node.exclusiveNeoantigens
#            neos_indel = self.node.fs_exclusiveNeoantigens
        neos = [neo for neo in neos if neo.kD <= kdthr]
#        neos_indel = [neo for neo in neos_indel if neo.kD <= kdthr]

        nqual = defaultdict(float)
        for nid in nq:
            nqual[nid] = nq[nid]

        neos = [[neo.id, -nqual[neo.id], neo] for neo in neos]
#        neos += [[neo.id, -nqual[neo.id], neo] for neo in neos+neos_indel]
        neos.sort(key=lambda neo: neo[1])
        data = []

        fit_components_columns = list(self.fitness_components.keys())
        fit_components_values = list(self.fitness_components.values())

        for elneo in neos:
            neo = elneo[2]
            quality_components = neo.qattributes.get_quality_components()
            quality_values = list(quality_components.values())
            quality_columns = list(quality_components.keys())

            gene = self.node.dmutations[neo.mid].gene

            el = [neo.id, sampleName, self.id, gene, neo.mtPeptide, neo.wtPeptide,
                  self.X, self.Y, neo.allele, neo.kD, neo.wtkD, self.fitness] + \
                [x for x in fit_components_values] + [nqual[neo.id]] + quality_values
            data.append(el)
        if len(data) > 0:
            data = pd.DataFrame(data)
            data.columns = ["neoantigen", "sample", "clone_number", "gene", "peptideMT", "peptideWT",
                            "X", "Y", "MHC_allele", "kDmt", "kDwt", "clone_fitness"] + \
                           ["fitness_component_"+x for x in fit_components_columns] + \
                           ["quality"] + quality_columns
            return data
        else:
            return None
#            el = [str(x) for x in el]
#            f.write("\t".join(el) + "\n")
#        f.close()

    def write_indel_neoantigens(self, sampleName="", exclusive=False, kdthr=500.):
        '''
        Writes basics indel neoantigen data to the output file.

        :param sampleName: str
            name the sample

        :param exclusive: bool
            whether to include only exclusive neoantigens.

        :param kdthr: float
            threshold on Kd of neoantigens to include (<=kdthr)

        :return pandas.DataFrame
        '''

        neos = self.node.fsneoantigens
        neos = [neo for neo in neos if neo.kD <= kdthr]
        neos = [[neo.id, neo] for neo in neos]
        neos.sort(key=lambda neo: neo[1].kD)
        data = []

        fit_components_columns = list(self.fitness_components.keys())
        fit_components_values = list(self.fitness_components.values())

        for elneo in neos:
            neo = elneo[1]

            gene = self.node.dmutations[neo.mid].gene
            el = [neo.id, sampleName, self.id, gene, neo.mtPeptide,
                  self.X, self.Y, neo.allele, neo.kD, self.fitness] + \
                 [x for x in fit_components_values]
            data.append(el)
        if len(data) > 0:
            data = pd.DataFrame(data)
            data.columns = ["neoantigen", "sample", "clone_number", "gene", "peptideMT",
                            "X", "Y", "MHC_allele", "kDmt", "clone_fitness"] + \
                           ["fitness_component_"+x for x in fit_components_columns]
            return data
        else:
            return None

    def write_mutations(self, opath, sampleName="", exclusive=False):
        '''
        Writes mutation data to the output file.

        :param opath: str
            output file path

        :param sampleName: str
            name the sample

        :param exclusive: bool
            whether to include only exclusive neoantigens.

        '''
        f = open(opath, 'a')
        muts = self.node.mutations
        if exclusive:
            muts = self.node.exclusiveMutations
        lines = map(lambda mut: "\t".join(map(str, [sampleName, self.id,
                                                    mut.id, self.X, self.Y,
                                                    mut.Taf, mut.Tcov, mut.Naf, mut.Ncov])), muts)
        for line in lines:
            f.write(line + "\n")
        f.close()

    def score_neoantigens(self, neo_qualities, tau=1, just_dominant=False):
        '''
        Score neoantigens in the node as log(self.Y) -neo.fitness*tau (just_dominant=False)
            or log(self.Y) -neo.fitness*tau if dominant else 0

        :param neo_qualities: dict: str->float
            dictionary mapping neoantigens to their qualities

        :param tau: float

        :param just_dominant: bool

        :return: dict: str -> float
            dictionary with scores for neoantigens
            the higher the score the better (impact on the tumor)
        '''

        dscores = {}
        if len(self.neoantigens) == 0:
            return dscores
        max_quality = max([neo_qualities[neo.id] for neo in self.neoantigens])

        for neo in self.neoantigens:
            neo_quality = neo_qualities[neo.id]
            #score = np.sqrt(self.Y*(max(0, self.Y- np.exp(np.log(self.Y) - neo_quality * tau))))
            score = self.Y*(max(0, self.Y- np.exp(np.log(self.Y) - neo_quality * tau)))
            if just_dominant:
                if neo_quality < max_quality:
                    score = 0
            dscores[neo.id] = score
        return dscores


    def toJSON(self):  # TODO  correct with SamleTree
        '''
        Creates json representation of the node.
        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''

        js = self.node.toJSON()
#        js['id'] = self.id
#        js['mutations'] = [str(mut.id) for mut in self.node.mutations]
#        js['neoantigens'] = [str(neo.id) for neo in self.node.neoantigens]
#        js['fitness'] = self.fitness
#        if len(self.node.children) > 0:
#            children = [cnode.toJSON() for cnode in self.children]
#            js['children'] = children
        js['X'] = self.X
        js['Y'] = self.Y
        return js


    def tree_node_to_sibyl(self):  # TODO  correct with SamleTree
        '''
        Creates json representation of the node.
        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''

        js = self.node.to_sibyl_JSON()
        # js['tree_id'] = tree_id
        # js['id'] = self.ids
    
        js['fitness'] = self.fitness
#        if len(self.node.children) > 0:
#            children = [cnode.toJSON() for cnode in self.children]
#            js['children'] = children
        js['X'] = self.X
        js['Y'] = self.Y
        # js['mutations'] = [str(mut.id) for mut in self.node.mutations]
        # js['neoantigens'] = [str(neo.id) for neo in self.node.neoantigens]
        return js
