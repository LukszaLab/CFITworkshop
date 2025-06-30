'''
Created on Mar 24, 2015

@author: mluksza
'''

from collections import defaultdict

import pandas as pd

from cfit.tree.node.GeneralNode import GeneralNode


class Node(GeneralNode):
    '''
    Class for representing tumor clones.

    Attributes:

        mutations: list
            list of all mutations in that clone, objects of class Mutation

        dmutations: dict
            str -> Mutation, maps mutation identifiers in that clone to objects of class Mutation

        exclusiveMutations: list
            list of mutations that originate in that clone only, objects of class Mutation

        cnvs: list
            list of cnvs identifiers, str

        exclusiveCnvs: list
             list of snvs that originate in that clone only, str

        neoantigens: list
            list of neoantigens in that clone, objects of class Neoantigen

        fsneoantigens: list
            list of neoantigens in that clone, objects of class FrameShiftNeoantigen

        exclusiveNeoantigens: list
            list of neoantigens that originate in that clone

        fsneoantigens: list
            list of neoantigens in that clone, objects of class Neoantigen

        fs_exclusiveNeoantigens: list
            list of neoantigens that originate in that clone

        mutation2neoantigens:  dict, defaultdict(list)
            maps neoantigens to mutations, mid -> list of neoantigens

        mutation2fsneoantigens:  dict, defaultdict(list)
            maps neoantigens to mutations, mid -> list of neoantigens


        __gene2mut_id: dict, defaultdict(list)
            maps genes to the list of mutation identifiers

    '''

    def __init__(self, nid):
        '''
        Constructor
        '''
        GeneralNode.__init__(self, nid)

        self.mutations = []  # list of all mutations in that clone, instances of class Mutation
        self.dmutations = {}  # dictionary, maps mid -> Mutation
        self.exclusiveMutations = []  # mutations that originate in that clone only

        self.cnvs = []
        self.exclusiveCnvs = []

        self.neoantigens = []
        self.exclusiveNeoantigens = []  # neoantigens that originate in that clone
        self.fsneoantigens = []
        self.fs_exclusiveNeoantigens = []  # frame shift neoantigens that originate in that clone
        self.mutation2neoantigens = defaultdict(list)  # maps neoantigens to mutations, mid -> list of neoantigens
        self.mutation2fsneoantigens = defaultdict(
            list)  # maps neoantigens to mutations, mid -> list of frame shift neoantigens

        self.__gene2mut_id = defaultdict(list)

    def n_nsyn(self):
        '''
        Returns number of non-synonymous mutations
        :return: int
        '''
        q = sum([(mut.is_nonsynonymous()) for mut in self.mutations])
        return q

    def n_syn(self):
        '''
        Returns number of synonymous mutations
        :return: int
        '''
        q = sum([(mut.is_synonymous()) for mut in self.mutations])
        return q

    def contains_mutated_gene(self, gene):
        '''
        Whether a given gene is mutated in that clone.
        :param gene: str
        :return: bool
        '''
        return int(len(self.__gene2mut_id[gene]) > 0)

    def add_mutation(self, mutation, check_unique=True):
        '''
        Adds a mutation to the node, and the nested (children) clones.

        :param mutation: cfit.tree.mutation.Mutation
            mutation object to be added

        :param check_unique: bool
            impose that the list of mutations is a set
        '''
        if check_unique:
            mids = [mut.id for mut in self.mutations]
            if mutation.id not in mids:
                self.mutations.append(mutation)
        else:
            self.mutations.append(mutation)
        self.dmutations[mutation.id] = mutation
#        self.logger("adding "+mutation.id+" to node "+str(self.id))
        self.__gene2mut_id[mutation.gene].append(mutation.id)
        for cnode in self.children:
            cnode.add_mutation(mutation, check_unique=check_unique)

    def add_cnv(self, cnvid):
        '''
        Add a cnv to the node and to the children nodes

        :param cnvid: str
        '''

        self.cnvs.append(cnvid)
        for cnode in self.children:
            cnode.add_cnv(cnvid)


    def add_neoantigen(self, neoantigen, recursive=True, check_unique=True):
        '''
        Add a neoantigen to the node
        :param neoantigen: cfit.tree.mutation.Neoantigen
            neoantigen object to be added

        :param recursive: bool
            if set to True, the neoantigen is added to the nested clones as well.

        :param check_unique: bool
            check that neoantigen list is a set

        :return:
        '''
        self.neoantigens.append(neoantigen)
        if check_unique:
            nset = set(self.neoantigens)
            nset.add(neoantigen)
            self.neoantigens = list(nset)
        if recursive:
            for cnode in self.children:
                cnode.add_neoantigen(neoantigen, recursive=recursive, check_unique=check_unique)
        self.mutation2neoantigens[neoantigen.mid].append(neoantigen)

    def add_frame_shift_neoantigen(self, neoantigen, recursive=True, check_unique=True):
        '''
        Add a neoantigen to the node
        :param neoantigen: cfit.tree.mutation.FrameShiftNeoantigen
            neoantigen object to be added

        :param recursive: bool
            if set to True, the neoantigen is added to the nested clones as well.

        :param check_unique: bool
            check that neoantigen list is a set

        :return:
        '''
        self.fsneoantigens.append(neoantigen)
        if check_unique:
            nset = set(self.fsneoantigens)
            nset.add(neoantigen)
            self.fsneoantigens = list(nset)
        if recursive:
            for cnode in self.children:
                cnode.add_frame_shift_neoantigen(neoantigen, recursive=recursive, check_unique=check_unique)
        self.mutation2fsneoantigens[neoantigen.mid].append(neoantigen)

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
            report the max presentation score over the neoantigens for mutation mid, or the sum of the scores, if false.

        :return: float
            the presentation score for the mutation
        '''

        v = [neo.presentation_score(kd0=kd0, strict=strict) for neo in self.mutation2neoantigens[mid]]
        if len(v) == 0:
            return 0.0
        if dominant:
            v = max(v)
        else:
            v = sum(v)
        return v

    def neo2mut_ratio(self, mode='all', nsyn=False, kd0=500., strict=True):
        '''
        Return the neoantigen-to-mutation ratio, the expected number of neoantigens per mutation in this clone.

        :param mode: str
            if set to "all" will set dominant parameter to False and will count all neoantigens per mutation
            otherwise reports only whether a mutation has a neoantigen.

        :param nsyn: bool
            include only non-synonymous mutations

        :param kd0: float
            threshold on neoantigen dissociation constant

        :param strict: bool
            whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :return: float
        '''

        mids = [mut.id for mut in self.mutations]
        if nsyn:
            mids = [mut.id for mut in self.mutations if mut.is_nonsynonymous()]

        dominant = True
        if mode == 'all':
            dominant = False
        num = sum([self.mutation_presentation_score(mid, kd0=kd0, strict=strict, dominant=dominant) for mid in mids])

        denom = len(mids)
        num = float(num)
        denom = float(denom)
        r = 0
        if denom > 0:
            r = num / denom
        return r

    def get_mutation_ids(self):
        '''
        Returns the list of mutation identifiers
        :return: list
        '''
        return [mut.id for mut in self.mutations]

    def get_neoantigen_ids(self):
        '''
        Returns the list of neoantigen identifiers
        :return: list
        '''
        return [neo.id for neo in self.neoantigens]

    def set_exclusive_mutations(self):
        '''
        Fills in the list of exclusive mutations and neoantigens at initialization.
        '''
        __id2mut = {}
        for mut in self.mutations:
            __id2mut[mut.id] = mut
        mids = set(self.get_mutation_ids()).difference(self.parent.get_mutation_ids())
        self.exclusiveMutations = [__id2mut[mid] for mid in mids]

        __id2neo = {}
        for neo in self.neoantigens:
            __id2neo[neo.id] = neo
        nids = set(self.get_neoantigen_ids()).difference(set(self.parent.get_neoantigen_ids()))
        self.exclusiveNeoantigens = [__id2neo[nid] for nid in nids]

    def set_exclusive_cnvs(self):
        '''
        Fills in the list of exclusive mutations and neoantigens at initialization.
        '''
        excnvids = set(self.cnvs).difference(self.parent.cnvs)
        self.exclusiveCnvs = excnvids

    def remove_mutations(self, muts):
        '''
        Remove mutations from the clone.

        :param muts: container (list, set)
            list of elements of class cfit.tree.mutation.Mutation

        '''
        smuts = set(self.mutations)
        self.mutations = list(smuts.difference(muts))
        smuts = set(self.exclusiveMutations)
        self.exclusiveMutations = list(smuts.difference(muts))

    def remove_neoantigens(self, neos):
        '''
        Remove neoantigens from the clone.
        :param neos: container (list, set)
            list of elements of class cfit.tree.mutation.Neoantigen
        '''
        sneos = set(self.neoantigens)
        self.neoantigens = list(sneos.difference(neos))
        sneos = set(self.exclusiveNeoantigens)
        self.exclusiveNeoantigens = list(sneos.difference(neos))

    def TMB(self):
        '''
        Returns the number of mutation in the clone.
        :return: int
        '''
        return len(self.mutations)

    def TMB_MHC(self, kd0=500, strict=True):
        '''
        Returns the number of mutations with a neoantigen.

        :param kd0: float
            threshold on neoantigen dissociation constant

        :param strict: bool
            whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :return: float
        '''

        q = sum(
            [self.mutation_presentation_score(mut.id, kd0=kd0, strict=strict, dominant=True) for mut in self.mutations])
        return q

    def neoantigen_load(self):
        '''
        Returns the number of neoantigens
        :param HLAW:
        :return: float
        '''
        nl = self.presentation_score(kd0=500, strict=True)
        return nl

    def presentation_score(self, kd0=50.0, strict=False):
        '''
        Aggregated presentation score over neoantigens in the clone.

        :param kd0: float
            threshold on neoantigen dissociation constant

        :param strict: bool
            whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :return: float
        '''
        pres_score = sum([neo.presentation_score(kd0=kd0, strict=strict) for neo in self.neoantigens])
        return pres_score

    def syn(self):
        '''
        Return the number of synonymous mutations

        :return: int
        '''
        s = sum(map(lambda mut: mut.is_synonymous(), self.mutations))
        return s


    def get_exclusive_mutations(self):
        '''
        Returns the set of mutations that originate in this clone.
        :return: list
            list of mutations, objects of class Mutation
        '''

        __id2mut = {}
        for mut in self.mutations:
            __id2mut[mut.id] = mut
        mids = set(self.get_mutation_ids()).difference(self.parent.get_mutation_ids())
        emuts = [__id2mut[mid] for mid in mids]
        return emuts

    def toVCF(self):
        '''
        Create a vcf-like DataFrame for all mutations in the node (clone)

        :return: pd.DataFrame
        '''

        def chr2num(chrom):
            try:
                chrnum = int(chrom)
            except:
                chrnum = 23
            return chrnum

        if len(self.mutations) == 0:
            return None
        cols = ["#CHROM", "POS", "ID", "REF", "ALT",
                "QUAL", "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR"]
        muts = self.mutations

        muts.sort(key=lambda mut: chr2num(mut.chrom) * 1e10 + int(mut.pos))
        dat = pd.DataFrame(map(lambda mut: mut.toVCF(), muts))
        dat.columns = cols
        return dat

    ####################
    # Sibyl
    ####################

    def toJSON(self):
        '''
        Creates json for Sibyl app.
        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''

        js = {}
        js['clone_id'] = self.id
        js['clone_mutations'] = [str(mut.id) for mut in self.exclusiveMutations]
#        js['neoantigens'] = [str(neo.id) for neo in self.neoantigens]
        if len(self.children) > 0:
            children = [cnode.toJSON() for cnode in self.children]
            js['children'] = children
        return js

    def to_sibyl_JSON(self):
        '''
        Creates json for Sibyl app.
        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''

        js = {}
        js['id'] = self.id
        # if self.id != self.parent.id:
            # js['p'] = self.parent.id
        # js['mutations'] = [str(mut.id) for mut in self.exclusiveMutations]
        # js['neoantigens'] = [str(neo.id) for neo in self.neoantigens]
        # if len(self.children) > 0:
            # children = [cnode.toJSON() for cnode in self.children]
            # js['children'] = children
        return js
    

    def get_tree_node(self):
        '''
        Creates json for Sibyl app.
        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''

        js = {}
        js['id'] = self.id
        if self.id != self.parent.id:
            js['p'] = self.parent.id
        js['mutations'] = [str(mut.id) for mut in self.exclusiveMutations]
        js['neoantigens'] = [str(neo.id) for neo in self.neoantigens]
        # if len(self.children) > 0:
            # children = [cnode.toJSON() for cnode in self.children]
            # js['children'] = children
        return js
