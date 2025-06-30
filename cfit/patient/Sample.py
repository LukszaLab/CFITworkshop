'''
Created on Mar 30, 2015

@author: mluksza
'''
import gzip
import json
import os
from collections import defaultdict
from random import shuffle

import numpy as np
import pandas as pd

from cfit.CoreObject import CoreObject
from cfit.tree.SampleTree import SampleTree
from cfit.tree.Tree import Tree
from cfit.tree.mutation.Mutation import Mutation
from cfit.tree.mutation.NonsynonymousMutation import NonsynonymousMutation
from cfit.tree.mutation.SampleMutation import SampleMutation
from cfit.tree.mutation.SynonymousMutation import SynonymousMutation
from cfit.tree.node.Node import Node
from cfit.util.Utils import Utils


class Sample(CoreObject):
    '''
    Class implements a single whole genome sequencing tumor sample.

    Attributes:
        name: str

        neoantigens: set
            set of neoantigens (objects of class cfit.tree.mutation.Neoantigen)

        fsneoantigens: set
            set of neoantigens (objects of class cfit.tree.mutation.Neoantigen)

        neoantigensMap: dict: str->Neoantigen
            dictionary mapping neoantigen id to the Neoantigen object

        neoantigenQualities: dict: str->float
            dictionary mapping neoantigen id to its computed quality in that sample

        neoantigenQualityComponents: dict: str->dict
            dictionary mapping neoantigen id to its quality components

        peptideMap: dict, defaultdict(list)
            maps neoantigen peptide id to the list of neoantigens

        mutations: dict: str -> SampleMutation
            Stores all mutations in the sample. Dictionary mapping mutation identifiers to SampleMutation objects

        mutation2neoantigens: dict, defaultdict(list)
            Dictionary mapping mutations (str, mutation id, <chr>_<pos>_<ref>_<alt>) to the list
            of neoantigen objects (Neoantigen)

        mutation2fsneoantigens: dict, defaultdict(list)
            Dictionary mapping mutations (str, mutation id, <chr>_<pos>_<ref>_<alt>) to the list
            of neoantigen objects (FrameShiftNeoantigen)

        affectedGenes: set
            set of genes that are in some way affected in the tumor

        synCount: int
            number of synonymous mutations

        nsynCount: int
            number of non-synonymous mutations

        trees: list
            the list of trees, objects of cfit.tree.SampleTree, sorted by their likelihood score. The trees represent
            possible clonal reconstructions.

        oneTree: cfit.tree.Tree:
            a tree representation of the sample with just one clone, ignoring mutation frequencies. Used for
            evaluating predictions that ignore tumor heterogeneity

        response: bool
            whether the patient is a responder to therapy

        survival: float
            survivral time of the patient

        dead: bool
            whether the patient is dead

        self.fitness: float
            predicted tumor fitness of the sample. Compute as a weigteed average over the clonal structures (trees).

        __tissue: str
            tumor tissue

        __q: float
            auxiliary attribute for storing sample quality statitstic.

        __normal_reads: int
            number of reads in the corresponding normal sample (if available)

        __tumor_reads: int
            number of reads in this tumor sample

        __gene2mut_id: dict, defaultdict(list)
            dictionary mapping gene names to the list of mutation identifiers (str) in this sample

        __sharing: float
            amount of shared volume with other sample(s) if multiple samples from the patient were aligned
            for clonal reconstruction

        __gene_expression: dict: str->float
            maps gene names (ENSG and descriptive names) to expression levels

        __mutation_expression: dict: str->float
            maps mutation identifiers to expression levels

    '''

    @staticmethod
    def scrambleListOfNeoantigens(neos, rpeptides=None):
        '''
        Randomization of neoantigens
        :param neos:
        :param rpeptides:
        :return:
        '''

        if rpeptides is None:
            newseqs = list(map(lambda neo: "", neos))
            for i in range(0, 9):
                chars = list(map(lambda neo: neo.mtPeptide[i], neos))
                shuffle(chars)
                for j in range(0, len(newseqs)):
                    newseqs[j] += chars[i]
        else:
            newseqs = rpeptides
            while len(newseqs) < len(neos):
                newseqs += rpeptides
        for i in range(0, len(neos)):
            neos[i].mtPeptide = newseqs[i]
            neos[i].wtPeptide = newseqs[i]
        return neos

    def __init__(self, name):
        '''
        Constructor method.

        :param name: str
            name of the sample
        '''
        self.name = name
        self.neoantigens = set()
        self.fsneoantigens = set()
        self.neoantigensMap = {}
        self.peptideMap = defaultdict(list)  # maps neoantigen peptide id to the list of neoantigen
        self.mutations = {}
        self.mutation2neoantigens = defaultdict(list)
        self.mutation2fsneoantigens = defaultdict(list)
        self.affectedGenes = set()
        self.synCount = 0
        self.nsynCount = 0
        self.trees = []
        self.oneTree = None
        self.timePoint = None
        self.response = False
        self.survival = 0.0
        self.dead = False
        self.fitness = 0.0
        self.__tissue = "Unknown"
        self.__q = 0.0
        self.__totalR = 0.0  # total number of aligned neoantigens
        self.__normal_reads = 0
        self.__tumor_reads = 0
        self.__gene2mut_id = defaultdict(list)
        self.__sharing = 0.0
        self.__gene_expression = None
        self.__mutation_expression = None
        self.neoantigenQualities = defaultdict(float)
        self.neoantigenQualityComponents = defaultdict(dict)

        self.vaccine_mutations = []

    @property
    def sharing(self):
        return self.__sharing

    @sharing.setter
    def sharing(self, sharing):
        self.__sharing = sharing

    @property
    def tissue(self):
        return self.__tissue

    @tissue.setter
    def tissue(self, tissue):
        self.__tissue = tissue

    @property
    def normal_reads(self):
        return self.__normal_reads

    @normal_reads.setter
    def normal_reads(self, normal_reads):
        self.__normal_reads = normal_reads

    @property
    def tumor_reads(self):
        return self.__tumor_reads

    @normal_reads.setter
    def tumor_reads(self, tumor_reads):
        self.__tumor_reads = tumor_reads

    @property
    def q(self):
        return self.__q

    @q.setter
    def q(self, q):
        self.__q = q

    def set_gene_expression(self, gene, expr_lev):
        '''
        Set gene expression level of gene
        :param gene: str
        :param expr_lev: float
        '''
        if self.__gene_expression is None:
            self.__gene_expression = defaultdict(lambda: None)
        self.__gene_expression[gene] = expr_lev

    def get_gene_expression(self, gene):
        if self.__gene_expression is None:
            return None
        return self.__gene_expression[gene]

    def set_mutation_expression_levels(self):
        '''
        Create dictionary mapping gene expression levels of mutations
        '''
        if self.__gene_expression is None:
            return
        self.__mutation_expression = defaultdict(float)
        for mut in self.mutations.values():
            if mut.is_coding():
                self.__mutation_expression[mut.ENSG] = self.__gene_expression[mut.ENSG]
                self.__mutation_expression[mut.gene] = self.__gene_expression[mut.gene]

    def __compute_gene2mut_id(self):
        '''
        Private method, computes __gene2mut_id auxiliary dictionary attribute that maps
        genes to mutation identifiers
        :return:
        '''
        for mid in self.mutations:
            mut = self.mutations[mid]
            self.__gene2mut_id[mut.gene].append(mid)

    def import_mutations_from_VCF_file(self, vcffile, mutations):
        '''
        Imports sample mutations from a vcf file. Sets up all the mutation related attributes.

        :param vcffile: str

        :param mutations: dict: str -> cfit.tree.mutation.Mutation
            dictionary mapping mutation identifiers to Mutation objects in the patient
        :return:
        '''
        f = open(vcffile)
        self.logger("importing " + vcffile)
        header = f.readline()

        columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                   "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR"]

        while header[:2] == "##":
            header = f.readline()
        header = header.strip().split("\t")
        line = f.readline()
        while line:
            tab = line.strip().split("\t")
            hdict = {}
            for i in range(0, len(header)):
                hdict[columns[i]] = tab[i]
            hdict["Sample"] = self.name
            if hdict["ID"][:3] == "chr":
                hdict["ID"] = hdict["ID"][3:]
            hdict["ID"] = hdict["ID"].upper()

            if hdict["ID"] in mutations:
                mutation = mutations[hdict["ID"]]
            else:
                infos = hdict["INFO"].split(",")
                infos_1 = [ann for ann in infos if "WARNING" not in ann.split("|")[-1]]
                if len(infos_1) > 0:
                    infos = infos_1
                info = infos[0]

                infos_1 = [info for info in infos if len(info.split("|")) > 1]
                # give preference to neoantigen missense mutations
                infos_1 = [info for info in infos_1 if "missense_variant" in info.split("|")[1]]
                #                infos_1 = [info for info in infos if "missense_variant" in info.split("|")[1]]
                if len(infos_1) > 0:
                    info_1 = infos_1[0]
                    info = info_1
                hdict["INFO"] = info
                mutation = Mutation()
                info = info.split("|")
                if len(info) <= 1:
                    line = f.readline()
                    continue
                if "synonymous_variant" in info[1]:
                    mutation = SynonymousMutation()
                elif "missense_variant" in info[1]:
                    mutation = NonsynonymousMutation()
                mutation.initialize_VCF(hdict)
                mutations[mutation.id] = mutation
            sampleMutation = SampleMutation(mutation)
            sampleMutation.initialize_VCF(hdict)
            self.add_mutation(sampleMutation)
            line = f.readline()
        f.close()
        self.__compute_gene2mut_id()
        return mutations

    def get_tumor_purity_scores(self, jsonpath, ord=0):
        '''
        Estimates tumor purity from PhyloWGS tree reconstruction. Purity is the CCF of the root of the tree.

        :param jsonpath: str
            patht to the PhyloWGS json file <sname>.summ.json.gz

        :param ord: int
            the number of the sample (if multiple samples were used in the reconstruction)

        :return: pd.DataFrame
            a table with tree statistics, including purity, over the family of trees.
            The columns are: Tree_num (tree rank, by score), llh (log-likelihood), poly (number of roots),
                            purity (estimated total CCF of the tree roots)
         '''

        if ".gz" in jsonpath:
            f = gzip.open(jsonpath)
        else:
            f = open(jsonpath)
        line = f.readline()
        f.close()
        js = json.loads(line)
        keys = list(js['trees'].keys())
        stats = []
        for key in keys:
            t = js["trees"][key]
            llh = t['llh']
            roots = t["structure"]["0"]
            poly = len(roots)
            purity = sum([t["populations"][str(root)]["cellular_prevalence"][ord] for root in roots])
            stats.append([key, llh, poly, purity])
        stats = pd.DataFrame(stats)
        stats.columns = ["Tree_num", "llh", "poly", "purity"]
        return stats

    def add_mutation(self, mut):
        '''
        Adds sample mutation to the sample. Updates mutations attribute of the sample.

        :param mut: cfit.tree.mutation.SampleMutation
        :return:
        '''
        self.mutations[mut.id] = mut

    def mutation_frequency_variance(self, filter_on_Taf=False):
        '''
        Compute frequency variance of mutations in the sample over trees in the sample

        :param filter_on_Taf: bool
            whether to exclude mutations not detected by reads in the sample

        :return: float
        '''
        cond = lambda mut: mut.Taf > 0 if filter_on_Taf else 1
        ldfreqs = []
        for tree in self.trees:
            dfreq = tree.get_mutation_frequencies()
            ldfreqs.append(dfreq)
        avev = np.mean(
            [np.var([dfreq[mid] for dfreq in ldfreqs]) for mid in self.mutations if cond(self.mutations[mid])])
        return avev

    def compute_adjusted_taf(self):
        '''

        '''
        max_taf = max([mut.Taf for mut in self.mutations.values()])
        for mut in self.mutations.values():
            mut.adjustedTaf = mut.Taf/max_taf

    def effective_sample_size(self, eps=0.0, filter_on_Taf=False):
        '''
        Compute the effective sample size (number of cells contributing to the heterogeneity of the
        sample)

        :param eps: float
            pseudocount for zero-variance samples

        :param filter_on_Taf: bool
            whether to exclude mutations not detected by reads in the sample

        :return: float
            estimate of effective sample size
        '''
        avev = self.mutation_frequency_variance(filter_on_Taf=filter_on_Taf)
        n = 1. / (avev + eps)
        return n

    def add_neoantigen(self, neoantigen):
        '''
        Adds a neoantigen to the sample. Checks if the mutation that generated the neoantigen is in the sample.

        :param neoantigen: cfit.tree.mutation.Neoantigen

        '''

        if neoantigen.mid in self.mutations.keys():
            self.neoantigens.add(neoantigen)
            self.neoantigensMap[neoantigen.id] = neoantigen
            self.peptideMap[neoantigen.peptide_id].append(neoantigen)
            self.mutation2neoantigens[neoantigen.mid].append(neoantigen)

    def add_frame_shift_neoantigen(self, neoantigen):
        '''
        Adds a neoantigen to the sample. Checks if the mutation that generated the neoantigen is in the sample.

        :param neoantigen: cfit.tree.mutation.FrameShiftNeoantigen

        '''

        if neoantigen.mid in self.mutations.keys():
            self.fsneoantigens.add(neoantigen)
            self.neoantigensMap[neoantigen.id] = neoantigen
            self.peptideMap[neoantigen.peptide_id].append(neoantigen)
            self.mutation2fsneoantigens[neoantigen.mid].append(neoantigen)

    #    def set_one_tree(self):
    #        '''
    #        Sets the single clone representation of the sample.
    #        :return:
    #        '''
    #        self.oneTree = Tree([None, None, None, None, ""], thesample=self)

    def get_tree_weights(self, beta=1.):
        '''
        Compute how the tree should be weigted, based on the tree likelihood scores.


        :param beta: float
            weight parameter. For beta=1 the weights are proportional to their likelihood.
            For beta=0 the weights are equal.
        :return: list
        '''

        weights = [tree.llh * beta for tree in self.trees]
        weights = Utils.log_norm(weights)
        weights = [np.exp(w) for w in weights]
        return weights

    def get_gene_frequency(self, gene, beta=1.0):
        '''
        Compute mutated gene frequency in the sample, by averaging over the trees.
        :param gene: str
            gene name
        :param beta: float
            weight parameter
        :return: float
        '''
        q = self.average_over_clones(Node.contains_mutated_gene, beta, gene=gene)
        return q

    def get_gene_mutation_frequency(self, gene, beta=1.0):
        '''
        Compute mutated gene frequency in the sample, by averaging over the trees.
        :param gene: str
            gene name
        :param beta: float
            weight parameter
        :return: float
        '''

        def tree_fun(tree):
            xx = [node.X for node in tree.nodes.values() if node.contains_mutated_gene(gene)]
            if len(xx) == 0:
                return 0
            else:
                return max(xx)

        q = self.average_tree_function(tree_fun, beta=beta)
        return q

    def get_mutation_frequency(self, mid, beta=1.0):
        '''
        Compute mutated gene frequency in the sample, by averaging over the trees.
        :param mid: str
            mutation identifier, <chr>_<pos>_<ref>_<alt>
        :param beta: float
            weight parameter
        :return: float
        '''

        nodefun = lambda node: (mid in [mut.id for mut in node.mutations])

        def tree_fun(tree):
            xx = [node.X for node in tree.nodes.values() if nodefun(node)]
            if len(xx) == 0:
                return 0
            else:
                return max(xx)

        q = self.average_tree_function(tree_fun, beta=beta)
        return q

    def get_trunk_volume(self, beta=1.0):
        '''
        The volume of the trunk clone, averaged over trees
        :param beta: float
            weight parameter
        :return: float
        '''
        q = self.average_tree_function(Tree.get_trunk_volume, beta)
        return q

    def get_trunk_fitness(self, node_fitness_fun, beta=1, **kwargs):
        '''
        The fitness of the trunk clone, averaged over trees.

        :param node_fitness_fun: function
            Node class method, to compute clone fitness

        :param beta: float
            weight parameter, for averaging over trees

        :param kwargs:

        :return: float
        '''

        kwargs['node_fitness_fun'] = node_fitness_fun
        q = self.average_tree_function(Tree.get_trunk_fitness, beta, **kwargs)
        return q

    def get_average_non_trunk_fitness(self, node_fitness_fun, beta=1, **kwargs):
        '''
        Complementary method, computes the average fitness of non-trunk clones.

        :param node_fitness_fun: function
            Node class method, to compute clone fitness

        :param beta: float
            weight parameter, for averaging over trees

        :param kwargs:

        :return: float
        '''
        kwargs['node_fitness_fun'] = node_fitness_fun
        q = self.average_tree_function(Tree.get_average_non_trunk_fitness, beta, **kwargs)
        return q

    def get_average_tree_height(self, beta=1, synonymous=True, xthr=0.0):
        '''
        Returns the tree height, computed as (the average) number of mutations in
        the deepest clone, averaging is done over the family of high-scoring
        reconstructed trees.

        :param beta: float
            llh weighting of trees

        :param synonymous: bool
            default=True, use only synonymous mutations

        :param xthr: float
            threshold on clone size

        :return: float
        '''

        weights = self.get_tree_weights(beta)
        hei = 0.0
        for wei, tree in zip(weights, self.trees):
            hei += tree.get_height(synonymous=synonymous, xthr=xthr) * wei
        return hei

    def entropy(self, beta, shared_only=False):
        '''
        Compute sample clonal heterogeneity, averaging ove the trees

        :param beta: float
            the tree weighting parameter

        :return: float
        '''
        q = self.average_tree_function(SampleTree.entropy, beta, shared_only=shared_only)
        return q

    def get_effective_number_of_clones(self, beta=1.0):
        '''
        exp(entropy), averaged over trees.

        :param beta:
            the tree weighting parameter

        :return: float
        '''
        ent = self.average_tree_function(SampleTree.get_effective_number_of_clones, beta=beta)
        return ent

    def get_evolved_effective_number_of_clones(self, beta, tau):
        '''
        The predicted effective number of clones (exp(entropy)), averaged over the trees.

        :param beta: float
            tree weighting parameter

        :param tau: float
            time parameter

        :return: float
        '''

        q = self.average_tree_function(SampleTree.get_evolved_effective_number_of_clones,
                                       beta=beta, tau=tau)
        return q

    def get_number_of_clones(self, beta=1.0):
        '''
        The number of clones in the sample, averaged over trees.

        :param beta: float
            tree weighting parameter

        :return: float
        '''
        q = self.average_tree_function(SampleTree.get_number_of_clones, beta=beta)
        return q

    '''
    Mutations and neoantigens
    '''

    def effective_TMB(self, beta=1.0):
        '''
        The effective number of all mutations, averaged over trees
        :param beta: float
            the tree weigting parameter

        :return: float
        '''
        q = self.average_over_clones(Node.TMB, beta)
        return q

    def effective_TMB_syn(self, beta=1.0):
        '''
        The effective number of synonymous mutations, averaged over trees
        :param beta: float
            the tree weigting parameter
        :return: float
        '''
        q = self.average_over_clones(Node.n_syn, beta)
        return q

    def effective_max_n_syn(self, beta=1.0):
        '''
        The molecular clock (age of tumor computed as the number of synonymous mutations in the deepest clones),
        averaged over trees

        :param beta: float
            the tree weigting parameter
        :return: float
        '''
        q = self.average_tree_function(Tree.max_syn, beta=beta)
        return q

    def effective_TMB_nsyn(self, beta=1.0):
        '''
        The effective number of non-synonymous mutations, averaged over trees
        :param beta: float
            the tree weigting parameter
        :return: float
        '''
        q = self.average_over_clones(Node.n_nsyn, beta)
        return q

    def TMB_syn(self, filter_on_Taf=True):
        '''
        The number of synonymous mutations

        :param filter_on_Taf: bool
            whether to condition on mutations with non-zero reads
        :return: float
            number of synonymous mutations
        '''
        cond = lambda mut: mut.Taf > 0 if filter_on_Taf else True
        q = sum([mut.is_synonymous() for mut in self.mutations.values() if cond(mut)])
        return q

    def TMB_nsyn(self, filter_on_Taf=True):
        '''
        The number of non-synonymous mutations

        :param filter_on_Taf: bool
            whether to condition on mutations with non-zero reads

        :return: float
            number of nonsynonymous mutations
        '''

        cond = lambda mut: mut.Taf > 0 if filter_on_Taf else True
        q = sum([mut.is_nonsynonymous() for mut in self.mutations.values() if cond(mut)])
        return q

    def TMB(self, filter_on_Taf=True):
        '''
        The number of mutations

        :param filter_on_Taf: bool
            whether to condition on mutations with non-zero reads

        :return: float
            number of nonsynonymous mutations
        '''

        cond = lambda mut: mut.Taf > 0 if filter_on_Taf else True
        q = sum([1. for mut in self.mutations.values() if cond(mut)])
        return q

    def mutation_presentation_score(self, mid, kd0=50., strict=True, dominant=True):
        '''

        Returns the number of presented neoantigens that are associated with the given mutation.

        :param mid: str
            mutation identifier, <chr>_<pos>_<ref>_<alt>

        :param kd0: float
            threshold on the dissociation constant of the neoantigen to be considered

        :param strict: bool
            whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :param dominant: bool
            how the neontigens are counted

        :return: float
        '''

        neos = self.mutation2neoantigens[mid]

        v = [neo.presentation_score(kd0=kd0, strict=strict) for neo in neos]
        if len(v) == 0:
            return 0.0
        if dominant:
            return max(v)
        else:
            return sum(v)

    def neo2mut_ratio(self, mode='all', nsyn=True, kd0=50., strict=True, filter_on_Taf=True):
        '''
        Computes the number of neoantigens per mutation in the sample, without averaging over
        the trees.

        :param mode: str
            if set to "all" will set dominant parameter to False and will count all neoantigens per mutation
            otherwise reports only whether a mutation has a neoantigen.

        :param nsyn: bool
            use non-synonymous mutations only

        :param kd0: float
            threshold on neoantigen Kd

        :param strict: bool
            whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :param filter_on_Taf: bool

        :return: float
        '''

        cond = lambda mut: mut.Taf > 0 if filter_on_Taf else True
        muts = [mut for mut in self.mutations.values() if cond(mut)]
        mids = [mut.id for mut in muts]
        if nsyn:
            mids = [mut.id for mut in self.mutations.values() if mut.is_nonsynonymous()]

        dominant = True
        if mode == 'all':
            dominant = False
        num = sum([self.mutation_presentation_score(mid, kd0=kd0, strict=strict, dominant=dominant) for mid in mids if
                   cond(self.mutations[mid])])

        denom = float(len(mids))
        num = float(num)

        r = 0
        if denom > 0:
            r = num / denom
        return r

    def effective_neo2mut_ratio(self, beta=1.0, mode='all', nsyn=False, kd0=50., strict=True):
        '''
        Number of neoantigens per mutation, averaged over the trees.

        :param beta: float
            tree weighting parameter

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

        '''

        q = self.average_over_clones(Node.neo2mut_ratio, beta, mode=mode, nsyn=nsyn, strict=strict, kd0=kd0)
        return q

    def presentation_score(self, kd0=50, strict=True, filter_on_Taf=True):
        '''
        Overall presentation score of the sample, equivalent to the neoantigen load. Computed without
        averaging over the trees.

        :param kd0: float
            threshold on neoantigen dissociation constant
        :param strict: bool
            How the Kd threshold should be applied - whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :param filter_on_Taf: bool
            Exclude mutations without reads in this sample

        :return: float
        '''

        cond = lambda mut: mut.Taf > 0 if filter_on_Taf else True
        pres_score = sum([neo.presentation_score(kd0=kd0, strict=strict) for neo in self.neoantigens if
                          cond(self.mutations[neo.mid])])
        return pres_score

    def effective_presentation_score(self, kd0=50, strict=False, beta=1.0):
        '''
        Overall effective presentation score of the sample, equivalent to the neoantigen load,
        taking into the clonal structures.

        :param kd0: float
            threshold on neoantigen dissociation constant

        :param strict: bool
            How the Kd threshold should be applied - whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :param beta: float
            three weighting parameter

        :return: float
        '''
        q = self.average_over_clones(Node.presentation_score, beta, kd0=kd0, strict=strict)
        return q

    def neoantigen_load(self, beta=1.0):
        '''
        Traditional effective neoantigen load computation by counting neoantigens.
        Redundant with presentation_score

        :param beta: float
            tree weighting parameter

        :return: float
        '''
        q = self.average_over_clones(Node.neoantigen_load, beta)
        return q

    #    def TMB(self, beta=1.0):
    #        '''
    #        Effective mutation load computation by counting mutations.
    #
    #        :param beta: float
    #            tree weighting parameter
    #
    #        :return: float
    #        '''
    #        q = self.average_over_clones(Node.TMB, beta)
    #        return q

    def TMB_MHC(self, kd0=500, strict=True, filter_on_Taf=True):
        '''
        Tumor mutational burden, including only mutations with a presentable neoantigen

        :param kd0: float
            threshold on neoantigen dissociation constant

        :param strict: bool
            whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :param filter_on_Taf: bool
            whether to exclude mutations not detected in the sample

        :return: float
            mutational load
        '''

        cond = lambda mut: mut.Taf > 0 if filter_on_Taf else True
        q = sum(
            [self.mutation_presentation_score(mut.id, kd0=kd0, strict=strict, dominant=True) for mut in
             self.mutations.values() if cond(mut)])
        return q

    def effective_TMB_MHC(self, beta=1.0, kd0=500, strict=True):
        '''
        Effective number of mutations with neoantigens, averaged over trees.

        :param beta: float
            tree weighting parameter

        :param kd0: float
            threshold on neoantigen dissociation constant

        :param strict: bool
            whether to use sharp presentation_score - neoantigen cutoff and count delta(neo.kd<=kd0) or a smooth one
            with kd0/(kd0+neo.kd) function

        :return: float
            effective mutational load, restricting t
        '''
        q = self.average_over_clones(Node.TMB_MHC, beta, kd0=kd0, strict=strict)
        return q

    '''
    Fitness
    '''

    def annotateP53Transactivity(self, p53transactivation):
        '''

        :param p53transactivation:
        :return:
        '''
        for tree in self.trees:
            tree.annotate_P53_transactivity(p53transactivation)
        self.oneTree.annotate_P53_transactivity(p53transactivation)

    '''
    Longitudinal analysis
    '''

    def get_shared_and_private(self, node_fitness_func):
        '''
        Report the amount of sharing between longitudinal samples, based on precomputed
        clone frequencies.

        :param node_fitness_func: function
            the statistic that is computed over the nodes
        :return: list
            statistics of shared and private clone frequencies,
            ave_transmitted_f, ave_private_f, ave_preserved_f, ave_lost_f
        '''

        ave_transmitted_f = self.average_over_clones(node_fitness_func, 1., shared=True, private=False, preserved=False,
                                                     lost=False)
        ave_private_f = self.average_over_clones(node_fitness_func, 1, private=True, shared=False, preserved=False,
                                                 lost=False)
        ave_preserved_f = self.average_over_clones(node_fitness_func, 1., preserved=True, lost=False, shared=False,
                                                   private=False)
        ave_lost_f = self.average_over_clones(node_fitness_func, 1, lost=True, shared=False, private=False,
                                              preserved=False)
        return ave_transmitted_f, ave_private_f, ave_preserved_f, ave_lost_f

    def get_shared_volume(self, beta=1.0):
        '''
        Shared volume between time points

        :param beta: float
            tree weighting parameter

        :return: float
        '''
        q = self.average_tree_function(Tree.get_shared_volume, beta=beta)
        return q

    def get_private_volume(self, beta=1.0):
        '''
        Private volume to this sample in this time point

        :param beta: float
            tree weighting parameter

        :return: float
        '''
        q = self.average_tree_function(Tree.get_private_volume, beta=beta)
        return q

    def P53report(self):
        '''
        Some reporting on P53 mutations
        :return:
        '''
        if len(self.trees) > 0:
            p53rep = self.trees[0].P53_report()
        else:
            p53rep = self.oneTree.P53_report()
        return p53rep

    def mutated_gene_Taf(self, gene):
        '''

        :param gene: str
        :return: float
        '''
        taf = 0.0
        mids = self.__gene2mut_id[gene]
        tafs = [self.mutations[mid].Taf for mid in mids]
        if len(tafs) > 0:
            taf = max(tafs)
        return taf

    def mutated_gene_Tcov(self, gene):
        '''

        :param gene: str
        :return: float
        '''
        taf = 0.0
        mids = self.__gene2mut_id[gene]
        tafs = [self.mutations[mid].Tcov for mid in mids]
        if len(tafs) > 0:
            taf = max(tafs)
        return taf

    '''
    Meta functions
    '''

    def average_over_clones(self, node_fun, beta=1.0, clonal=True,
                            shared=False, private=False, preserved=False, lost=False, **kwargs):
        '''
        Averages the function over the designated part of the tumor. The respective clone frequencies
        are precomputed.

        :param node_fun: function
            Node class method

        :param beta: float
            tree weigting parameter

        :param clonal: bool
            1 - use clonality, 0 - don't

        :param shared: bool
            use shared clones

        :param private: bool
            use private (new) clones

        :param preserved: bool
            use preserved clones

        :param lost: bool
            use clones that were lost

        :param kwargs: dict
            parameters to the node_fun method

        :return: float
        '''
        trees = self.trees if clonal else [self.oneTree]
        if clonal:
            weights = self.get_tree_weights(beta)
        else:
            weights = [1.]

        q = [w * tree.average_over_nodes(node_fun, shared=shared, private=private, preserved=preserved, lost=lost,
                                         **kwargs) for (w, tree) in zip(weights, trees)]
        q = sum(q)
        return q

    def average_flux_over_clones(self, node_fun, beta=1.0, **kwargs):
        '''
        Averages the function over the designated part of the tumor. The respective clone frequencies
        are precomputed.

        :param node_fun: function
            Node class method

        :param beta: float
            tree weigting parameter

        :param kwargs: dict
            parameters to the node_fun method

        :return: float
        '''
        trees = self.trees
        weights = self.get_tree_weights(beta)

        q = [w * tree.average_flux_over_nodes(node_fun, **kwargs) for (w, tree) in zip(weights, trees)]
        q = sum(q)
        return q


    def average_tree_function(self, tree_fun, beta=1.0, **kwargs):
        '''
        Averages over trees

        :param tree_fun: function
            cfit.tree.Tree method

        :param beta:
            tree weighting parameter

        :param kwargs: dict
            parameters to tree_fun

        :return: float
        '''
        weights = self.get_tree_weights(beta)
        q = [w * tree_fun(tree, **kwargs) for (w, tree) in zip(weights, self.trees)]
        q = sum(q)
        return q

    def average_ntau_over_trees(self, node_fitness_fun, tau, beta=1.0, **kwargs):
        '''
        Computes n(tau) average over the trees

        :param node_fitness_fun: function
            cfit.tree.node.Node class method to compute fitness

        :param tau: float
            the time parameter

        :param beta: float
            the tree weighting parameter

        :param kwargs: dict
            parmaters to the node_fitness_fun function

        :return: float
        '''
        weights = self.get_tree_weights(beta)
        q = [w * tree.ntau(node_fitness_fun, tau, **kwargs) for (w, tree) in zip(weights, self.trees)]
        q = sum(q)
        return q

    #    def meta_function(self, function, **kwargs):
    #        for tree in self.trees:
    #            tree.meta_function(function=function, **kwargs)

    def meta_tree_function(self, tree_function, **kwargs):
        for tree in self.trees:
            tree_function(tree, **kwargs)

    '''
    Statistics
    '''

    def get_stats(self, beta=1.0):
        '''

        :param beta: float

        :return: list
        '''
        mload = self.TMB(beta=beta)
        nload = self.neoantigen_load(beta=beta)
        ent = np.exp(self.entropy(beta=beta))
        nclones = self.average_tree_function(lambda tree: len(tree.nodes) - 1)
        # stats.columns = ["sample", "m", "n", "ave_m", "ave_n", "nclones", "exp_entropy"]
        return [self.name, len(self.mutations), len(self.neoantigens), mload, nload, nclones, ent]

    def average_llh(self, beta=1.0):
        '''
        The average log-likelihood of the trees
        :param beta:
            the tree weighting parameter
        :return: float
        '''
        weights = self.get_tree_weights(beta)
        q = [w * tree.llh for (w, tree) in zip(weights, self.trees)]
        q = sum(q)
        return q

    def get_intratumor_fitness_std(self, beta=1):
        '''
        Fitness standard deviation on the trees. (old)
        :param beta: float
            the tree weighting parameter
        :return: float
        '''
        weights = self.get_tree_weights(beta)
        std = sum(list(map(lambda el: el[0] * np.sqrt(el[1].variance), zip(weights, self.trees))))
        return std

    def get_intratumor_fitness_variance(self, beta=1):
        '''
        Fitness variance on the trees. (old)
        :param beta: float
            the tree weighting parameter
        :return: float
        '''

        weights = self.get_tree_weights(beta)
        std = sum(list(map(lambda el: el[0] * el[1].variance, zip(weights, self.trees))))
        return std

    '''
    Patient classification statistics
    '''

    def compute_population_fitness(self, clonal=True, beta=1, tau=0.001):
        '''
        Compute the average tumor population fitness;

        :param clonal: bool
            if True use heterogenous model predictions
        :param beta: float
            the tree weighting parameter
        :param tau: float
            the time parameter
        :return: float
        '''
        if clonal:
            f = self.average_tree_function(SampleTree.get_population_fitness, beta=beta, tau=tau)
        else:
            f = self.oneTree.get_population_fitness(tau)
        self.__q = f
        return f

    def compute_predicted_population_size(self, clonal, beta, tau):
        '''
        Compute n(tau), averaged over trees.

        :param clonal: bool
            if True use heterogenous model predictions

        :param beta: float
            the tree weighting parameter

        :param tau: float
            the time parameter

        :return: float
        '''

        if clonal:
            s = self.average_tree_function(SampleTree.get_predicted_population_size, beta=beta, tau=tau)
        else:
            s = self.oneTree.get_predicted_population_size(tau)
        self.__q = s
        return s

    def compute_predicted_population_size_derivative(self, clonal, beta, tau):
        '''
        Compute n(tau) derivative over tau, averaged over trees

        :param clonal: bool
            if True use heterogenous model predictions

        :param beta: float
            the tree weighting parameter

        :param tau: float
            the time parameter

        :return: float
        '''
        s = 0
        if clonal:
            s = self.average_tree_function(SampleTree.get_predicted_population_size_derivative,
                                           beta=beta, tau=tau)
        else:
            s = self.oneTree.get_predicted_population_size_derivative(tau)
        self.__q = s
        return s

    def compute_predicted_integrated_population_size(self, clonal, beta, tau):
        '''

        :param clonal: bool
            if True use heterogenous model predictions

        :param beta: float
            the tree weighting parameter

        :param tau: float
            the time parameter

        :return: float
        '''
        if clonal:
            s = self.average_tree_function(SampleTree.get_predicted_integrated_population_size,
                                           beta=beta, tau=tau)
        else:
            s = self.oneTree.get_population_fitness(tau)
        self.__q = s
        return s

    '''
    Other
    '''

    def set_tree_self_copies(self):
        '''
        Copy trees
        :return:
        '''
        for tree in self.trees:
            tree.set_self_copy()

    '''
    Output
    '''

    def write_alignments(self, of):
        '''
        Old
        :param of: file handle, open for writing
        :return:
        '''
        for neo in self.neoantigens:
            for ae in neo.alignedEpitopes.values():
                mutpos = 1
                for i in range(0, 9):
                    if neo.wtPeptide[i] != neo.mtPeptide[i]:
                        mutpos = i + 1
                line = [neo.mid, self.name, neo.id, ae.epitopeName, ae.mutScoreSW, ae.wtScoreSW, ae.mutSW, ae.wtSW,
                        neo.kD, neo.wtkD, mutpos]
                line = "\t".join(list(map(str, line)))
                of.write(line + "\n")

    def write_neoantigens(self, opath, exclusive=False, clonal=1):
        '''


        :param opath: str
            output file path

        :param exclusive: bool
            report only exclusive neoantigens

        :param clonal: bool
            if true reports clonal structures of the top scoring tree

        :return:
        '''
#        self.logger("Writing neoantigens for sample "+self.name)
        if clonal == 1:
#            self.logger("number of trees: "+str(len(self.trees)))
#            for i, tree in enumerate(self.trees):
#                self.logger("\t"+str(i)+": "+str(len(tree.nodes)))
#                nodes = tree.nodes.values()
#                self.logger([len(node.node.dmutations) for node in nodes])

            self.trees[0].write_neoantigens(opath, sampleName=self.name, exclusive=exclusive, nq=self.neoantigenQualities)
        else:
            self.oneTree.write_neoantigens(opath, sampleName=self.name, exclusive=exclusive, nq=self.neoantigenQualities)

    def write_indel_neoantigens(self, opath, exclusive=False, clonal=1):
        '''
        Writes basic indel neoantigen data to a file

        :param opath: str
            output file path

        :param exclusive: bool
            report only exclusive neoantigens

        :param clonal: bool
            if true reports clonal structures of the top scoring tree

        :return:
        '''
        if clonal == 1:
            self.trees[0].write_indel_neoantigens(opath, sampleName=self.name, exclusive=exclusive)
        else:
            self.oneTree.write_indel_neoantigens(opath, sampleName=self.name, exclusive=exclusive)

    def rank_neoantigens(self, opath, writeheader=True, clonal=1):
        '''


        :param opath: str
            output file path

        :param exclusive: bool
            report only exclusive neoantigens

        :param writeheader: bool
            include header

        :param clonal: bool
            if true reports clonal structures of the top scoring tree

        :return:
        '''

        if clonal == 1:
            self.trees[0].write_neoantigens(opath, sampleName=self.name, exclusive=exclusive, writeheader=writeheader)
        else:
            self.oneTree.write_neoantigens(opath, sampleName=self.name, exclusive=exclusive, writeheader=writeheader)

    def get_mutation_data(self, exclusive=True, clonal=True):
        '''
        Returns a DataFrame with data about mutations in the sample,
        grouped by clones.

        :param exclusive: bool
            if True will report each mutation only once,
            in the clone where it originates
        :param clonal: bool
            if True/1 will report mutation data grouped by clones
            of the highest scoring tree, otherwise will assume homogeneity

        :return: pd.DataFrame
        '''

        mut2CCF = self.get_mutation_frequencies(beta=1, exclusive=False)

        mdata = [[mid, self.mutations[mid].gene, mut2CCF[mid],
                  self.mutations[mid].Taf, self.mutations[mid].Naf,
                  self.mutations[mid].Tcov, self.mutations[mid].Ncov
                  ] for mid in mut2CCF.keys()]
        mdata = pd.DataFrame(mdata)
        mdata.columns = ["Mutation", "Gene", "CCF", "Taf", "Naf", "Tcov", "Ncov"]

    #        if clonal:
#            mdata = self.trees[0].get_mutation_data(sampleName=self.name, exclusive=exclusive)
#        else:
#            mdata = self.oneTree.write_mutations(sampleName=self.name, exclusive=exclusive)
#        mdata["Taf"] = [self.mutations[mid].Taf for mid in mdata.Mutation]
#        mdata["Naf"] = [self.mutations[mid].Naf for mid in mdata.Mutation]
#        mdata["Tcov"] = [self.mutations[mid].Tcov for mid in mdata.Mutation]
#        mdata["Ncov"] = [self.mutations[mid].Ncov for mid in mdata.Mutation]
        return mdata

    def write_mutations(self, opath, exclusive=True, clonal=True):
        '''
        Writes mutations from a sample to a file

        :param opath: str
            output file path

        :param exclusive: bool
            if True will report each mutation only once,
            in the clone where it originates

        :param clonal: bool
            if True/1 will report mutation data grouped by clones
            of the highest scoring tree, otherwise will assume homogeneity
        :return:
        '''

        self.logger(self.name + " writing mutations")
        mdata = self.get_mutation_data(exclusive, clonal)
        mdata.to_csv(opath, index=False, sep="\t")

    def get_tree_clone(self, tree_id, clone_id):
        '''

        Return clone frequency in a given tree

        :param tree_id: int
            tree index
        :param clone_id: int
            clone index
        :return: float, float
            clone frequency, tree weight

        '''

        node = self.trees[tree_id].nodes[clone_id]
        return node

    def get_tree_clone_frequency(self, tree_id, clone_id, exclusive=False):
        '''

        :param tree_id: int
        :param clone_id: int
        :param exclusive: bool
        :return: float
        '''

        node = self.get_tree_clone(tree_id, clone_id)
        freq = node.Y if exclusive else node.X
        return freq

    def get_mutation_frequencies(self, beta=1.0, exclusive=False, nonsynonymous_only=False, kd_threshhold=None):
        '''
        Create a dictionary of averaged mutation frequencies

        :param beta: float
            tree weighting parameter

        :param exclusive: bool
            whether to report Y or X.

        :param kd_threshhold: float
            threshold on neoantigen kd to include mutations with neoantigesns

        :return: dict (defaultdict(float)
            mutation identifiers mapped to the averaged frequencies (over trees)
        '''

        weights = self.get_tree_weights(beta)
        mut2CCF = {}#defaultdict(float)
        for w, tree in zip(weights, self.trees):
            for nid in tree.nodes:
                samplenode = tree.nodes[nid]
                node = samplenode.node
                for mut in node.exclusiveMutations:
                    if nonsynonymous_only and not mut.is_nonsynonymous():
                        continue

                    if kd_threshhold is not None:
                        neos = self.mutation2neoantigens[mut.id]
                        if len(neos) == 0:
                            continue
                        kds = [neo.kD for neo in neos]
                        wtkds = [neo.wtkD for neo in neos]
                        if min(kds) > kd_threshhold and min(wtkds) > kd_threshhold:
                            continue
                    if mut.id not in mut2CCF:
                        mut2CCF[mut.id] = 0
                    mut2CCF[mut.id] += w * (samplenode.Y if exclusive else samplenode.X)
        return mut2CCF

    def get_mutation_fitness(self, beta=1.0):
        '''
        Create a dictionary of averaged mutation fitness

        :param beta: float
            tree weighting parameter

        :param exclusive: bool
            whether to report Y or X.

        :return: dict: str->float, dict: str->float
            mutation identifiers mapped to the averaged fitness, relative fitness,  (over trees)
        '''

        weights = self.get_tree_weights(beta)
        mut2fitness = defaultdict(float)
        mut2relfitness = defaultdict(float)
        for w, tree in zip(weights, self.trees):
            for nid in tree.nodes:
                samplenode = tree.nodes[nid]
                node = samplenode.node
                for mut in node.exclusiveMutations:
                    mut2fitness[mut.id] += w * samplenode.fitness
                    mut2relfitness[mut.id] += w * samplenode.rfitness

        return mut2fitness, mut2relfitness

    def write_CCFs(self, opath, beta=1.0):
        '''
        Writes inclusive clone frequencies to a file

        :param opath: str
            output file path

        :param beta: float
            tree weighting parameter

        :return:
        '''

        mut2CCF = self.get_mutation_frequencies(beta=beta)
        muttab = [[mid, mut2CCF[mid]] for mid in mut2CCF.keys()]
        muttab = pd.DataFrame(muttab)
        muttab.columns = ["mutation", "CCF"]
        muttab.to_csv(opath, sep="\t", index=False)

    def toVCF(self, outdir, trunkal=False):
        '''
        Write vcf files of clone nodes to files in outdir.
        The files are written to the outdir folder and
        are named as <samplename>_<cloneid>.vcf

        :param outdir: str
            output directory

        :param trunkal: bool
            If True will output 2 files only, one for clone 1 (trunk)
            and one for the remaining clones
        :return:
        '''

        if len(self.trees) > 0:
            besttree = self.trees[0]
        else:
            besttree = self.oneTree
        vcfjs = besttree.toVCF()
        if trunkal:
            vcf = vcfjs[1]
            if vcf is not None:
                filename = os.path.join(outdir, self.name + "_clonal.vcf")
                vcf.to_csv(filename, sep="\t", index=False)
            cloneids = list(vcfjs.keys())
            cloneids = filter(lambda cid: cid > 1, cloneids)
            if len(cloneids) > 0:
                vcf = None
                for cid in cloneids:
                    if vcf is None:
                        vcf = vcfjs[cid]
                    else:
                        vcf = vcf.append(vcfjs[cid])
                if vcf is not None:
                    vcf = vcf.drop_duplicates()
                    filename = os.path.join(outdir, self.name + "_subclonal.vcf")
                    vcf.to_csv(filename, sep="\t", index=False)
        else:
            for cloneid in vcfjs:
                vcf = vcfjs[cloneid]
                if vcf is not None:
                    filename = os.path.join(outdir, self.name + "_" + str(cloneid) + ".vcf")
                    vcf.to_csv(filename, sep="\t", index=False)

    def score_neoantigens(self, beta=1, tau=1, just_dominant=False, clonal=True):
        '''
        Score neoantigens in the trees and their nodes as log(node.X) -neo.fitness*tau (just_dominant=False)
            or log(node.X)+ neo.fitness*tau if dominant else 0
        Average the score over trees

        :param beta: float
            tree weighting parameter

        :param tau: float

        :param just_dominant: bool

        :param clonal: bool
            True - use reconstructed clonal structures, False - homogenous structure

        :return: pd.DataFrame
            table with scores for neoantigens
        '''

        if clonal:
            weights = self.get_tree_weights(beta)
            trees = self.trees
        else:
            weights = [1.]
            trees = [self.oneTree]
        dscores_major = defaultdict(float)
        dscores_minor = defaultdict(float)

        dX = defaultdict(float)
        dY = defaultdict(float)
        dclone_fitness = defaultdict(float)

        #print(f"scoring sample {sample.name}")
        for (w, tree) in zip(weights, trees):
            dscores1 = tree.score_neoantigens(neo_qualities=self.neoantigenQualities,
                                              tau=tau, just_dominant=just_dominant)
            for node in tree.nodes.values():
                nids = set([neo.id for neo in node.exclusiveNeoantigens])
                for nid in nids:
                    dscores_major[nid] += w * dscores1[nid]

            for node in tree.nodes.values():
                nids = set([neo.id for neo in node.exclusiveNeoantigens])
                for nid in nids:
                    dX[nid] += w * node.X
                    dY[nid] += w * node.Y
                    dclone_fitness[nid] = w * node.fitness

        if just_dominant:
            # prepare secondary ranking
            for (w, tree) in zip(weights, self.trees):
                dscores1 = tree.score_neoantigens(neo_qualities=self.neoantigenQualities,
                                                  tau=tau, just_dominant=False)
                for node in tree.nodes.values():
                    nids = set([neo.id for neo in node.exclusiveNeoantigens])
                    for nid in nids:
                        dscores_major[nid] += w * dscores1[nid]
#                for nid in dscores1:
#                    dscores_minor[nid] += w * dscores1[nid]
        dscores = {}
        for nid in dscores_major:
            dscores[nid] = (dscores_major[nid] + 0.0001 * dscores_minor[nid]) / 1.0001

        neo_order = [nid for nid in dscores]
        neo_order.sort(key=lambda nid: -dscores[nid])
        columns1 = ["Neoantigen", "Sample", "Gene", "peptideMT", "peptideWT", "kdMT", "kdWT", "Quality"]
        columns2 = ["X", "Y", "Clone_fitness", "Predicted_impact"]

        lines = []
        columns = columns1 + columns2
        for nid in neo_order:
            ninfo = self.neoantigen_info(nid)
            neo = self.neoantigensMap[nid]
            components = list(neo.qattributes.get_quality_components().keys())
            columns = columns1 + components +columns2
            ninfo["X"] = dX[nid]
            ninfo["Y"] = dY[nid]
            ninfo["Clone_fitness"] = dclone_fitness[nid]
            ninfo["Predicted_impact"] = dscores_major[nid]
            line = [ninfo[col] for col in columns]
            lines.append(line)
        lines = pd.DataFrame(lines, columns=columns)
        return lines

    def neoantigen_info(self, nid):
        '''

        :param nid: str
        :return: dict
        '''

        neo = self.neoantigensMap[nid]
        dinfo = {"Neoantigen": neo.id,
                 "Sample": self.name,
                 "Gene": self.mutations[neo.mid].gene,
                 "peptideMT": neo.mtPeptide,
                 "peptideWT": neo.wtPeptide,
                 "kdMT": neo.kD,
                 "kdWT": neo.wtkD,
                 }
        components = neo.qattributes.get_quality_components()
        for comp in components:
            dinfo[comp] = components[comp]
        dinfo['Quality'] = self.neoantigenQualities[neo.id]
        return dinfo

    def write_trees(self, odir):
        '''

        :param odir: str

        '''

        if not os.path.exists(odir):
            os.mkdir(odir)
        for i, tree in enumerate(self.trees):
            ofile = os.path.join(odir, "tree_" + self.name + "_" + str(i) + ".txt")
            tree.write(ofile)


    def toJSON(self):
        '''
        Creates json for Sibyl app.
        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''

        js = {}
        js['id'] = self.name
        js["sample_trees"] = [tree.toJSON() for tree in self.trees]

        #js["trees"] = jstrees
#        jstree = self.trees[0].toJSON()
#        js['besttree'] = jstree
#        jsneo = {}
#        for neo in self.neoantigens:
#            jsneo[neo.id] = neo.toJSON()
#        js['neoantigens'] = jsneo
#        jsmut = {}
#        for mutid in self.mutations:
#            jsmut[mutid] = self.mutations[mutid].toJSON()
#        js['mutations'] = jsmut
        return js

    def to_sibyl_JSON(self):
        '''
        Creates json for Sibyl app.
        :return: dict
            dictionary that represents the tree and can be written in a json file.
        '''

        js = {}
        js['id'] = self.name
        # js["sample_trees"] = [tree.toJSON() for tree in self.trees]

        # js["trees"] = jstrees
        jstree = self.trees[0].toJSON()
        js['besttree'] = jstree
        jsneo = {}
        for neo in self.neoantigens:
            jsneo[neo.id] = neo.toJSON()
        js['neoantigens'] = jsneo
        jsmut = {}
        for mutid in self.mutations:
           jsmut[mutid] = self.mutations[mutid].toJSON()
        js['mutations'] = jsmut
        return js
    
    def tree_nodes_to_sibyl(self):
        return [{'tree_id': id, 'sample_tree_nodes': tree.tree_nodes_to_sibyl()} for id, tree in enumerate(self.trees)]
          
        # print(self.name, '----------------------')
        # res = [{'tree_id': id, 'sample_tree_nodes': tree.tree_nodes_to_sibyl()} for id, tree in enumerate(self.trees)]
        # # res = self.trees[0].tree_nodes_to_sibyl(0)
        # print(res)
        # return res
    
