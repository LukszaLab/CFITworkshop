'''
Created on Jan 1, 2020

@author: mluksza
'''

import os
from collections import defaultdict

import numpy as np

from cfit.CoreObject import CoreObject
from cfit.tree.node.SampleNode import SampleNode


class TimePoint(CoreObject):
    '''
    Class implementing a time point of the tumor in a patient, example primary / met, pre / post therapy

    Attributes:

        name: str
            time point name, eg. Pre, Post, Primary

        __samples: dict
            maps sample names to cfit.patient.Sample objects

    '''

    def __init__(self, name=None):
        '''
        Constructor

        :param name: str
            name of the timepoint
        '''
        self.name = name
        self.__samples = {}

    @property
    def samples(self):
        return self.__samples

    @samples.setter
    def samples(self, samples):
        self.__samples = samples

    @property
    def tissues(self):
        return [sample.tissue for sample in self.__samples.values()]

    @property
    def sampleIDs(self):
        return ",".join([sample.name for sample in self.__samples.values()])

    def all_weighted_nodes(self, beta=1):
        '''

        :param beta: float
        :return: list
            (SampleNode, int, float)
            (node, tree index, weight)
        '''
        nodes = []
        for sample in self.samples.values():
            weights = sample.get_tree_weights(beta=beta)
            for tind, (tree, weight) in enumerate(zip(sample.trees, weights)):
                for node in tree.nodes.values():
                    nodes.append((node, tind, weight))

#        for (node, tind, weight) in nodes:
            #try:
#            self.logger(self.name+"\t"+str([(tind, node.id, node.cY, weight)]))
            #except:
            #    self.logger(self.name+"\tERROR"+str([(tind, node.id, weight)]))

        return nodes

    def add_sample(self, sample):
        '''
        Add sample to the time point.

        :param sample: cfit.patient.Sample
        '''
        self.__samples[sample.name] = sample

    def add_neoantigen(self, neo):
        '''
        Adds neoantigen to the samples in the time point

        :param neo: cfit.tree.mutation.Neoantigen

        '''
        for sname in self.__samples:
            sample = self.__samples[sname]
            sample.add_neoantigen(neo)

    def add_frame_shift_neoantigen(self, neo):
        '''
        Adds neoantigen to the samples in the time point

        :param neo: cfit.tree.mutation.FrameShiftNeoantigen

        '''
        for sname in self.__samples:
            sample = self.__samples[sname]
            sample.add_frame_shift_neoantigen(neo)

    def trees(self, num=0, just_one=False, sample_name=None):
        '''
        Returns the list of trees ranked at number "num" from all samples. All samples have the same tree topology,
        but different frequencies of clones.

        :param num: int
            index (rank) of the tree. if -1 return all trees

        :param just_one: bool
            just one sample from the time point

        :param sample_name: str
            sample name (optional)

        :return: list
            list of cfit.tree.SampleTree objects

        '''
        trees = []
        samples = self.__samples
        if just_one:
            samples = [list(self.__samples.keys())[0]]
        if sample_name is None:
            for sname in samples:
                sample = self.__samples[sname]
                if num >= 0:
                    trees.append(sample.trees[num])
                else:
                    for tree in sample.trees:
                        trees.append(tree)
        else:
            sample = self.__samples[sample_name]
            if num >= 0:
                trees.append(sample.trees[num])
            else:
                for tree in sample.trees:
                    trees.append(tree)
        return trees

    #    def set_mutation_node_index(self):
    #        '''
    #        Sets the mutation-node index, for faster access to node mutation content
    #        '''
    #        alltrees = self.trees(num=-1)
    #        for tree in alltrees:
    #            tree.set_mutation_node_index()

    def fitness(self, beta=1, shared=False, private=False):
        '''

        :param shared: bool
            use shared clones

        :param private: bool
            use private clones

        :return: float
        '''

        avefs = [sample.average_over_clones(SampleNode.get_fitness, beta=beta, shared=shared, private=private) for
                 sample in self.samples.values()]
        return np.mean(avefs)

    def fitness_flux(self, beta=1):
        '''

        :param shared: bool
            use shared clones

        :param private: bool
            use private clones

        :return: float
        '''

        avefs = [sample.average_flux_over_clones(SampleNode.get_fitness, beta=beta) for
                 sample in self.samples.values()]
        return np.mean(avefs)

    def get_tree_clone_frequency(self, tree_id, clone_id, exclusive=False):
        '''

        Return effective clone frequency (averaged over samples)

        :param tree_id: int
            tree index
        :param clone_id: int
            clone index
        :param exclusive: bool
        :return: float
        '''
        nodes = [sample.get_tree_clone(tree_id, clone_id) for sample in self.samples.values()]
        avefs = [node.Y if exclusive else node.X for node in nodes]
        return np.mean(avefs)

    def get_tree_clone_fitness(self, tree_id, clone_id, absolute=True):
        '''

        :param tree_id: int

        :param clone_id: int

        :param absolute: bool
        :return:
        '''
        nodes = [sample.get_tree_clone(tree_id, clone_id) for sample in self.samples.values()]
        avefs = [node.fitness if absolute else node.rfitness for node in nodes]
        return np.mean(avefs)

    def compute_fitness(self, params, include_tau=True):
        '''
        Compute fitness and rfitness attributes for all nodes in all trees
        :param params: dict
            {"tau": None/float,
            "weights":
                {<component_name>: str -> <component_value>: None/float}
            }
            dictionary defining the set of parameters to optimize over or to fix
        '''

        #        self.logger(params)
        trees = self.trees(num=-1)
        for tree in trees:
            nodes = tree.nodes.values()
            for node in nodes:
                node.fitness = sum(
                    [node.fitness_components[name] * params["weights"][name] for name in params["weights"]])
                if include_tau:
                    node.fitness *= params["tau"]
                #self.logger(node.fitness)
            avef = sum([node.fitness * node.Y for node in nodes])
            for node in nodes:
                node.rfitness = node.fitness - avef

    def number_of_new_clones(self):
        '''

        :return: float
        '''

        tree_fun_n_private_clones = lambda tree: sum([(node.privateY > 0) for node in tree.nodes.values()])
        n = [sample.average_tree_function(tree_fun_n_private_clones, beta=1) for sample in self.samples.values()]
        return np.mean(n)

    def average_over_clones(self, node_fun, beta=1.0,
                            shared=False, private=False, **kwargs):
        '''
        Averages the function over the designated part of the tumor. The respective clone frequencies
        are precomputed.

        :param node_fun: function
            Node class method

        :param beta: float
            tree weigting parameter

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

        avefs = [sample.average_over_clones(node_fun, beta=beta, shared=shared, private=private, **kwargs)
                 for sample in self.samples.values()]
        return np.mean(avefs)

    def set_tree_self_copies(self, minimal=True):
        '''
        Copies the trees in all samples - sets Tree.copy attribute in the trees.

        :param minimal: bool

        '''
        trees = self.trees(-1)
        for tree in trees:
            if minimal:
                tree.set_minimal_self_copy()
            else:
                tree.set_self_copy()

    #        for sname in self.__samples:
    #            sample = self.__samples[sname]
    #            sample.set_tree_self_copies()

    def get_mutation_frequencies(self, beta=1.0, exclusive=False, nonsynonymous_only=False, kd_threshhold=None):
        '''
        Get the list of averaged mutation frequencies in the time point

        :param beta: float
            tree weighting parameter

        :param exclusive: bool
            whether to report Y or X.

        :param kd_threshhold: float
            threshold on neoantigen kd. Only mutations with neoantigens of kd below the threshold will be included
        :return: dict
            dictionary mapping mutation identifiers to their averaged frequencies.
        '''

        mids = list(set([mid for sample in self.__samples.values() for mid in list(sample.mutations.keys())]))
        mut2CCF = {}
        for mid in mids:
            mut2CCF[mid] = 0.0
        for sample in self.__samples.values():
            mc = sample.get_mutation_frequencies(beta=beta, exclusive=exclusive,
                                                 nonsynonymous_only=nonsynonymous_only,
                                                 kd_threshhold=kd_threshhold)
            for mid in mc:
                mut2CCF[mid] += mc[mid]
        n = len(self.__samples)
        for mid in mids:
            mut2CCF[mid] /= n

        return mut2CCF

    def average_over_samples(self, sample_fun, **kwargs):
        '''
        Meta function to average arbitrary function over samples

        :param sample_fun: function
            Sample class method or function that takes a Sample class object as the first argument.

        :param kwargs: dict
            Parameters passed to samle_fun

        :return: float

        '''
        ave = np.mean([sample_fun(sample, **kwargs) for sample in self.samples.values()])
        return ave

    def average_over_sample_trees(self, tree_fun, beta=1, **kwargs):
        '''

        :param tree_fun: function

        :param kwargs: dict

        :return: float
        '''

        ave = np.mean([sample.average_tree_function(tree_fun, beta=beta, **kwargs) for sample in self.samples.values()])
        return ave

    def write_trees(self, odir):
        '''

        :param odir: str

        '''

        if not os.path.exists(odir):
            os.mkdir(odir)
        tpodir = os.path.join(odir, self.name)
        for sample in self.samples.values():
            sample.write_trees(tpodir)

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
        gmut2CCF0 = defaultdict(list)
        for sample in self.samples.values():
            mut2CCF = sample.get_mutation_frequencies(beta, exclusive, nonsynonymous_only, kd_threshhold)
            for mid in mut2CCF:
                gmut2CCF0[mid].append(mut2CCF[mid])
        gmut2CCF = defaultdict(float)
        for mid in gmut2CCF0:
            gmut2CCF[mid] = np.mean(gmut2CCF0[mid])
        return gmut2CCF

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

        gmut2fitness0 = defaultdict(list)
        gmut2relfitness0 = defaultdict(list)

        for sample in self.samples.values():
            mut2fitness, mut2relfitness = sample.get_mutation_fitness(beta)
            for mid in mut2fitness:
                gmut2fitness0[mid].append(mut2fitness[mid])
                gmut2relfitness0[mid].append(mut2relfitness[mid])

        gmut2fitness = defaultdict(float)
        gmut2relfitness = defaultdict(float)
        for mid in gmut2fitness0:
            gmut2fitness[mid] = np.mean(gmut2fitness0[mid])
            gmut2relfitness[mid] = np.mean(gmut2relfitness0[mid])
        return gmut2fitness, gmut2relfitness


    def toJSON(self):
        js = {"name": self.name}
        js["samples"] = [sample.toJSON() for sample in self.samples.values()]
        return js
    
    # def get_samples_sibyl(self, time_point_key):
    #     res = [{ 'time_point_id': time_point_key, 'sample_id': sample.name, 'nodes': sample.tree_nodes_to_sibyl()} for sample in self.__samples.values()]
    #     print('--------------')
    #     print(res)
    #     return res
    
    
    def get_samples_sibyl(self, time_point_key):
        res = []
        for sample in self.__samples.values():
            for samples in sample.tree_nodes_to_sibyl():
                # if 'sample_tree_nodes' in samples:
                nodes = samples['sample_tree_nodes']  # Extract sample_tree_nodes
                tree_id = samples['tree_id']
                sibyl_data = {
                    'time_point_id': time_point_key,
                    'sample_id': sample.name,
                    'tree_id': tree_id,
                    'nodes': nodes  # Assign nodes to 'nodes' key in the result
                }
                res.append(sibyl_data)

        # print('--------------')
        # print(res)
        return res
    
    
    
    
    # def get_samples_sibyl(self, time_point_key):
    #     res = []
    #     for sample in self.__samples.values():
    #         for tree_id, sample_tree_nodes in enumerate(sample.tree_nodes_to_sibyl(), 1):
    #             if 'nodes' in sample_tree_nodes:
    #                 for node in sample_tree_nodes['nodes']:
    #                     sibyl_data = {
    #                         'time_point_id': time_point_key,
    #                         'sample_id': sample.name,
    #                         'tree_id': tree_id,
    #                         'nodes': node
    #                     }
    #                     res.append(sibyl_data)
    #     print('--------------')
    #     print(res)
    #     return res
