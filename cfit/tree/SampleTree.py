'''
Created on July 31, 2021

@author: mluksza
'''
from collections import defaultdict

import numpy as np
import pandas as pd

from cfit.CoreObject import CoreObject
from cfit.tree.node.Node import Node
from cfit.tree.node.SampleNode import SampleNode
from cfit.util.Utils import Utils


class SampleTree(CoreObject):
    '''
    Class for representing tree properties in a given sample.

    Attributes:
        tree: cfit.tree.Tree
            tree topology

        copy: cfit.tree.SampleTree
            copied object, used eg. for predictions

        nodes dict: int -> cfit.tree.node.SampleNode
            dictionary mapping nodes on the tree to sample-specific nodes SampleNode

        neoantigenQualities: dict: str->float
            dictionary mapping neoantigen id to its computed quality in that sample
    '''

    def __init__(self, tree, param=None, ord=-1, format='phylowgs'):
        '''
        :params tree: cfit.tree.Tree

        :param param: arbitrary
            object storing frequency information of clones in trees

        :param ord: int
            sample number

        :param format: str
            handling phylowgs/pairtree

        '''

        self.nodes = {}
        self.purity = 1
        ####################################
        # sample specific stuff
        ####################################
        # purity
        self.tree = tree
        self.processed = False

        if format.lower() == 'phylowgs':
            self.__init_phylowgs(tree, param, ord)
        elif format.lower() == 'pairtree':
            tid, npz = param
            self.__init_pairtree(tree, tid, npz, ord)
        self.neoantigenQualities = defaultdict(float)


    def __init_phylowgs(self, tree, jstree=None, ord=-1):
        '''
        Constructor

        :param tree: cfit.tree.Tree

        :param jstree: dict (json)
            from phylowgs output, file storing frequencies of clones

        :param ord: int
            Sample number for multiple-sample reconstructions

        '''

        if jstree is not None:
            self.tree = tree
            roots = jstree["structure"]["0"]
            self.purity = sum([jstree["populations"][str(root)]["cellular_prevalence"][ord] for root in roots])

            # Frequencies:
            # Fill in inclusive (X) and exclusive (Y) clone frequencies.
            # Inclusive frequencies are synonymous to the cellular cancer fraction,
            # they include the frequencies of nested clones.
            # Exclusive frequencies sum up to 1, they exclude nested clones.
            # 1. Find sample purity and normalize CCFs so that the the root frequency is 1.

            frequencies = {}
            populations = jstree['populations']
            for node in populations:
                cid = int(node)
                frequencies[cid] = populations[node]['cellular_prevalence'][ord]  # [num]

            for node in populations:
                cid = int(node)
                tree_node = self.tree.nodes[cid]
                samnode = SampleNode(tree_node)
                frequencies[cid] /= self.purity
                if cid == 0:
                    frequencies[cid] = 1.0
                samnode.X = frequencies[cid]
                self.nodes[cid] = samnode

            # exclusive frequencies: Y_i = X_i - sum_{c \in children(i)} X_c
            for node in self.tree.nodes.values():
                nid = node.id
                samnode = self.nodes[nid]
                csamnodes = [self.nodes[cnode.id] for cnode in node.children]
                samnode.Y = samnode.X - sum([cnode.X for cnode in csamnodes])
        
        else:  # one node tree
#            self.tree = tree.one_node_tree()
            for nid in self.tree.nodes:
                self.nodes[nid] = SampleNode(self.tree.nodes[nid])
            self.nodes[0].X = 1.0
            self.nodes[1].X = 1.0
            self.nodes[0].Y = 0.0
            self.nodes[1].Y = 1.0

    @property
    def mutation_node_index(self):
        return self.tree.mutation_node_index

    @property
    def variance(self):
        return self.tree.variance

    @property
    def F0(self):
        return self.tree.F0

    @property
    def llh(self):
        return self.tree.llh

    @property
    def root(self):
        root = self.nodes[self.tree.root.id]
        return root


    def __init_pairtree(self, tree, tid=0, npz=None, ord=-1):
        #replace tid and npz params with one parameter that just contains clone frequencies of tree
        '''
        Constructor

        :param tree: cfit.tree.Tree

        :param tid: int
            index of tree from pairtree output

        :param npz: npz file
            from pairtree output, contains all info about reconstructed trees
        
        :param ord: int
            sample number for multiple-sample reconstructions 
        '''
        self.tree = tree
        self.nodes = {}
        self.purity = 1
        if not npz is None:
            trees_data = np.load(npz)
            if tid >= len(trees_data['struct']):
                return
            parent_vector = trees_data['struct'][tid]
            n = len(parent_vector)
            root_children = []
            for i in range(n):
                if parent_vector[i] == 0:
                    root_children += [i+1]
            phi = trees_data['phi'][tid]
            freqs = {k: phi[k][ord] for k in range(n+1)}  # population frequencies provided with pairtree output are inclusive
            self.purity = sum(freqs[j] for j in root_children)
            
           
            for cid in range(n+1):
                tree_node = self.tree.nodes[cid]
                samnode = SampleNode(tree_node)
                freqs[cid] /= self.purity # why do we do this??
                if cid == 0:
                    freqs[cid] = 1.0
                samnode.X = freqs[cid]
                self.nodes[cid] = samnode
            
            # exclusive frequencies: Y_i = X_i - sum_{X_c for c in children(i)}
            for node in self.tree.nodes.values():
                nid = node.id
                samnode = self.nodes[nid]
                csamnodes = [self.nodes[cnode.id] for cnode in node.children]
                samnode.Y = samnode.X - sum([cnode.X for cnode in csamnodes])
                #print(f"Sample node {nid}: Y = {samnode.Y}, X = {samnode.X}")

        else:  # one node tree
#            self.tree = tree.one_node_tree()
            for nid in self.tree.nodes:
                self.nodes[nid] = SampleNode(self.tree.nodes[nid])
            self.nodes[0].X = 1.0
            self.nodes[1].X = 1.0
            self.nodes[0].Y = 0.0
            self.nodes[1].Y = 1.0

            

    def one_node_tree(self, thepatient):
        '''
        Creates a tree with a singleton clone.

        :param thepatient: cfit.patient.Sample
            The sample for which the tree is reconstructed.
        '''

        self.tree = thepatient.trees[0].one_one_tree()
        self.purity = self.tree.purity
        self.nodes[0].X = 1.0
        self.nodes[1].X = 1.0
        self.nodes[0].Y = 0.0
        self.nodes[1].Y = 1.0

    def path_to_root(self, nid):
        '''
        Return path of nodes to the root from node nid.

        :param nid: int
            node identifier

        :return: list
            list of SampleNode objects
        '''

        node = self.nodes[nid]
        rpath = node.node.path_to_root()
        srpath = [self.nodes[node.id] for node in rpath]
        return srpath

    '''
    Frequencies and averaging
    '''

    def meta_beta_function(self, node_fun=None, **kwargs):
        '''
        A meta method that applies function to each node and averages the output by the exclusive clone frequencies.
        Redundant with self.average_over_nodes(function, False, False, False, False, **kwargs)

        :param node_fun: function
            Node class method to be applied on the node, should return a numeric value

        :param kwargs: dict
            parameters passed to the function

        :return: float
            The averaged statistics over the nodes on the tree
        '''
        q = sum([node_fun(node.node, **kwargs) * node.Y for node in self.nodes.values()])
        return q

    def average_over_nodes(self, node_fun, shared=False, private=False, preserved=False, lost=False, **kwargs):
        '''
        A meta function to compute the average over nodes in the tree

        :param node_fun: function
            Node class method that computes the quantity to be averaged over nodes

        :param shared: bool
            restrict to the nodes that are shared between time points (annotation precomputed)

        :param private: bool
            restrict to the nodes that are private to this tree over time points (annotation precomputed)

        :param preserved: bool
            restrict to the nodes that are preserved in this tree over time points (annotation precomputed)

        :param lost: bool
            restrict to the nodes that are lost (marginal) in this tree over time points (annotation precomputed)

        :param kwargs:
            additional parameters passed to node_fun

        :return: float
            the averaged statistic over the relevant nodes on the tree.

        '''

        def get_y(node):
            if shared:
                return node.sharedY
            elif private:
                return node.privateY  # gained (new) clones
            elif preserved:
                return node.preservedY
            elif lost:
                return node.lostY
            else:
                return node.Y

        values = [get_y(node) * node_fun(node, **kwargs) for node in self.nodes.values()]
        ave = sum(values)
        return ave

    def set_primary_ancestors(self):
        '''
        Set a map of ancestors from the primary tumors.
        The map is store on .anc attribute of self.nodes (of class SampleNode)
        '''

        tnodes = self.tree.get_pre_order() # list of nodes of class Node, in pre-order
        for tnode in tnodes:
            nid = tnode.id
            node = self.nodes[nid]
            if nid == self.root.id or node.sharedY > 0:
                node.anc = nid
            else:
                node.anc = self.nodes[node.node.parent.id].anc

    def average_flux_over_nodes(self, node_fun, **kwargs):
        '''
        A meta function to compute the average over nodes in the tree

        :param node_fun: function
            Node class method that computes the quantity to be averaged over nodes

        :param kwargs:
            additional parameters passed to node_fun

        :return: float
            the averaged statistic over the relevant nodes on the tree.

        '''


#        def ancestor(node):
#            anc = node.node.parent.id
#            ancnode = self.nodes[anc]
#            while ancnode.privateY > 0 and self.root.id != anc: #stop condition: anc is the identifier of the primary tumor ancestor
#                anc = ancnode.node.parent.id
#                ancnode = self.nodes[anc]
#            return ancnode
#        pairs = [(ancestor(node).id, node) for node in self.nodes.values()]

        self.set_primary_ancestors()
        pairs = [(node.anc, node) for node in self.nodes.values()]

        xrec = defaultdict(float)

        for el in pairs:
            xrec[el[0]] += el[1].Y
#        values = [node.Y * ((node_fun(node, **kwargs) - node_fun(ancestor(node), **kwargs)) if node.privateY>0 else 0) for node in self.nodes.values()]

        values = [node.Y * ((node_fun(node, **kwargs) - node_fun(self.nodes[node.anc], **kwargs)) if node.privateY>0 else 0) for node in self.nodes.values()]
        phi2 = sum(values)

        return phi2


    def get_shared_volume(self):
        '''

        :return: float
            the total volume of clones shared across time points (annotation precomputed)
        '''
        return sum([node.Y for node in self.nodes.values() if node.sharedY > 0])

    def get_private_volume(self):
        '''

        :return: float
            the total volume of clones that are private to this tree across time points (annotation precomputed)
        '''
        return sum([node.Y for node in self.nodes.values() if node.privateY > 0])

    def get_trunk_volume(self):
        '''

        :return: float
            total volume of the trunk clones (children of clone 0)
        '''
        cids = self.tree.root.children
        cnodes = [self.nodes[cid] for cid in cids]
        trunk_Y = sum([clone.Y for clone in cnodes])
        return trunk_Y

    def get_trunk_fitness(self, node_fitness_fun, **kwargs):
        '''
        Fitness of the trunk clones, in most cases it is a single clone. If multiple, the fitness is averaged.
        :param node_fitness_fun: function
            the fitness function, method of Node class

        :param kwargs:
            parameters to node_fitness_fun

        :return: float
            the averaged fitness of the trunk clones
        '''
        trunk_Y = self.get_trunk_volume()
        if trunk_Y == 0:  # not sure how this would happen, but just in case
            return 0
        trunk_F = sum([clone.Y * node_fitness_fun(clone, **kwargs) for clone in self.root.children]) / trunk_Y
        return trunk_F

    def get_average_non_trunk_fitness(self, node_fitness_fun, **kwargs):
        '''
        Method complementary to get_trunk_fitness, computes the average fitness of the non-trunk clones.

        :param node_fitness_fun: function
            the fitness function, method of Node

        :param kwargs:
            parameters to node_fitness_fun

        :return: float
            the averaged fitness of the non-trunk clones
        '''

        other_clones = [clone for clone in self.nodes.values() if
                        clone.node.parent.id != self.tree.root.id and clone.node.id != self.tree.root.id]
        rest_Y = sum([clone.Y for clone in other_clones])
        if rest_Y == 0:
            return 0
        rest_F = sum([clone.Y * node_fitness_fun(clone, **kwargs) for clone in other_clones]) / rest_Y
        return rest_F

    def get_clone_Ys(self):
        '''

        :return: dict
            Dictionary mapping clone identifiers to their exclusive frequencies
        '''
        Ys = {}
        for node in self.nodes.values():
            Ys[node.id] = node.Y
        return Ys

    def get_clone_Xs(self):
        '''

        :return: dict
            Dictionary mapping clone identifiers to their inclusive frequencies
        '''
        Xs = {}
        for node in self.nodes.values():
            Xs[node.id] = node.X
        return Xs

    ##############################################
    # Predictions
    ##############################################

    def predicted_tree(self, node_fitness_fun, tau, **kwargs):
        '''
        Creates a tree, that is a copy of self, but with clone frequencies predicted with the given fitness model function

        :param node_fitness_fun: function
            fitness model function, operating on a node

        :param tau: float
            time for the predictions

        :param kwargs:
            parameters passed to node_fitness_fun

        :return: cfit.tree.Tree
            a new tree, with clone frequencies as given by the prediction with the fitness model
        '''

        oy = {}
        for nid in self.nodes:
            #            tt.nodes[nid] = copy.deepcopy(self.nodes[nid])
            oy[nid] = self.nodes[nid].Y
        Z = 0
        lv = []
        for node, node0 in zip(self.copy.nodes.values(), self.nodes.values()):
            ly = np.log(node0.Y) + node_fitness_fun(node0, **kwargs) * tau
            lv.append(ly)
            node.lY = ly
        logZ = Utils.log_sum(lv)
        for node in self.copy.nodes.values():
            node.Y = np.exp(node.lY - logZ)

        self.copy.set_X_based_on_Y()
        return self.copy

    def relative_entropy(self, other_tree):
        '''
        Computes the relative entropy for the clone size distribution.

        :param other_tree: cfit.tree.Tree
            tree of the same topology, differs by clone sizes only (Y)

        :return: float
        '''

        re = 0
        for nid in self.nodes:
            node1 = self.nodes[nid]
            node2 = other_tree.nodes[nid]
            if node1.Y > 0:
                re += node1.Y * np.log(node1.Y / max(1e-5, node2.Y))
        return re

    def __set_node_X_based_on_Y(self, nid):
        '''

        :param nid: int
            node identifier
        :return:
        '''

        node = self.nodes[nid]
        if len(node.children) == 0:
            node.X = node.Y
        else:
            for cnode in node.children:
                self.__set_node_X_based_on_Y(cnode.id)
            csnodes = [self.nodes[cnode.id] for cnode in node.children]
            node.X = node.Y + sum([cnode.X for cnode in csnodes])

    def set_X_based_on_Y(self):
        '''
        Computes X as a sum of Y of self and the daughter clones.
        '''
        self.__set_node_X_based_on_Y(self.tree.root.id)

    def mutational_load(self):
        '''
        Returns the average mutational load over the nodes on the tree.
        :return: float
        '''
        wl = self.average_over_nodes(Node.TMB)
        return wl

    def neoantigen_load(self):
        '''
        Returns the average neoantigen load over the nodes on the tree.

        :return: float
        '''
        wl = self.average_over_nodes(Node.neoantigen_load)
        return wl

    def neo2mut_ratio(self):
        '''
        Returns the number of neoantigens per mutation, averaged over the nodes on the tree.
        :return: float
        '''
        wl = self.average_over_nodes(Node.neo2mut_ratio)
        return wl

    def syn(self):
        '''
        Returns the number of synonymous mutations, averaged over the nodes on the tree.
        :return: float
        '''
        wl = self.average_over_nodes(Node.syn)
        return wl

    def max_syn(self):
        '''
        Returns the height of the tree in the number of synonymous mutations
        :return: float
        '''
        wl = max([node.syn() for node in self.nodes.values()])
        return wl

    def non_clonal_mutational_load(self):
        '''

        :return: float
        '''
        wl = sum([node.Y * len(node.node.exclusiveMutations) for node in self.nodes.values()])
        return wl

    def non_clonal_neoantigen_load(self, HLAW, wt=False):
        '''

        :param HLAW: cfit.fitness.HLAweights

        :param wt: bool

        :return: float
        '''
        wl = sum([sum([neo.load(HLAW, wt) for neo in node.node.exclusiveNeoantigens]) for node in self.nodes.values()])
        return wl

    '''
    Heterogeneity measures
    '''

    def entropy(self, shared_only=False):
        '''
        Heterogeneity measure, the entropy of the clone size distribution
        :return: float
        '''
        ys = [node.cY if shared_only else node.Y for node in self.nodes.values()]
        ys = [y for y in ys if y > 0]
        entropy = -sum([y * np.log(y) for y in ys])
        return entropy

    def simpson_index(self):
        '''
        Heterogeneity measure
        :return: float
        '''

        ys = [node.Y * node.Y for node in self.nodes.values()]
        ind = sum(ys)
        return ind

    def get_effective_number_of_clones(self):
        '''
        exp(entropy) over the clone size distribution
        :return: float
            effective number of clobes
        '''
        e = self.entropy()
        nc = np.exp(e)
        return nc

    def get_evolved_effective_number_of_clones(self, tau):
        '''
        Returns the effective number of clones in the tumor evolved according to the fitness model.
        :param tau: float
            time parameter for predictions
        :return:  float
        '''
        fil = [Utils.zlog(node.Y) + tau * node.fitness for node in self.nodes.values()]
        lZ = Utils.log_sum(fil)
        y = [np.exp(el - lZ) for el in fil]
        y = [el for el in y if el > 0]
        e = sum([-yy * np.log(yy) for yy in y])
        e = np.exp(e)
        return e

    def get_number_of_clones(self):
        '''

        :return: int
            the number of distinct nodes in the tree
        '''
        return len(list(self.nodes.values())) - 1


    def get_shared_CCF(self, mid1, mid2):
        '''
        Get the fraction of shared frequency of two mutations

        :param mid1: str
            identifier of mutation 1

        :param mid2: str
            identifier of mutation 2

        :return: float
            shared frequency in the tumor
        '''

        shared_ccf = 0
        for node in self.nodes.values():
            mids = set([mut.id for mut in node.mutations])
            if mid1 in mids and mid2 in mids:
                shared_ccf += node.Y

        return shared_ccf

    def get_CCF(self, mid):
        '''
        Get the CCF of mutation mid

        :param mid: str
            identifier of the mutation

        :return: float
            CCF in the tumor
        '''

        ccf = 0
        for node in self.nodes.values():
            mids = set([mut.id for mut in node.mutations])
            if mid in mids:
                ccf += node.Y

        return ccf


    '''
    Patient classification criteria
    '''

    def get_predicted_population_size(self, tau):
        '''
        Predicted n(tau)

        :param tau: float
            The time parameter

        :return: float
            n(tau) value
        '''
        fil = [Utils.zlog(node.Y) + tau * node.fitness for node in self.nodes.values()]
        fi = Utils.log_sum(fil)
        return fi

    def get_population_fitness(self, tau=0):
        '''
        average evolved fitness at time tau.
        :param tau: float
            The time parameter

        :return: float
            average evolved fitness at time tau.
        '''
        nodes = [node for node in self.nodes.values() if node.Y > 0]
        v = [tau * (node.fitness) + np.log(node.Y) for node in nodes]
        lZ = Utils.log_sum(v)
        fil = sum([node.fitness * np.exp(tau * (node.fitness) + np.log(node.Y) - lZ) for node in nodes])
        return fil

    def get_predicted_population_size_derivative(self, tau):
        '''
        Compute d n(tau) / d tau

        :param tau: float
            The time parameter

        :return: float
            d n(tau) / d tau
        '''
        nodes = [node for node in self.nodes.values() if node.Y > 0]
        fil = sum([node.fitness * node.Y * np.exp(tau * (node.fitness)) for node in nodes])
        return fil

    def get_predicted_integrated_population_size(self, tau):
        '''

        :param tau: float

        :return: float
        '''
        nodes = list(self.nodes.values())
        _fun = lambda node: Utils.zlog(node.Y) - np.log(abs(node.fitness)) + np.log(1 - np.exp(tau * node.fitness))
        fil = [_fun(node) for node in nodes]
        fi = Utils.log_sum(fil)
        return fi

    '''
    Output
    '''

    def write_fitness(self, outfile):
        '''
        Writes clone frequencies and fitness of all nodes to the output file.

        :param outfile: str
            path to the output file

        '''
        wl = [[node.id, node.X, node.Y, node.fitness] for node in self.nodes.values()]
        of = open(outfile, 'w')
        wl.sort(key=lambda el: el[0])
        for el in wl:
            of.write("\t".join(list(map(str, el))) + "\n")
        of.close()

    def write_neoantigens(self, opath, sampleName, exclusive=False, nq={}):
        '''

        :param opath: str
            output file

        :param sampleName: str
            name of the sample

        :param exclusive: bool
            whether to output neoantigens only along with the clones where they originate (exclusive==True),
            or to include all nested neoantigens, including these that originated in the ancestral clones.

        :param nq: dict
            dictionary mapping neoantigen identifiers to qualities

        '''

        nodes = list(self.nodes.values())
        nodes.sort(key=lambda node: -node.X)

        datasets = [node.write_neoantigens(sampleName=sampleName, exclusive=exclusive, nq=nq) for node in nodes]
        datasets = [data for data in datasets if data is not None]
        if len(datasets) > 0:
            data = pd.concat(datasets)
            data.to_csv(opath, sep="\t", index=False)

    def write_indel_neoantigens(self, opath, sampleName, exclusive=False):
        '''
        Prints basic indel neoantigen data to a file

        :param opath: str
            output file

        :param sampleName: str
            name of the sample

        :param exclusive: bool
            whether to output neoantigens only along with the clones where they originate (exclusive==True),
            or to include all nested neoantigens, including these that originated in the ancestral clones.

        '''

        nodes = list(self.nodes.values())
        nodes.sort(key=lambda node: -node.X)

        datasets = [node.write_indel_neoantigens(sampleName=sampleName, exclusive=exclusive) for node in nodes]
        datasets = [data for data in datasets if data is not None]
        if len(datasets) > 0:
            data = pd.concat(datasets)
            data.to_csv(opath, sep="\t", index=False)

    def get_mutation_data(self, sampleName, exclusive=True):
        '''
        Returns a DataFrame with data about mutations in the sample,
        grouped by clonality.

        :param sampleName: str
            name of the sample

        :param exclusive: bool
            if True will report each mutation only once, in the clone where it originates

        :return: pd.DataFrame
            data frame with data about the mutations
        '''

        nodes = list(self.nodes.values())
        nodes.sort(key=lambda node: -node.X)
        mdata = None
        for node in nodes:
            nmdata = node.get_mutation_data(sampleName, exclusive=exclusive)
            if nmdata is not None:
                if mdata is None:
                    mdata = nmdata
                else:
                    mdata = pd.concat([mdata, nmdata])
        return mdata

    def write_mutations(self, opath, sampleName, exclusive=True, writeheader=True):
        '''

        :param opath: str
            the output file path

        :param sampleName: str
            sample name
        :param exclusive: bool
            if True will report each mutation only once, in the clone where it originates

        :param writeheader: bool
            whether to include header in the output.

        :return:
        '''
        header = ["Sample", "CloneNumber", "mutation", "gene",
                  "X", "Y", "Taf", "Tcov", "Naf", "Ncov"]
        if writeheader:
            f = open(opath, 'w')
            f.write("\t".join(header) + "\n")
            f.close()
        nodes = list(self.nodes.values())
        nodes.sort(key=lambda node: -node.X)
        muts = set()
        for node in nodes:
            for mut in node.exclusiveMutations:
                muts.add(mut.id)
            node.write_mutations(opath, sampleName, exclusive=exclusive)

    def write_clone_fitness(self, opath):
        '''

        :param opath: str
        :return:
        '''
        of = open(opath, 'w')
        header = ["Clone", "Y", "IM", "DA"]
        of.write("\t".join(header) + "\n")
        for node in self.nodes.values():
            line = "\t".join(list(map(str, [node.id, node.Y, node.fitness, 0]))) + "\n"
            of.write(line)
        of.close()

    def compute_predicted_frequencies(self, tau):
        '''
        Computes the predicted clone frequencies, according to the fitness values set by node.fitness attribute

        :param tau: float
            the time parameter

        :return:
            Modifies the tree in place. Set the node.Y2 attribute to the predicted exclusive frequencies.
            Sets node.X2 values based on node.Y2 values.
        '''

        Z = 0
        for node in self.nodes.values():
            node.Y2 = node.Y * np.exp(node.fitness * tau)
            Z += node.Y2
        for node in self.nodes.values():
            node.Y2 /= Z
        self.root.set_X2()

    def get_mutation_frequencies(self):
        '''
        Returns a dictionary that maps node identifiers to their inclusive frequencies (cancer cellular fraction)

        :return: dict
        '''

        freqs = defaultdict(float)
        for node in self.nodes.values():
            for mut in node.exclusiveMutations:
                mid = mut.id
                freqs[mid] = node.X
        return freqs

    def setF0(self):
        '''
        Sets the average clone fitness

        :return:
            Modifies the tree in place. Sets self.F0 attribute
        '''
        self.F0 = sum([node.Y * node.fitness for node in self.nodes.values()])

    def get_height(self, synonymous=True, xthr=0.0):
        '''
        Estimates the height of the tree, by including the depth of clones of CCF larger the then xthr value.
        These empty clones may be included if the tree reconstruction is done on multiple samples.
        :param synonymous: bool
            use synonymous mutations only

        :param xthr: float
            frequency threshold

        :return: float
            the depth of the tree, counted in the number of nested mutations
        '''

        nodes = [node for node in self.nodes.values() if node.X >= xthr]
        if len(nodes) == 0:
            return 0
        if synonymous:
            s1 = self.nodes[1].node.TMB_syn()
            height = max([node.node.TMB_syn() for node in nodes])
        else:
            s1 = len(self.nodes[1].node.mutations)
            height = max([len(node.mutations) for node in nodes])
        height -= s1
        return height

    '''
    Fitness related methods
    '''

    def set_predicted_frequencies(self, tau, node_fitness_fun=None, **kwargs):
        '''

        :param tau: float

        :param node_fitness_fun: function

        :param kwargs: dict

        :return:
        '''

        values = [[node, np.log(node.Y) + node_fitness_fun(node, **kwargs) * tau] for node in
                  self.nodes.values() if node.Y > 0.0]
        z = Utils.log_sum([x[1] for x in values])
        for node, val in values:
            node.predictedY = np.exp(val - z)

    def ntau(self, node_fitness_fun, tau, shared=False, **kwargs):
        '''
        Compute general n(tau), the predicted cancer population size for a given fitness function
        :param node_fitness_fun: function
            the fitness function, method of Node

        :param tau: float
            the projected time parameter

        :param shared: bool
            use only the clones shared across time points (annotation precomputed)

        :param kwargs: dict
            parameters passed to node_fitness_fun

        :return: float
            n(tau) value
        '''
        values = [np.log(node.sharedY if shared else node.Y) + node_fitness_fun(node, **kwargs) * tau for node in
                  self.nodes.values() if node.Y > 0.0]
        ntau = np.exp(Utils.log_sum(values))
        return ntau

    ######################################
    #   Tree modifications, node removing
    ######################################

    #    def remove_nodes(self, nids):
    #        '''
    #        Removes a given node from the tree
    #        :param nids: list
    #            list of integer identifiers of nodes to be removed
    #        '''
    #
    #        for nid in nids:
    #            self.remove_node(nid)
    #
    #        #correct frequencies
    #        ZX = [node.X for node in self.root.children]
    #        for node in self.nodes.values():
    #            node.X /= ZX
    #
    #        for node in self.nodes.values():
    #            node.Y = node.X - sum([cnode.X for cnode in node.children])

    def get_parent(self, node):
        '''
        Return parent node of node

        :param node: cfit.tree.node.SampleNode
        :return: cfit.tree.node.SampleNode
        '''
        pnode = self.nodes[node.node.parent.id]
        return pnode

    def soft_remove_nodes(self, nids, include_nested=True, overwrite=True):
        '''
        Removes a given node from the tree by setting their frequencies to 0.
        Sets cY and cX frequency attributes of nodes.

        :param nids: list
            list of integer identifiers of nodes to be removed

        :param include_nested: bool
            whether to include frequency of the remove clones in the parent clones

        :param overwrite: bool
            whether to overwrite Y and X attributes of nodes as well.

        :return float
            returns sum of frequencies, should be 1 or 0, if all clades were removed
        '''
        self.processed = True
        innids = defaultdict(lambda: False)
        for nid in nids:
            innids[nid] = True
        for node in self.nodes.values():
#            self.logger("setting node cY:"+str(node.id)+" "+str(node.Y))
#            node.logger("set")
            node.cY = node.Y
            node.cX = node.X

        for nid in nids:
            node = self.nodes[nid]
            if include_nested:
                pnode = self.get_parent(node)
                dropit = False
                while innids[pnode.id] and pnode.id != self.root.id:
                    nid = pnode.id
                    pnode = self.get_parent(pnode)
                    if nid == pnode.id:
                        dropit = True
                        break
                if pnode.id != self.root.id and not dropit:
                    pnode.cY += node.cY
            node.cY = 0.0
            node.cX = 0.0
        Z = sum([node.cY for node in self.nodes.values()])
        # if not include_nested:
        if 0 < Z < 1:
            for node in self.nodes.values():
                node.cY /= Z
            for tnode in self.tree.get_post_order():
                node = self.nodes[tnode.id]
                cnodes = [self.nodes[tcnode.id] for tcnode in tnode.children]
                node.cX = node.cY + sum([cnode.cX for cnode in cnodes])
        Z = sum([node.cY for node in self.nodes.values()])
        if overwrite:
            for node in self.nodes.values():
                node.Y = node.cY
                node.X = node.cX
        return Z

    def remove_new_clones(self):
        '''
        Remove clones that are not present in the other time point
        '''

        nids = [node.id for node in self.nodes.values() if node.privateY > 0]
        self.soft_remove_nodes(nids)

    def set_minimal_self_copy(self):
        '''
        Copies self into attribute self.copy. The only copied attributes of nodes are X, Y and fitness
        '''
        #        self.copy = Tree(params=None, thesample=self.sample)
        self.copy = SampleTree(self.tree)
        for node in self.nodes.values():
            self.copy.nodes[node.id] = SampleNode(node.node)
            self.copy.nodes[node.id].X = node.X
            self.copy.nodes[node.id].Y = node.Y
            self.copy.nodes[node.id].fitness = node.fitness
            self.copy.nodes[node.id].rfitness = node.rfitness

    def score_neoantigens(self, neo_qualities, tau=1, just_dominant=False):
        '''
        Score neoantigens in the nodes as log(node.X) -neo.fitness*tau (just_dominant=False)
            or log(node.X) -neo.fitness*tau if dominant else 0

        :param neo_qualities: dict: str->float
            dictionary mapping neoantigens to their qualities

        :param tau: float

        :param just_dominant: bool

        :return: dict: str -> float
            dictionary with scores for neoantigens
            the lower the score the better (impact on the tumor)

        '''
        dscores = defaultdict(float)
        for node in self.nodes.values():
            dscores1 = node.score_neoantigens(neo_qualities=neo_qualities, tau=tau,
                                              just_dominant=just_dominant)
            for nid in dscores1:
                dscores[nid] += dscores1[nid]
        return dscores


    def score_neoantigens_global(self, neo_qualities, beta=1):
        '''

        Compute recognition probabilities using a model of immunodominance on the whole tumor level

        Q'(n) = (1-Pr^exhaustion(t(n)))*Pr^visible(clone_size(n))*Q(n)

        Immunodominance: only top peptides in the tumor are recognized, but that depends on the tumor size, frequency and t-cell exhaustion
        Pr^reg(n) = exp(beta Q'(n))/\sum_n' exp(beta Q'(n'))

        '''


        dscores = {}
        scores = []
        for node in self.nodes.values():
            scores += [(neo.id, neo_qualities[neo.id], node.X) for neo in self.exclusiveNeoantigens]
        Z = sum([np.exp(el[2]*beta) for el in scores])

        for el in scores:
            nid = el[0]
            q = el[1]
            dscores[nid] = q*np.exp(beta*q/Z)
        for node in self.nodes.values():
            node.set_neoantigen_immunogenicities(dscores)




    def write(self, ofile):
        nodes = []
        for node in self.nodes.values():
            nodes.append([node.id, node.X, node.Y, node.fitness])
        nodes = pd.DataFrame(nodes)
        nodes.columns = ["Node", "X", "Y", "fitness"]
        nodes.to_csv(ofile, sep="\t", index=False)


    ##############
    # write
    ##############

    def toJSON(self):
        jstree = self.tree.toJSON()
        nodes = [jstree["topology"]]
        while nodes:
            node = nodes[0]
            nodes = nodes[1:]
            if "children" in node:
                for child in node["children"]:
                    nodes.append(child)
            clone = self.nodes[node["clone_id"]]
            node["X"] = clone.X
            node["x"] = clone.Y
            if self.processed:
#                node["tilde_X"] = clone.cX
                node["tilde_x"] = clone.cY
#                self.logger("tree processed")
#            else:
#                self.logger("tree not processed")
            node["new_x"] = clone.privateY
#            node["shared_x"] = clone.sharedY
#            node["shared_clone"] = clone.shared_clone
        return jstree
    
    
    def tree_nodes_to_sibyl(self):
        # print('tree, '+str(self.root.id), 'nodes: '+str(len(self.nodes.values())))
        js = list(map(lambda node: node.tree_node_to_sibyl(), self.nodes.values()))
        # self.nodes.values()

        return js
        # jstree = self.tree.toJSON()
#         nodes = [jstree["topology"]]
#         while nodes:
#             node = nodes[0]
#             nodes = nodes[1:]
#             if "children" in node:
#                 for child in node["children"]:
#                     nodes.append(child)
#             clone = self.nodes[node["clone_id"]]
#             node["X"] = clone.X
#             node["x"] = clone.Y
#             if self.processed:
# #                node["tilde_X"] = clone.cX
#                 node["tilde_x"] = clone.cY
# #                self.logger("tree processed")
# #            else:
# #                self.logger("tree not processed")
#             node["new_x"] = clone.privateY
# #            node["shared_x"] = clone.sharedY
# #            node["shared_clone"] = clone.shared_clone
#         return jstree
