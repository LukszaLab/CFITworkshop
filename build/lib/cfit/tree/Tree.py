'''
Created on Mar 24, 2015

@author: mluksza
'''

import copy
from collections import defaultdict

import numpy as np
import pandas as pd
import json

from cfit.tree.GeneralTree import GeneralTree
from cfit.tree.node.Node import Node


class Tree(GeneralTree):
    '''
    Class for clonal structure representation.

    Attributes:
        __variance: float

        __F0: float
            the average fitness over clones of the tree

        __llh: float

        mutation_node_index: dictionary, str->int->bool
            set characteristic function, whether a given mutation (mutation.id) belongs
            to a given node (node.id)

    '''

    def __init__(self, params, thepatient, nodeclass=Node, format='phylowgs'):
        GeneralTree.__init__(self)
        self.__variance = 0.0
        self.__F0 = 0.0
        self.__llh = 0.0
        self.__nodeclass = nodeclass

        if len(params):
                if format == 'phylowgs':
                    self.__init_phylowgs(params, thepatient, nodeclass)
                elif format == "pairtree":
                    self.__init_pairtree(params, thepatient, nodeclass)
        else:
            self.__init_one_node_tree(thepatient=thepatient)

    def __init_one_node_tree(self, thepatient):
        '''
        Initialize a mock one node tree, all mutations assigned to that node

        :param thepatient: cfit.patient.Patient
            patient for whom the tree is reconstructed

        '''
        self.llh = 1.
        self.add_node(0, 0)
        self.add_node(1, 0)
        self.connect_nodes(1, 0)

        mids = list(thepatient.mutations.keys())
        mids.sort()
        for mid in mids:
            mutation = thepatient.mutations[mid]
            self.nodes[1].add_mutation(mutation, check_unique=False)
#        self.logger("Tree of patient "+thepatient.name)
#        self.logger(self.nodes[1].dmutations)
        # add neoantigens (if filled in)
        try:  # if mid is pre thepatient.mutation2neoantigens
            neos = thepatient.mutation2neoantigens[mid]
            for neoantigen in neos:
                self.nodes[1].add_neoantigen(neoantigen, recursive=True, check_unique=False)
        except:
            pass


    def __init_pairtree(self, params=None, thepatient=None, nodeclass=Node):
        '''
        Constructor method

        :param params: list of Pairtree related parameters: 

            resultnpz: str
                .npz file output from Pairtree algorithm

            ord: int
                tree index for multiple-trees reconstructed by Pairtree
            
            sid2mutname: dict
                Dictionary mapping Pairtree mutation identifiers s<num> to global mutation identifier <chrom>_<pos>_<ref>_<alt>/
        
        :param thepatient: cfit.patient.Patient
                patient for which the tree is reconstructed

        :param nodeclass: cfit.tree.node.Node
                a node class for the tree
        '''

        resultnpz = None
        if not params is None:
            [resultnpz, ord, sid2mutname] = params
        if resultnpz is not None:
            self.init_from_pairtree(thepatient,
                                    resultnpz, ord, sid2mutname)
        else:
            self.__init_one_node_tree(thepatient)

        self.copy = None
        self.mutation_node_index = defaultdict(lambda: defaultdict(bool))

    def __init_phylowgs(self, params=None, thepatient=None, nodeclass=Node):
        '''
        Constructor method

        :param params: list of PhyloWGS related parameters: [jstree, jsmutass, jsssms, thepatient, sid2mutname, ord]
        where

            jstree:  str
                json file, <sample_name>.summ.json.gz (or old summ_<sample_name>.json.gz)

            jsmutass: str
                json file <name>.mutass.zip (or mutass_<name>.zip)
                        best tree mutation assignment, dictionary of mutation assignments
            jsssms: str
                json file

            sid2mutname: dict
                Dictionary mapping PhyloWGS mutation identifiers s<num> to global mutation identifier <chrom>_<pos>_<ref>_<alt>/
            ord: int
                sample number for multiple-sample reconstructions
        :param thepatient: cfit.patient.Patient
                patient for which the tree is reconstructed

        :param nodeclass: cfit.tree.node.Node
            a node class for the tree
        '''

        jstree = None
        if not params is None:
            [jstree, jsmutass, jsssms, sid2mutname] = params
        if jstree is not None:
            self.init_from_phylowgs(thepatient,
                                    jstree, jsmutass, jsssms, sid2mutname)
        else: # there is no reconstructed phylogeny - create single node tree
            self.__init_one_node_tree(thepatient)

        self.copy = None
        self.mutation_node_index = defaultdict(lambda: defaultdict(bool))

    @property
    def variance(self):
        return self.__variance

    @variance.setter
    def variance(self, variance):
        self.__variance = variance

    @property
    def F0(self):
        return self.__F0

    @F0.setter
    def F0(self, F0):
        self.__F0 = F0

    @property
    def llh(self):
        return self.__llh

    @llh.setter
    def llh(self, llh):
        self.__llh = llh

    def init_from_phylowgs(self, thepatient,
                           jstree, jsmutass, jsssms, sid2mutname=None):
        '''
        Initializes trees from PhyloWGS files

        :param thepatient: cfit.patient.Patient
            Patient for which the tree is reconstructed

        :param jstree:  str
            tree dictionary read from the json file with clone frequencies,
            <sample_name>.summ.json.gz  (or old summ_<sample_name>.json.gz)

        :param jsmutass: str
            zip file with the family of trees, each file contains mutation assignment to clones.
            best tree mutation assignment, dictionary of mutation assignments
            json file <sample_name>.mutass.zip (or old mutass_<name>.zip)

        :param jsssms: str
            mutation data: identifiers, reads etc
            json file <sample_name>.muts.json.gz (or old muts_<sample_name>.json.gz)

        :param sid2mutname: dict
            Dictionary mapping PhyloWGS mutation identifiers s<num> to global mutation identifier <chrom>_<pos>_<ref>_<alt>
            optionoal - the mapping should me provided in the muts_*.json.gz file

        '''

        # structure and frequencies
        structure = jstree['structure']

        # loglikelihood
        llh = float(jstree["llh"])
        if np.isnan(llh):
            llh = -1e10
        if np.isneginf(llh):
            llh = -1e10
        self.llh = llh

        # set the root clone 0
        self.add_node(0, 0)

        # add other clones to the tree
        for parent in structure:
            children = structure[parent]
            for child in children:
                self.add_node(child, parent)

        # connect clones to their parent clones
        for parent in structure:
            children = structure[parent]
            for child in children:
                self.connect_nodes(int(child), int(parent))

        # mapping mutations from phylowgs identifiers
        for key in jsmutass: #key: "0", "1",
            cid = int(key) #cid: 0, 1
            # mutation info
            sids = jsmutass[key]['ssms'] # mutations from ssm file
            cnvids = jsmutass[key]['cnvs'] # cnvs from cnv file
            for sid in sids:
                # sid: s<num>
                # sid = str(sid)
                try:
                    mid = jsssms[sid]['name']
                except:
                    mid = sid2mutname[sid]
                mid = "_".join(mid.split("_")[:4])  # <chrom>_<pos>_<ref>_<alt>
                # ## test to use mutation ids without sample name
                (mid, checked) = thepatient.check_mutation(mid)
                if not checked:
                    self.logger(str(mid) + " is missing in " + thepatient.name)
                    continue
                mutation = thepatient.mutations[mid]
                self.nodes[cid].add_mutation(mutation, check_unique=False)
                # add neoantigens (if filled in)
                try:  # if mid is pre thepatient.mutation2neoantigens
                    neos = thepatient.mutation2neoantigens[mid]
                    for neoantigen in neos:
                        self.nodes[cid].add_neoantigen(neoantigen, recursive=True, check_unique=False)
                except:
                    pass
            for cnvid in cnvids:
                self.nodes[cid].add_cnv(cnvid)

        for node in self.nodes.values():
            node.mutations = list(set(node.mutations))
            node.neoantigens = list(set(node.neoantigens))
        self.set_exclusive_mutations()
        self.set_exclusive_cnvs()


    def init_from_pairtree(self, thepatient, 
                           resultnpz, ord, sid2mutname):
        '''
        Initializes trees from Pairtree files

        :param thepatient: cfit.patient.Patient
            Patient for which the tree is reconstructed

        :param resultnpz: str
            file path to .npz file outputted by pairtree, with all neseccary clone and tree info
           
        :param ord: int
            tree index of desired tree in .npz results file

        :param sid2mutname: dict
            Dictionary mapping Pairtree mutation identifiers s<num> to global mutation identifier <chrom>_<pos>_<ref>_<alt>
        '''

        data = np.load(resultnpz)
        
        #if ord >= len(data['struct']):
        #return 
        
        self.add_node(0,0)
        
        parent_vector = data['struct'][ord]
        n_clones = len(parent_vector)

        self.llh = float(data['llh'][ord])

        for i in range(n_clones):
            cid = i+1
            pid = int(parent_vector[i])
            self.add_node(cid, pid)
        
        for i in range(n_clones):
            cid = i+1
            pid = int(parent_vector[i])
            self.connect_nodes(cid, pid)


        # read cluster info from .npx file
        #    to create cluster2muts dict, where key = cluster number and value = mutations exclusive to cluster
        clusters = json.loads(data["clusters.json"])
        cluster2muts = dict()
        for i in range(1, n_clones + 1):
            cluster2muts[i] = clusters[i-1]

        for nid in self.nodes.keys(): # assuming self.nodes is list of clones 
            if nid == 0:
                continue
            for sid in cluster2muts[nid]:
                self.logger("mapping mutations")
                mid = sid2mutname[sid]
                self.logger("sid: "+str(sid)+" mid: "+str(mid))
                mutation = thepatient.mutations[mid]
                self.nodes[nid].add_mutation(mutation, check_unique=False)
                try:
                    neos = thepatient.mutation2antigens[mid]
                    for neoantigen in neos:
                        self.nodes[nid].add_neoantigen(neoantigen, recursive=True, check_unique=False)
                except:
                    pass
        
        for node in self.nodes.values():
            node.mutations = list(set(node.mutations))
            node.neoantigens = list(set(node.neoantigens))

        self.set_exclusive_mutations()
        self.set_exclusive_cnvs()

    def one_node_tree(self, thepatient):
        '''
        Creates a homogenous (single clone) tree version of the tree, with all mutations in the single clone
        :return: cfit.tree.Tree

        '''
        otree = Tree([], thepatient)
        onode = otree.nodes[1]
        mutations = set()
        neoantigens = set()
        for node in self.nodes.values():
            mutations = mutations.union(node.mutations)
            neoantigens = neoantigens.union(node.neoantigens)
        onode.mutations = list(mutations)
        onode.neoantigens = list(neoantigens)
        otree.set_exclusive_mutations()
        return otree

    def set_self_copy(self):
        '''
        Copies self into attribute self.copy
        '''

        if self.copy is not None:
            self.copy = None
        tt = copy.deepcopy(self)
        self.copy = tt

    def set_exclusive_mutations(self):
        '''
        Assigns mutations to nodes where they originate.
        '''
        for node in self.nodes.values():
            node.set_exclusive_mutations()

    def set_exclusive_cnvs(self):
        '''
        Assigns cnvs to nodes where they originate.
        '''
        for node in self.nodes.values():
            node.set_exclusive_cnvs()

    def distribute_neoantigens_to_clones(self, thepatient):
        '''
        Assigns neoantigens to respective clones on the tree

        :param thepatient: Patient
            The patient for which the tree is reconstructed

        '''
        for node in self.nodes.values():
            for mut in node.mutations:
                mid = mut.id
                if mid in thepatient.mutation2neoantigens:
                    neos = thepatient.mutation2neoantigens[mid]
                    for neoantigen in neos:
                        node.add_neoantigen(neoantigen, recursive=False)
                if mid in thepatient.mutation2fsneoantigens:
                    neos = thepatient.mutation2fsneoantigens[mid]
                    for neoantigen in neos:
                        node.add_frame_shift_neoantigen(neoantigen, recursive=False)

    #    def one_node_tree(self, thepatient):
    #        '''
    #        Creates a tree with a singleton clone.
    #
    #        :param thepatient: cfit.patient.Patient
    #            The patient for which the tree is reconstructed.
    #
    #        '''
    #        GeneralTree.__init__(self)
    #        self.llh = 1.
    #        self.add_node(0, 0)
    #        self.add_node(1, 0)
    #        self.connect_nodes(1, 0)
    #        if False:
    #            for node in self.nodes.values():
    ##                for mutation in node.exclusiveMutations:
    #                    self.nodes[1].add_mutation(mutation, check_unique=False)
    #                for neoantigen in node.exclusiveNeoantigens:
    #                    self.nodes[1].add_neoantigen(neoantigen, check_unique=False, recursive=True)
    #
    #    #        for mid in thesample.mutations:
    #    #            mutation = thesample.mutations[mid]
    #    #            self.nodes[1].add_mutation(mutation, check_unique=False)
    #    #            if mid in thesample.mutation2neoantigens:
    #    #                neos = thesample.mutation2neoantigens[mid]
    #    #                for neoantigen in neos:
    #    #                    self.nodes[1].add_neoantigen(neoantigen, check_unique=False, recursive=True)
    #            self.nodes[0].set_exclusive_mutations()
    #            self.nodes[1].set_exclusive_mutations()

    def set_root(self, rnode):
        '''
        Sets the root on the tree to rnode

        :param baseNode: cfit.tree.node.Node
            the node object that is set to be the root of the tree.

        '''
        for cnode in rnode.children:
            cnode.parent = rnode
        self.nodes[rnode.id] = rnode
        self.root = rnode
        self.root.parent = self.root

    def add_node(self, nid, pnid):
        '''
        Add new Node to the tree with nid identifier, as a child to node pid.

        :param nid: int
            New node identifier

        :param pnid: int
            parent node identifier

        '''
        if (nid not in self.nodes.keys()):
            node = self.__nodeclass(nid)
            self.nodes[nid] = node
        node = self.nodes[nid]
        if nid == pnid or nid == 0:
            self.set_root(node)
        self.nodes[nid] = node

    def get_mutations(self):
        '''
        Return all mutations on the tree.
        :return: set
            The set of mutations stored in all nodes of the tree
        '''
        mutations = set()
        for node in self.nodes.values():
            mutations = mutations.union(node.exclusiveMutations)
        return mutations

    def get_neoantigens(self):
        '''
        Returns all neoantigens on the tree.

        :return: set
            The set of neoantigens stored in all nodes of the tree
        '''
        neoantigens = set()
        for node in self.nodes.values():
            neoantigens = neoantigens.union(node.exclusiveNeoantigens)
        return neoantigens

    def set_mutation_node_index(self):
        '''
        Sets the mutation-node index, for faster access to node mutation content
        '''
        self.mutation_node_index = defaultdict(lambda: defaultdict(bool))
        for node in self.nodes.values():
            for mut in node.mutations:
                self.mutation_node_index[mut.id][node.id] = True

    def toVCF(self):
        '''
        Returns a dictionary of vcf-like data frames for each node (clone) in the tumor. The dictionary
        maps clone id to the DataFrame

        :return: dict
            dictionary mapping node identifiers to vcf formatted list of attributes
        '''
        js = {}
        for nid in self.nodes:
            node = self.nodes[nid]
            vcf = node.toVCF()
            js[nid] = vcf
        return js

    def get_stats(self):
        '''
        Reports the basic statistics: of clones - number of mutations, new mutations, neoantigens and CCF in a data frame.

        :return: pd.DataFrame
        '''

        stats = []
        for nid in self.nodes:
            node = self.nodes[nid]
            stats.append(
                [nid, node.TMB(), len(node.get_exclusive_mutations()), node.neoantigen_load(), node.X])
        stats = pd.DataFrame(stats)
        stats.columns = ["Node", "SSMs", "New_SSMs", "Neoantigens", "CCF"]
        return stats

    def get_mutations_by_node(self, exclusive=False):
        '''

        :param exclusive: bool
            whether to report the exclusive or inclusive list of mutations.

        :return: dict (collection.defaultdict(list))
            dictionary mapping node identifiers to the list of mutation identifiers of
            the form <chr>_<pos>_<ref>_<alt>
        '''

        lmuts = defaultdict(list)
        for nid in self.nodes:
            node = self.nodes[nid]
            if exclusive:
                mids = [mut.id for mut in node.get_exclusive_mutations()]
            else:
                mids = [mut.id for mut in node.mutations]
            lmuts[nid] = mids
        return lmuts

    def toJSON(self):
        '''
        Creates json for Sibyl app.

        :return: dict
            dictionary that represents the tree and can be written in a json file.

        '''
        js = {}
        jstree = self.root.toJSON()
        js['topology'] = jstree
        js['score'] = self.llh
        return js

    def get_tree_nodes(self):
        # print('tree, '+str(self.root.id), 'nodes: '+str(len(self.nodes.values())))
        js = list(map(lambda node: node.get_tree_node(), self.nodes.values()))
        # self.nodes.values()
    
        return js


    def to_Cardelino(self):
        '''
        Creates a table with cardelino output.

        return: pandas.DataFrame
        '''

        mids = set()
        dmutations = defaultdict(lambda: defaultdict(int))
        for nid in self.nodes:
            node = self.nodes[nid]
            for mutation in node.mutations:
                dmutations[mutation.id][nid] = 1
                mids.add(mutation.id)

        tab = pd.DataFrame([[dmutations[mid][nid] for nid in self.nodes] for mid in mids], columns=["clone"+str(i) for i in self.nodes])
        tab.index = list(mids)
        return tab


