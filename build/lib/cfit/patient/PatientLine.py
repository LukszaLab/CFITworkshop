'''
Created on Oct 17, 2018

@author: mluksza
'''

import gzip
import json
import math
import os
import zipfile
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.optimize import LinearConstraint
from scipy.optimize import basinhopping
from scipy.stats import spearmanr
from scipy.stats import wasserstein_distance

from cfit.patient.Patient import Patient
from cfit.patient.Sample import Sample
from cfit.patient.TimePoint import TimePoint
from cfit.tree.SampleTree import SampleTree
from cfit.tree.Tree import Tree
from cfit.tree.mutation.FrameShiftNeoantigen import FrameShiftNeoantigen
from cfit.tree.mutation.Neoantigen import Neoantigen
from cfit.tree.node.Node import Node
from cfit.util.Utils import Utils


class PatientLine(Patient):
    '''
    Class implementing longitudinally organized samples in a tumor from a patient.

    Attributes:

        __pname: str
            name

        __orderedTimePoints: list
            list of str, the names of time points, ordered by time

        timePoints: dict #defaultdict(lambda: cfit.patient.TimePoint)
            dictionary mapping time point names (str) to time points objects

        samples : dict: str->cfit.patient.Sample
            maps sampleid to sample object

        time_pairs: list
            to be reviewed



    '''

    def __init__(self, params):
        '''
        Constructor method
        '''
        Patient.__init__(self, params)
        self.__pname = self.name
        # time ordered samples
        self.__orderedTimePoints = []
        self.timePoints = {}  # defaultdict(lambda: None)  # time point name -> TimePoint
        self.sampleMap = {}
        self.time_pairs = []
        # maps time points to TimePoints

    @property
    def pname(self):
        return self.__pname

    @pname.setter
    def pname(self, pname):
        self.__pname = pname

    @property
    def orderedTimePoints(self):
        return self.__orderedTimePoints

    @property
    def all_samples(self):
        samples = []
        for tp in self.timePoints:
            samples += list(self.timePoints[tp].samples.values())
        return samples

    @property
    def all_trees(self):
        trees = []
        for sample in self.all_samples:
            for tree in sample.trees:
                trees.append(tree)
        return trees

    @property
    def all_nodes(self):
        nodes = []
        for sample in self.all_samples:
            for tree in sample.trees:
                for node in tree.nodes.values():
                    nodes.append(node)
        return nodes

    def get_time_pairs_by_prefix(self, pref1, pref2):
        '''
        Organizes the list of longitudinal pairs of time points, sets the time_pairs attribute

        :param pref1: str
            example "Pre", "Prim"

        :param pref2: str
            example "Met", "Post"

        :return: list
            of [str, str]
        '''

        tpoints1 = [tp for tp in self.timePoints if tp.startswith(pref1)]
        tpoints2 = [tp for tp in self.timePoints if tp.startswith(pref2)]
        time_pairs = [[tp1, tp2] for tp1 in tpoints1 for tp2 in tpoints2]
        return time_pairs

    def time_point_similarities(self, num):
        '''
        Compute relative entropies between all time point pairs.

        :param num: int
            which tree index to use. S
            TBD: the method should be rewritten to do proper averaging over all trees.
        :return: dict
            time point name 1 (str) -> time point name 2 (str) -> float
        '''
        ltrees = {}
        for tp in self.timePoints:
            tpoint = self.timePoints[tp]
            ltrees[tp] = tpoint.trees(num)
        sims = {}
        eps = 1e-6
        for tp1 in ltrees:
            trees1 = ltrees[tp1]
            sims[tp1] = {}
            for tp2 in ltrees:
                sims[tp1][tp2] = []
                trees2 = ltrees[tp2]
                for tree1 in trees1:
                    for tree2 in trees2:
                        kl1 = 0
                        kl2 = 0
                        for nid in tree1.nodes:
                            node1 = tree1.nodes[nid]
                            node2 = tree2.nodes[nid]
                            if node1.Y > 0:
                                kl1 += node1.Y * np.log(node1.Y / max(node2.Y, eps))
                            if node2.Y > 0:
                                kl2 += node2.Y * np.log(node2.Y / max(node1.Y, eps))
                        sims[tp1][tp2].append([kl1, kl2])
        return sims

    def analyze_sample_purity(self, mapping_params, sample_order, treedir, odir):
        '''

        :param mapping_params: dict
            json dictionary describing the structure of samples, time points and patients.

        :param sample_orderr: list
            list of sample names (str), in the order they are reffered to in phylowgs
            eg. ["S1", "S2"]

        :param treedir: str
            path to the tree directory

        :param odir: str
        :return:
        '''
        for sample_data in mapping_params['samples']:
            self.logger(sample_data, 5)
            [patientname, samplename, snum, vcffile, timepoint, tissue, normal_reads, tumor_reads] = sample_data
            fname = vcffile.split("_full")[0]
            sample = Sample(snum)
            sample.tumor_reads = tumor_reads
            sample.normal_reads = normal_reads
            sample.tissue = tissue
            jsonpath = os.path.join(treedir, self.name + ".summ.json.gz")
            self.logger("Tree file :" + jsonpath, 5)
            if not os.path.exists(jsonpath):
                continue
            ord = sample_order.index(snum)
            purity_stats = sample.get_tumor_purity_scores(jsonpath, ord=ord)
            opath = os.path.join(odir, fname + ".txt")  # self.name+"_"+snum)
            purity_stats = purity_stats.sort_values('llh', ascending=False)
            purity_stats.to_csv(opath, sep="\t", index=False)

    def create_samples(self, mapping_params, sample_order, vcfdir, treedir, ntrees,
                       tree_format='phylowgs'):
        '''

        Creates all data structures to represent samples, mutations and their phylogenies.


        :param mapping_params: dict
            json dictionary describing patient sample structure, time points etc.

        :param sample_order: list
            ordered list of sample numbers, e.g. ["S1", "S2",...]

        :param vcfdir: str
            path to the directory with vcf files

        :param treedir: str
            path to the directory with tree files of that patient

        :param ntrees: int
            the number of top scoring trees to be used.

        :param tree_format: str
            currently 'phylowgs' or 'pairtree'

        '''

        # 1. create samples, without trees
        for sample_data in mapping_params['samples']:
            [patientname, samplename, snum, vcffile, timepoint, tissue, normal_reads, tumor_reads] = sample_data
            if "_full_pan_sample" in vcffile:
                sampleID = vcffile.split("_full_pan_sample")[0]
            else:
                sampleID = samplename
            self.logger("Creating sample " + patientname + " " + samplename)
            sample = Sample(samplename)
            sample.tumor_reads = tumor_reads
            sample.normal_reads = normal_reads
            sample.tissue = tissue
            sample.sampleID = sampleID
            self.add_sample(sample, timepoint)
            self.mutations = sample.import_mutations_from_VCF_file(os.path.join(vcfdir, vcffile), self.mutations)

        # 2. Read trees and clone frequencies
#        self.logger("patient mutations:")
#        self.logger(self.mutations)
        if tree_format == 'phylowgs':
            self.import_phylowgs_sample_tree(mapping_params, sample_order, treedir, ntrees)

        elif tree_format == 'pairtree':
            self.import_pairtree_sample_tree(mapping_params, sample_order, treedir, ntrees)

    def import_phylowgs_sample_tree(self, mapping_params, sample_order, treedir, ntrees):
        '''
        Import tree and frequencies from PhyloWGS output files

        :param mapping_params: dict
            json dictionary describing patient sample structure, time points etc.

        :param sample_order: list
            ordered list of sample numbers, e.g. ["S1", "S2",...]

        :param treedir: str
            path to the directory with tree PhyloWGS files

        :param ntrees: int
            the number of top scoring trees to be used.


        '''

        # import trees
        jsonpath = os.path.join(treedir, self.name + ".summ.json.gz")
        self.logger("jsonpath = " + jsonpath + " " + str(os.path.exists(jsonpath)))
        if not os.path.exists(jsonpath):
            jsonpath = os.path.join(treedir, "summ_" + self.name + ".json.gz")
        js, bestnums = self.__read_phylowgs_clonal_structures(jsonpath, ntrees)
        self.logger("read clonal structures for patient "+self.name+" "+str(len(self.trees)))
        self.oneTree = self.trees[0].one_node_tree(self)
        # Fill in sample trees
        for sample_data in mapping_params['samples']:
            [patientname, samplename, snum, vcffile, timepoint, tissue, normal_reads, tumor_reads] = sample_data
            if "_full_pan_sample" in vcffile:
                sampleID = vcffile.split("_full_pan_sample")[0]
            else:
                sampleID = samplename
            sample = self.sampleMap[sampleID]
            sample.oneTree = SampleTree(self.trees[0])

            ord = sample_order.index(snum)
            if js is not None:
                for i, bestnum in enumerate(bestnums):
                    jstree = js['trees'][bestnum]
                    sampletree = SampleTree(self.trees[i], param=jstree, ord=ord)
                    sample.trees.append(sampletree)
            else:
                for i in range(ntrees):
                    sampletree = SampleTree(self.oneTree)
                    sample.trees.append(sampletree)


    def import_pairtree_sample_tree(self, mapping_params, sample_order, treedir, ntrees):
        '''

        Import tree and frequencies from Pairtree output files

        :param mapping_params: dict
            json dictionary describing patient sample structure, time points etc.

        :param sample_order: list
            ordered list of sample numbers, e.g. ["S1", "S2",...]

        :param treedir: str
            path to the directory with tree Pairtree files for patient 

        :param ntrees: int
            the number of top scoring trees to be used.

        '''
        npzpath = os.path.join(treedir, self.name+"_results.npz")
        self.logger("npzpath = " + npzpath + " " + str(os.path.exists(npzpath)))

        ssmpath = os.path.join(treedir, self.name+".ssm")
        self.logger("ssmpath = " + ssmpath + " " + str(os.path.exists(ssmpath)))
        
        sid2mutname = defaultdict(lambda: "")
        mdata = pd.read_csv(ssmpath, sep='\t')
        
        for (sid, mid) in zip(mdata.id, mdata.name):
            sid2mutname[sid] = mid

        num_trees = len(np.load(npzpath)['struct'])
        
        self.logger("Importing trees from " + npzpath)
        self.logger("Importing top " + str(ntrees) + " of " + str(num_trees) + " trees...")

        for i in range(min(ntrees, num_trees)):
            params = [npzpath, i, sid2mutname]
            #self.logger(str(i))
            #self.logger(params)
            t = Tree(params, thepatient=self, nodeclass=Node, format='pairtree')
            self.trees.append(t)


        # This line creates a single node tree with all mutations, needed for
        # some model comparisons when we disregard tumor heterogeneity
        self.oneTree = self.trees[0].one_node_tree(self)

        # 2. Write code that creates corresponding SampleTree objects for each sample of the patient

        for sample_data in mapping_params['samples']: #iterates over sample names
            [patientname, samplename, snum, vcffile, timepoint, tissue, normal_reads, tumor_reads] = sample_data
            if "_full_pan_sample" in vcffile:
                sampleID = vcffile.split("_full_pan_sample")[0]
            else:
                sampleID = samplename
            sample = self.sampleMap[sampleID]
            sample.oneTree = SampleTree(self.trees[0])
            sample.timePoint = timepoint

            ord = sample_order.index(snum)
            for i in range(min(ntrees, num_trees)):
                param = [i, npzpath]
                samptree = SampleTree(self.trees[i], param=param, ord=ord, format="pairtree")
                sample.trees.append(samptree)

    def add_neoantigens(self, neofile, kd_thr=2000, remove_MTandX=False, ns=None):
        '''

        :param neofile: str
            path to neoantigen file, neodir/neoantigens_<patient_name>.txt. The format of the file is:

            ID	MUTATION_ID	Sample	MT_Peptide	WT_Peptide	Allele	MT_Score	WT_Score
            1_1148420_G_A_9	1_1148420_G_A	5-LTS	TATQDTVCC	TATQDTVCR	HLA-A02:01	27395.0	28728.0

        :param kd_thr: float or None
            threshold on Kd of neoantigens

        :param remove_MTandX: bool

        :param ns: list
            list of acceptable peptide lengths

        '''
        kd_thr2 = 1e8 if kd_thr is None else kd_thr
        self.logger("Adding neoantigens from file " + neofile)
        neos = pd.read_csv(neofile, sep="\t")
        if ns is not None:
            nsset = set(ns)
            neos = neos[[len(pep) in nsset for pep in neos.MT_Peptide]]

        neoid = {}
        for line in neos.itertuples():
            if line.MT_Score > kd_thr2:
                continue
            if remove_MTandX:
                if line.PEPTIDE_ID.split("_")[0].isalpha():  # remove MT and X
                    continue
            allele = line.Allele
            allele = allele.replace("HLA-", "").replace(":", "")
            try:
                nid = line.PEPTIDE_ID + "_" + allele
                peptide_id = line.PEPTIDE_ID
            except AttributeError:
                nid = line.ID + "_" + allele
                peptide_id = line.ID
            params = [nid, peptide_id, line.MUTATION_ID, None, line.WT_Peptide, line.MT_Peptide, line.Allele,
                      line.WT_Score, line.MT_Score, None, 1]
            params[0] = nid
            pep_allele_id = peptide_id + "_" + line.Allele
            if pep_allele_id in neoid:
                neo = neoid[pep_allele_id]
            else:
                neo = Neoantigen(params)
                neoid[pep_allele_id] = neo
            for tp in self.timePoints:
                tpoint = self.timePoints[tp]
                tpoint.add_neoantigen(neo)
            self.neoantigens[neo.id] = neo
            self.mutation2neoantigens[neo.mid].append(neo)
        self.logger("done.")

    def add_frame_shift_neoantigens(self, neofile, kd_thr=2000, remove_MTandX=False, ns=None):
        '''
        Add frameshift neoantigens

        :param neofile: str
            path to neoantigen file, neodir/neoantigens_<patient_name>.txt. The format of the file is:

            ID	MUTATION_ID	Sample	MT_Peptide	WT_Peptide	Allele	MT_Score	WT_Score
            1_1148420_G_A_9	1_1148420_G_A	5-LTS	TATQDTVCC	-	HLA-A02:01	27395.0	-1

        :param kd_thr: float
            threshold on Kd of neoantigens

        :param ns: list
            list of acceptable peptide lengths
        '''
        self.logger("Adding fs neoantigens from file " + neofile)
        kd_thr2 = 1e8 if kd_thr is None else kd_thr
        neos = pd.read_csv(neofile, sep="\t")
        if ns is not None:
            nsset = set(ns)
            neos = neos[[len(pep) in nsset for pep in neos.MT_Peptide]]

        neoid = {}
        for line in neos.itertuples():
            if line.MT_Score > kd_thr2:
                continue
            if remove_MTandX:
                if line.PEPTIDE_ID.split("_")[0].isalpha():  # remove MT and X
                    continue
            allele = line.Allele
            allele = allele.replace("HLA-", "").replace(":", "")
            try:
                nid = line.PEPTIDE_ID + "_" + allele
                peptide_id = line.PEPTIDE_ID
            except AttributeError:
                nid = line.ID + "_" + allele
                peptide_id = line.ID
            params = [nid, peptide_id, line.MUTATION_ID, None, line.WT_Peptide, line.MT_Peptide, line.Allele,
                      line.WT_Score, line.MT_Score, None, 1]
            params[0] = nid
            pep_allele_id = peptide_id + "_" + line.Allele
            if pep_allele_id in neoid:
                neo = neoid[pep_allele_id]
            else:
                neo = FrameShiftNeoantigen(params)
                neoid[pep_allele_id] = neo
            for tp in self.timePoints:
                tpoint = self.timePoints[tp]
                tpoint.add_frame_shift_neoantigen(neo)
            self.fsneoantigens[neo.id] = neo
            self.mutation2fsneoantigens[neo.mid].append(neo)
        self.logger("done.")

    def set_exclusive_mutations(self):
        '''
        Sets exclusive mutations in clones in the trees in samples.
        '''
        self.logger(self.name + ": setting exclusive mutations.")
        for tree in self.trees:
            tree.set_exclusive_mutations()
        self.oneTree.set_exclusive_mutations()
        self.logger("done.")

    def distribute_neoantigens_to_clones(self):
        '''
        Distributes neoantigens in over clones of trees.
        '''
        self.logger(self.name + ": distributing neonatigens over clones")
        for tree in self.trees:
            tree.distribute_neoantigens_to_clones(self)
        self.oneTree.distribute_neoantigens_to_clones(self)
        self.logger("done.")

    def add_sample(self, sample, timePoint="T1"):
        '''
        Add sample to the designated time point.

        :param sample: cfit.patient.Sample
        :param timePoint: str
        '''
        self.sampleMap[sample.sampleID] = sample
        if timePoint not in self.__orderedTimePoints:
            self.__orderedTimePoints.append(timePoint)
        Patient.add_sample(self, sample)
        if timePoint not in self.timePoints:
            self.timePoints[timePoint] = TimePoint(timePoint)
        self.timePoints[timePoint].add_sample(sample)

    def get_samples(self, simple=False):
        '''
        Returns a list of samples with the name of time point.

        :return: list
        '''

        lsamples = []
        for tp in self.__orderedTimePoints:
            if simple:
                for sample in self.timePoints[tp].samples.values():
                    lsamples.append(sample)
            else:
                lsamples.append([tp, list(self.timePoints[tp].samples.values())])
        return lsamples

    def get_mutation(self, mid):
        '''

        :param mid: str
            mutation identified

        :return: Mutation
        '''
        sample0 = self.all_samples[0]
        return sample0.mutations[mid]

    def get_mutation_frequencies(self, beta=1.0, exclusive=False, by_sample=False,
                                 nonsynonymous_only=False, kd_threshhold=None):
        '''
        Get mutation frequencies in the patient, across the time points

        :param beta: float
            tree weighting parameter

        :param exclusive: bool
            whether to report Y or X.

        :param by_sample: bool
            whether to report frequencies for each sample, or averaged by time point

        :param nonsynonymous_only: bool
            whether to include only nonsynonymous mutations

        :param kd_threshhold: float
            whether to include only mutations with neoantigens of kd < kd_threshold
        :return: pd.DataFrame
        '''
        dCCFs = defaultdict(lambda: defaultdict(float))
        mids = set()
        if by_sample:
            cnames = [sample.sampleID for sample in self.all_samples]
            for sample in self.all_samples:
                mut2CCF = sample.get_mutation_frequencies(beta=beta, exclusive=exclusive,
                                                          nonsynonymous_only=nonsynonymous_only,
                                                          kd_threshhold=kd_threshhold)
                for mid in mut2CCF:
                    dCCFs[mid][sample.sampleID] = mut2CCF[mid]
                    mids.add(mid)
        else:
            cnames = self.__orderedTimePoints
            for tp in self.timePoints:
                tpoint = self.timePoints[tp]
                mut2CCF = tpoint.get_mutation_frequencies(beta=beta, exclusive=exclusive,
                                                          nonsynonymous_only=nonsynonymous_only,
                                                          kd_threshhold=kd_threshhold)
                for mid in mut2CCF:
                    dCCFs[mid][tp] = mut2CCF[mid]
                    mids.add(mid)

        mids = list(mids)
        rank = defaultdict(float)
        for mid in mids:
            for cname in cnames:
                rank[mid] += dCCFs[mid][cname]
        mids.sort(key=lambda mid: -rank[mid])

        tab = [[mid] + [dCCFs[mid][colname] for colname in cnames] for mid in mids]
        tab = pd.DataFrame(tab)
        tab.columns = ["Mutation_ID"] + cnames
        return tab

    def clone_trajectories(self, treenum=0, inclusive=True):
        '''
        Reports clone frequency trajectories for the given tree over the time points.

        :param treenum: int
            index (ranking) of the tree to use

        :param inclusive: bool
            whether to use inclusive or exclusive trajectories

        :return: pd.DataFrame
            format - columns are time point names, rows are clones.
        '''

        freqs = defaultdict(lambda: defaultdict(float))
        timepoints = list(self.timePoints)
        nids = set()
        for timepoint in timepoints:
            sample = self.timePoints[timepoint]
            tree = sample.trees[treenum]
            for nid in tree.nodes:
                node = tree.nodes[nid]
                if inclusive:
                    freqs[timepoint][nid] = node.X
                else:
                    freqs[timepoint][nid] = node.Y
                nids.add(nid)

        nids = list(nids)
        nids.sort()
        freqtab = [[nid] + [freqs[tp][nid] for tp in timepoints] for nid in nids]
        freqtab = pd.DataFrame(freqtab)
        freqtab.columns = ["Clone"] + timepoints
        return freqtab

    def clone_neoantigen_load_trajectories(self, treenum=0):
        '''

        Reports effective neoantigen load trajectories for the given tree over the time points.

        :param treenum: int
            index (ranking) of the tree to use

        :return: pd.DataFrame
            format - columns are time point names, rows are clones.
        '''
        freqs = defaultdict(lambda: defaultdict(float))
        timepoints = list(self.timePoints)
        nids = set()
        for timepoint in timepoints:
            sample = self.timePoints[timepoint]
            tree = sample.trees[treenum]
            for nid in tree.nodes:
                node = tree.nodes[nid]
                freqs[timepoint][nid] = node.neoantigen_load()
                nids.add(nid)
        nids = list(nids)
        nids.sort()
        freqtab = [[nid] + [freqs[tp][nid] for tp in timepoints] for nid in nids]
        freqtab = pd.DataFrame(freqtab)
        freqtab.columns = ["Clone"] + timepoints
        return freqtab

    def clone_fitness_trajectories(self, treenum=0):
        '''

        Reports the average fitness trajectories for the given tree over the time points.

        :param treenum: int
            index (rank) of the tree to use

        :return: pd.DataFrame
            format - columns are time point names, rows are clones.
        '''

        ffs = defaultdict(lambda: defaultdict(float))
        rffs = defaultdict(lambda: defaultdict(float))
        timepoints = list(self.timePoints)
        nids = set()
        for timepoint in timepoints:
            sample = self.timePoints[timepoint]
            tree = sample.trees[treenum]
            for nid in tree.nodes:
                node = tree.nodes[nid]
                ffs[timepoint][nid] = node.fitness
                rffs[timepoint][nid] = node.rfitness
                nids.add(nid)
        nids = list(nids)
        nids.sort()

        freqtab = [[nid] + [ffs[tp][nid] for tp in timepoints] for nid in nids]
        freqtab = pd.DataFrame(freqtab)
        freqtab.columns = ["Clone"] + timepoints

        rfreqtab = [[nid] + [rffs[tp][nid] for tp in timepoints] for nid in nids]
        rfreqtab = pd.DataFrame(rfreqtab)
        rfreqtab.columns = ["Clone"] + timepoints

        return freqtab, rfreqtab

    ##############################################
    # Predictions statistics
    ##############################################

    def compute_observed_fitness_flux(self, treenum=0, timepoints=None, randomize=False):
        '''

        :param treenum: int
            index (rank) of the tree to use

        :param timepoints: list

        :param randomize: bool
        :return:
        '''
        import random
        if timepoints is None:
            timepoints = list(self.timePoints)
        ltrees = [self.timePoints[t].trees(treenum) for t in timepoints]
        nids = ltrees[0][0].nodes.keys()

        nids2 = [_nid for _nid in nids]
        if randomize:
            random.shuffle(nids2)

        aves = defaultdict(lambda: defaultdict(list))
        for fpoint in timepoints:
            for ypoint in timepoints:
                for sample in self.timePoints[fpoint].samples.values():
                    ftree = sample.trees[treenum]
                    fnodes = [ftree.nodes[_nid] for _nid in nids2]
                    for sample2 in self.timePoints[ypoint].samples.values():
                        ytree = sample2.trees[treenum]
                        ynodes = [ytree.nodes[_nid] for _nid in nids]
                        aves[fpoint][ypoint].append(sum([fn.nfitness * yn.Y for (yn, fn) in zip(ynodes, fnodes)]))
                aves[fpoint][ypoint] = np.mean(aves[fpoint][ypoint])
                # self.logger("ave fitness " + fpoint + " " + ypoint + " " + str(aves[fpoint][ypoint]), 5)

        dff = defaultdict(lambda: [0.0, 0.0])
        dkl = defaultdict(lambda: [0.0, 0.0])

        if len(self.time_pairs) == 0:
            tpairs = zip(timepoints[:-1], timepoints[1:])
        else:
            tpairs = self.time_pairs

        self.logger("patient " + self.name, 5)
        self.logger(tpairs, 5)

        for (t1, t2) in tpairs:
            ldff = []
            lkl = []
            # iterates over all samples from the time point 1
            for samp1 in self.timePoints[t1].samples.values():
                tree1 = samp1.trees[treenum]
                # iterates over all samples from the time point 2
                for samp2 in self.timePoints[t2].samples.values():
                    dff_x = 0.0
                    dff_y = 0.0
                    kl1 = 0.0
                    kl2 = 0.0
                    tree2 = samp2.trees[treenum]

                    ys1 = [tree1.nodes[nid].Y + 1e-5 for nid in nids]
                    z = sum(ys1)
                    ys1 = [y / z for y in ys1]

                    ys2 = [tree2.nodes[nid].Y + 1e-5 for nid in nids2]
                    z = sum(ys2)
                    ys2 = [y / z for y in ys2]

                    for nid, nid2, y1, y2 in zip(nids, nids2, ys1, ys2):
                        ynodet1 = tree1.nodes[nid]
                        ynodet2 = tree2.nodes[nid]
                        kl1 += y1 * np.log(y1 / y2)
                        kl2 += y2 * np.log(y2 / y1)
                        nfitnessY = np.log(y2 / y1)
                        nfitnessX = np.log(ynodet2.X / max(1e-5, ynodet1.X))
                        dff_x += nfitnessY * (y2 - y1)
                        dff_y += nfitnessX * (y2 - y1)
                        self.logger(
                            "node " + " ".join([str(x) for x in [nid, t1, t2, nfitnessX, nfitnessY, y1, y2, y2 - y1]]),
                            5)
                    ldff.append([dff_y, dff_x])
                    lkl.append([kl1, kl2])

            ldff = [np.mean([x[0] for x in ldff]), np.mean([x[1] for x in ldff])]
            lkl = [np.mean([x[0] for x in lkl]), np.mean([x[1] for x in lkl])]
            dff[t1 + "-" + t2] = ldff
            dkl[t1 + "-" + t2] = lkl

        return dff, dkl, tpairs

    #    def compute_KL_between_time_points(self, tp1, tp2, eps):
    #        self.prepare_posterior(tp1, tp2, eps)
    #        tpoint1 = self.timePoints[tp1]
    #        tpoint2 = self.timePoints[tp2]
    #        trees1 = tpoint1.trees(n=-1)
    #        trees2 = tpoint2.trees(n=-1)
    #        tree_pairs = zip(trees1, trees2)

    def compute_KL(self, treenum=0, timepoints=None, randomize=False):
        '''
        Kullback-Leibler divergence between the time points

        :param treenum: int

        :param timepoints: list

        :param randomize: bool

        :return: dict
        '''
        import random
        if timepoints is None:
            timepoints = list(self.timePoints)
        ltrees = [self.timePoints[t].trees(treenum) for t in timepoints]
        nids = ltrees[0][0].nodes.keys()

        nids2 = [_nid for _nid in nids]
        if randomize:
            random.shuffle(nids2)
        # trees = {}
        dkl = defaultdict(float)
        tpairs = zip(timepoints[:-1], timepoints[1:])
        for (t1, t2) in tpairs:
            kls = []
            for samp1 in self.timePoints[t1].samples.values():
                tree1 = samp1.trees[treenum]
                for samp2 in self.timePoints[t2].samples.values():
                    kl = 0.0
                    tree2 = samp2.trees[treenum]
                    for nid, nid2 in zip(nids, nids2):
                        ynodet1 = tree1.nodes[nid]
                        ynodet2 = tree2.nodes[nid]

                        if ynodet2.Y > 0:
                            kl += ynodet2.Y * np.log(ynodet2.Y / max(1e-3, ynodet1.Y))
                        kls.append(kl)
            kls = np.mean(kls)
            dkl[t1] = kls
        return dkl

    def set_tree_self_copies(self):
        self.logger("Copying trees in patient " + self.pname, 5)
        for tpoint in self.timePoints.values():
            tpoint.set_tree_self_copies()
        self.logger("Done", 5)

    def get_tree_pairs(self, tp1, tp2, beta=1.0):
        '''
        Returns paired trees between time points

        :param tp1: str

        :param tp2: str

        :param beta: float

        :return: list
            (cfit.tree.SampleTree objects, cfit.tree.SampleTree objects, float)
        '''
        tpoint1 = self.timePoints[tp1]
        tpoint2 = self.timePoints[tp2]
        trees1 = tpoint1.trees(num=-1, just_one=True)
        trees2 = tpoint2.trees(num=-1, just_one=True)
        weights = np.array([tree.llh for tree in trees1])
        zw = Utils.log_sum(beta * weights)
        weights = np.exp(weights - zw)
        pairs = list(zip(trees1, trees2, weights))
        return pairs

    def observed_distance(self, tp1, tp2, eps=0.0, beta=1., normalize=False):
        '''

        :param tp1: str
            "from" time point

        :param tp2: str
            "to" time point

        :param eps: float
            threshold on frequency for clones in tp1 to be declared absent

        :param beta: float
            tree weight parameter

        :param normalize: bool
            whether to normalize the distance by entropy

        :return list
            [kl1, kl2, wd], float
            KL from tree1 to tree2, from tree2 to tree1 and EMD
        '''

        self.set_tree_self_copies()

        timePoint1 = self.timePoints[tp1]
        timePoint2 = self.timePoints[tp2]
        tissues1 = "_".join(timePoint1.tissues)
        tissues2 = "_".join(timePoint2.tissues)
        tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)
        res = []
        thr = 1e-10
        for tree1, tree2, w in tree_pairs:
            new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
            tree1.copy.soft_remove_nodes(new_clones)
            tree2.copy.soft_remove_nodes(new_clones)

            Y1 = [node.Y for node in tree1.copy.nodes.values()]
            Y2 = [node.Y for node in tree2.copy.nodes.values()]
            # p = [(y1, y2) for (y1, y2) in zip(Y1, Y2) if min(y1, y2) < 0.0001]
            # print(p)
            nn = len([x for x in zip(Y1, Y2) if sum(x) > 0])  # select non-empty clones
            norm = np.log(nn)
            # KL
            kl1 = sum([y1 * np.log(y1 / y2) for (y1, y2) in zip(Y1, Y2) if abs(y1) > thr and abs(y2) > thr and y1 != 0])
            kl2 = sum([y2 * np.log(y2 / y1) for (y1, y2) in zip(Y1, Y2) if abs(y1) > thr and abs(y2) > thr and y2 != 0])
            if normalize:
                kl1 /= norm
                kl2 /= norm

            # WD
            vals1 = [node.id for node in tree1.nodes.values()]
            vals2 = [node.id for node in tree2.nodes.values()]
            w1 = [node.Y for node in tree1.nodes.values()]
            w2 = [node.Y for node in tree2.nodes.values()]
            wd = wasserstein_distance(vals1, vals2, w1, w2)

            res.append([w, kl1, kl2, wd])

        kl1 = sum([x[0] * x[1] for x in res])
        kl2 = sum([x[0] * x[2] for x in res])
        wd = sum([x[0] * x[3] for x in res])

        return kl1, kl2, wd

    def distance(self, Yh, Y2, measure):
        '''
        Compute distance between clone frequncies

        :param Yh: list
            list of floats representing the probability distribution (prediction)
        :param Y2: list
            list of floats representing the probability distribution

        :param measure: str

        '''
        if measure == 'KL':
            dist = sum([y2 * np.log(y2 / yh) for (yh, y2) in zip(Yh, Y2) if y2 > 0])
        elif measure == "euclid":
            dist = np.sqrt(sum((Yh - Y2) ** 2))
        elif measure == "pearson":
            dist = np.corrcoef(Yh, Y2)[0, 1]
        elif measure == "spearman":
            dist, _ = spearmanr(Yh, Y2)
        elif measure == "EMD":
            vals = list(range(len(Yh)))
            dist = wasserstein_distance(vals, vals, Yh, Y2)
        return dist

    def log_distance(self, lYh, Y2, measure):
        '''

        '''
        Yh = [np.exp(lyh) for lyh in lYh]
        if measure == 'KL':
            dist = sum([y2 * (np.log(y2) - lyh) for (lyh, y2) in zip(lYh, Y2) if y2 > 0])
        elif measure == "euclid":
            dist = np.sqrt(sum((Yh - Y2) ** 2))
        elif measure == "pearson":
            dist = np.corrcoef(Yh, Y2)[0, 1]
        elif measure == "spearman":
            dist, _ = spearmanr(Yh, Y2)
        elif measure == "EMD":
            vals = list(range(len(Yh)))
            dist = wasserstein_distance(vals, vals, Yh, Y2)
        return dist

    def prediction_distance(self, tp1, tp2, eps, taus, beta=1., include_nested=True):
        '''

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param eps: float
            threshold on cluster size

        :param taus: list
            of floats

        :param beta: float

        :param include_nested: bool

        :return: tuple (dict, dict, float, float, str, str


        '''

        timePoint1 = self.timePoints[tp1]
        timePoint2 = self.timePoints[tp2]
        tissues1 = "_".join(timePoint1.tissues)
        tissues2 = "_".join(timePoint2.tissues)

        tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)

        dres = {}
        dres_obs = {}
        measures = ["KL", "euclid", "pearson", "spearman", "EMD"]
        for measure in measures:
            dres[measure] = np.array([0. for _ in taus])
            dres_obs[measure] = np.array([0. for _ in taus])
        #        self.logger("Predictions for patient " + self.name)
        Z1, Z2 = 0.0, 0.0
        for tree1, tree2, w in tree_pairs:
            new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
            z1 = tree1.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
            z2 = tree2.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
            Z1 += w * z1
            Z2 += w * z2
            # Y1 = np.array([node.Y for node in tree1.copy.nodes.values() if node.id != tree1.copy.root.id])
            nids = list(tree1.nodes.keys())
            nodes1 = [tree1.nodes[nid] for nid in nids]
            nodes2 = [tree2.nodes[nid] for nid in nids]
            nodes = [ns for ns in zip(nodes1, nodes2) if ns[0].cY > 0 and ns[0].id != tree1.root.id]
            nodes1 = [n1 for (n1, n2) in nodes]
            nodes2 = [n2 for (n1, n2) in nodes]
            Y1 = np.array([node.cY for node in nodes1])
            Y2 = np.array([node.cY for node in nodes2])
            lY1 = np.log(Y1)
            #            Y2 = np.array([node.cY for node in tree2.nodes.values() if node.id != tree2.root.id])
            dists = {}
            obs_dists = {}
            for measure in measures:
                dists[measure] = []
                obs_dists[measure] = []
            for tau in taus:
                node_logfit = lambda node: -math.inf if node.cY == 0 else np.log(node.cY) + node.fitness * tau
                #                Yh =[node_logfit(node) for node in tree1.nodes.values() if node.id != tree1.root.id]
                #                fit = [node.fitness for node in nodes1]
                Yh0 = [node_logfit(node) for node in nodes1]
                lYh = Utils.log_norm(Yh0)
                #                Yh = np.array([np.exp(yh) for yh in Utils.log_norm(Yh0)])
                #                nids = [node.id for node in nodes1]
                for measure in dres:
                    dists[measure].append(self.log_distance(lYh, Y2, measure))
                    obs_dists[measure].append(self.log_distance(lY1, Y2, measure))
            for measure in measures:
                dres[measure] += w * np.array(dists[measure])
                dres_obs[measure] += w * np.array(obs_dists[measure])
        # tmetrics = pd.DataFrame([pmetrics[metric] for metric in pmetrics])
        # tmetrics.columns = list(pmetrics.keys())

        #        klhs = [sum([x[0] * x[4 + i] for x in res]) for i in range(len(taus))]
        #        sites1 = list(sites1)
        #        sites1.sort()
        #        sites1 = ",".join(sites1)
        #        sites2 = list(sites2)
        #        sites2.sort()
        #        sites2 = ",".join(sites2)

        return dres, dres_obs, Z1, Z2, tissues1, tissues2

    def distribution_distance(self, tp1, tp2, eps, weights, tau,
                              reference_time_point=2,
                              beta=1., include_nested=True, measure="KL", tree_pairs=None):
        '''

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param eps: float
            threshold on cluster size

        :param weights: dict: str -> float
            dictionary mapping fitness components to their weights

        :param tau: float

        :param reference_time_point: int
            1 -> tp1, 2->tp2

        :param beta: float

        :param include_nested: bool

        :param measure: str
            distance measure to use, defaults to Kullback Leibler divergence

        :param tree_pairs: list
            list of paired trees from time points tp1 and tp2, already initialized, with shared clones marked
            [(tree1, tree2, llh),...]
        :return: float
            distance between the predicted distribution and the reference distribution
        '''

        if tree_pairs is None:
            tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)
            for tree1, tree2, w in tree_pairs:
                new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
                tree1.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
                tree2.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)

        di = 0
        for tree1, tree2, w in tree_pairs:
            nids = list(tree1.nodes.keys())
            nodes1 = [tree1.nodes[nid] for nid in nids]
            nodes2 = [tree2.nodes[nid] for nid in nids]
            nodes2 = [nodes1, nodes2][reference_time_point - 1]
            for node in nodes1:
                node.f = node.fitness
                node.f = sum([node.fitness_components[name] * weights[name] for name in weights])

            nodes = [ns for ns in zip(nodes1, nodes2) if ns[0].cY > 0 and ns[0].id != tree1.root.id]
            nodes1 = [n1 for (n1, n2) in nodes]
            nodes2 = [n2 for (n1, n2) in nodes]

            Y2 = np.array([node.cY for node in nodes2])
            ys = [node.cY for node in nodes1]
            fs = [node.f for node in nodes1]
            Yh0 = [-math.inf if y == 0 else np.log(y) + f * tau for y, f in zip(ys, fs)]
            lYh = Utils.log_norm(Yh0)
            di1 = self.log_distance(lYh, Y2, measure)
            di += w * di1
        return di

    def initialize_paired_trees(self, tp1, tp2, eps, beta=1., include_nested=True):
        '''

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param eps: float
            threshold on cluster frequency to define new clones

        :param beta: float
            weight parameter for trees

        :param include_nested: bool

        :return: list
            [(tree1, tree2, llh), ...]
        '''

        tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)

        tree_pairs2 = []
        for i, (tree1, tree2, w) in enumerate(tree_pairs):
            self.logger("tree pair "+str(i)+" weight "+str(w))
            new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
            self.logger("Removing new clones from tree1")
            tree1.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
            for node in tree1.nodes.values():
                node.logger("\tnode "+str(node.id)+" "+str(node.cY)+" "+str(node))

            self.logger("Removing new clones from tree2")
            tree2.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
            for node in tree2.nodes.values():
                node.logger("\tnode "+str(node.id)+" "+str(node.cY)+" "+str(node))
            tree_pairs2.append((tree1, tree2, w))
        return tree_pairs2

    def identify_new_clones(self, tree1, tree2, eps=0.03):
        '''
        Identify the clones that are new in tree2 and not shared with tree1

        :param tree1: SampleTree
            eg. first time point tree

        :param tree2: SampleTree
            eg. second time point tree

        :param eps: float
            threshold on frequency in tree1

        :return: list
            list on identifiers of clones detected as new in tree2
        '''
        pnodes = zip(tree1.nodes.values(), tree2.nodes.values())
        new_clones = [node1.id for (node1, node2) in pnodes if
                      node1.X < eps and node2.X >= eps and node2.X > node1.X]
        return new_clones


    def example(self, tp1, tp2, eps, params, beta=1., include_nested=True, tree_pairs=None):
        '''

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param eps: float
            threshold on cluster size

        :param params: dict
            {"tau": None/float,
            "weights":
                {<component_name>: str -> <component_value>: None/float}
            }
            dictionary defining the set of parameters to optimize over or to fix

        :param beta: float

        :param include_nested: bool

        :param tree_pairs: list
            if the tree pairs have been already initialized

        :return: tuple (params: dict, optimized_distance: float)
            dict: {"tau": str -> float, "weights": {<component_name>: str -> <component_weight> :float}}
            the dictionary of parameters includes both optimized and fixed components

        '''

        if tree_pairs is None:
            tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)

            for tree1, tree2, w in tree_pairs:
                new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
                new_clones = set([node.id for node in tree1.nodes.values() if node.X < eps]).union(new_clones)
                tree1.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
                tree2.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)

            # node_logfit = lambda node, tau: -math.inf if node.cY == 0 else np.log(node.cY) + node.fitness * tau
        (tree1, tree2, w) = tree_pairs[0]
        weights = params["weights"]
        nodes1 = tree1.nodes.values()
        nodes2 = [tree2.nodes[node.id] for node in nodes1]
        fs = [sum([node1.fitness_components[name] * weights[name] for name in weights]) for node1 in nodes1]
        avef = sum([node1.Y * f for (node1, f) in zip(nodes1, fs)])

        dat = [[node1.Y, node2.Y, node1.X, node2.X, f, f - avef] for (node1, node2, f) in zip(nodes1, nodes2, fs)]
        dat = pd.DataFrame(dat)
        dat.columns = ["Y1", "Y2", "X1", "X2", "F", "nF"]
        return dat

    def optimized_prediction_distance(self, tp1, tp2, eps,
                                      params,
                                      beta=1., include_nested=True, tree_pairs=None,
                                      inparams0=None):
        '''
        Find fitness model weights that minimize the KL-distance between prediction from time point 1 and
        clone composition of time point 2.

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param eps: float
            threshold on clone size in time point 1. Smaller clones are removed with assumption they are new in
            time point 2.

        :param params: dict
            {"tau": None/float,
            "weights":
                {<component_name>: str -> <component_value>: None/float}
            }
            dictionary defining the set of parameters to optimize over or to fix

        :param beta: float

        :param include_nested: bool

        :param tree_pairs: list
            if the tree pairs have been already initialized

        :return: tuple (params: dict, optimized_distance: float)
            dict: {"tau": str -> float, "weights": {<component_name>: str -> <component_weight> :float}}
            the dictionary of parameters includes both optimized and fixed components

        '''

        if tree_pairs is None:

            tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)

            for tree1, tree2, w in tree_pairs:
                new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
                tree1.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
                tree2.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)

        def distance_fun(x, *args):
            '''

            :param x: tuple of floats
                (weights), (tau, weights)
                parameters to be optimized, can include tau, proceeded by components weights

            :param args: list
                (tree_pairs: list,
                distance_measure: str,
                tau: float if fixed or None if tau be optimized,
                optimized_fitness_components: list of str
                    names of the components to be optimized
                fixed_components: dict: str->float
                    names of the fixed components
                )
                ref_time_point: 1/2 - if 1, the starting time point, if 2, the target time point is used as the reference for distance
            :return:
            '''

            tree_pairs = args[0]
            measure = args[1]
            tau = args[2]
            optimized_fitness_components = args[3]
            fixed_components = args[4]

            ind = 0
            if tau is None:
                tau = x[0]
                ind = 1

            weights = {}
            nlen = len(optimized_fitness_components)
            for i, comp in enumerate(optimized_fitness_components):
                if i == nlen - 1:
                    weights[comp] = 1 - sum(weights.values())
                else:
                    weights[comp] = x[ind + i]

            for comp in fixed_components:
                weights[comp] = fixed_components[comp]

            di = 0
            isigma = 0.3
            #            if "immune" in weights:
            #                pen = (weights["immune"]*isigma)**2
            #                di += pen

            for tree1, tree2, w in tree_pairs:
                nids = list(tree1.nodes.keys())
                nodes1 = [tree1.nodes[nid] for nid in nids]
                nodes2 = [tree2.nodes[nid] for nid in nids]
                for node in nodes1:
                    node.f = node.fitness
                    node.f = sum([node.fitness_components[name] * weights[name] for name in weights])
                nodes = [ns for ns in zip(nodes1, nodes2) if ns[0].cY > 0 and ns[0].id != tree1.root.id]
                nodes1 = [n1 for (n1, n2) in nodes]
                nodes2 = [n2 for (n1, n2) in nodes]
                Y2 = np.array([node.cY for node in nodes2])
                ys = [node.cY for node in nodes1]
                fs = [node.f for node in nodes1]
                Yh0 = [-math.inf if y == 0 else np.log(y) + f * tau for y, f in zip(ys, fs)]
                lYh = Utils.log_norm(Yh0)
                di1 = self.log_distance(lYh, Y2, measure)
                di += w * di1
            #            self.logger("f(" + ",".join([str(pa) for pa in x])+") = "+str(di))
            return di

        def kl0(tree_pairs, measure):
            di = 0
            for tree1, tree2, w in tree_pairs:
                nids = list(tree1.nodes.keys())
                nodes1 = [tree1.nodes[nid] for nid in nids]
                nodes2 = [tree2.nodes[nid] for nid in nids]
                nodes = [ns for ns in zip(nodes1, nodes2) if ns[0].cY > 0 and ns[0].id != tree1.root.id]
                nodes1 = [n1 for (n1, n2) in nodes]
                nodes2 = [n2 for (n1, n2) in nodes]
                Y2 = np.array([node.cY for node in nodes2])
                Yh = [node.cY for node in nodes1]
                #                Yh0 = [-math.inf if y == 0 else np.log(y) for y in ys]
                #                lYh = Utils.log_norm(Yh0)
                di1 = self.distance(Yh, Y2, measure)
                di += w * di1
            return di

        # initialize parameters

        zerodi = kl0(tree_pairs, measure='KL')

        x0 = []
        xnull0 = []
        optimized_components = []
        optimized_fitness_components = []
        lb = []
        ub = []
        diag = []

        if params["tau"] is None:
            optimized_components = ["tau"]

            lb.append(0.)
            ub.append(np.infty)
            diag.append(1.)

            if inparams0 is None:
                x0.append(0)

            elif inparams0["tau"] is not None:
                x0.append(inparams0["tau"])

            xnull0.append(0)

        fixed_components = {}
        tot = 1

        for comp in params["weights"]:
            if params["weights"][comp] is None:
                optimized_components.append(comp)
                optimized_fitness_components.append(comp)
                # sum_constr.append(1)

            else:
                fixed_components[comp] = params["weights"][comp]
                tot -= fixed_components[comp]
        for _ in optimized_fitness_components[:-1]:
            diag.append(1.)
            lb.append(0)
            ub.append(1)
            x0.append(tot / len(optimized_fitness_components))
            xnull0.append(tot / len(optimized_fitness_components))

        x0 = np.array(x0)
        xnull0 = np.array(xnull0)
        linearConstraints = []
        if len(diag) > 0:
            diag = np.diag(diag)
            linear_constraint1 = LinearConstraint(diag, lb, ub)
            linearConstraints.append(linear_constraint1)
        #        self.logger(params)
        #        self.logger("linear constraint: " + str(diag) + " " + str(list(zip(lb, ub))))
        minkwargs = {"args": (tree_pairs, 'KL', params["tau"],
                              optimized_fitness_components, fixed_components),
                     "constraints": linearConstraints}

        solved = False
        while not solved:
            #            self.logger("Optimizing with (tau-" + str(params["tau"] is None) + " " + str(optimized_fitness_components)+")")
            #            self.logger(x0)
            val = basinhopping(distance_fun, x0, niter=500,
                               minimizer_kwargs=minkwargs)
            ind = 0
            x0 = []
            opttau = None
            if params["tau"] is None:
                opttau = abs(val["x"][0])
                x0.append(abs(val["x"][0]))
                ind = 1

            weights = {}
            for i, comp in enumerate(optimized_fitness_components[:-1]):
                weights[comp] = abs(val["x"][ind + i])
                x0.append(weights[comp])
            if len(optimized_fitness_components) > 0:
                weights[optimized_fitness_components[-1]] = 1. - sum(weights.values())

            for comp in fixed_components:
                weights[comp] = fixed_components[comp]
            x0 = np.array(x0)
            optdi = distance_fun(x0, *(tree_pairs, 'KL',
                                       params["tau"], optimized_fitness_components, fixed_components))
            # self.logger("optdi: " + str(optdi) + " zerodi: " + str(zerodi))
            if optdi > zerodi and opttau is not None:
                if abs(optdi - zerodi) > 1e-5:
                    x0 = xnull0
                    # self.logger("optimization not solved, repeating")
                    print(val)
                else:
                    solved = True
            else:
                solved = True

        optparams = {"tau": opttau, "weights": weights}

        #        self.logger(optparams)
        return optparams, optdi

    def fitness_flux(self, tp1, tp2, beta=1., include_nested=True, tree_pairs=None):
        '''
        Compute frequency projections, save them on SampleNode.predictedY and .predictedY attributes

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param eps: float
            threshold on cluster size

        :param params: dict
            {"tau": None/float,
            "weights":
                {<component_name>: str -> <component_value>: None/float}
            }
            dictionary defining the set of parameters to optimize over or to fix

        :param beta: float

        :param include_nested: bool

        :param tree_pairs: list
            if the tree pairs have been already initialized

        '''

        if tree_pairs is None:
            tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)

        #weights = params["weights"]
        #tau = params["tau"]
        phis = []
        for tree1, tree2, w in tree_pairs:
            tree2.set_primary_ancestors()
            ancnids = set([node.anc for node in tree2.nodes.values()])
            nids = list(tree1.nodes.keys())
            nodes1 = [tree1.nodes[nid] for nid in nids]
            nodes2 = [tree2.nodes[nid] for nid in nids]

            for node in nodes1:
                node.ancY = 0.
            ancnodes1 = [tree1.nodes[nid] for nid in ancnids]
            Z = sum(node.Y for node in ancnodes1)
            self.logger("Z="+str(Z))
            if Z>0:
                for node in ancnodes1:
                    node.ancY = node.Y/Z

#            for node in nodes1:
#                node.f = sum([node.fitness_components[name] * weights[name] for name in weights])
#            for node in nodes2:
#                node.f = sum([node.fitness_components[name] * weights[name] for name in weights])
            nodes = list(zip(nodes1, nodes2))
            prim_nodes = [(node1, node2) for (node1, node2) in nodes if node2.sharedY>0 or node1.id == tree1.root.id]
            anc2_Y_rec = defaultdict(float)
            for node2 in nodes2:
                anc2_Y_rec[node2.anc] += node2.Y
            Phi1 = sum([node1.fitness*(anc2_Y_rec[node1.id] - node1.Y) for (node1, node2) in prim_nodes])
            Phi1p = sum([node1.fitness*(anc2_Y_rec[node1.id] - node1.ancY) for (node1, node2) in prim_nodes])
            Phi2 = sum([(node2.fitness-tree1.nodes[node2.anc].fitness)*node2.Y for node2 in nodes2])
            phis.append((Phi1, Phi1p, Phi2, w))
        Phi1 = sum([el[0]*el[3] for el in phis])
        Phi1p = sum([el[1]*el[3] for el in phis])
        Phi2 = sum([el[2]*el[3] for el in phis])
        return Phi1, Phi1p, Phi2



    def predict(self, tp1, tp2, eps, params, beta=1., include_nested=True, tree_pairs=None):
        '''
        Compute frequency projections, save them on SampleNode.predictedY and .predictedY attributes

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param eps: float
            threshold on cluster size

        :param params: dict
            {"tau": None/float,
            "weights":
                {<component_name>: str -> <component_value>: None/float}
            }
            dictionary defining the set of parameters to optimize over or to fix

        :param beta: float

        :param include_nested: bool

        :param tree_pairs: list
            if the tree pairs have been already initialized

        '''

        if tree_pairs is None:
#            tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)

        #            for tree1, tree2, w in tree_pairs:
        #                new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
        #                tree1.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
        #                tree2.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
            tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)

            for tree1, tree2, w in tree_pairs:
                new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
                tree1.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)
                tree2.soft_remove_nodes(new_clones, include_nested=include_nested, overwrite=False)


        weights = params["weights"]
        tau = params["tau"]
        for tree1, tree2, w in tree_pairs:

            # for node in tree2.nodes.values():
            #    node.predictedY = 0

            nids = list(tree1.nodes.keys())
            nodes1 = [tree1.nodes[nid] for nid in nids]
            nodes2 = [tree2.nodes[nid] for nid in nids]
            for node in nodes1:
                # node.f = node.fitness
                node.f = sum([node.fitness_components[name] * weights[name] for name in weights])
            nodes = list(zip(nodes1, nodes2))
            nodes1 = [n1 for (n1, n2) in nodes]
            nodes2 = [n2 for (n1, n2) in nodes]
            # ys = [node.cY for node in nodes1]
            #            ys = [node.Y for node in nodes1]
            #            fs = [node.f for node in nodes1]
            for node1, node2 in nodes:
                node2.yh0 = -math.inf if node1.Y == 0 else np.log(node1.Y) + node1.f * tau
            lZ = Utils.log_sum([node2.yh0 for node2 in nodes2])
            for node2 in nodes2:
                node2.predictedY = np.exp(node2.yh0 - lZ)
            #            Yh0 = [-math.inf if y == 0 else np.log(y) + f * tau for y, f in zip(ys, fs)]
            #            lYh = Utils.log_norm(Yh0)
            #
            #            for lyh, node2 in zip(lYh, nodes2):
            #                node2.predictedY = np.exp(lyh)

            for tnode in tree2.tree.get_post_order():
                node = tree2.nodes[tnode.id]
                cnodes = [tree2.nodes[tcnode.id] for tcnode in tnode.children]
                node.predictedX = node.predictedY + sum([cnode.predictedX for cnode in cnodes])

    def compute_fitness(self, params):
        '''
        Compute fitness and rfitness attributes for all nodes of all trees of all samples
        :param params: dict
            {"tau": None/float,
            "weights":
                {<component_name>: str -> <component_value>: None/float}
            }
            dictionary defining the set of parameters to optimize over or to fix
        '''

        trees = self.all_trees
        for tree in trees:
            nodes = tree.nodes.values()
            for node in nodes:
                node.fitness = sum(
                    [node.fitness_components[name] * params.weights[name] * params["tau"] for name in params.weights])
            avef = [node.fitness * node.Y for node in nodes]

            for node in nodes:
                node.rfitness = node.fitness - avef

    def paired_clones(self, tp1, tp2):
        tpoint1 = self.timePoints[tp1]
        tpoint2 = self.timePoints[tp2]
        tree1 = tpoint1.trees(num=0)[0]
        tree2 = tpoint2.trees(num=0)[0]
        pclones = []
        for nid in tree1.nodes.keys():
            node1 = tree1.nodes[nid]
            node2 = tree2.nodes[nid]
            pclones.append([node1, node2])
        return pclones

    def compute_predicted_fitness_flux(self, treenum=0, timepoints=None, randomize=False):
        '''

        :param treenum: int

        :param timepoints: list

        :param randomize: bool

        :return:
        '''
        import random
        if timepoints is None:
            timepoints = list(self.timePoints)
        ltrees = [self.timePoints[t].trees(treenum) for t in timepoints]
        nids = ltrees[0][0].nodes.keys()

        nids2 = [_nid for _nid in nids]
        if randomize:
            random.shuffle(nids2)
        # trees = {}
        aves = defaultdict(lambda: defaultdict(list))
        for fpoint in timepoints:
            # trees[tpoint1] = self.timePoints[tpoint1].trees[treenum]
            for ypoint in timepoints:
                for sample in self.timePoints[fpoint].samples.values():
                    ftree = sample.trees[treenum]
                    fnodes = [ftree.nodes[_nid] for _nid in nids2]
                    for sample2 in self.timePoints[ypoint].samples.values():
                        ytree = sample2.trees[treenum]
                        ynodes = [ytree.nodes[_nid] for _nid in nids]
                        aves[fpoint][ypoint].append(sum([fn.nfitness * yn.Y for (yn, fn) in zip(ynodes, fnodes)]))
                aves[fpoint][ypoint] = np.mean(aves[fpoint][ypoint])

        dff = defaultdict(lambda: [0.0, 0.0])
        if len(self.time_pairs) == 0:
            tpairs = zip(timepoints[:-1], timepoints[1:])
        else:
            tpairs = self.time_pairs

        clone_flux = []
        for (t1, t2) in tpairs:
            ldff = []
            for samp1 in self.timePoints[t1].samples.values():
                tree1 = samp1.trees[treenum]
                for samp2 in self.timePoints[t2].samples.values():
                    dff1 = [0, 0, 0, 0]
                    tree2 = samp2.trees[treenum]
                    for nid, nid2 in zip(nids, nids2):
                        ynodet1 = tree1.nodes[nid]
                        ynodet2 = tree2.nodes[nid]
                        fnodet1 = tree1.nodes[nid2]
                        fnodet2 = tree2.nodes[nid2]
                        div = 1 / max(1e-03, self.PFS)
                        ff1 = (fnodet1.nfitness - aves[t1][t1]) * (ynodet2.Y - ynodet1.Y)
                        ff2 = (fnodet2.nfitness - aves[t2][t1]) * (ynodet2.Y - ynodet1.Y)
                        dff1 = [dff1[0] + ff1, dff1[1] + ff2]  # , dff1[2] + normff1, dff1[3] + normff2]
                        ldff.append(dff1)
                        clone_flux.append(
                            [self.name, samp1.name, samp2.name, nid, ynodet1.nfitness, ynodet1.Y, ynodet2.Y])
            ldff = [np.mean([x[0] for x in ldff]), np.mean([x[1] for x in ldff])]
            dff[t1 + "-" + t2] = ldff

        return dff, clone_flux, tpairs

    '''
    Longitudinal analysis
    '''

    def predict_tree(self, tp1, tp2, tau):
        '''
        Computes distance between the predicted and observed trees

        :param tp1: str
            baseline time point

        :param tp2: str
            predictions time point

        :param tau: float
            time parameter


        :return: list
            [relative entropy to timepoint1, relative entropy to the true tree]

        '''
        tpoint1 = self.timePoints[tp1]
        tpoint2 = self.timePoints[tp2]
        res = []
        for sample in tpoint1.samples.values():
            trees = sample.trees
            predicted_trees = [tree.predicted_tree(Node.get_fitness, tau) for tree in trees]
            weights = sample.get_tree_weights(beta=1.)
            for sample2 in tpoint2.samples.values():
                trees2 = sample2.trees
                re = sum([tree2.relative_entropy(ptree) * w for (tree, tree2, ptree, w) in
                          zip(trees, trees2, predicted_trees, weights)])
                re0 = sum([tree2.relative_entropy(tree) * w for (tree, tree2, ptree, w) in
                           zip(trees, trees2, predicted_trees, weights)])
                self.logger("Relative entropy " + self.name + " " + tp1 + " " + tp2, 5)
                res.append([re, re0])

        res = [np.mean([x[0] for x in res]), np.mean([x[1] for x in res])]
        return res

    def nclones(self, beta=1.):
        '''
        Return the number of clones

        :param beta: float
            tree weighting parameter

        :return: float
        '''

        tp = list(self.timePoints.values())[0]
        trees = tp.trees(num=-1, just_one=True)
        weights = [tree.llh for tree in trees]
        zw = Utils.log_sum([beta * w for w in weights])
        weights = [np.exp(w - zw) for w in weights]
        ns = sum([(len(tree.nodes.values()) - 1) * w for tree, w in zip(trees, weights)])
        return ns

    def mark_shared_clones(self, tp1, tp2, eps):
        '''

        :param tp1: str
            from (source tumor)

        :param tp2: str
            to (target tumor)

        :param eps: float
            threshold on the clone in the source to be considered

        :return:
            The node of trees in samples in tp2 get an attribute, sharedY
        '''

        tpoint1 = self.timePoints[tp1]
        tpoint2 = self.timePoints[tp2]

        ntrees = len([*tpoint1.samples.values()][0].trees)
        self.logger("marking shared clones " + tp1 + " -> " + tp2 + ", number of tress=" + str(ntrees))
        # prepare the weights for trees
        for sample in tpoint2.samples.values():
            for tree in sample.trees:
                tree.source_Ys = defaultdict(dict)
                for node in tree.nodes.values():
                    node.privateY = 0.0
                    node.sharedY = 0.0
                    node.shared_clone = 0

        for i in range(ntrees):  # for each tree (these are shared with all other samples)
            # time point 1 data - get average clone frequencies
            self.logger(self.name + " tree " + str(i))
            ave_Ys = defaultdict(float)
            ave_Ys2 = defaultdict(float)

            lYs = defaultdict(list)
            for sample in tpoint1.samples.values():
                Ys = sample.trees[i].get_clone_Ys()
                for cid in Ys:
                    lYs[cid].append(Ys[cid])
            for cid in lYs:
                ave_Ys[cid] = np.mean(lYs[cid])

            lYs = defaultdict(list)
            for sample in tpoint2.samples.values():
                Ys = sample.trees[i].get_clone_Ys()
                for cid in Ys:
                    lYs[cid].append(Ys[cid])
            for cid in lYs:
                ave_Ys2[cid] = np.mean(lYs[cid])

            source_cids = set([cid for cid in ave_Ys if ave_Ys[cid] >= eps and ave_Ys2[cid] >= eps])
            private_cids = set([cid for cid in ave_Ys if ave_Ys[cid] < eps and ave_Ys2[cid] > eps])
            self.logger(self.name + ", tree " + str(i) + ": private clones in time point " + str(tp2) + ":")
            self.logger([(cid, ave_Ys[cid], ave_Ys2[cid]) for cid in private_cids])
            for sample1 in tpoint1.samples.values():
                tt1 = sample1.trees[i]
                for node in tt1.nodes.values():
                    node.shared_clone = 0
                for nid in source_cids:
                    tt1.nodes[nid].shared_clone = 1

            for sample2 in tpoint2.samples.values():
                tt2 = sample2.trees[i]
                z_shared = sum([tt2.nodes[cid].Y for cid in source_cids])
                z_private = sum([tt2.nodes[cid].Y for cid in private_cids])
                for node in tt2.nodes.values():
                    node.privateY = 0.0
                    node.sharedY = 0.0
                    node.shared_clone = 0

                for node in tt2.nodes.values():
                    node.shared_clone = 0
                    if node.id in source_cids:
                        node.shared_clone = 1
                        node.privateY = 0.0
                        if z_shared > eps:
                            node.sharedY = node.Y / z_shared
                    elif node in private_cids:
                        node.sharedY = 0.0
                        if z_private > eps:
                            node.privateY = node.Y / z_private
        return True

    def average_function(self, tp1, tp2, node_function, beta=1., use_shared=False, use_new=False, eps=0.0, **kwargs):
        '''
        A meta function - the average function of time point 1 in time point 2, computed
        as the the weighted value of the function of tp1 with weights being the exclusive frequencies from tp2:
        score(tp1|tp2) = \sum_\alpha Y_\alpha(tp_2) * function_\alpha(tp1)

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param node_function: function
            function applied on a node

        :param beta: float
            tree weighting parameter

        :param use_shared: bool
            whether to use only shared clones

        :param use_new: bool
            whether to use only new clones in tp2

        :param kwargs: dict
            positional parameters passed to node_function

        :return: (float, float)
            - average of nodes in tree2 with freqs from tree1
            - average of nodes in tree1 with freqs from tree2
        '''

        tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)

        if use_shared:
            self.logger("Setting copies in " + self.name, 5)
            self.set_tree_self_copies()

        scores = []
        for tree1, tree2, w in tree_pairs:
            if use_shared:
                # self.logger("removing soft", 5)
                new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
                tree1.copy.soft_remove_nodes(new_clones, include_nested=True)
                # self.logger("done", 5)
                tree2.copy.soft_remove_nodes(new_clones, include_nested=True)
                nodes1 = tree1.copy.nodes.values()
                nodes2 = tree2.copy.nodes.values()
            elif use_new:
                # self.logger("removing soft", 5)
                new_clones = self.identify_new_clones(tree1, tree2, eps=eps)
                nodes1 = [tree1.nodes[nid] for nid in new_clones]
                nodes2 = [tree2.nodes[nid] for nid in new_clones]
            else:
                nodes1 = tree1.nodes.values()
                nodes2 = tree2.nodes.values()
            score1 = sum([node_function(node1, **kwargs) * node1.Y for (node1, node2) in zip(nodes1, nodes2)])
            score2 = sum([node_function(node1, **kwargs) * node2.Y for (node1, node2) in zip(nodes1, nodes2)])
            z1 = sum([node.Y for node in nodes1])
            z2 = sum([node.Y for node in nodes2])
            if z1 > 0:
                score1 /= z1
            if z2 > 0:
                score2 /= z2
            scores.append([w, score1, score2])
        score1 = sum([s[0] * s[1] for s in scores])
        score2 = sum([s[0] * s[2] for s in scores])
        return score1, score2

    def get_clone_info(self, tp1, tp2, beta=1., gene_list=[], mutation_list=[], in_exclusive=True):
        '''
        Get info about clones
        :param tp1: str

        :param tp2: str

        :param beta: float
            tree weighting parameter

        :param gene_list: list
            list of gene names to be checked - clones with a mutation in these genes will be marked with "1"

        :param mutation_list: list
            list of mutation names (<chr>_<position>_<ref>_<alt>) to be checked

        :param in_exclusive: bool
            report only new mutations

        :return: pd.DataFrame
        '''
        tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)
        cinfo = []
        for tree1, tree2, w in tree_pairs:
            cids = list(tree1.nodes.keys())
            avef = tree1.average_over_nodes(node_fun=lambda node: node.fitness)
            aveml = tree1.average_over_nodes(node_fun=lambda node: len(node.mutations))
            avenl = tree1.average_over_nodes(node_fun=lambda node: len(node.neoantigens))
            for cid in cids:
                node1 = tree1.nodes[cid]
                node2 = tree2.nodes[cid]
                mload = len(node1.mutations)
                nload = len(node2.neoantigens)
                if in_exclusive:
                    driver_genes = [
                        sum([mut.gene == gene for mut in node1.exclusiveMutations if mut.substitution != ""]) for
                        gene in gene_list]
                    driver_mutations = [sum([mut.id == mid for mut in node1.exclusiveMutations]) for
                                        mid in mutation_list]
                else:
                    driver_genes = [sum([mut.gene == gene for mut in node1.mutations if mut.substitution != ""])
                                    for
                                    gene in gene_list]
                    driver_mutations = [sum([mut.id == mid for mut in node1.mutations]) for
                                        mid in mutation_list]

                line = [cid, self.name, tp1 + "->" + tp2, self.PFS, self.OS, w, node1.X, node2.X, node1.Y, node2.Y,
                        mload, mload - aveml, nload, nload - avenl, node1.fitness,
                        node1.fitness - avef] + driver_genes + driver_mutations
                cinfo.append(line)
        cinfo = pd.DataFrame(cinfo)
        cinfo.columns = ["Clone", "Patient", "TimePoints", "PFS", "OS", "Weight", "X1", "X2", "Y1", "Y2", "Mload",
                         "aveMload",
                         "Nload", "aveNload", "F", "aveF"] + gene_list + mutation_list
        return cinfo

    def trunkality(self, tp1, tp2, beta=1., use_shared=False, eps=0.0):

        '''
        Measure for the "trunkality" of time point 1 in time point 2, computed
        as the the weighted CCF of tp1 with weights being the exclusive frequencies from tp2:
        trunkality(tp1|tp2) = \sum_\alpha Y_\alpha(tp_2) * X_\alpha(tp1)

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param beta: float
            tree weighting parameter

        :param use_shared: bool
            whether to use only shared clones

        :return: (float, float)
            the trunkality measure (tp1->tp2, tp2->tp1)
        '''

        node_function = lambda node: node.X
        score1, score2 = self.average_function(tp1, tp2, node_function, beta=beta, use_shared=use_shared, eps=eps)
        return score1, score2

    def mark_clone_sharing(self, tp1, tp2, eps, dir=1):
        '''
        Fills in sharedY and privateY attributes of nodes in the trees between time points tp1 and tp2.
        sharedY sum up to 1 over clones, so to privateY.

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param eps: float
            threshold on cluster size to decide it exists

        :param dir: int
            1 - Pre->Post, sets in sharedY and privateY in tp2, 0: Post -> Pre, sets in preservedY and lostY in tp1

        :return: bool
            The nodes of trees in samples in tp1/tp2 get the attributes
        '''

        if dir:
            tpoint1 = self.timePoints[tp1]
            tpoint2 = self.timePoints[tp2]
        else:
            tpoint1 = self.timePoints[tp2]
            tpoint2 = self.timePoints[tp1]
        ntrees = len([*tpoint1.samples.values()][0].trees)
        if dir:
            self.logger("marking shared and private clones " + tp1 + " -> " + tp2 + " ntrees=" + str(ntrees))
        else:
            self.logger("marking lost and preserved clones" + tp1 + " -> " + tp2 + " ntrees=" + str(ntrees))
        # prepare the node frequencies - initialization
        for sample in tpoint2.samples.values():  # iterating over samples in the second time point
            for tree in sample.trees:  # iterating over the trees in the second time point sample
                tree.source_Ys = defaultdict(dict)
                for node in tree.nodes.values():
                    if dir:
                        node.sharedY = 0.0
                        node.privateY = 0.0
                    else:
                        node.preservedY = 0.0
                        node.lostY = 0.0

        dYs = defaultdict(lambda: defaultdict(float))
        dYs2 = defaultdict(lambda: defaultdict(float))
        dXs = defaultdict(lambda: defaultdict(float))
        dXs2 = defaultdict(lambda: defaultdict(float))
        # clone frequencies in time point 1 - averaged over samples if
        # more than 1 exists.

        for i in range(ntrees):  # for each tree (the tree topology is commone for all other samples)
            # time point 1 data - compute average clone frequencies
            dYs[i] = defaultdict(float)
            dYs2[i] = defaultdict(float)
            dXs[i] = defaultdict(float)
            dXs2[i] = defaultdict(float)

            lYs = defaultdict(list)
            lXs = defaultdict(list)
            for sample in tpoint1.samples.values():  # iterate over samples in this time point (usually just 1)
                Ys = sample.trees[i].get_clone_Ys()
                Xs = sample.trees[i].get_clone_Xs()
                for cid in Ys:
                    lYs[cid].append(Ys[cid])
                    lXs[cid].append(Xs[cid])
            for cid in lYs:
                dYs[i][cid] = np.mean(lYs[cid])
                dXs[i][cid] = np.mean(lXs[cid])
            lYs = defaultdict(list)
            lXs = defaultdict(list)
            for sample in tpoint2.samples.values():  # iterate over samples in this time point (usually just 1)
                Ys = sample.trees[i].get_clone_Ys()
                Xs = sample.trees[i].get_clone_Xs()
                for cid in Ys:
                    lYs[cid].append(Ys[cid])
                    lXs[cid].append(Xs[cid])
            for cid in lYs:
                dYs2[i][cid] = np.mean(lYs[cid])
                dXs2[i][cid] = np.mean(lXs[cid])

        for i in range(ntrees):  # for each tree (these are shared with all other samples)
            # time point 1 data - get average clone frequencies
            ave_Ys = dYs[i]

            source_cids = set(
                [cid for cid in dXs[i] if dXs[i][cid] >= eps and dYs[i][cid] > 0])  # clones shared with tp1
            private_cids = set([cid for cid in dYs[i] if
                                dYs[i][cid] < eps and dYs2[i][cid] >= eps and dYs2[i][cid] > 3 * dYs[i][
                                    cid]])  # clones not shared with tp1, new in tp2

            #            source_cids = set([cid for cid in dYs[i] if dYs[i][cid] >= eps and dYs2[i][cid]>=eps])  # clones shared with tp1
            #            private_cids = set([cid for cid in dYs[i] if dYs[i][cid] < eps and dYs2[i][cid]>=eps and dYs2[i][cid]>3*dYs[i][cid]])  # clones not shared with tp1, new in tp2
            #            self.logger("New clones:")
            #            self.logger([(cid, dYs[i][cid], dYs2[i][cid]) for cid in private_cids])
            for sample2 in tpoint2.samples.values():  # iterate over all samples in tp2, set private and shared node frequencies
                tt2 = sample2.trees[i]
                z_shared = sum([tt2.nodes[cid].Y for cid in source_cids])
                z_private = sum([tt2.nodes[cid].Y for cid in private_cids])
                # self.logger(self.name + " private volume in " + sample2.name +" = "+str(z_private))
                for node in tt2.nodes.values():
                    if dir:
                        node.sharedY = 0.0
                        node.privateY = 0.0
                    else:
                        node.lostY = 0.0
                        node.preservedY = 0.0

                for node in tt2.nodes.values():
                    if node.id in source_cids:
                        if dir:
                            node.privateY = 0.0
                            node.sharedY = node.Y / z_shared if z_shared > eps else 0.0
                        else:
                            node.lostY = 0.0
                            node.preservedY = node.Y / z_shared if z_shared > eps else 0.0
                    elif node.id in private_cids:
                        if dir:
                            node.sharedY = 0.0
                            node.privateY = node.Y / z_private if z_private > eps else 0.0
                        else:
                            node.preservedY = 0.0
                            node.lostY = node.Y / z_private if z_private > eps else 0.0
        return True

    def average_transmitted_clones(self, tp1, tp2, eps, node_function, beta=1, **kwargs):
        '''
        Average node_function over the transmitted clones.

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param eps: float
            threshold on cluster size to decide it exists

        :param node_function: function
            Node class function

        :param beta: float
            tree weighting function

        :param kwargs: dict
            named arguments to node_function

        :return: float
        '''

        tpoint1 = self.timePoints[tp1]
        tpoint2 = self.timePoints[tp2]
        ntrees = len([*tpoint1.samples.values()][0].trees)

        def process_tree(i):
            def get_val_and_z(sample):
                tweight = sample.tweights[i]
                ave_Ys = sample.ave_Ys[i]
                transmitted_clones = [clone for clone in sample.trees[i].nodes.values() if ave_Ys[clone.id] > eps]
                z = sum([clone.Y for clone in transmitted_clones])
                val = sum([node_function(clone, **kwargs) * clone.Y for clone in transmitted_clones]) / z
                return val * tweight, z * tweight

            vz_list = [get_val_and_z(sample) for sample in tpoint2.samples.values()]
            v = np.mean([vz[0] for vz in vz_list])
            z = np.mean([vz[1] for vz in vz_list])
            return v, z

        for sample in tpoint2.samples.values():
            sample.tweights = sample.get_tree_weights(beta=beta)
            sample.ave_Ys = defaultdict(dict)

        for i in range(ntrees):
            # time point 1 data - get average clone frequencies
            ave_Ys = defaultdict(float)
            lYs = defaultdict(list)
            for sample in tpoint1.samples.values():
                Ys = sample.trees[i].get_clone_Ys()
                for cid in Ys:
                    lYs[cid].append(Ys[cid])
            for cid in lYs:
                ave_Ys[cid] = np.mean(lYs[cid])
            for sample2 in tpoint2.samples.values():
                sample2.ave_Ys[i] = ave_Ys

        # compute statistics of transmitted clone in time point 2

        vz_list = [process_tree(i) for i in range(ntrees)]
        v = sum([vz[0] for vz in vz_list])
        z = sum([vz[1] for vz in vz_list])
        return v, z

    '''
    Meta functions
    '''

    def meta_tree_function(self, tree_function, **kwargs):
        '''
        Applies the function to all trees in all samples

        :param tree_function: function
            the tree object function to be applied

        :param kwargs: dict
            named parameters of the tree_function
        '''

        for tp in self.timePoints:
            tpoint = self.timePoints[tp]
            for sample in tpoint.samples.values():
                sample.meta_tree_function(tree_function=tree_function, **kwargs)

    '''
    Output
    '''

    def writeJSON(self, outdir, prefix=""):
        '''
        Writes json trees of samples to a folder.

        :param outdir: str
            output folder where json files will be written to
        :return:
        '''
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        for tp in self.timePoints:
            tpoint = self.timePoints[tp]

            for sample in tpoint.samples.values():
                sid = self.name + "|" + sample.name + "|" + self.cohort + "|" + tp + "|" + sample.tissue
                sid = sid.replace("_", "")
                sid = prefix + sid
                jpath = os.path.join(outdir, "tree_" + sid + ".json")
                js = self.toJSON(PFS=True)

                jssamples = {sid: sample.toJSON()}
                js['samples'] = jssamples
                js['id'] = sid
                with open(jpath, 'w') as of:
                    json.dump(js, of, indent=True)

    def __read_phylowgs_clonal_structures(self, jsonpath, ntrees):
        '''
        Import clonal structures for the top scoring "ntrees" trees.

        :param jsonpath: str
            path to the tree summary file, <sname>.summ.json.gz

        :param ntrees: int
            number of highest scoring trees to import

        :return:
            creates a cfit.tree.Tree object and updates the trees attribute.
        '''

        self.logger("read clonal structure " + jsonpath)
        jsondir = os.path.dirname(jsonpath)
        sname = os.path.basename(jsonpath)
        sname = sname.split(".json")[0]
        sname = sname.replace("summ_", "").replace(".summ", "")

        self.logger("Importing tree from " + jsonpath)
        self.logger(sname)
        if not os.path.exists(jsonpath):
            t = Tree([], thepatient=self) #creates one clone tree with all mutations in one clone of frequency 1
            self.trees.append(t)
            js, bestnums = None, None
        else:
            if ".gz" in jsonpath:
                f = gzip.open(jsonpath)
            else:
                f = open(jsonpath)
            line = f.readline()
            f.close()
            #            self.logger("loading line from "+jsonpath)
            js = json.loads(line)
            #            self.logger("done.")
            keys = js['trees'].keys()
            skeys = sorted(keys, key=lambda key: -js['trees'][key]['llh'])

            path03 = os.path.join(jsondir, "muts_" + sname + ".json.gz")
            if not os.path.exists(path03):
                path03 = os.path.join(jsondir, sname + ".muts.json.gz")

            #            self.logger("opening zipped file " + path03)
            with gzip.open(path03, 'rb') as zf:
                js_muts = json.load(zf)
            #            self.logger("done.")
            js_ssms = js_muts['ssms']
            mpath = os.path.join(jsondir, "mutass_" + sname + ".zip")
            if not os.path.exists(mpath):
                mpath = os.path.join(jsondir, sname + ".mutass.zip")

            self.logger("Importing top " + str(ntrees) + " of " + str(len(skeys)) + " trees...", 4)

            sid2mutname = defaultdict(lambda: "")
            pathssm = os.path.join(jsondir, sname + "_ssm.txt")
            if os.path.exists(pathssm):
                # read the _ssm mutation file, input to phylowgs, to recover mutation names
                # and mapping of mutation identifiers. This is not necessary if the phylowgs has
                # annotation of mutations (which is the preferable solution)

                pathssm = os.path.join(jsondir, sname + "_ssm.txt")
                mdata = pd.read_csv(pathssm, sep="\t")
                for (sid, mid) in zip(mdata.id, mdata.gene):
                    sid2mutname[sid] = mid

            bestnums = []
            for i in range(0, min(ntrees, len(skeys))):
                # self.logger("\t" + str(i+1))

                bestnum = skeys[i]
                bestnums.append(bestnum)
                js_tree = js['trees'][bestnum]

                # assignment of mutations to clones
                with zipfile.ZipFile(mpath) as zf:
                    data = zf.read(str(bestnum) + '.json')
                    js_muts = json.loads(data)
                    js_mutass = js_muts['mut_assignments']

                params = [js_tree, js_mutass, js_ssms, sid2mutname]
                t = Tree(params, thepatient=self)
                self.trees.append(t)
        return js, bestnums

    def write_trees(self, odir):
        '''

        :param odir: str

        '''
        if not os.path.exists(odir):
            os.mkdir(odir)

        dtrees = {}
        ntrees = 0
        for tpoint in self.timePoints.values():
            dtrees[tpoint.name] = tpoint.trees(num=-1)
            ntrees = len(dtrees[tpoint.name])
        self.logger(self.name + ", ntrees=" + str(ntrees))
        for i in range(ntrees):
            tinfo = []
            mutinfo = []
            ord = self.orderedTimePoints[:]
            ord.reverse()
            for tp in ord:
                tree = dtrees[tp][i]
                nums = list(tree.nodes.keys())
                nums.sort()
                tinfo.append([tp] + [tree.nodes[nid].X for nid in nums])
            tinfo.append(["n_syn"] + [tree.nodes[nid].node.n_syn() for nid in nums])
            tinfo.append(["n_nsyn"] + [tree.nodes[nid].node.n_nsyn() for nid in nums])
            tinfo.append(["TMB"] + [tree.nodes[nid].node.TMB() for nid in nums])
            tinfo.append(["F"] + [tree.nodes[nid].fitness for nid in nums])
            tinfo.append(["Relative_F"] + [tree.nodes[nid].rfitness for nid in nums])
            tinfo = pd.DataFrame(tinfo)
            #           columns = ["Name"] + [str(nid) for nid in nums]
            tinfo = tinfo.transpose()
            # tinfo.columns =
            tinfo.to_csv(os.path.join(odir, "tree_" + str(i + 1) + "_" + self.name + ".txt"), sep="\t", index=False,
                         header=False)

    def clone_frequency_predictions(self, tp1, tp2, beta=1.0):
        '''

        :param tp1: str
            name of time point 1

        :param tp2: str
            name of time point 2

        :param beta: float
            tree weighting parameter

        :return: pandas.DataFrame, dict
            table with frequency trajectories, prediction and fitness of all clones
            fraction of clones predicted to grow/decline, fraction of the volume of the tumor correctly predicted
        '''

#        tpoint1 = self.timePoints[tp1]
#        tpoint2 = self.timePoints[tp2]


        tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)
        dnodes1 = defaultdict(lambda: defaultdict(list))
        dnodes2 = defaultdict(lambda: defaultdict(list))
        dweight = {}

        for tind, (tree1, tree2, w) in enumerate(tree_pairs):
            new_clones = self.identify_new_clones(tree1, tree2, eps=0.03)
            tree1.soft_remove_nodes(new_clones, include_nested=False, overwrite=False)
            tree2.soft_remove_nodes(new_clones, include_nested=False, overwrite=False)
            dweight[tind] = w
            for node in tree1.nodes.values():
                dnodes1[tind][node.id].append(node)
            for node in tree2.nodes.values():
                dnodes2[tind][node.id].append(node)


#        tpairs = self.initialize_paired_trees(tp1, tp2, 0.03, beta=1., include_nested=False)

#        nodes1 = tpoint1.all_weighted_nodes(beta=beta)  # list of (Node, tree_id, tree_weight)
#        dnodes1 = defaultdict(lambda: defaultdict(list))  # tree_index -> clone id -> tuple (clone, weight)
#
#        dweight = {}
#        for node1, tind1, w1 in nodes1:
#            dnodes1[tind1][node1.id].append(node1)  # over samples in the same time point
#            dweight[tind1] = w1#

#        nodes2 = tpoint2.all_weighted_nodes(beta=beta)
#        dnodes2 = defaultdict(lambda: defaultdict(list))  # tree_index -> clone id -> tuple (clone, weight)
#        for node2, tind2, w2 in nodes2:
#            dnodes2[tind2][node2.id].append(node2)

        clone_data = []
        tindices = list(dnodes1.keys())
        tindices.sort()
        for tind in tindices:
            ncount = 0
            for nid in dnodes1[tind]:
                ncount += 1
                X1 = np.mean([node.X for node in dnodes1[tind][nid]])
                Y1 = np.mean([node.Y for node in dnodes1[tind][nid]])
                X2 = np.mean([node.X for node in dnodes2[tind][nid]])
                Y2 = np.mean([node.Y for node in dnodes2[tind][nid]])

                cY1 = np.mean([node.cY for node in
                               dnodes1[tind][nid]])  # frequency only over clones shared between the time points

                cY2 = np.mean([node.cY for node in
                               dnodes2[tind][nid]])  # frequency only over clones shared between the time points
                Xh = np.mean([node.predictedX for node in dnodes2[tind][nid]])
                Yh = np.mean([node.predictedY for node in dnodes2[tind][nid]])

                fit = np.mean([node.fitness for node in dnodes1[tind][nid]])
                clone_data.append([nid, X1, X2, Y1, Y2, cY1, cY2, Xh, Yh, fit, tind, dweight[tind]])

        clone_data = pd.DataFrame(clone_data)
        clone_data.columns = ["Node", "X1", "X2", "Y1", "Y2", "cY1", "cY2", "Xh", "Yh",
                              "Fitness", "Tree", "Weight"]
        return clone_data

    def write_mutation_predictions(self, tp1, tp2, odir, beta=1.0):
        '''

        :param tp1: str

        :param tp2: str

        :param odir: str

        :param beta: float

        '''

        if not os.path.exists(odir):
            os.mkdir(odir)

#        tpoint1 = self.timePoints[tp1]
#        tpoint2 = self.timePoints[tp2]

#        nodes1 = tpoint1.all_weighted_nodes(beta=beta)  # list of (Node, tree_id, tree_weight)
#        dnodes1 = defaultdict(lambda: defaultdict(list))  # tree_index -> clone id -> list of (clone, weight)

#        dweight = {}
#        for node1, tind1, w1 in nodes1:
#            dnodes1[tind1][node1.id].append(node1)  # over samples in the same time point
#            dweight[tind1] = w1

#        nodes2 = tpoint2.all_weighted_nodes(beta=beta)
#        dnodes2 = defaultdict(lambda: defaultdict(list))  # tree_index -> clone id -> list of (clone, weight)
#        for node2, tind2, w2 in nodes2:
#            dnodes2[tind2][node2.id].append(node2)

        tree_pairs = self.get_tree_pairs(tp1, tp2, beta=beta)
        dnodes1 = defaultdict(lambda: defaultdict(list))
        dnodes2 = defaultdict(lambda: defaultdict(list))
        dweight = {}

        for tind, (tree1, tree2, w) in enumerate(tree_pairs):
            new_clones = self.identify_new_clones(tree1, tree2, eps=0.03)
            tree1.soft_remove_nodes(new_clones, include_nested=False, overwrite=False)
            tree2.soft_remove_nodes(new_clones, include_nested=False, overwrite=False)
            dweight[tind] = w
            for node in tree1.nodes.values():
                dnodes1[tind][node.id].append(node)
            for node in tree2.nodes.values():
                dnodes2[tind][node.id].append(node)



        clone_data = []
        for tind in dnodes1:
            ncount = 0
            for nid in dnodes1[tind]:
                ncount += 1
                X1 = np.mean([node.X for node in dnodes1[tind][nid]])
                Y1 = np.mean([node.Y for node in dnodes1[tind][nid]])
                X2 = np.mean([node.X for node in dnodes2[tind][nid]])
                Y2 = np.mean([node.Y for node in dnodes2[tind][nid]])
                cY1 = np.mean([node.cY for node in dnodes1[tind][nid]])
                cY2 = np.mean([node.cY for node in dnodes2[tind][nid]])
                Xh = np.mean([node.predictedX for node in dnodes2[tind][nid]])
                Yh = np.mean([node.predictedY for node in dnodes2[tind][nid]])
                fit = np.mean([node.fitness for node in dnodes1[tind][nid]])
                rfit = np.mean([node.rfitness for node in dnodes1[tind][nid]])
                rfit2 = np.log(Xh / X1) if Xh > 0 else -10.
                sample_node = dnodes1[tind][nid][0]
                for mutation in sample_node.node.mutations:
                    clone_data.append([mutation.id, int(mutation.is_synonymous()), int(mutation.is_nonsynonymous()),
                                       nid, X1, X2, Y1, Y2, cY1, cY2, Xh, Yh, fit, rfit, rfit2, tind, dweight[tind]])

        clone_data = pd.DataFrame(clone_data)
        clone_data.columns = ["Mutation", "Synonymous", "Nonsynonymous", "Node", "X1", "X2", "Y1", "Y2",
                              "cY1", "cY2", "Xh", "Yh",
                              "Fitness", "Rel_fitness", "Rel_fitness2", "Tree", "Weight"]
        clone_data.to_csv(os.path.join(odir, "Mutation_clone_data_tree_" + self.name + "_" + tp1 + "_" + tp2 + ".txt"),
                          sep="\t", index=False)

    def write_mutation_fitness(self, tp1, tp2, odir, beta=1.0):
        '''

        :param tp1: str
        :param tp2: str
        :param odir: str
        :param beta: float
        :return:
        '''

        if not os.path.exists(odir):
            os.mkdir(odir)

        mids = list(self.mutations.keys())
        is_syn = [int(self.mutations[mid].is_synonymous()) for mid in mids]
        is_nsyn = [int(self.mutations[mid].is_nonsynonymous()) for mid in mids]

        tpoint1 = self.timePoints[tp1]
        mut2CCF = tpoint1.get_mutation_frequencies(beta=beta, exclusive=False)
        mut2fitness, mut2relfitness = tpoint1.get_mutation_fitness(beta=beta)
        freq1 = [mut2CCF[mid] for mid in mids]
        fit1 = [mut2fitness[mid] for mid in mids]
        relfit1 = [mut2relfitness[mid] for mid in mids]

        tpoint2 = self.timePoints[tp2]
        mut2CCF = tpoint2.get_mutation_frequencies(beta=beta, exclusive=False)
        mut2fitness, mut2relfitness = tpoint2.get_mutation_fitness(beta=beta)
        freq2 = [mut2CCF[mid] for mid in mids]
        fit2 = [mut2fitness[mid] for mid in mids]
        relfit2 = [mut2relfitness[mid] for mid in mids]

        data = list(zip(mids, is_syn, is_nsyn, freq1, freq2, relfit1, relfit2, fit1, fit2))
        data = pd.DataFrame(data)
        data.columns = ["Mutation", "Synonymous", "Nonsynonymous",
                        "X1", "X2",
                        "Rel_fitness1", "Rel_fitness2",
                        "Fitness1", "Fitness2"]

        data.to_csv(os.path.join(odir, "mutation_fitness_" + self.name + "_" + tp1 + "_" + tp2 + ".txt"), sep="\t",
                    index=False,
                    header=True)

    def write_neoantigen_fitness(self, tp1, tp2, odir, beta=1.0):
        '''

        :param tp1: str
        :param tp2: str
        :param odir: str
        :param beta: float
        :return:
        '''

        if not os.path.exists(odir):
            os.mkdir(odir)
        tpoint1 = self.timePoints[tp1]

        mids = list(self.mutations.keys())
        mut2CCF1 = tpoint1.get_mutation_frequencies(beta=beta, exclusive=False)
        mut2fitness1, mut2relfitness1 = tpoint1.get_mutation_fitness(beta=beta)
        # freq1 = [mut2CCF[mid] for mid in mids]
        # fit1 = [mut2fitness[mid] for mid in mids]
        # relfit1 = [mut2relfitness[mid] for mid in mids]

        tpoint2 = self.timePoints[tp2]
        mut2CCF2 = tpoint2.get_mutation_frequencies(beta=beta, exclusive=False)
        mut2fitness2, mut2relfitness2 = tpoint2.get_mutation_fitness(beta=beta)
        # freq2 = [mut2CCF[mid] for mid in mids]
        # fit2 = [mut2fitness[mid] for mid in mids]
        # relfit2 = [mut2relfitness[mid] for mid in mids]

        data = []
        for mid in mids:
            neos = self.mutation2neoantigens[mid]
            for neo in neos:
                line = [neo.id, mid, neo.quality, mut2CCF1[mid], mut2CCF2[mid],
                        mut2relfitness1[mid], mut2relfitness2[mid],
                        mut2fitness1[mid], mut2fitness2[mid]]
                data.append(line)
        data = pd.DataFrame(data)
        data.columns = ["Neoantigen", "Mutation", "Quality",
                        "X1", "X2",
                        "Rel_fitness1", "Rel_fitness2",
                        "Fitness1", "Fitness2"]

        data.to_csv(os.path.join(odir, "neoantigen_fitness_" + self.name + "_" + tp1 + "_" + tp2 + ".txt"), sep="\t",
                    index=False,
                    header=True)

    def toJSON(self):
        js = super(PatientLine, self).toJSON()
        js["time_point_order"] = self.orderedTimePoints
        js["time_points"] = [tp.toJSON() for tp in self.timePoints.values()]

        return js
    
    # def print_attributes(self):
    #     print("Attributes of the object:", self.name)
    #     for key, value in self.__dict__.items():
    #         if isinstance(value, (int, float, str, bool)):
    #             print(f"{key}: {value}\n")
    #     print(f"samples: {len(self.samples)}")
                    
    def get_patient_to_sibyl(self):
       
        js = {
            'id': self.name,
            'OS': self.OS,
            'PFS': self.PFS,
            'dead': bool(self.dead),
            'response': self.response,
            'cohort': self.cohort,
            'type': self.type,
            'HLA_genes': self.HLAS,
            # 'topology': self.trees[0].topology()
            
            # 'topologies': [{'id': id, 'topology': tree.topology()} for id, tree in enumerate(self.trees)]
            
            # 'mutations': [mut.toJSON() for mut in self.mutations.values()],
            # 'neoantigens': [neo.toJSON() for neo in self.neoantigens.values()]
        }
        
       
        # print_attributes()
        #js["mutations"] = [mut.toJSON() for mut in self.mutations.values()]
        #js["neoantigens"] = [neo.toJSON() for neo in self.neoantigens.values()]
        # jssamples = {}
        # for sample in self.samples:
        #    jssamples[sample.name] = sample.to_sibyl_JSON()
        # js['samples'] = jssamples
        return js

#         jstrees = [tree.toJSON()  for tree in self.trees]
#         js["trees"] = jstrees
#         js['HLA_genes'] = self.HLAS
#         js["mutations"] = [mut.toJSON() for mut in self.mutations.values()]
#         js["neoantigens"] = [neo.toJSON() for neo in self.neoantigens.values()]
#         return js
#           id: { type: String, required: true },
#   HLA_genes: { type: Array, required: true }, //Patient.HLAS
#   OS: { type: Number, required: true }, //Patient.OS
#   dead: { type: Boolean, required: true }, //Patient.dead
#   response: { type: Boolean, required: true }, // Patient.response
#   samples: { type: Object, required: true }, //PatientLine.samples

    def get_mutations_to_sibyl(self):
        # return { 'patient_id': self.name, 'mutations': [mut.toJSON() for mut in self.mutations.values()] }
        return [mut.toJSON() for mut in self.mutations.values()]
    def get_neoantigens_to_sibyl(self):
        # return { 'patient_id': self.name, 'neoantigens': [mut.toJSON() for mut in self.mutations.values()] }
        return [mut.toJSON() for mut in self.mutations.values()]
    
               
    def get_tree_nodes_to_sibyl(self):
        return [{'id': id, 'nodes': tree.get_tree_nodes()} for id, tree in enumerate(self.trees)]
    
    def get_sample_tree_nodes_to_sibyl(self):
        time_points = []
        for time_point_key, time_point in self.timePoints.items():
            time_points.extend(time_point.get_samples_sibyl(time_point_key))
        # time_points = [time_point.get_samples_sibyl(time_point_key) for time_point_key, time_point in self.timePoints.items()]
        # print(self.name, time_points)
        return time_points

    # def samples_to_sibyl_JSON(self):
    #     # js = {
    #     #     'id': self.name,
    #     #     'OS': self.OS,
    #     #     'PFS': self.PFS,
    #     #     'dead': bool(self.dead),
    #     #     'response': self.response,
    #     #     'cohort': self.cohort,
    #     #     'type': self.type,
    #     #     'HLA_genes': self.HLAS,
    #     #     # 'mutations': [mut.toJSON() for mut in self.mutations.values()],
    #     #     # 'neoantigens': [neo.toJSON() for neo in self.neoantigens.values()]
    #     # }
       
    #     # print_attributes()
    #     #js["mutations"] = [mut.toJSON() for mut in self.mutations.values()]
    #     #js["neoantigens"] = [neo.toJSON() for neo in self.neoantigens.values()]
    #     jssamples = {}
    #     for sample in self.samples:
    #        jssamples[sample.name] = sample.trees_to_sibyl_JSON()
    #     # js['samples'] = jssamples
    #     return jssamples
    
    def sample_tree_nodes_to_sibyl_JSON(self):
        # js = {
        #     'id': self.name,
        #     'OS': self.OS,
        #     'PFS': self.PFS,
        #     'dead': bool(self.dead),
        #     'response': self.response,
        #     'cohort': self.cohort,
        #     'type': self.type,
        #     'HLA_genes': self.HLAS,
        #     # 'mutations': [mut.toJSON() for mut in self.mutations.values()],
        #     # 'neoantigens': [neo.toJSON() for neo in self.neoantigens.values()]
        # }
       
        # print_attributes()
        #js["mutations"] = [mut.toJSON() for mut in self.mutations.values()]
        #js["neoantigens"] = [neo.toJSON() for neo in self.neoantigens.values()]
        jstrees = {}
        for sample in self.samples:
            print('sample, '+sample.name)
            jstrees[sample.name] = sample.trees_to_sibyl_JSON()
        # js['samples'] = jssamples
        return jstrees

