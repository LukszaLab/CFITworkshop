'''
Created on Oct 27, 2015

@author: mluksza
'''

import contextlib
import copy
import io
import json
import math
import os
import random
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from pandas import DataFrame
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
from scipy.stats import mannwhitneyu
from sklearn import metrics

from cfit.CoreObject import CoreObject
from cfit.fitness.FitnessModel import FitnessModel
from cfit.fitness.HLAweights import HLAweights
from cfit.fitness.HWTHLAweights import HWTHLAweights
from cfit.patient.PatientLine import PatientLine
from cfit.patient.Sample import Sample
from cfit.util.Utils import Utils
from cfit.tree.node.SampleNode import SampleNode

import matplotlib.pyplot as plt

@contextlib.contextmanager
def nostdout():
    save_stdout = sys.stdout
    sys.stdout = io.BytesIO()
    yield
    sys.stdout = save_stdout

class AnalysisArgs(CoreObject):
    '''
    If not using command line, create the argument object
    '''
    def __init__(self):
        self.kd_thr = 500
        self.dir = ""
        self.ns="9"
        self.ntrees = 10
        self.tree_format = "phylowgs"
        self.config = os.path.join(self.dir, "config.json")
        self.mapping = os.path.join(self.dir, "mapping.json")
        self.PDAC = False
        self.IO = False
        self.netMHC = "40"
        self.ep_dist_model_name = "all_tcr_all_combos_model"

class Analysis(CoreObject):
    '''
    The class that organizes data import and fitness model computation on a cohort of patients.

        patients: dict: str -> Patient
            all patients in cohort, mapping of patient names to objects of class Patient / PatientLine. A patient can have multiple samples

        patientMask: defaultdict(lambda: 1)
            dictionary to exclude some patients from analysis

        all_patients: list

        samplename2patient: dict
            maps samples to patients

        patient2samples: defaultdict(lambda: defaultdict(lambda: ""))
            maps patients to time point to sample name

        clonal: int
            1 - with clonality, 0 - ignore clonality

        self.npos: str
            1-9H - removes neoantigens with non-hydrophobic residue of wildtype peptide on position 2 or 9
            1-9 - no epitope filtering

        model:  str
            model name, for reporting

        quantile: float
            splitting fraction for survival analysis, default is median

        charAbsTime: float
            characteristic time tau, currently ignored

        criterion: str
            eg. "POPULATION_SIZE", patient classification criterion

        ntrees: int
            number of top trees from a sample over which predictions are averaged

        Rdict: defaultdict(float)

        outdir: str
            default output directory

        patient_names:
            x

        survival: defaultdict(lambda: [0., 1.])
            dictionary mapping survival and status

        kdMT_threshold: float
            threshold on neoantigens


    Global variables - to be moved from Utils class
            MATCH_ALLELES = False
            KS
            AS
            TAUS
            method
            a
            k
            KDnormalize = 1
            HYDROPHOBIC_RESIDUES = "AILMFWYV"
            WEIRD_RESIDUES = "CGP"
            AGGR = "MAX"

            AGGRnum = float("inf")

    '''

    MHC_VERSION = "34"  # default MHC calling algorithm

    def __init__(self):
        self.patients = {}  # all patients in cohort, objects of class Patient. A patient can have multiple samples
        self.patientMask = defaultdict(lambda: 1)  # dictionary to exclude some patients from analysis
        self.all_patients = []
        self.samplename2patient = {}  # maps samples to patients
        self.patient2samples = defaultdict(
            lambda: defaultdict(lambda: ""))  # maps patients to time point to sample name

        self.clonal = 1  # 1 - with clonality, 0 - ignore clonality
        # self.positions = {1, 2, 3, 4, 5, 6, 7, 8, 9}  # neoantigen positions
        self.npos = "1-9"  # 1-9H - removes neoantigens with non-hydrophobic residue of wildtype peptide on position 2 or 9
        # 1-9 - no epitope filtering
        self.model = ""
        self.quantile = 0.5  # splitting fraction for survival analysis, default is median
        self.charAbsTime = 1.0  # characteristic time tau, currently ignored
        self.criterion = "POPULATION_SIZE"  # patient classification

        self.ntrees = 5  # number of top trees from a sample over which predictions are averaged
        self.Rdict = defaultdict(lambda: 0.0)
        self.outdir = ""
        self.patient_names = None
        self.survival = defaultdict(lambda: [0., 1.])
        self.kdMT_threshold = Utils.INF
        self.Fmodel = FitnessModel()

    @property
    def samples(self):
        dsamples = {}
        for patient in self.patients.values():
            lsamples = patient.get_samples(simple=True)
            for sample in lsamples:
                dsamples[sample.name] = sample
        return dsamples

    def initialize_config(self, config, mapping, args=None, dir=None, kd_thr=500, ns=None,
                          tree_format='phylowgs'):
        '''
        Initializes all the data


        :param config: dict
            dictionary with the configuration variables, pointing to the data folders
            'vcf_dir', 'ssm_dir', 'tree_dir'

        :param mapping: dict
            dictionary describing the cohort structure

        :param args: dict
            argparse object with arguments.
            args.dir :str
            args.kd_thr:float
            args.ns: str
            args.tree_format: str
        :param dir: str
            path to the data directory

        :param kd_thr: float or None
            threshold on kd for neoantigens

        :param ns: list
            lengths of neoantigens

        :param tree_format: str
            tree format: phylowgs, pairtree...
        '''

        verb_level = 4
        vcfdir = config['vcf_dir']
        treedir = config['tree_dir']
        if "ssm_dir" in config:
            ssmdir = config['ssm_dir']
        else:
            ssmdir = None

        if args is not None:
            dir = args.dir

            if "-" in args.ns:
                ns = [int(x) for x in  args.ns.split("-")[:2]]
                ns = list(range(ns[0], ns[1]+1))
            else:
                ns = [int(x) for x in (args.ns).split(",")]
            kd_thr = args.kd_thr
            tree_format = args.tree_format

        # Ininitializes the cohort data
        neodir = os.path.join(dir, config['neo_dir'], self.MHC_VERSION)
        if not os.path.exists(neodir):
            neodir = os.path.join(dir, config['neo_dir'])

        hlafile = os.path.join(dir, "HLA", "HLA_calls.txt")
        dhlas = defaultdict(list)
        if os.path.exists(hlafile):
            hlatab = pd.read_csv(hlafile, sep="\t", header=None)
            for row in hlatab.itertuples():
                dhlas[row[1]] = row[2].split(",")

        for patdata in mapping:
            '''
            Creates list of patients, and their samples
            '''
            name = patdata['name']
            if 'pname' in patdata:
                pname = patdata['pname']
            else:
                pname = name
            pat_tree_dir = os.path.join(dir, treedir, name)
            pat_vcf_dir = os.path.join(dir, vcfdir, pname)

            if not os.path.exists(
                    pat_vcf_dir):  # if there is single vcf file for the patient, it doesn't require a directory
                pat_vcf_dir = os.path.join(dir, vcfdir)

            pat_type = patdata['type']
            patient = PatientLine([name, None, None, None])
            patient.set_HLAs(dhlas[name])
            self.logger("Setting hlas for "+str(name)+" "+str(dhlas[name]))
            patient.cohort = patdata['cohort']
            patient.pname = pname
            patient.type = pat_type
            try:
                patient.OS = patdata['OS']
            except KeyError:
                patient.OS = 0
            try:
                patient.PFS = patdata['PFS']
            except KeyError:
                patient.PFS = 0

            patient.survival = patient.OS
            patient.dead = patdata["dead"]
            if "PFS_status" in patdata:
                patient.PFS_status = patdata["PFS_status"]
            else:
                patient.PFS_status = patdata["dead"]
            if ssmdir is not None:
                with open(os.path.join(dir, ssmdir, name, 'params.json'), 'r') as fp:
                    tree_params = json.load(fp)
                    sample_order = tree_params["samples"]
            else:
                sample_order = [el[2] for el in patdata["samples"]]

            patient.create_samples(patdata, sample_order, pat_vcf_dir, pat_tree_dir, self.ntrees,
                                   tree_format=tree_format)

            for sample in patient.samples:
                sample.HLAS = patient.HLAS

            self.logger("Patient " + patient.name + " has " + str(len(patient.samples)) + " samples.", verb_level)
            self.patients[name] = patient

            # Import neoantigens
            neofile = os.path.join(neodir, "neoantigens_" + pname + ".txt")
            if os.path.exists(neofile):
                patient.add_neoantigens(neofile, kd_thr=kd_thr, ns=ns)

            neofile = os.path.join(neodir, "neoantigens_other_" + pname + ".txt")
            if os.path.exists(neofile):
                patient.add_frame_shift_neoantigens(neofile, kd_thr=kd_thr, ns=ns)
            patient.distribute_neoantigens_to_clones()
            patient.set_exclusive_mutations()

    def set_MHC_version(self, version):
        '''
        Set the netmhc version

        :param version: str
            version of netMHC, 3.4 (34), 4.0 (40)
        '''
        if version in ["3.4", "34"]:
            self.MHC_VERSION = "34"
        elif version in ["4.0", "40"]:
            self.MHC_VERSION = "40"
        elif version in ["pan4.1", "pan41"]:
            self.MHC_VERSION = "pan41"
        else:
            raise Exception("netMHC version not handled")


    def set_fitness_model_component(self, fcomponent, name, weight):
        '''

        :param fcomponent: cfit.fitness.CloneFitness

        :param name: str

        :param weight: float

        '''

        self.Fmodel.add_component(fcomponent, name, weight)


    def reset_fitness_model_components(self):
        '''
        Removes fitness model components
        '''
        self.Fmodel.reset_fitness_model_components()


    def set_neantigen_quality_model(self, Qmodel, **kwargs):
        '''
        Initializes the model to compute neoantigen quality. Fills in precomputed neoantigen qattributes, if neeeded.

        :param NQclass: cfit.fitness.NeoantigenQuality

        :param kwargs: dict
            arguments to the class constructor
        '''

        Qmodel.initialize_neoantigens(self, **kwargs)
        self.Qmodel = Qmodel

    ####################################################################################################################
    # Methods for auxiliary access to objects
    ####################################################################################################################


    def get_patients(self):
        '''

        :return: list
            list of patient objects
        '''
        return list(self.patients.values())


    def get_samples(self):
        '''
        List of all samples in the cohort

        :return: list
            list of sample objects

        '''
        samples = [sample for patient in self.patients.values() for sample in patient.get_samples(simple=True)]

        return samples


    def get_neoantigens(self):
        '''
        List of all neoantigens in the patients of the cohort

        :return: list
        '''
        neos = [neo for patient in self.patients.values() for neo in patient.neoantigens.values()]
        #        neos = [neo for patient in self.patients.values() for sample in patient.samples for neo in sample.neoantigens]
        return neos


    def get_neoantigen_sample_pairs(self):
        '''
        List of all neoantigens-sample pairs in the samples of the cohort

        :return: list
        '''
        neosamples = [(neo, sample) for patient in self.patients.values() for sample in patient.samples for neo in
                      sample.neoantigens]
        return neosamples


    def get_nodes(self):
        '''
        List of all clones and samples of the cohort

        :return: list
            [(SampleNode, SampleTree, Sample),...]
        '''
        nodes = [(node, tree, sample) for patient in self.patients.values() for sample in patient.samples for tree in
                 sample.trees for
                 node in tree.nodes.values()]
        onodes = [(node, sample.oneTree, sample) for patient in self.patients.values() for sample in patient.samples for node in
                  sample.oneTree.nodes.values()]
        return nodes + onodes


    def get_mutations(self):
        '''
        List of all mutations in the samples of the cohort

        :return: list
        '''
        lsamples = self.get_samples()
        mutations = [mut for sample in lsamples for mut in sample.mutations.values()]
        return mutations


    def get_trees(self):
        '''
        List of all trees in the samples of the cohort

        :return: list
        '''

        lsamples = self.get_samples()
        trees = [tree for sample in lsamples for tree in sample.trees]
        return trees

    ####################################################################################################################
    # Methods for fitness
    ####################################################################################################################


    def compute_neoantigen_qualities(self, **kwargs):
        '''
        Computes neoantigen qualities of all neoantigens, setting up neo.quality attributes and sample.neoantigenQualities

        :param args: list

        :param kwargs: dict
        '''
        #        neos = self.get_neoantigens()
        neo_pats = [(neo, patient) for patient in self.patients.values() for neo in patient.neoantigens.values()]
        for (neo, pat) in neo_pats:
            quality = self.Qmodel.compute_quality(neo, **kwargs)
            for sample in pat.get_samples(simple=True):
                sample.neoantigenQualities[neo.id] = quality
                for tree in sample.trees:
                    tree.neoantigenQualities[neo.id] = quality


    def compute_neoantigen_sample_qualities(self, components=False, **kwargs):
        '''
        Computes neoantigen qualities of all neoantigens, in the context of each sample (may depend  on gene expression).
        sets sample.neoantigenQualitities attribute

        :param components: bool

        :param kwargs: dict
        '''
        neosamples = self.get_neoantigen_sample_pairs()
        for (neo, sample) in neosamples:
            quality = self.Qmodel.compute_neoantigen_sample_quality(neo, sample, **kwargs)
            sample.neoantigenQualities[neo.id] = quality
            sample.neoantigenQualities[neo.id] = max(0, quality)
            for tree in sample.trees:
                tree.neoantigenQualities[neo.id] = max(0, quality)
            if components:
                sample.neoantigenQualityComponents[neo.id] = self.Qmodel.compute_neoantigen_sample_quality_components(
                    neo, sample, **kwargs)


    def compute_neoantigen_sample_qualities_function(self, quality_function=lambda n, s: -n.kD, **kwargs):
        '''
        Computes neoantigen qualities of all neoantigens, in the context of each sample (may depend  on gene expression).
        sets sample.neoantigenQualitities attribute

        :param quality_function: function

        :param kwargs: dict
        '''
        neosamples = self.get_neoantigen_sample_pairs()
        for (neo, sample) in neosamples:
            quality = quality_function(neo, sample, **kwargs)
            sample.neoantigenQualities[neo.id] = quality
            sample.neoantigenQualities[neo.id] = max(0, quality)
            for tree in sample.trees:
                tree.neoantigenQualities[neo.id] = max(0, quality)


    def compute_node_fitness(self, recompute_components=True):
        '''
        Computes fitness of all nodes using the fitness model set in the self.Fmodel attribute.
        After setting node.fitness attributes, sets the relative fitness in each tree, node.rfitness.

        :param recompute_components: bool
            whether to recompute the components (set to False if only weights are changed)
        '''

        # absolute fitness
        nodes = self.get_nodes()
        for (node, sampleTree, sample) in nodes:
            self.Fmodel.compute_node_fitness(node, sampleTree, sample, recompute_components=recompute_components)

        # relative fitness
        trees = self.get_trees()
        for tree in trees:
            avef = tree.average_over_nodes(lambda node: node.fitness)
            for node in tree.nodes.values():
                node.rfitness = node.fitness - avef


    def standardize_node_fitness(self):
        '''
        Standardizes precomputed fitness of all nodes (as store in node.fitness):
        Computes the standard deviation of fitness over all nodes

        After setting node.fitness attributes, sets the relative fitness in each tree, node.rfitness.

        '''
        # absolute fitness
        nodes = self.get_nodes()
        self.Fmodel.normalize_fitness(nodes)
        # relative fitness
        trees = self.get_trees()
        for tree in trees:
            avef = tree.average_over_nodes(lambda node: node.fitness)
            for node in tree.nodes.values():
                node.rfitness = node.fitness - avef


    def set_time_pairs_by_prefix(self, pref1, pref2):
        '''

        :param pref1: str
            name of the time point, eg. "primary", "pre", "prim"
        :param pref2: str
            name of the second time point, eg. "met", "post"

        '''
        for patient in self.patients.values():
            patient.set_time_pairs_by_prefix(pref1, pref2)

    def init_npos(self, npos):
        self.npos = npos
        if npos == "1-9H":
            self.HLAW = HWTHLAweights(npos)
        else:
            self.HLAW = HLAweights(npos)


    def write_mutation_frequencies(self, outdir, beta=1.0, exclusive=False, by_sample=False):
        '''

        :param outdir: str
            path to the output directory

        :param beta: float
            tree weighting parameter

        :param exclusive: bool
            whether to report Y or X.

        :param by_sample: bool
            whether to report frequencies for each sample, or averaged by time point

        :return:
        '''
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        pref = "Y" if exclusive else "X"
        pref = pref + "_sample" if by_sample else pref
        for patient in self.patients.values():
            mccfs = patient.get_mutation_frequencies(beta=beta, exclusive=exclusive, by_sample=by_sample)
            mccfs.to_csv(os.path.join(outdir, pref + "_" + patient.name + ".txt"), sep="\t", index=False)


    def mutated_genes_CCF(self):
        '''
        Creates a table with frequency data on mutated genes
        :return: pd.DataFrame
        '''
        all_samples = []
        for pat in self.patients.values():
            ltps = pat.get_samples()
            for [tp, samples] in ltps:
                all_samples += [(pat.name, tp, sample) for sample in samples]
        genes = []
        for el in all_samples:
            (pname, tpname, sample) = el
            sgenes = list(set([mut.gene for mut in sample.mutations.values()]))
            genes += sgenes
        genes = list(set(genes))
        genes.sort()
        gene_tab = []
        columns = ["gene"] + [pname + "_" + tpname + "_" + sample.name for pname, tpname, sample in all_samples]
        for gene in genes:
            row = [sample.get_gene_frequency(gene) for _, _, sample in all_samples]
            gene_tab.append([gene] + row)
        gene_tab = pd.DataFrame(gene_tab)
        gene_tab.columns = columns
        ave = gene_tab.mean(axis=1)
        gene_tab["ave"] = ave
        gene_tab = gene_tab.sort_values("ave", ascending=False).iloc[:, :-1]
        return gene_tab

    def write_mutations(self, outdir, exclusive=True, clonal=True):
        '''
        Write mutations of patient samples to files.
        :param outdir: str
            path to the output directory

        :param exclusive: bool

        :param clonal: bool

        '''

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        for pname in self.patients:
            patient = self.patients[pname]
            pdir = os.path.join(outdir, pname)
            patient.write_mutations(pdir, exclusive=exclusive,
                                    clonal=clonal)

    def write_CCFs(self, outdir):
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        pnames = list(self.patients.keys())
        pnames.sort()
        for pname in pnames:
            patient = self.patients[pname]
            for sample in patient.samples:
                opath = os.path.join(outdir, "CCF_" + pname + "_" + sample.name + ".txt")
                sample.write_CCFs(opath)


    def write_sample_statistics(self, outdir, xthr=0.0):
        '''

        :param outdir: str
            path to the output directory

        :param xthr: float

        '''
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        stats = []
        pnames = list(self.patients.keys())
        pnames.sort()
        for pname in pnames:
            patient = self.patients[pname]
            for sample in patient.samples:
                line = [sample.name,
                        sample.get_effective_number_of_clones(),
                        sample.get_number_of_clones(),
                        len(sample.mutations),
                        -sample.TMB(),
                        len(sample.neoantigens),
                        -sample.neoantigen_load(),
                        sample.q,
                        sample.neo2mut_ratio(),
                        sample.get_average_tree_height(synonymous=True, xthr=xthr)]
                stats.append(line)
        stats = pd.DataFrame(stats)
        stats.columns = ["sample", "effective_number_of_clones",
                         "number_of_clones", "mutational_load",
                         "average_mutational_load", "neoantigen_load",
                         "average_neoantigen_load", "average_fitness",
                         "neo_to_mut_ratio",
                         "mol_clock"]
        stats.to_csv(os.path.join(outdir, "sample_statistics.txt"), sep="\t", index=False)


    def write_sibyl(self, outdir, prefix=""):
        '''
        Write Sibyl app output to json files

        :param outdir: str

        :param prefix: str

        '''

        self.logger("writing Sibyl...")
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        for pname in self.patients:
            patient = self.patients[pname]
            pdir = os.path.join(outdir, pname)
            patient.writeJSON(pdir, prefix=prefix)
        self.logger("Done with writing Sibyl...")


    def write_neoantigen_fitness(self, outdir, exclusive=True, longitudinal=False,
                                 clonal=True):
        '''

        :param outdir: str
            output directory

        :param exclusive: bool
            whether to include neoantigens only in the clone their originate

        :param longitudinal: bool

        '''
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if longitudinal:
            for pname in self.patients:
                patient = self.patients[pname]
                for tp in patient.orderedTimePoints:
                    tpoint = patient.timePoints[tp]
                    include_snumber = (len(tpoint.samples) > 1)
                    for i, sname in enumerate(tpoint.samples):
                        sample = tpoint.samples[sname]
                        opath = os.path.join(outdir, "nf_" + sname + ".txt")
                        sname0 = sample.name
                        sample.name = pname + "_" + tp
                        if include_snumber:
                            sample.name = sample.name + "_" + sname
                        sample.write_neoantigens(opath, exclusive=exclusive, clonal=clonal)
                        opath_indel = os.path.join(outdir, "nf_indel_" + sname + ".txt")
                        sample.write_indel_neoantigens(opath_indel, exclusive=exclusive, clonal=clonal)

                        sample.name = sname0
        else:
            for sname in self.samples:
                sample = self.samples[sname]
                opath = os.path.join(outdir, "nf_" + sname + ".txt")
                sample.write_neoantigens(opath, exclusive=exclusive, clonal=clonal)
                opath_indel = os.path.join(outdir, "nf_indel_" + sname + ".txt")
                sample.write_indel_neoantigens(opath_indel, exclusive=exclusive, clonal=clonal)


    def sample_statistics(self, kd_thr=500):
        '''

        :param kdthr: float
            threshold on neoantigen Kd
        '''
        kdthr = 1e10 if kd_thr is None else kd_thr
        lstats = []
        for sname in self.samples:
            sample = self.samples[sname]

            n_nsyn = sample.TMB_nsyn() # absolute number of missense mutation
            n_syn = sample.TMB_syn() # absolute number of synonymous mutations
            n_MHC = sample.TMB_MHC() # absolute number of mutations with a MHC-presented neoantigens

            n_nsyn_eff = sample.effective_TMB_nsyn(beta=1) # number of missense mutations, averaged over clones
            n_syn_eff = sample.effective_TMB_syn(beta=1) # number of synonymous mutations, averaged over clones
            n_MHC_eff = sample.effective_TMB_MHC(beta=1) # number of mutations with a MHC-presented neoantigens,
                                                         # averaged over clones

            tmb = sample.TMB() # absolute number of all mutations
            tnb = sample.presentation_score(kd0=kdthr, strict=True, filter_on_Taf=False)
            # absolute number of neoantigens
            tmb_eff = sample.effective_TMB() # number of all mutations, averaged over clones
            tnb_eff = sample.effective_presentation_score(kd0=kdthr, strict=True, beta=1)
            # number of neoantigens, averaged over clones

            entr = sample.entropy(beta=1) # tumor heterogeneity
            nclones = sample.get_number_of_clones()
            neff_clones = sample.get_effective_number_of_clones(beta=1)

            stats = [sname] + [tmb, n_nsyn, n_syn, n_MHC, tnb] + [tmb_eff, n_nsyn_eff, n_syn_eff, n_MHC_eff, tnb_eff]
            stats += [nclones, neff_clones, entr]
            lstats.append(stats)
            
        lstats = pd.DataFrame(lstats)
            #dstats[sname] = stats
        colnames = ["Sample", "TMB", "TMB_nsyn", "TMB_syn", "TMB_MHC", "TNB"] + \
                   ["TMB_eff", "TMB_nsyn_eff", "TMB_syn_eff", "TMB_MHC_eff", "TNB_eff"] + \
                   ["nclones", "neff_clones", "entropy"]
        lstats.columns = colnames
        return lstats

    def fitness_sample_statistics(self):
        '''

        :param kdthr: float
            threshold on neoantigen Kd
        '''

        lstats = []
        for sname in self.samples:
            sample = self.samples[sname]
            ntau = sample.q
            ave_fitness = sample.average_over_clones(SampleNode.get_fitness, beta=1, clonal=self.clonal)

            stats = [sname] + [ntau, ave_fitness]
            lstats.append(stats)
        lstats = pd.DataFrame(lstats)
        #dstats[sname] = stats
        colnames = ["Sample", "ntau", "ave_fitness"]
        lstats.columns = colnames
        return lstats


    def reset_patient_masking(self):
        for pname in self.all_patients:
            self.patientMask[pname] = 1


    def mask_patient(self, pname):
        '''
        Set to ignore the given patient
        :param pname: str
            patient to mask
        :return:
        '''
        self.patientMask[pname] = 0


    def logrank_Pvalue(self, merge_samples=False, OS=True, quantile=0.5, thrval=None):
        '''

        Computes logrank score pvalue for separation of patients based on the precomputed
        patients scores.


        :param merge_samples: bool
            will average ntau of patients with names: patientname_<num>.
            There can only be one space in the name (backward compatibility for ICGC)

        :param OS: bool

        :param quantile: float

        :param thrval: float

        :return: list
            [p-value, score, quantile value of the score, sum(groups), n-sum(groups), n]
        '''

        def _check(df):
            '''
            Check if the survival curves are in the right "order"
            - if the predicted responders have better prognosis than
            the predicted non-responders

            :param df: dict
                dictionary with keys:

            :return: int
            '''

            dat = list(zip(df['durations'], df['events'], df['groups']))
            times = defaultdict(list)
            events = defaultdict(list)
            kmf = {}
            for (t, e, g) in dat:
                times[g].append(t)
                events[g].append(e)
            destimates = {}
            for g in times:
                km = KaplanMeierFitter()
                with nostdout():
                    km.fit(times[g], events[g], label=str(g))
                kmf[g] = km
                destimates[g] = getattr(km, "survival_function_")
            new_index = np.concatenate((destimates[1].index, destimates[0].index))
            new_index = np.unique(new_index)
            dvals = {}
            for g in destimates:
                destimates[g] = destimates[g].reindex(new_index, method='ffill')
                dvals[g] = [line[0] for line in destimates[g].itertuples(index=False)]
            # we want the predicted responder curve to be above the predicted non-responder
            if thrval is None:
                med = np.quantile([el[0] - el[1] for el in zip(dvals[1], dvals[0])], quantile)
            else:
                med = thrval
            #            med = np.median([el[0] - el[1] for el in zip(dvals[1], dvals[0])])
            return (1 if med >= 0 else -1)

        pats = list(self.patients.values())
        if OS:
            pats = [pat for pat in pats if not math.isnan(pat.OS)]
        else:
            pats = [pat for pat in pats if not math.isnan(pat.PFS)]
        #pats = [pat for pat in pats]
        qs = [pat.q for pat in pats]
        durations = [pat.OS if OS else pat.PFS for pat in pats]
        events = [pat.dead for pat in pats]
        if thrval is None:
            medq = np.quantile(qs, quantile)
        else:
            medq = thrval
        # group 1: predicted responders
        # group 0: no predicted response
        groups = [1 if pat.q <= medq else 0 for pat in pats]
        n = len(pats)

        if merge_samples:
            stem2patients = defaultdict(list)
            for pat in pats:
                stem = pat.name.split("_")
                if len(stem) == 2:
                    stem = "_".join(stem[-1:])
                else:
                    stem = pat.name
                stem2patients[stem].append(pat)
            stems = list(stem2patients.keys())
            qs = [np.mean([pat.q for pat in stem2patients[stem]]) for stem in stems]
            durations = [np.mean([pat.OS if OS else pat.PFS for pat in stem2patients[stem]]) for stem in stems]
            events = [int(np.mean([pat.OS if OS else pat.PFS for pat in stem2patients[stem]])) for stem in stems]
            if thrval is None:
                medq = np.quantile(qs, quantile)
            else:
                medq = thrval
            # medq = np.median(qs)
            groups = [1 if q <= medq else 0 for q in qs]
            n = len(stems)

        df = pd.DataFrame({
            'durations': durations,
            'groups': groups,  # could be strings too
            'events': events, })
        try:
            check = _check(df)
        except:
            check = 1

        results = multivariate_logrank_test(df['durations'], df['groups'], df['events'])
        rpval = results.p_value
        rtest = results.test_statistic * check
        return [rpval, rtest, medq, sum(groups), n - sum(groups), n]


    def survival_AUC(self, OS=True, quants=np.arange(0.1, 1., 0.1)):
        '''

        Computes logrank score pvalue for separation of patients based on the precomputed
        patients fitness scores.

        :param OS: bool
            whether to use OS or PFS in survival analysis

        :param quants: list
            list of quantiles at which to test the survival separation

        :return: list, list
            list of log-rank scores at quantiles, and list of log10 p-values
        '''

        def _check(df):
            '''
            Check if the survival curves are in the right "order"
            - if the predicted responders have better prognosis than
            the predicted non-responders

            :param df: dict
                dictionary with keys:

            :return: int
            '''

            dat = list(zip(df['durations'], df['events'], df['groups']))
            times = defaultdict(list)
            events = defaultdict(list)
            kmf = {}
            for (t, e, g) in dat:
                times[g].append(t)
                events[g].append(e)
            destimates = {}

            for g in times:
                km = KaplanMeierFitter()
                with nostdout():
                    km.fit(times[g], events[g], label=str(g))
                kmf[g] = km
                destimates[g] = getattr(km, "survival_function_")
            try:
                new_index = np.concatenate((destimates[1].index, destimates[0].index))
                new_index = np.unique(new_index)
                dvals = {}
                for g in destimates:
                    destimates[g] = destimates[g].reindex(new_index, method='ffill')
                    dvals[g] = [line[0] for line in destimates[g].itertuples(index=False)]
                diff = [el[0] - el[1] for el in zip(dvals[1], dvals[0])]
                med = np.median(diff)
            except:
                return 0
            # we want the predicted responder curve to be above the predicted non-responder
            return (1 if med >= 0 else -1)

        pats = list(self.patients.values())
        pats.sort(key=lambda pat: pat.q)
        qvals = [pat.q for pat in pats]
        durations = [pat.OS if OS else pat.PFS for pat in pats]
        events = [pat.dead for pat in pats]
        rtests = []
        rpvals = []
        n = len(pats)

        for quant in quants:
            thrval = np.quantile(qvals, quant)
            groups = [1 if q <= thrval else 0 for q in qvals]
            df = pd.DataFrame({
                'durations': durations,
                'groups': groups,  # could be strings too
                'events': events, })
            check = _check(df)
            check = max(0., check)
            results = multivariate_logrank_test(df['durations'], df['groups'], df['events'])
            rtest = results.test_statistic * check
            rpval = np.log10(max(1e-10, results.p_value)) if check == 1 else 0.
            rtests.append(rtest)
            rpvals.append(rpval)

        return rtests, rpvals


    ##############################################
    # Validation
    ##############################################

    def classify_survival(self, beta=1., tau=None,
                          scoreFunction=None, patient_names=None,
                          outdir=None, merge_samples=False,
                          criterion='POPULATION_SIZE',
                          compute_score=True, OS=True, PFS=False, quantile=0.5, thrval=None,
                          **kwargs):
        '''
        Perform survival analysis on the cohort.

        :param beta: float
            tree weighting parameter

        :param tau: float
            time parameter

        :param scoreFunction: function

        :param patient_names: list

        :param outdir: str
            output directory path (optional)

        :param quantile: float
            survival curves separation threshold, defaults to median (0.5)

        :param thrval: float
            separation threshold for patient.q values

        :param kwargs: dict

        :return: list
            [logrankscore, pvalue, n]
        '''

        if tau is None:
            tau = 1.  # self.charAbsTime
        if patient_names is None:
            patient_names = self.patients.keys()
        self.monthLimit = 1000

        self.criterion = criterion

        if compute_score:
            if scoreFunction is None:
                [patient.compute_score(self.criterion, self.clonal, beta, tau) for patient in self.patients.values()]
            else:
                for patient in self.patients.values():
                    patient.meta_beta_function(function=scoreFunction, beta=beta, **kwargs)
        AUC = []
        n1 = -1
        n2 = -1
        if OS:
            resOS = self.logrank_Pvalue(merge_samples=merge_samples, OS=True, quantile=quantile, thrval=thrval)
#            thrval = resOS[2]
            logrankscore = resOS[1]
            pvalue = resOS[0]
            n1 = resOS[3]
            n2 = resOS[4]
            # [rpval, rtest, medq, sum(groups), n - sum(groups), n]
            #            AUC = self.survival_AUC(OS=True, quants=AUC_quantiles)
            AUC = ([0], [0])
        else:
#            thrval = 0
            logrankscore = 0
            pvalue = 1

        if PFS:
            resPFS = self.logrank_Pvalue(merge_samples=merge_samples, OS=False, quantile=quantile, thrval=thrval)
            logrankscorePFS = resPFS[1]
            pvaluePFS = resPFS[0]
            n1 = resPFS[3]
            n2 = resPFS[4]

            #            AUC = self.survival_AUC(OS=False, quants=AUC_quantiles)
            AUC = ([0], [0])
        else:
            logrankscorePFS = 0
            pvaluePFS = 1

        return {"score_OS": logrankscore, "pval_OS": pvalue,
                "score_PFS": logrankscorePFS, "pval_PFS": pvaluePFS,
                #"thrval": thrval,
                "AUC": AUC, "n1": n1, "n2": n2}


    def classify_AUC(self, patient_class,
                     beta=1., tau=None, scoreFunction=None,
                     criterion='POPULATION_SIZE', compute_score=True, ofile=None,
                     **kwargs):
        '''
        Perform classification and AUC

        :param patient_class: dict

            classification of patients str->int (1 -> negative, 2 -> positive (eg.responders, LTS))

        :param beta: float
            tree weighting parameter

        :param tau: float
            time parameter

        :param scoreFunction: function

        :param patient_names: list

        :param kwargs: dict

        :return: list

        '''

        if tau is None:
            tau = 1.  # self.charAbsTime

        self.criterion = criterion

        if compute_score:
            if scoreFunction is None:
                [patient.compute_score(self.criterion, self.clonal, beta, tau) for patient in
                 self.patients.values()]
            else:
                for patient in self.patients.values():
                    patient.meta_beta_function(function=scoreFunction, beta=beta, **kwargs)

        patients = list(self.patients.values())

        labels = list(set(patient_class.values()))
        labels.sort()
        label1 = labels[0]
        label2 = labels[1]
        x1 = [pat.q for pat in patients if patient_class[pat.name] == label1]
        x2 = [pat.q for pat in patients if patient_class[pat.name] == label2]
        mwpval = mannwhitneyu(x1, x2, alternative='two-sided').pvalue

        y = np.array([patient_class[pat.name] for pat in patients])
        pred = np.array([-patient.q for patient in patients])

        fpr1, tpr1, thresholds1 = metrics.roc_curve(y, pred, pos_label=label1)
        fpr2, tpr2, thresholds2 = metrics.roc_curve(y, pred, pos_label=label2)

        auc1 = metrics.auc(fpr1, tpr1)
        auc2 = metrics.auc(fpr2, tpr2)

        if auc1 > auc2:
            fpr, tpr, thresholds, auc = fpr1, tpr1, thresholds1, auc1
        else:
            fpr, tpr, thresholds, auc = fpr2, tpr2, thresholds2, auc2

        if ofile is not None:
            plt.plot(fpr,tpr)
            plt.ylabel('True Positive Rate')
            plt.xlabel('False Positive Rate')
            plt.title("AUC="+str(round(auc,2)))
            plt.savefig(ofile)
            plt.close()
        return {"FPR": fpr, "TPR": tpr,
                "thresholds": thresholds, "AUC": auc,
                "MW_pval": mwpval, "mean_NR": np.mean(x1), "mean_R": np.mean(x2)}


    def randomize_labels(self):
        '''
        Shuffles patient labels for significance analysis
        '''
        patient_labels = {}
        patients = self.patients.values()
        labels = []
        for pat in patients:
            patient_labels[pat.name] = (pat.OS, pat.PFS, pat.dead, pat.cohort)
            labels.append([pat.OS, pat.PFS, pat.dead, pat.cohort])
        random.shuffle(labels)

        for pat, lab in zip(patients, labels):
            pat.OS = lab[0]
            pat.PFS = lab[1]
            pat.dead = lab[2]
            pat.cohort = lab[3]

        self.real_labels = patient_labels


    def reset_labels(self):
        '''
        Resets patient labels after randomization
        '''
        for pat in self.patients.values():
            pat.OS = self.real_labels[pat.name][0]
            pat.PFS = self.real_labels[pat.name][1]
            pat.dead = self.real_labels[pat.name][2]
            pat.cohort = self.real_labels[pat.name][3]

    ##############################################
    # Model fitting
    ##############################################


    def init_tree_pairs(self, tp_pref1, tp_pref2, eps=0.03, beta=1, include_nested=True):
        '''

        :param tp_pref1: str
            prefix of the name of the first time point

        :param tp_pref2: str
            prefix of the name of the second time point

        :param eps: float
            frequency threshold for new clones

        :param beta: float
            weight parameter for trees

        :param include_nested: bool
            whether to include volume of new clones into the clones present in the first time point

        '''
        for pat in self.patients.values():
            self.logger("pat " + pat.name)
            pat.tree_pairs = {}
            tppairs = pat.get_time_pairs_by_prefix(tp_pref1, tp_pref2)
            pat.time_point_pairs = tppairs
            for tp1, tp2 in tppairs:
                self.logger("pat " + pat.name + " " + tp1 + " " + tp2)
                tpairs = pat.initialize_paired_trees(tp1, tp2, eps=eps, beta=beta, include_nested=include_nested)
                pat.tree_pairs[(tp1, tp2)] = tpairs

        tpatients = [(pat, tp1, tp2, pat.name + "_" + tp1 + "_" + tp2) for pat in self.patients.values() for
                     (tp1, tp2) in pat.tree_pairs.keys()]
        self.tpatients = tpatients

    ##############################################
    # Plots
    ##############################################

    def plot_ROC(self, patient_class, ofile):
        '''
        :param patient_class: dict

        :param ofile: str

        '''
        patients = list(self.patients.values())

        y = np.array([patient_class[pat.name] for pat in patients])
        pred = np.array([-patient.q for patient in patients])
        fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=2)


        auc = metrics.auc(fpr, tpr)
        plt.plot(fpr,tpr)
        plt.ylabel('True Positive Rate')
        plt.xlabel('False Positive Rate')
        plt.title("AUC="+str(round(auc,2)))
        plt.savefig(ofile)
        plt.close()


    def plot_survival(self, outfile="", quantile=0.5, prism=False, pval=None, OS=True, thrval=None, show=False):
        '''

        :param outfile: str
            the pdf file to output the plot to

        :param quantile: float
            split fraction, defaults to median, belongs to [0,1]
        :param prism: bool
            whether to do prism output

        :param pval: float

        :param OS: bool

        :thrval: float

        '''

        def deNone(x):
            if x is None:
                return ""
            return x

        import matplotlib.pyplot as plt
        #        if os.environ.get('DISPLAY', '') == '':
        #            self.logger('no display found. Using non-interactive Agg backend')
        #            matplotlib.use('Agg')

        def assign_sample(sample, medfit):
            if sample.q <= medfit:
                return "negative"
            else:
                return "positive"

        patients = [patient for patient in self.patients.values() if self.patientMask[patient.name]]
        if OS:
            patients = [pat for pat in patients if not math.isnan(pat.OS)]
        else:
            patients = [pat for pat in patients if not math.isnan(pat.PFS)]
        if thrval is None:
            medfit = np.quantile([patient.q for patient in patients], quantile)
        else:
            medfit = thrval

        odir = os.path.dirname(outfile)
        self.monthLimit = 10000
        #        try:
        if True:
            T = [min(patient.OS if OS else patient.PFS, self.monthLimit) for patient in patients]
            E = [patient.dead for patient in patients]
            #self.logger(T)

            group = [assign_sample(patient, medfit) for patient in patients]
            pnames = [patient.name for patient in patients]
            if prism:
                patients1 = [el for el in zip(pnames, T, E, group) if el[3] == 'negative']
                patients2 = [el for el in zip(pnames, T, E, group) if el[3] == 'positive']
                prismfile = os.path.join(odir, "PRISM_" + self.model + ".txt")
                of = open(prismfile, 'w')
                for pat in patients1:
                    [pname, t, dead, _] = pat
                    if dead:
                        dead = 1
                    else:
                        dead = 0
                    of.write(pname + "\t" + str(t) + "\t" + str(dead) + "\t\n")
                for pat in patients2:
                    [pname, t, dead, _] = pat
                    if dead:
                        dead = 1
                    else:
                        dead = 0
                    of.write(pname + "\t" + str(t) + "\t" + "\t" + str(dead) + "\n")
                of.close()
            dat = {}
            dat["T"] = T
            dat["E"] = E
            dat["group"] = group
            df = DataFrame(data=dat)
            T = df['T']
            E = df['E']
            groups = df['group']
            ix = (groups == 'negative')
            patvalues = patients
            group = [assign_sample(patient, medfit) for patient in patvalues]

            n = len([el for el in group if el == "negative"])
            p = len([el for el in group if el == "positive"])

            dat = {}
            dat["T"] = T
            dat["E"] = E
            dat["group"] = group
            df = DataFrame(data=dat)
            kmf = KaplanMeierFitter()
            kmfp = KaplanMeierFitter()
            kmfn = KaplanMeierFitter()
            c = 0
            neg = []
            pos = []
            for i in ix:
                pname = patients[c].name
                if i:
                    neg.append(pname)
                else:
                    pos.append(pname)
                c += 1
            kmfp.fit(T[~ix], E[~ix], label='positive')
            kmfn.fit(T[ix], E[ix], label='negative')
            pestimate = getattr(kmfp, "survival_function_")
            nestimate = getattr(kmfn, "survival_function_")
            new_index = np.concatenate((nestimate.index, pestimate.index))
            new_index = np.unique(new_index)
            kmf.fit(T[~ix], E[~ix], label=r'$\log n(\tau) > $' + str(round(medfit, 3)) + " (n=" + str(len(pos)) + ")")
            ax = kmf.plot(show_censors=True)
            title = "q = "+str(quantile)
            if not pval is None:
                title = ", p-value = " + str(pval)

            ax.set_title(title)
            kmf.fit(T[ix], E[ix], label=r'$\log n(\tau) \leq$' + str(round(medfit, 3)) + " (n=" + str(len(neg)) + ")")
            kmf.plot(ax=ax, show_censors=True)
            plt.ylim((0., 1.))
            if show:
                plt.show()
            plt.savefig(outfile)

            plt.close()

    #        except:
    #            self.logger("Exception in plotting", 0)

    ##############################################
    # Other output
    ##############################################


    def toVCF(self, outdir, trunkal=False):
        '''
        Creates VCF files for each clone in each sample.
        Writes the files to outdir folder. The vcf files are
        called after sample name and clone identifiers,
        <samplename>_<clone_id>.vcf


        :param outdir: str

        :param trunkal: bool
            if True it will only create 2 vcf files, to separate
            trunkal from non-trunkal mutations. Otherwise a separate vcf file is created for each clone.

        '''

        for sample in self.get_samples():
            sample.toVCF(outdir, trunkal=trunkal)

    ############################################
    # formats
    ############################################

    def toJSON(self):
        js = [patient.toJSON() for patient in self.patients.values()]
        return js
