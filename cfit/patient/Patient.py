'''
Created on Feb 21, 2017

@author: mluksza
'''

import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd

from cfit.CoreObject import CoreObject


class Patient(CoreObject):
    '''
    Class representing a patient, with possibly many samples.

    Attributes:
        name: str
            patient identifier

        survival: float
            survival time

        __OS: float:
            survival time

        __PFS: float
            progression-free survival

        dead: bool
            OS status of the patient

        PFS_status: bool
            PFS status of the patient

        response: bool
            whether the patient is a responder

        __samples: list
            list of samples of that patient (instances of class cfit.patient.Sample)

        __q: float
            score due to fitness (eg. n(tau))

        HLAS: list
            list of patient MHC alleles (str)

        mutPres: defaultdict(lambda: defaultdict(lambda: False))
            dictionary that maps whether a given mutation (by mutation identifier) is present in a given sample (by sample name)
            mid -> sample_name -> bool

        __cohort: str
            name of the cohort

        __type: str
            type of patients tumor (eg. Metachronous/Synchronous)

        neoantigens: dict
            maps neo identifier (str) to Neoantigen object

        fsneoantigens: dict
            maps neo identifier (str) to Neoantigen object

        mutations: dict
            maps mutation identifier (str) to Mutation object

        mutation2neoantigens: dict: str -> list
            maps mutation identifier (str) to the list of Neoantigen objects

        mutation2fsneoantigens: dict: str -> list

        trees: list of cfit.tree.Tree
            list of trees

        oneTree: cfit.tree.Tree
            homogenous structure version

    '''

    def __init__(self, params):
        '''
        Constructor
        '''
        [name, survival, dead, classification] = params
        self.name = name  # patient identifier
        self.survival = survival  # survival time
        self.__OS = survival
        self.__PFS = survival
        self.dead = dead  # status
        self.response = (classification == "response")
        self.__samples = []  # list of samples of that patient (instances of class Sample)
        self.__q = 0  # score due to fitness (eg. n(tau))
        self.HLAS = []
        # dictionary mut_id (chrom_pos_ref_alt) -> sample_name -> bool
        self.mutPres = defaultdict(lambda: defaultdict(lambda: False))
        self.__cohort = None
        self.__type = None
        self.neoantigens = dict()
        self.mutations = dict()
        self.mutation2neoantigens = defaultdict(list)

        self.fsneoantigens = dict()
        self.mutation2fsneoantigens = defaultdict(list)

        self.trees = []
        self.vaccine_mutations = []
    @property
    def OS(self):
        return self.__OS

    @OS.setter
    def OS(self, OS):
        self.__OS = OS

    @property
    def PFS(self):
        return self.__PFS

    @PFS.setter
    def PFS(self, PFS):
        self.__PFS = PFS

    @property
    def q(self):
        return self.__q

    @q.setter
    def q(self, q):
        self.__q = q

    @property
    def samples(self):
        return self.__samples

    @samples.setter
    def samples(self, samples):
        self.__samples = samples

    @property
    def cohort(self):
        return self.__cohort

    @cohort.setter
    def cohort(self, cohort):
        self.__cohort = cohort

    @property
    def type(self):
        return self.__type

    @type.setter
    def type(self, type):
        self.__type = type

    def rename(self, newname):
        self.name = newname

    def short_sample_name(self, fullsamplename):
        '''
        Returns short name of a sample. Getting obsolete.

        :param fullsamplename: str
        :return: str
        '''
        return fullsamplename.replace(self.name, "")

    def add_sample(self, sample, **kwargs):
        '''

        :param sample: cfit.patient.Sample
            sample object to be added
        :param kwargs: dict
            not used in this class

        :return:
        '''
        sample.response = self.response
        self.__samples.append(sample)
        sname = self.short_sample_name(sample.name)
        for mutid in sample.mutations:
            self.mutPres[mutid][sname] = True

    def remove_sample(self, sname):
        '''
        Remove a sample of a given name
        :param sname: str
        '''
        self.__samples = [sample for sample in self.__samples if sample.name != sname]
        for mutid in self.mutPres:
            if sname in self.mutPres[mutid]:
                del self.mutPres[mutid][sname]

    def mark_shared_clones(self, tp1, tp2, eps):
        '''

        :param tp1: str
        :param tp2: str
        :param eps: float
        :return:
        '''
        self.logger("Method not implemented / relevant.")
        return False

    def mark_lost_clones(self, tp1, tp2, eps):
        '''

        :param tp1: str
        :param tp2: str
        :param eps: float
        :return:
        '''
        self.logger("Method not implemented / relevant.")
        return False

    def get_samples(self, simple=False):
        '''
        Returns a list of samples with the name of time point.
        :return: list
        '''
        if simple:
            return [sample for sample in self.samples]
        else:
            return [['TP1', sample] for sample in self.samples]

    def compute_score(self, criterion, clonal, beta, tau):
        '''
        Computes the score which will be the biomarker based on which samples will classified in survival analysis.

        :param criterion: str
            which patient scoring criterion to use:
                n(tau) - POPULATION_SIZE
                <f> - FITNESS
                d n(tau)/d tau - SPEED

        :param clonal: bool
            if True tumor heterogeneity and trees are taken into account

        :param beta: float
            tree weighting parameter

        :param tau: float
            time parameter

        :return:
            sets patient's attribute __q.
        '''
        self.__q = 0.0
        if criterion == "FITNESS" or tau == 0:
            qs = [sample.compute_population_fitness(clonal, beta, tau) for sample in self.samples]
        elif criterion == "POPULATION_SIZE":
            qs = [sample.compute_predicted_population_size(clonal, beta, tau) for sample in self.samples]
        elif criterion == "SPEED":
            qs = [sample.compute_predicted_population_size_derivative(clonal, beta, tau) for sample in self.samples]
        if len(self.samples) > 0:
            self.__q = np.mean(qs)

    def check_mutation(self, mid):
        '''
        Checks if a given mutation belongs to the list of mutations. Allows for shorter alt alleles
        :param mid: str
            <chrom>_<pos>_<ref>_<alt>

        :return: (str, bool)
            (mutation_id, whether it belongs here)
            mutation_id may be corrected (phylowgs seems to be cutting alt allele to one character)
        '''
        if mid in self.mutations:
            return (mid, True)
        else:
            for mid0 in self.mutations:
                if mid in mid0:
                    return (mid0, True)
        return (mid, False)

    def meta_beta_function(self, function, beta=1.0, **kwargs):
        '''
        Evaluates an arbitrary function over trees in samples in a patient.

        :param function: function

        :param beta: float

        :param kwargs: dict
        '''
        q = np.mean([sample.meta_beta_function(function, beta=beta, **kwargs) for sample in self.samples])
#        q = float(sum(
#            list(map(lambda sample: sample.meta_beta_function(function, beta=beta, **kwargs), self.samples)))) / len(
#            self.samples)
        self.__q = q

    def compute_predicted_population_size(self, clonal, beta, tau):
        '''
        n(tau)

        :param clonal: bool
            if True tumor heterogeneity and trees are taken into account

        :param beta: float
            tree weighting parameter

        :param tau: float
            time parameter

        :return:
            sets patient's attribute __q.

        '''
        self.__q = 0.0
        map(lambda sample: sample.compute_predicted_population_size(clonal, beta, tau), self.samples)
        if len(self.samples) > 0:
            self.__q = float(sum(list(map(lambda sample: sample.q, self.samples)))) / len(self.samples)

    def entropy(self, beta=1.0):
        '''
        Tree based entropy of the clone sizes, averaged over samples.

        :param beta: float
            Tree weighting parameter

        :return: float
        '''
        if len(self.samples) == 0:
            return 0
        vs = [sample.entropy(beta=beta) for sample in self.samples]
        return sum(vs) / len(vs)

    def get_effective_number_of_clones(self, beta):
        '''

        :param beta: float
            tree weighting parameter
        :return: float
        '''
        try:
            v = float(sum(list(map(lambda sample: sample.get_effective_number_of_clones(beta), self.samples))))
            v /= len(self.samples)
            return v
        except:
            return 0.

    def gather_HLAs_from_neoantigens(self):
        '''
        sets HLA alleles
        :return:
        '''
        alleles = list(set([neo.allele for sample in self.samples for neo in sample.neoantigens]))
        alleles.sort()
        self.HLAS = alleles

    def set_HLAs(self, allele_list):
        '''
        :param allele_list: list
            list of str
        '''
        self.HLAS = allele_list

    def get_mutation_data(self, exclusive=True, clonal=True):
        '''
        Returns a DataFrame for each sample,
        with data about mutations in the sample,
        grouped by clonality of the highest scoring trees.

        :param exclusive: bool
            if True will report each mutation only once,
            in the clone where it originates

        :param clonal: bool
            if True/1 will report mutation data grouped by clones
            of the highest scoring tree, otherwise will
            assume homogeneity

        :return: pd.DataFrame
        '''

        dmdatas = {}
        snames = [self.short_sample_name(sn.name) for sn in self.samples]
        for sample in self.samples:
            # get mutation data for every sample
            sname = (sample.name).replace(self.name, "")
            mdata = sample.get_mutation_data(exclusive=exclusive, clonal=clonal)
            pres = []
            # annotate mutations with their presence in other samples
            for line in mdata.itertuples():
                pline = [int(self.mutPres[line.Mutation][sn]) for sn in snames]
                pres.append(pline)

            pres = pd.DataFrame(pres)
            pres.columns = snames
            for col in snames:
                mdata[col] = list(pres[col])
            dmdatas[sname] = mdata
        return dmdatas

    def syn(self):
        '''
        Returns a DataFrame for each sample,
        with synonymous mutation counts
        '''
        dmdatas = {}
        for sample in self.samples:
            # get mutation data for every sample
            sname = (sample.name).replace(self.name, "")
            dmdatas[sname] = [sample.TMB_syn(), sample.effective_TMB_syn()]
        return dmdatas

    def set_mutation_node_index(self):
        '''
        Sets the mutation-node index, for faster access to node mutation content
        '''
        for tree in self.trees:
            tree.set_mutation_node_index()

    def write_mutations(self, outdir, *args, **kwargs):
        '''
        Write patient samples' mutations to files
        Parameters:
            - outdir - output folder
        '''
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        dmdatas = self.get_mutation_data(args, kwargs)
        for sname in dmdatas:
            path = os.path.join(outdir, self.name + sname + ".txt")
            mdata = dmdatas[sname]
            mdata.to_csv(path, index=False, sep="\t")

    def set_vaccine_mutations(self, mutation_ids):
        '''
        :param mutation_ids: list
            list if mutation identifiers (str)
        '''
        self.vaccine_mutations = mutation_ids
        for sample in self.__samples:
            sample.vaccine_mutations = mutation_ids

    def topTreeStats(self):
        sample = self.samples[0]
        if len(sample.trees) > 0:
            tree = sample.trees[0]
        else:
            tree = sample.oneTree
        stats = tree.get_stats()
        return stats

    def writeJSON(self, outdir, prefix=""):
        '''
        Writes json trees of samples to a folder.
        Parameters:
             - outdir : output folder where json files will be
                         written to
        '''
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        for sample in self.samples:
            sname = prefix + sample.name
            jpath = os.path.join(outdir, "tree_" + sname + ".json")
            js = self.toJSON()
            jssamples = {sname: sample.toJSON()}
            js['samples'] = jssamples
            js['id'] = sample.name
            with open(jpath, 'w') as of:
                json.dump(js, of, indent=True)


    def toJSON(self, include_response=False):
        js = {'id': self.name,
              'OS': self.OS,
              "PFS": self.PFS,
              'status': self.dead,
              'cohort': self.cohort}
        if include_response:
            js['response'] = self.response
#        jssamples = {}
#        for sample in self.samples:
#            jssamples[sample.name] = sample.toJSON()
#        js['samples'] = jssamples
        jstrees = [tree.toJSON()  for tree in self.trees]
        js["trees"] = jstrees
        js['HLA_genes'] = self.HLAS
        js["mutations"] = [mut.toJSON() for mut in self.mutations.values()]
        js["neoantigens"] = [neo.toJSON() for neo in self.neoantigens.values()]
        return js
