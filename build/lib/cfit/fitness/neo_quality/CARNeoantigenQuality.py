#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 09:57:18 2020

@author: zacharysethna
"""

import numpy as np

from cfit.CoreObject import CoreObject
from cfit.fitness.neo_quality.ARNeoantigenQuality import ARNeoantigenQuality
from cfit.fitness.neo_quality.EpitopeDistance import EpitopeDistance
from cfit.util.Utils import Utils
from cfit.fitness.neo_quality.QualityAttributes import QualityAttributes


class CAR_attributes(QualityAttributes):
    '''
    A class for neoantigen.qattributes, for storing ((1-w)*logA+w*logC)xR specific pre-computed data.
    '''

    def __init__(self, a=26):
        self.GbSum = -a
        self.epitopeAlignments = []
        self.alignedEpitopes = {}
        self.maxEpitope = None

        self.mt_iedb_minD = IEDBminD()
        self.wt_iedb_minD = IEDBminD()

        self.A = 0
        self.R = 0
        self.D = 0
        self.R_d = 0

    def get_quality_components(self):
        di = {"A": self.A, "R": self.R, "D": self.D}
        return di

    def out(self):
        epitopename = "None"
        epitopeseq = "None"
        eid = -1
        iedbid = -1
        shift = 0
        mtfrom = -1
        mtto = -1
        mtscore = 0
        mtaln = "None"
        species = "None"
        aln = ["-", "-"]
        if self.maxepitope is not None:
            epitopename = self.maxepitope[0].epitopeName
            epitopeseq = self.maxepitope[0].epitopeSeq
            eid = self.maxepitope[0].epitope.id
            iedbid = self.maxepitope[0].epitope.iedbid
            shift = self.maxepitope[0].epitope.shift
            species = self.maxepitope[0].epitope.species
            mtfrom = self.maxepitope[0].mtfrom
            mtto = self.maxepitope[0].mtto
            mtaln = self.maxepitope[0].mtaln
            mtscore = self.maxepitope[0].mutScoreSW
            # aln = AlignedEpitope.align(neo.mtPeptide, epitopeseq)[:2]

        return {"D": self.D, "A": self.A, "R": self.R, "epitopename": epitopename, "epitopeseq": epitopeseq,
                "eid": eid, "iedbid": iedbid, "shift": shift, "mtfrom": mtfrom, "mtto": mtto, "mtscore": mtscore,
                "mtaln": mtaln, "species": species, "aln0": aln[0], "aln1": aln[1]}


class IEDBminD(CoreObject):
    def __init__(self):
        self.D = np.inf
        self.min_dist_pep = ''
        self.min_dist_pep_id = None


class CARNeoantigenQuality(ARNeoantigenQuality, EpitopeDistance):
    '''
    Implement (w*logC + (1-w)logA)xR neoantigen quality

    Attributes:

        alndir: str
            directory with the precomputed neoantigen alignments for each sample

        epitopes: dict
            dictionary mapping IEDB epitope identifiers to cfit.fitness.Epitope objects

        eids: list
            list of epitope identifiers

        epitopemask: dict
            dictionary mapping IEDB epitope identifiers to float values, allowing for masking epitopes from
            contributing to neoantigen quality (if the mapped value is 0)

        zeromask: dict

        wt: bool
    '''

    L = 1.  # concentration
    M = 1.  # mutant peptide concentration
    W = 1.  # wildtype peptide concentration
    meanD = 11.5
    stdD = 2
    WEPS = 0.0
    WTCAP = Utils.INF

    def __init__(self, iedbfasta="", alndir="", matrix="blosum62",
                 ep_dist_model_name="all_tcr_all_combos_model"):
        '''

        :param iedbfasta: str
            path to the IEDB fasta file

        :param alndir: str
            path to the alignment directory, neoantigens/*

        :param matrix: str

        :param ep_dist_model_name: str

        '''
        ARNeoantigenQuality.__init__(self, iedbfasta=iedbfasta, alndir=alndir, matrix=matrix)
        EpitopeDistance.__init__(self, model_name=ep_dist_model_name)
        if ep_dist_model_name is not None:
            self.set_model(ep_dist_model_name)

        self.wt = False
        self.d_scaling = lambda x: x

    def compute_quality(self, neo, a=26, k=1, wd=1.0, kdthr=None, modelName=None,
                        include_R=True):
        '''

        :param neo: cfit.tree.mutation.Neoantigen

        :param a: float

        :param k: float

        :param wd: float
            weight of D term

        :param kdthr: float
            hard threshold on the dissociation constant

        :param include_R: bool
            whether to include R

        :return: float
            recognition_potential

        '''

        #if len(neo.mtPeptide) != 9:
        #    return 0

        neo.qattributes.GbSum = -a

        A = self.compute_A(neo, KDnormalize=1)
        D = self.epitope_dist(neo.wtPeptide, neo.mtPeptide)
        R = 1.
        if include_R:
            R = self.compute_recognition_probability(neo, a, k)

        neo.qattributes.D = D
        neo.qattributes.A = A
        neo.qattributes.R = R
        neo.qattributes.neotype_prob = self.compute_neotype_prob(neo)
        neo.qattributes.mut_pres_prob = self.compute_mut_presentation(neo)

        quality = ((1 - wd) * np.log(A) + wd * D) * R #* 1/neo.kD

        if kdthr is not None:
            if neo.kD > kdthr:
                quality = 0.0

        neo.quality = quality
        return quality


    def compute_neoantigen_sample_quality(self, neo, sample, a=26, k=1, w=0.5, kdthr=None,
                                 include_R=True, just_R=False):
        '''

        Compute quality of a neoantigen in a sample.

        :param neo: cfit.tree.mutation.Neoantigen

        :param sample: cfit.patient.Sample

        :param a: float
            shift parameter of the R term

        :param k: float
            slope component of the R term

        :param w: float
            relative weight of logC term [ w*log C + (1-w)*log A ]

        :param kdthr: float
            hard threshold on the dissociation constant

        :param include_R: bool
            whether to include R: if false (*)*R model is used

        :param just_R: bool
            whether to include only R: if true R model is used

        :return: float
            neoantigen quality

        '''

#        if len(neo.mtPeptide) != 9:
#            return 0
        neo.qattributes.GbSum = -a
        neo_concentration = sample.mutations[neo.mid].expression

        D = self.epitope_dist(neo.wtPeptide, neo.mtPeptide)
        A = self.compute_A(neo, KDnormalize=1, neo_concetration=neo_concentration)

        R = 1.
        if include_R:
            R = self.compute_recognition_probability(neo, a, k)

        neo.qattributes.D = D
        neo.qattributes.A = A
        neo.qattributes.R = R
        neo.qattributes.neotype_prob = self.compute_neotype_prob(neo)
        neo.qattributes.mut_pres_prob = self.compute_mut_presentation(neo)


        if just_R:
            quality = R
        else:
            quality = ((1 - w) * np.log(A) + w * D) * R

        if kdthr is not None:
            if neo.kD > kdthr:
                quality = 0.0
        return quality

    def compute_neoantigen_sample_quality_components(self, neo, sample, a=26, k=1, w=0.5,
                                                     include_R=True, just_R=False):
        '''

        :param neo: cfit.tree.mutation.Neoantigen

        :param sample: cfit.patient.Sample

        :param a: float

        :param k: float

        :param w: float
            relative weight of logC term [ w*logC + (1-w)*logA ]

        :param kdthr: float
            hard threshold on the dissociation constant

        :param include_R: bool
            whether to include R

        :return: float
            neoantigen quality

        '''

        neo_concentration = sample.mutations[neo.mid].expression

        D = self.epitope_dist(neo.wtPeptide, neo.mtPeptide)
        A = self.compute_A(neo, KDnormalize=1, neo_concetration=neo_concentration)
        R = 1.
        if include_R or just_R:
            R = self.compute_recognition_probability(neo, a, k)

        return {"D": D, "logA": np.log(A), "R": R}

    def get_parameters(self):
        p = super(CARNeoantigenQuality, self).get_parameters()
        return p

    def compute_mut_presentation(self, neo, mtKd_thresh=500, mtKd_n=1):
        return 1 / (1 + (neo.kD / mtKd_thresh) ** mtKd_n)

    def compute_neotype_prob(self, neo, wtKd_thresh=500, wtKd_n=1, D_thresh=4, D_n=2):

        wtKd = neo.wtkD
        D = self.epitope_dist(neo.wtPeptide, neo.mtPeptide)
        P_t1 = 1 / (1 + (wtKd_thresh / wtKd) ** wtKd_n)
        P_t2 = 1 / (1 + np.exp(D_n * (D_thresh - D)))

        return {'P_t1': P_t1, 'P_t2': P_t2}

    def initialize_neoantigens(self, anl, recompute_scores=True):
        '''
        Initializes neoantigen quality related data for the neoantigens in the cohort

        :param anl: cfit.util.Analysis

        :param recompute_scores: bool
            whether to recompute the alignment scores provided in the alignment files. If yes, the
            precomputed alignments serve only as a list of candidates.
        '''
        patients = list(anl.patients.values())
        for patient in patients:
            for sample in patient.samples:
                for neo in sample.neoantigens:
                    neo.qattributes = CAR_attributes()
        self.import_patient_data_alignments(self.alndir, patients, recompute_scores=recompute_scores)


    def initialize_on_patients(self, patients, recompute_scores=True):
        '''
        Initialization

        :param patients: list
            list of cfit.patient.Patient objects

        :param recompute_scores: bool
            whether to recompute the alignment scores provided in the alignment files. If yes, the
            precomputed alignments serve only as a list of candidates.

        '''
        super(CARNeoantigenQuality).initialize_on_patients(patients, recompute_scores=recompute_scores)
        for patient in patients:
            for sample in patient.samples:
                for neo in sample.neoantigens:
                    neo.qattributres = CAR_attributes()


    def compute_D(self, neo):
        '''
        Compute Distance between MT and WT peptides
        '''
        return self.d_scaling(self.epitope_dist(neo.wtPeptide, neo.mtPeptide))
