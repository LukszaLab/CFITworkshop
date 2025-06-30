import json
import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML

import cfit

from cfit.fitness.neo_quality.AlignedEpitope import AlignedEpitope
from cfit.fitness.neo_quality.Epitope import Epitope
from cfit.fitness.neo_quality.EpitopeMask import EpitopeMask
from cfit.fitness.neo_quality.NeoantigenQuality import NeoantigenQuality
from cfit.util.Utils import Utils
from cfit.fitness.neo_quality.QualityAttributes import QualityAttributes

class AR_attributes(QualityAttributes):
    '''
    A class for neoantigen.qattributes, for storing AxR specific pre-computed data.
    '''

    def __init__(self, a=26):
        self.GbSum = -a
        self.epitopeAlignments = []
        self.alignedEpitopes = {}
        self.maxEpitope = None

        self.A = 0
        self.R = 0

    def get_quality_components(self):
        di = {"A": self.A, "R": self.R}
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

        return {"A": self.A, "R": self.R, "epitopename": epitopename, "epitopeseq": epitopeseq,
                "eid": eid, "iedbid": iedbid, "shift": shift, "mtfrom": mtfrom, "mtto": mtto, "mtscore": mtscore,
                "mtaln": mtaln, "species": species, "aln0": aln[0], "aln1": aln[1]}


class ARNeoantigenQuality(NeoantigenQuality):
    '''
    Implement AxR neoantigen quality

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

        a: float

        k: float

        include_A: bool

        include_R: bool

        kdthr: float

        KDnormalize: float

        wt: bool
    '''

    L = 1.  # concentration
    M = 1.  # mutant peptide concentration
    W = 1.  # wildtype peptide concentration
    WEPS = 0.0
    WTCAP = Utils.INF

    def __init__(self, iedbfasta="", alndir="", matrix="blosum62"):
        '''

        :param iedbfasta: str

        :param alndir: str

        :param matrix: str
        '''
        '''

        :param iedbfasta: str
            path to the IEDB fasta file

        :param alndir: str
            path to the alignment directory, neoantigens/*

        '''
        super(ARNeoantigenQuality, self).__init__()
        self.alndir = alndir
        self.epitopes = {}
        self.eids = set()
        self.epitopemask = EpitopeMask(["empty", {}])
        self.zeromask = EpitopeMask(["empty", {}])
        self.__read_epitope_data_from_fasta_file(fastafile=iedbfasta)
        self.wt = False

        self.a = None
        self.k = None
        self.include_A = True
        self.include_R = True
        self.kdthr = Utils.INF
        self.KDnormalize = 1.

        smatrixfile = os.path.join(os.path.dirname(cfit.__file__), "data", "matrices", matrix + ".json")
        with open(smatrixfile) as f:
            submatrix = json.load(f)
            sm = {}
            for el in submatrix:
                ref = el[0]
                alt = el[-1]
                sm[(ref, alt)] = submatrix[el]
        self.substmatrix = sm

    def get_parameters(self):
        '''
        Dictionary of parameters specific to the model
        :return: dict
            str -> value
            name of the parameter -> value
        '''

        p = {}
        p["a"] = self.a
        p["k"] = self.k
        p["include_A"] = self.include_A
        p["include_R"] = self.include_R
        p["kdthr"] = self.kdthr
        p["KDnormalize"] = self.KDnormalize
        model = "noname"
        if self.include_A and self.include_R:
            if self.KDnormalize == 1:
                model = "A"
            else:
                model = "Kd"
            model += "R"
        elif self.include_A:
            if self.KDnormalize == 1:
                model = "A"
            else:
                model = "Kd"
        elif self.include_R:
            model = "R"
        p["model_name"] = model
        return p

    def compute_quality(self, neo, a=26, k=1, HLAW=None, include_A=True, include_R=True, kdthr=None, KDnormalize=1.):
        '''

        :param neo: cfit.tree.mutation.Neoantigen

        :param a: float

        :param k: float

        :param include_A: bool
            whether to evaluate A

        :param include_R: bool
            whether to evaluate R

        :param kdthr: float
            hard threshold on the dissociation constant

        :return: float
            recognition_potential - AxR

        '''

        self.a = a
        self.k = k
        self.include_A = include_A
        self.include_R = include_R
        self.kdthr = kdthr
        self.KDnormalize = KDnormalize

        neo.qattributes.GbSum = -a
        A = self.compute_A(neo, KDnormalize=self.KDnormalize)
        # A = neo.wtkD / neo.kD
        if not include_A:
            A = 1.0
        R = self.compute_recognition_probability(neo, a, k)
        if not include_R:
            R = 1.0
        neo.qattributes.A = A
        neo.qattributes.R = R
        quality = A * R
        if HLAW is not None:
            W = HLAW.get_weight(
                neo)  # masks some neoantigens, if set (eg. position 2&9 with hydrophobic residue for dpos="1-9H"
        else:
            W = self.W
        # was used for IO datasets)
        quality *= W
        if kdthr is not None:
            if neo.kD > kdthr:
                quality = 0.0

        neo.quality = quality
        return quality

    def compute_neoantigen_sample_quality(self, neo, sample,
                                          a=26, k=1, HLAW=None, include_A=True, include_R=True, kdthr=None,
                                          KDnormalize=1.):
        '''

        Computes quality of the neoantigen, modifies sample.neoantigenQualities attribute

        :param neo: cfit.tree.mutation.Neoantigen

        :param sample: cfit.patient.Sample

        :param a: float

        :param k: float

        :param include_A: bool
            whether to evaluate A

        :param include_R: bool
            whether to evaluate R

        :param kdthr: float
            hard threshold on the dissociation constant

        :return: float
            recognition_potential - AxR
        '''

        self.a = a
        self.k = k
        self.include_A = include_A
        self.include_R = include_R
        self.kdthr = kdthr
        self.KDnormalize = KDnormalize

        neo.qattributes.GbSum = -a
        A = self.compute_A(neo, KDnormalize=self.KDnormalize, neo_concetration=sample.mutations[neo.mid].expression)
        if not include_A:
            A = 1.0
        R = self.compute_recognition_probability(neo, a, k)
        if not include_R:
            R = 1.0
        neo.qattributes.A = A
        neo.qattributes.R = R
        quality = A * R
        if HLAW is not None:
            W = HLAW.get_weight(
                neo)  # masks some neoantigens, if set (eg. position 2&9 with hydrophobic residue for dpos="1-9H"
        else:
            W = self.W
        # was used for IO datasets)
        quality *= W
        if kdthr is not None:
            if neo.kD > kdthr:
                quality = 0.0

        return quality

    def initialize_on_patients(self, patients, recompute_scores=True):
        '''
        Initialization

        :param patients: list
            list of cfit.patient.Patient objects

        :param recompute_scores: bool
            whether to recompute the alignment scores provided in the alignment files. If yes, the
            precomputed alignments serve only as a list of candidates.

        '''
        for patient in patients:
            for sample in patient.samples:
                for neo in sample.neoantigens:
                    neo.qattributes = AR_attributes()
        self.import_patient_data_alignments(self.alndir, patients, recompute_scores=recompute_scores)

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
                    neo.qattributes = AR_attributes()
        self.import_patient_data_alignments(self.alndir, patients, recompute_scores=recompute_scores)

    def correct_neoantigen_WTkd(self, neo):
        '''
        Correct the wildtype dissociation constant with the epsilon correction
        :param neo: cfit.tree.mutation.Neoantigen

        :return: float
            corrected wtKd
        '''
        kd = neo.wtkD
        prb = self.W / (self.W + kd)
        pru = kd / (self.W + kd)
        eps = self.WEPS
        prb += eps
        pru += eps
        z = 1 + 2 * eps
        prb /= z
        pru /= z
        return pru / prb

    def compute_A(self, neo, KDnormalize=1, neo_concetration=1.):
        '''

        Computes the "amplitude", depending on the external setting of Utils.KDnormalize:
        Default - KDNormalize = 1

        :param neo: cfit.tree.mutation.Neoantigen

        :param KDnormalize: float
            how to compute A

        :param neo_concetration: float
            relative concentration of the neoantigen
        :return: float
        '''

        if KDnormalize == 0:
            A = 1.
        elif KDnormalize == 1:  # pmb/pmu *pwu/pwb, DEFAULT
            A = neo_concetration / neo.kD * self.correct_neoantigen_WTkd(neo)
        elif KDnormalize == 1.5:  # pmb * pwu/pwb
            A = neo_concetration / (neo.kD + neo_concetration) * self.correct_neoantigen_WTkd(neo)
        elif KDnormalize == 0.5:
            A = np.log(self.M * self.correct_neoantigen_WTkd(neo) / (neo.kD))
        elif KDnormalize == 0.1:
            A = neo_concetration / (neo.kD + neo_concetration) * np.log(self.correct_neoantigen_WTkd(neo))
        elif KDnormalize == 2:  # tolerance only T
            A = self.correct_neoantigen_WTkd(neo)
        elif KDnormalize == 3:  # mutant peptide binding
            A = neo_concetration / (neo.kD + neo_concetration)
        elif KDnormalize == 4:  # mutant only M
            A = 1 / self.kD
        elif KDnormalize == 5:
            A = neo_concetration / (neo_concetration + neo.kD) * neo.wtkD / (neo.wtkD + self.W)
        return A

    def compute_recognition_probability(self, neo, a, k, additional_score=0):
        '''
        Compute R value for the given neoantigen.

        :param neo: cfit.tree.mutation.Neoantigen
            neoantigen object

        :param a: float
            sigmoid function shift

        :param k: float
            sigmoid function slope

        :param additional_score: float

        :return: float
            TCR recognition probability
        '''
        neo.qattributes.maxepitope = None
        if len(neo.qattributes.alignedEpitopes.values()) == 0:
            return 0.0
        bindProb = 1.0
        if a >= 0:
            # if kd0>0:
            #           #     tolscore = self.epitope_dist(neo.wtPeptide, neo.mtPeptide)*kd0/(kd0+neo.wtkD)*weight_tol
            #    tolscore = self.epitope_dist(neo.wtPeptide, neo.mtPeptide)*np.exp(-neo.wtkD/kd0)*weight_tol
            # else:
            #    tolscore = self.epitope_dist(neo.wtPeptide, neo.mtPeptide) * weight_tol

            tolscore = additional_score
            bindProb = 0.0
            if len(neo.qattributes.alignedEpitopes) > 0:
                ascores = list(map(lambda aepitope: [aepitope,
                                                     aepitope.get_score(self.wt)
                                                     * self.epitopemask.mask_value(aepitope.get_epitope_id())
                                                     ],
                                   neo.qattributes.alignedEpitopes.values()))
                ascores.sort(key=lambda _el: -_el[1])
                v = [-k * (a - score[1] - tolscore) for score in ascores]
                lgb = Utils.log_sum(v)
                neo.qattributes.GbSum = -lgb
                lZ = Utils.log_sum2(0, lgb)
                bindProb = np.exp(lgb - lZ)
                neo.qattributes.maxepitope = ascores[0]
                neo.qattributes.maxepitope[1] = bindProb
        return bindProb

    '''
    Private methods 
    '''

    def __read_epitope_data_from_fasta_file(self, fastafile):
        '''
        Initializes the object and imports epitope sequences and data from the IEDB.fasta file

        :param fastafile: str
            IEDB.fasta file path

        '''

        f = open(fastafile)
        seqs = SeqIO.parse(f, "fasta")
        vepitopes = []
        for seq in seqs:
            eid = self.__get_epitope_id(seq.description)
            line = str(seq.seq) + "\t" + seq.description
            epitope = Epitope(line, eid)
            vepitopes.append(epitope)
        f.close()
        self.epitopes = {}
        for epitope in vepitopes:
            self.epitopes[epitope.id] = epitope
            self.eids.add(epitope.id)
            self.epitopemask.set_mask_value(epitope.id, 1)
            self.zeromask.set_mask_value(epitope.id, 1)
        self.epitopemask.set_epitopes(self.epitopes)
        self.zeromask.set_epitopes(self.epitopes)

    def __get_epitope_id(self, desc):
        '''
        Get IEDB epitope identifier

        :param desc: str
            extracts epitope identifier from description of the fasta header

        :return: int
            the identifier
        '''
        try:
            tab = desc.split("|")
            eid = int(tab[0])
        except AttributeError:
            eid = desc
        return eid

    def __find_epitope(self, desc):
        '''
        Returns epitope with a matching epitope identifier as the fasta file header

        :param desc: str
            fasta file entry header

        :return: cfit.fitness.Epitope
            the mapped epitope
        '''
        eid = self.__get_epitope_id(desc)
        return self.epitopes[eid]

    def __add_aligned_epitope(self, neo, epitope, mt_score=0, wt_score=0, recompute_scores=True):
        '''
        Creates an aligned epitope and updates neoantigen object

        :param neo: cfit.tree.mutation.Neoantigen
            neoantigen that is being aligned

        :param epitope: cfit.fitness.Epitope
            the IEDB epitope

        :return: cfit.fitness.AlignedEpitope
            aligned epitope.
            It also modifies in place the neoantigen object, by appending setting qattributes attribute.
        '''

        aepitope = AlignedEpitope(epitope, neo, substmatrix=self.substmatrix,
                                  mt_score=mt_score, wt_score=wt_score,
                                  recompute_scores=recompute_scores)
        if epitope.id in neo.qattributes.alignedEpitopes:
            # records best alignment with a given epitope for a givent neoantigen
            if aepitope.mutScore > neo.qattributes.alignedEpitopes[epitope.id].mutScore:
                neo.qattributes.alignedEpitopes[epitope.id] = aepitope
            else:
                return None
        else:
            neo.qattributes.alignedEpitopes[epitope.id] = aepitope
        return aepitope

    def __read_sample_blast_alignments_xml_format(self, sample, path):
        '''
        Reads precomputed Blast alignments in XML format (old)

        :param sample: cfit.patient.Sample
            The sample for which the alignments are imported

        :param path: str
            path to the precomputed alignments, in the xml format

        :return: list
            list of aligned epitopes, cfit.fitness.AlignedEpitope. It modifies neoantigen objects of the sample,
            with add_aligned_epitope method.
        '''

        f = open(path)
        blast_records = NCBIXML.parse(f)
        nc = {}
        nc["WT"] = 0
        nc["MUT"] = 0
        alignedEpitopes = list()

        try:
            for brecord in blast_records:
                desc = str(brecord.query).split("|")
                typ = desc[1]
                nid = int(desc[2])
                if nid not in sample.neoantigensMap:
                    continue
                neo = sample.neoantigensMap[nid]
                for alignment in brecord.alignments:
                    desc = alignment.hit_def
                    epitope = self.__find_epitope(desc)
                    if '*' in neo.wtPeptide or "-1" in neo.wtPeptide or '*' in neo.mtPeptide or '*' in str(epitope.seq):
                        continue
                    aepitope = self.__add_aligned_epitope(neo, epitope)
                    if not aepitope is None:
                        alignedEpitopes.append(aepitope)
                        nc[typ] += 1
                        for hsp in alignment.hsps:
                            if typ == "WT":
                                aepitope.set_wildtype_score(hsp.score, [str(hsp.query), str(hsp.sbjct), str(hsp.match)])
                            else:
                                aepitope.set_mutant_score(hsp.score, [str(hsp.query), str(hsp.sbjct), str(hsp.match)])
        except ValueError:
            pass
        f.close()
        return alignedEpitopes

    def __read_sample_blast_table_format(self, sample, path, recompute_scores=True):
        '''
        Prepares neoantigens in the sample: reads precomputed Blast alignments
        in table format with columns
        ID        Epitope_ID        Epitope        Score_MT        Score_WT
        ID - peptide id, <chr>_<pos>_<ref>_<alt>_<mutated_position>

        :param sample: cfit.patient.Sample
            The sample for which the alignments are imported

        :param path: str
            path to the precomputed alignments, in the xml format

        :param recompute_scores: bool
            whether to recompute the alignment scores provided in the alignment files. If yes, the
            precomputed alignments serve only as a list of candidates.

        '''

        tab = pd.read_csv(path, sep="\t")
        # alignedEpitopes = []
        for line in tab.itertuples(index=False):
            nid = line.ID  # peptide_id
            if nid not in sample.peptideMap:  # sample.neoantigensMap:
                # self.logger("nid " + str(nid) + " not in the neoantigenMap!", 0)
                # self.logger(sample.name)
                continue
            neos = sample.peptideMap[line.ID]  # there can be multiple neoantigens, for different HLA alleles
            epitope = self.__find_epitope(line.Epitope_ID)
            mt_score = line.Score_MT
            wt_score = line.Score_WT
            for neo in neos:
                if '*' in neo.wtPeptide or "-1" in neo.wtPeptide or '*' in neo.mtPeptide or '*' in str(epitope.seq):
                    continue
                aepitope = self.__add_aligned_epitope(neo, epitope, mt_score=mt_score, wt_score=wt_score,
                                                      recompute_scores=recompute_scores)
                # if aepitope is not None:
                #    alignedEpitopes.append(aepitope)

    #        return alignedEpitopes

    def import_patient_data_alignments(self, alndir, patients, recompute_scores=True):
        '''

        :param alndir: str
            directory with neoantigen-epitope alignments (Neoantigens/)

        :param patients: list
            list of cfit.patient.Patient objects

        :return:
            Initializes the alignements of the object
        '''

        for patient in patients:
            pname = patient.pname
            alpath = os.path.join(alndir, "neoantigens_" + pname + "_iedb.xml")  # obsolete format
            if not os.path.exists(alpath):
                alpath = os.path.join(alndir, "alignment_" + pname + ".txt")
            for sample in patient.samples:
                if ".xml" in alpath:
                    self.__read_sample_blast_alignments_xml_format(sample, alpath)
                else:
                    self.__read_sample_blast_table_format(sample, alpath, recompute_scores=recompute_scores)
