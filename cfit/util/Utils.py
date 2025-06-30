'''
Created on Mar 30, 2015

@author: mluksza
'''

import numpy as np
import argparse
from cfit.util.Log import Log



class Utils(object):
    '''
    classdocs
    '''
    MATCH_ALLELES = False
    #    KS=[0.414214, 0.534511, 0.668179, 0.820679, 1.,1.2185, 1.49661, 1.87087, 2.41421, 4.9, 3.29655, 5.02734,10.1531]
    #    AS=range(15,35)
    KS = [0.414214, 0.534511, 0.668179, 0.820679, 1., 1.2185, 1.49661, 1.87087, 2.41421, 4.9, 3.29655, 5.02734, 10.1531]
    #    AS=range(15,40)
    AS = list(np.linspace(15, 40, (40 - 15 + 1) * 2 - 1))

    TAUS = [0.000001]
    tt = range(1, 10)
    TAUS += list(map(lambda i: float(i) / 1000., tt))
    TAUS += list(map(lambda i: float(i) / 100., tt))
    TAUS += list(map(lambda i: float(i) / 100, range(11, 20)))
    TAUS += list(map(lambda i: float(i) / 10., tt))
    TAUS += list(map(lambda i: float(i), tt))
    TAUS += list(map(lambda i: float(i) * 10, tt))

    method = "SIG"  # EXP|LIN
    INF = float("inf")
    a = 23.0
    k = 1.0

    HYDROPHOBIC_RESIDUES = "AILMFWYV"
    WEIRD_RESIDUES = "CGP"

    AGGR = "MAX"
    KDnormalize = 1
    AGGRnum = float("inf")
    PTAB = None

    @staticmethod
    def make_cohort_parser():
        Log.logger("make_parser")
        parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument("-d", "--dir", help="data directory")
        parser.add_argument("-config", default=None, help="config.json file path")
        parser.add_argument("-mapping", default=None, help="mapping.json file path")
        parser.add_argument("-PDAC", action='store_true', help="use PDAC default parameters")
        parser.add_argument("-IO", action='store_true', help="use IO default parameters")
        parser.add_argument("-PFS", action='store_true', help="use PFS instead of OS")
        parser.add_argument("-netMHC", default="40", help="netMHC version, 34/40/pan41")
        parser.add_argument("-kd_thr", type=float, default=None)
        parser.add_argument("-ntrees", type=int, default=5)
        parser.add_argument("-ep_dist_model_name", default="all_tcr_all_combos_model")
        parser.add_argument('-ns', default='9', help="peptide lengths, e.g '9', '8,9', '8-14")
        parser.add_argument("-tf", "--tree_format", default='phylowgs', help="format of trees, currently handling 'phylowgs' and 'pairtree'")
        return parser

    @staticmethod
    def fill_up_response(mapping, OS=None, PFS=None):
        if OS is not None:
            for el in mapping:
                if el["OS"] < OS:
                    el["response"] = "SD"
                else:
                    el["response"] = "CR"
        elif PFS is not None:
            for el in mapping:
                if el["PFS"] < PFS:
                    el["response"] = "SD"
                else:
                    el["response"] = "CR"
        return mapping

    @staticmethod
    def set_ptab_params(ptab, pvalue_thr, AUC=False, use_max=False, quantile=None, OS=True, PFS=False):

        #if no_w:  # only A
        #    ptab = ptab[ptab.w == 0]
        #if just_w:  # only D
        #    ptab = ptab[ptab.w == 1]
        #if no_DG:  # no driver genes
        #    ptab = ptab[ptab.sigma == 1]
        #if just_DG:
        #    ptab = ptab[ptab.sigma == 0]

        cols = list(ptab.columns)
        cparams = cols[5:cols.index('criterion')]

        if AUC:
            scores = ptab.AUC
            pvals = ptab.MW_pval


        #            ws = aucs*500
        #            ws = np.array(Utils.log_norm(ws))
        #            ws = np.exp(ws)
        #            q = np.quantile(ws, 0.9995)
        #            ind = [i for i, w in enumerate(ws) if w >= q]
        #            relevant_scores = [s for (s, w) in zip(aucs, ws) if w >= q]
        #            ws = [w for w in ws if w >=q]
        #            logrank_score = sum([s * w for s, w in zip(relevant_scores, ws)])
        else:
            if quantile is not None:
                ptab = ptab[ptab["Q"]==quantile]
            if OS:
                scores = ptab.score_OS
                pvals = ptab.pvalue_OS
            elif PFS:
                scores = ptab.score_PFS
                pvals = ptab.pvalue_PFS
        #        ws = [s if s>0 and p <= pvalue_thr else -Utils.INF for (s,p) in zip(scores, pvals)]
        Log.logger("ptab shape = "+str(ptab.shape))
        ws = [s for (s, p) in zip(scores, pvals) if s > 0 and p <= pvalue_thr]
        if use_max or len(ws) == 0:
            maxscore = max(scores)
            ws = [1]
            Log.logger("max score=" + str(maxscore) + " len(scores)=" + str(len(scores)))
            ind = [[i for i, (s, p) in enumerate(zip(scores, pvals)) if s >= maxscore][0]]
            relevant_scores = [[s for (s, p) in zip(scores, pvals) if s >= maxscore][0]]
        else:
            ind = [i for i, (s, p) in enumerate(zip(scores, pvals)) if s > 0 and p <= pvalue_thr]
            relevant_scores = [s for (s, p) in zip(scores, pvals) if s > 0 and p <= pvalue_thr]
            ws = np.array(Utils.log_norm(ws))
            ws = np.exp(ws)

        opt_score = sum([s * w for s, w in zip(relevant_scores, ws)])
        Log.logger('optimal score = ' + str(opt_score))
        # avearage over the relevant log rank score landscape
        optvalues = np.dot(ws, ptab.loc[:, cparams].iloc[ind, :])

        dopt = {}
        Log.logger('opt parameters: ')
        for par, val in zip(cparams, optvalues):
            dopt[par] = val
            Log.logger("\t"+par+"="+str(val))
        Utils.params = {"weights": ws,
                        "params": cparams,
                        "values": ptab.loc[:, cparams].iloc[ind, :],
                        "optimal_params": dopt,
                        "optimal_score": opt_score}

    @staticmethod
    def residueChangeClass(res1, res2):
        code = ""
        if res1 in Utils.HYDROPHOBIC_RESIDUES:
            code += "H"
        elif res2 in Utils.WEIRD_RESIDUES:
            code += "W"
        else:
            code += "N"
        if res2 in Utils.HYDROPHOBIC_RESIDUES:
            code += "H"
        elif res2 in Utils.WEIRD_RESIDUES:
            code += "W"
        else:
            code += "N"
        return code

    @staticmethod
    def log_sum2(v1, v2):
        '''

        '''
        ma = np.maximum(v1, v2)
        if ma == -Utils.INF:
            return -Utils.INF
        return ma + np.log(np.exp(v1 - ma) + np.exp(v2 - ma))

    @staticmethod
    def log_sum(v):
        '''

        :param v: np.array
        '''
        v = np.array(v)
        ma = np.max(v)
        if ma == -Utils.INF:
            return -Utils.INF
        res = np.log(sum(np.exp(v-ma)))+ma
        return res

    @staticmethod
    def log_norm(v):
        '''
        :param v: list
        :return: list
        '''
        logZ = Utils.log_sum(v)
        
        return list(map(lambda x: x - logZ, v))

    @staticmethod
    def fagsig(x):
        if Utils.a == -100:
            return 1
        e0 = float(Utils.k)
        d0 = float(Utils.a)
        v1 = (x - d0) / (e0)
        ma = np.maximum(-v1, 0)
        res = ma + np.log(np.exp(0 - ma) + np.exp(-v1 - ma))
        v = np.exp(-res)
        return v

    @staticmethod
    def fagexp(x):
        d0 = float(Utils.a)
        v = -(x / (0.5 * d0))
        return -np.exp(v)

    @staticmethod
    def faglin(x):
        return max((x - Utils.a) / 30, 0)

    @staticmethod
    def fagllin(x):
        return x

    @staticmethod
    def fag(x):
        if Utils.method == "SIG":
            return Utils.fagsig(x)
        elif Utils.method == "EXP":
            return Utils.fagexp(x)
        elif Utils.method == "LIN":
            return Utils.faglin(x)
        else:
            return Utils.fagllin(x)

    @staticmethod
    def zlog(v):
        if v <= 0:
            return -Utils.INF
        return np.log(v)

    @staticmethod
    def writelineHeader(pof, modelParamNames=["a", "k", "tau"], othercolumns=None):
        header = ["model", "positions", "clonality", "aggr", "mask"] + modelParamNames + ["criterion"]
        header += ["score_OS", "pvalue_OS", "BIC_OS", "AIC_OS",
                   "score_PFS", "pvalue_PFS", "BIC_PFS", "AIC_PFS",
                   "AUC_score", "AUC_pvalue",
                   "quantile", "neps", "N"]
        if othercolumns is not None:
            header += othercolumns
        header = "\t".join(header) + "\n"
        pof.write(header)

    @staticmethod
    def writeline(anl, name2, res, modelParamValues, eps, model, pof, patient_names=None, otherparams=None):
        '''

        :param anl: Analysis

        :param name2: str

        :param res: list

        :param modelParamValues: list

        :param eps: float

        :param model: str

        :param pof: file handle

        :param patient_names:

        :param otherparams:

        :return:
        '''
        if patient_names is None:
            patient_names = anl.patients.keys()
        params = [anl.npos, anl.clonal, Utils.AGGR, name2] + modelParamValues + [anl.criterion]
        try:
            AUC_score = sum(res["AUC"][0])
            AUC_pval = sum(res["AUC"][1])
        except:
            AUC_score = 0.0
            AUC_pval = 0.0

#        n = modelParamValues[-2]
        n = len(patient_names)
        k = len(modelParamValues)
        #k = modelParamValues[-1]


        BIC_OS = k * np.log(n) - 2 * res["score_OS"]
        BIC_PFS = k * np.log(n) - 2 * res["score_PFS"]

        AIC_OS = 2 * k - 2 * res["score_OS"]
        AIC_PFS = 2 * k - 2 * res["score_PFS"]

        line = [model] + params + [res["score_OS"], res["pval_OS"], BIC_OS, AIC_OS] \
               + [res["score_PFS"], res["pval_PFS"], BIC_PFS, AIC_PFS, AUC_score, AUC_pval] \
               + [anl.quantile, eps, len(patient_names)]
        if otherparams is not None:
            line += otherparams
        line = "\t".join(list(map(lambda el: str(el), line))) + "\n"
        pof.write(line)

    @staticmethod
    def writelineHeader_AUC(pof, modelParamNames=["a", "k", "tau"], othercolumns=None):
        '''

        :param pof: file handle open for writing

        :param modelParamNames: list

        :param othercolumns: list

        '''
        header = ["model", "positions", "clonality", "aggr", "mask"] + modelParamNames + ["criterion"]
        header += ["AUC", "MW_pval", "mean_NR", "mean_R", "neps", "N"]
        if othercolumns is not None:
            header += othercolumns
        header = "\t".join(header) + "\n"
        pof.write(header)

    @staticmethod
    def writeline_AUC(pr, name2, res, modelParamValues, eps, model, pof, patient_names=None, otherparams=None):
        '''

        :param pr:

        :param name2: str

        :param res: list

        :param modelParamValues: dict

        :param eps: float

        :param model: name

        :param pof: file handle

        :param patient_names: list

        :param otherparams: list

        '''
        if patient_names is None:
            patient_names = pr.patients.keys()
        params = [pr.npos, pr.clonal, Utils.AGGR, name2] + modelParamValues + [pr.criterion]
        line = [model] + params + [res["AUC"], res["MW_pval"], res["mean_NR"], res["mean_R"]] \
               + [eps, len(patient_names)]
        if otherparams is not None:
            line += otherparams
        line = "\t".join(list(map(lambda el: str(el), line))) + "\n"
        pof.write(line)

