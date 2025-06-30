import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cfit.fitness.ImmuneCloneFitness import ImmuneCloneFitness
from cfit.fitness.DGCloneFitness import DGCloneFitness
from cfit.fitness.neo_quality.CARNeoantigenQuality import CARNeoantigenQuality
from cfit.util.Analysis import Analysis
from cfit.util.Log import Log
from cfit.util.Utils import Utils
from cfit.tree.mutation.Neoantigen import Neoantigen

if __name__ == "__main__":

    '''
    Run as python compute_ntau_CAR.py -d <data_directory>
    '''

    parser = Utils.make_cohort_parser()
#    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    #NQ model params - fixed
#    parser.add_argument("-quantile", type=float, default=0.5, help="quantile split for survival analysis")
    parser.add_argument("-quantile_ptab", type=float, default=0.5, help="quantile split for survival analysis")
    parser.add_argument("-a", type=float, default=None)
    parser.add_argument("-k", type=float, default=None)
    parser.add_argument("-sigmaI", type=float, default=1.)
    parser.add_argument("-sigmaP", type=float, default=0)
    parser.add_argument("-w", type=float, default=None)
    #or from training

    parser.add_argument("-ptab", type=str, default=None, help="path to the pvalues.txt file to choose optimal parameters")
    parser.add_argument("-pval_thr", type=float, default=0.01, help="threshold on significance in ptab")
    parser.add_argument("-exclude_R", action='store_true', help="")
    parser.add_argument("-just_R", action='store_true', help="")
    parser.add_argument("-ICGC", action='store_true')
    parser.add_argument("-just_immune", action='store_true')

    parser.add_argument("-o", "--odir", default=None)
    args = parser.parse_args()

    include_R = not args.exclude_R
    just_R = args.just_R

    if args.odir is None:
        odir = os.path.join(args.dir, "ntau_CAR_"+args.netMHC)
    else:
        odir = args.odir

    if not os.path.exists(odir):
        os.mkdir(odir)

    suffix = "_PFS" if args.PFS else "_OS"
    aggrfun = max
    name = ""
    sigmaI = args.sigmaI
    sigmaP = args.sigmaP

    if args.ptab is not None:
        Log.logger("Using optimized parameters from "+args.ptab)
        if "sum" in args.ptab:
            aggrfun = sum
        ptab = pd.read_csv(args.ptab, sep="\t")
        if "score_OS" in list(ptab.columns):
            AUC=False
#            suffix = "_logrank"
        elif "AUC" in list(ptab.columns):
            AUC = True
#            suffix = "_AUC"

        name = os.path.basename(args.ptab).replace("pvalues_","").replace("AUC_","").replace(".txt", "").replace("logrank_","")
        print(ptab.columns)
        ptab = ptab[ptab.Q == args.quantile_ptab]
        Utils.set_ptab_params(ptab, args.pval_thr, AUC=AUC, quantile=args.quantile_ptab, OS=(not args.PFS), PFS=args.PFS)
        a = Utils.params["optimal_params"]["a"]
        k = Utils.params["optimal_params"]["k"]
        w = Utils.params["optimal_params"]["w"]
        tau = Utils.params["optimal_params"]["tau"]
#        Q = Utils.params["optimal_params"]["Q"]
    elif args.ICGC:
            (a, k, tau, w) = (33.0, 1.0, 3.94421, 0.16)
            sigmaI = 1
            sigmaP = 0
    elif args.a is not None:
        a = args.a
        k = args.k
        sigmaI = args.sigmaI
        sigmaP = args.sigmaP
        w = args.w
    else:
        a = 22.9
        k = 1
        w = 0.22
        if not args.just_immune:
            sigmaI = 1.39
            sigmaP = 4.68
        tau = sigmaI + sigmaP
        sigma = sigmaI/tau

    modelParamValues = [a, k, tau, w]
    modelParamNames = ["a", "k", "tau", "w"]
    params = pd.DataFrame([modelParamValues])
    params.columns = modelParamNames
    params.to_csv(os.path.join(odir, "parameters_"+name+suffix+".txt"), sep="\t", index=False)

    anl = Analysis()
    anl.set_MHC_version(args.netMHC)

    if args.mapping is None:
        mappingfile = os.path.join(args.dir, "mapping.json")
    else:
        mappingfile = args.mapping
    with open(mappingfile) as f:
        mappingjs = json.load(f)

    mappingjs = Utils.fill_up_response(mappingjs, OS=36)

    if args.config is None:
        configfile = os.path.join(args.dir, "config.json")
    else:
        configfile = args.config
    with open(configfile) as f:
        configjs = json.load(f)


    vcfdir = os.path.join(args.dir, configjs["vcf_dir"])
    alndir = os.path.join(args.dir, configjs["aln_dir"])
    iedbfasta = None
    if "iedb_file" in configjs:
        iedbfasta = os.path.join(args.dir, configjs["iedb_file"])
    else:
        iedbfasta = os.path.join(alndir, "enemy.fasta")

    Log.logger("Initializing...")
    Log.logger("iedbfasta = "+iedbfasta)

    anl.clonal = 1  # whether to take the trees into account in the analysis, 1 - yes
    anl.ntrees = args.ntrees  # number of top scoring trees

    anl.initialize_config(configjs, mappingjs, args=args)
    anl.model = "CAR"

    Qmodel = CARNeoantigenQuality(alndir=alndir, iedbfasta=iedbfasta,
                                  ep_dist_model_name=args.ep_dist_model_name
                                )

    anl.set_neantigen_quality_model(Qmodel)
    # set fitness components
    fitnessModelComp1 = ImmuneCloneFitness(aggrfun=aggrfun)
    #fitnessModelComp2 = DGCloneFitness()

    anl.set_fitness_model_component(fitnessModelComp1, "immune", 1)
    #anl.set_fitness_model_component(fitnessModelComp2, "dg", 1.-sigma)


    Utils.AGGR = "MAX"  # legacy of CancerIO, needed only for output, used in Utils.writeline

    # compute neoantigen qualities
    anl.compute_neoantigen_sample_qualities(a=a, k=k, w=w,
                                            kdthr=args.kd_thr,
                                            include_R=include_R, just_R=just_R)
    # compute fitness of all clones
    anl.compute_node_fitness(recompute_components=True)

    # write neoantigens with their quality componenst and frequencies from the top tree
    anl.write_neoantigen_fitness(os.path.join(odir, "Neoantigen_qualities_"+name+suffix),
                                 exclusive=True,
##                                 longitudinal=False,
                                 longitudinal=True,
                                 clonal=True)

    anl.write_mutations(os.path.join(odir, "Mutations"))
    anl.write_CCFs(os.path.join(odir, "Mutations"))


    otdir = os.path.join(odir, "Trees")
    if not os.path.exists(otdir):
        os.mkdir(otdir)
    for pat in anl.patients.values():
        pjson = pat.toJSON()
        with open(os.path.join(otdir, "trees_"+pat.name+".json"), 'w') as of:
            json.dump(pjson, of, indent=True)


    # survival analysis
    pof_surv = open(os.path.join(odir, "pvalues_survival_"+name+suffix+".txt"), 'w')
    Utils.writelineHeader(pof_surv, modelParamNames=modelParamNames)

#    for quantile in np.arange(0.1, 1.,0.1):
#        quantile = round(quantile,1)
    quantile = 0.5
    res_surv = anl.classify_survival(beta=1., tau=tau, outdir=None, OS=True, PFS=True, quantile=quantile)
    modelParamValues = [a, k, tau, w, quantile]
    Utils.writeline(anl, "CAR", res_surv, modelParamValues, Neoantigen.WEPS,
                        "CAR", pof=pof_surv)
    anl.plot_survival(os.path.join(odir, "survival_plot_OS_q_"+str(quantile)+"_"+name+suffix+".pdf"), OS=True, pval=res_surv["pval_OS"], quantile=quantile)
    anl.plot_survival(os.path.join(odir, "survival_plot_PFS_q_"+str(quantile)+"_"+name+suffix+".pdf"), OS=False, pval=res_surv["pval_PFS"],quantile=quantile)
    pof_surv.close()
#    anl.plot_survival(os.path.join(odir, "survival_plot_OS_"+name+suffix+".pdf"), OS=True, pval=res_surv["pval_OS"], quantile=args.quantile)
#    anl.plot_survival(os.path.join(odir, "survival_plot_PFS_"+name+suffix+".pdf"), OS=False, pval=res_surv["pval_PFS"],quantile=args.quantile)

    # AUC/ROC analysis
    perform_AUC = len(set([el["response"] for el in mappingjs])) > 1

    if perform_AUC:
        patient_class = {}
        for el in mappingjs:
            patient_class[el["name"]] = 2 if el["response"] in ["CR", "PR"] else 1

        res_AUC = anl.classify_AUC(patient_class, beta=1, tau=tau, plot=True, ofile=os.path.join(odir, "ROC"+suffix+".pdf"))
        pof_AUC = open(os.path.join(odir, "pvalues_AUC"+suffix+".txt"), 'w')
        Utils.writelineHeader_AUC(pof_AUC, modelParamNames=modelParamNames)
        Utils.writeline_AUC(anl, "CAR", res_AUC, modelParamValues, Neoantigen.WEPS,
                            "CAR", pof=pof_AUC)
        pof_AUC.close()

    sample_statistics = anl.sample_statistics(kd_thr=args.kd_thr)
    fitness_statistics = anl.fitness_sample_statistics()
    sample_statistics.to_csv(os.path.join(odir, "sample_statistics.txt"), sep="\t", index=False)
    fitness_statistics.to_csv(os.path.join(odir, "fitness_sample_statistics.txt"), sep="\t", index=False)

    ntaus = [np.exp(pat.q) for pat in anl.patients.values()]
    ntaus.sort()
    plt.hist(ntaus,bins=np.arange(0,1.01,0.05))
    plt.savefig(os.path.join(odir, "ntau_hist.pdf"))
    plt.close()

    plt.scatter(range(len(ntaus)), ntaus)
    plt.savefig(os.path.join(odir, "ntau_all.pdf"))
    plt.close()

    avefit = list(fitness_statistics.ave_fitness)
    avefit.sort()
    plt.hist(avefit)
    plt.savefig(os.path.join(odir, "avefit_hist.pdf"))
    plt.close()

    plt.scatter(range(len(avefit)), avefit)
    plt.savefig(os.path.join(odir, "avefit_all.pdf"))
    plt.close()

