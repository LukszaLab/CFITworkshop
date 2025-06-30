import json
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from cfit.fitness.HLAweights import HLAweights
from cfit.fitness.HWTHLAweights import HWTHLAweights
from cfit.fitness.ImmuneCloneFitness import ImmuneCloneFitness
from cfit.fitness.neo_quality.ARNeoantigenQuality import ARNeoantigenQuality
from cfit.tree.mutation.Neoantigen import Neoantigen
from cfit.util.Analysis import Analysis
from cfit.util.Log import Log
from cfit.util.Utils import Utils
from cfit.plot.PlotTreeAll import PlotTreeAll
if __name__ == "__main__":

    '''
    
    Run as python compute_ntau_AR.py -d <data_directory> -p -PDAC

    dir=$HOME/Dropbox/Research/10_PDAC/Modeling/Data_orig_fs
    cdir=$HOME/Workspace/CFIT/bin/
    python $cdir/compute_ntau_AR.py -d $dir -PDAC -kd_thr 500 -ns 9
    
    '''

    parser = Utils.make_cohort_parser()
    parser.add_argument("-quantile_ptab", type=float, default=0.5, help="quantile split for survival analysis")
    parser.add_argument("-ptab", type=str, default=None, help="path to the pvalues.txt file to choose optimal parameters")
    parser.add_argument("-pval_thr", type=float, default=0.01, help="threshold on significance in ptab")
    parser.add_argument("-a", type=float, default=None)
    parser.add_argument("-k", type=float, default=None)
    parser.add_argument("-tau", type=float, default=None)
    parser.add_argument("-OS_based_response", type=float, default=None, help="if set use as a threshold for response")
    parser.add_argument("-q", "--quantile", default=0.5, type=float)
    parser.add_argument("-o", "--odir", default=None)

    args = parser.parse_args()
    name = ""

    odir = os.path.join(args.dir, "ntau_AR_"+args.netMHC)
    if not os.path.exists(odir):
        os.mkdir(odir)

    a = args.a
    k = args.k
    tau = args.tau

    npos = "1-9"
    Neoantigen.WEPS = 0.0

    suffix = ""
    aggrfun = max
    if args.PDAC:
        Log.logger("Using PDAC parameters from 2017")
        if a is None:
            a = 26.
        if k is None:
            k = 1.
        Neoantigen.WEPS = 0.0
        if tau is None:
            tau = 0.04
        npos = "1-9"
    elif args.IO:
        Log.logger("Using IO parameters from 2017")
        if a is None:
            a = 26.
        if k is None:
            k = 4.86936
        Neoantigen.WEPS = 0.0003
        if tau is None:
            tau = 0.09
        npos = "1-9H"
    elif args.ptab is not None:
        Neoantigen.WEPS = 0.0
        npos = "1-9"
        Log.logger("Using optimized parameters from "+args.ptab)
        ptab = pd.read_csv(args.ptab, sep="\t")
        AUC = False
        if "score_OS" in list(ptab.columns):
            AUC=False
#            suffix = "_logrank"
#            if args.ptab_Q is not None:
#                if "Q" in ptab.columns:
#                    ptab = ptab[ptab.Q==args.ptab_Q]
        elif "AUC" in list(ptab.columns):
            AUC = True
#            suffix = "_AUC"
            if args.ptab_T is not None:
                if "T" in ptab.columns:
                    ptab = ptab[ptab.T==args.ptab_T]

        name = os.path.basename(args.ptab).replace("pvalues_","").replace("AUC_","").replace(".txt", "").replace("logrank_","")
        ptab = ptab[ptab.Q == args.quantile_ptab]
        print(ptab.shape)
        Utils.set_ptab_params(ptab, args.pval_thr, AUC=AUC, quantile=args.quantile_ptab, OS=(not args.PFS), PFS=args.PFS)

        a = Utils.params["optimal_params"]["a"]
        k = Utils.params["optimal_params"]["k"]
        tau = Utils.params["optimal_params"]["tau"]

    suffix = "_PFS" if args.PFS else "_OS"
    aggrfun = max

    modelParamValues = [a, k, tau]
    modelParamNames = ["a", "k", "tau"]
    params = pd.DataFrame([modelParamValues])
    params.columns = modelParamNames
    params.to_csv(os.path.join(odir, "parameters_"+suffix+".txt"), sep="\t", index=False)

    anl = Analysis()
    anl.set_MHC_version(args.netMHC)
    if args.mapping is None:
        mappingfile = os.path.join(args.dir, "mapping.json")
    else:
        mappingfile = args.mapping
    with open(mappingfile) as f:
        mappingjs = json.load(f)

    if args.OS_based_response is not None:
        mappingjs = Utils.fill_up_response(mapping=mappingjs, OS=args.OS_based_response)

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

    anl.model = "AR"
    anl.npos = npos
    if npos == "1-9H":
        HLAW = HWTHLAweights(npos)
    else:
        HLAW = HLAweights(npos)

    Qmodel = ARNeoantigenQuality(alndir=alndir, iedbfasta=iedbfasta)
    
    anl.set_neantigen_quality_model(Qmodel)
    # set fitness components
    fitnessModelComp = ImmuneCloneFitness(aggrfun=aggrfun)
    Utils.AGGR = "MAX"  # needed only for output, used in Utils.writeline

    anl.set_fitness_model_component(fitnessModelComp, "immune", 1.)

    # perform evolutionary predictions:
    # 1. compute neoantigen qualities
    anl.compute_neoantigen_sample_qualities(a=a, k=k, HLAW=HLAW, include_A=True, include_R=True,
                                            kdthr=args.kd_thr, KDnormalize=1.)

    # 2. compute fitness of all clones
    anl.compute_node_fitness(recompute_components=True)

    # write neoantigens with their quality componenst and frequencies from the top tree
    anl.write_neoantigen_fitness(os.path.join(odir, "Neontigen_qualities_"+suffix), exclusive=True, longitudinal=False,
                                 clonal=True)

    #anl.write_mutations(os.path.join(odir, "Mutations"))
    #anl.write_CCFs(os.path.join(odir, "Mutations"))

    otdir = os.path.join(odir, "Trees")
    if not os.path.exists(otdir):
        os.mkdir(otdir)
    for pat in anl.patients.values():
        pjson = pat.toJSON()
        with open(os.path.join(otdir, "trees_"+pat.name+".json"), 'w') as of:
            json.dump(pjson, of, indent=True)

    # Survival analysis
    pof_surv = open(os.path.join(odir, "pvalues_survival_"+suffix+".txt"), 'w')
    Utils.writelineHeader(pof_surv, modelParamNames=modelParamNames)

    quantile = 0.5
    res_surv = anl.classify_survival(beta=1., tau=tau, outdir=None, OS=True, PFS=True, quantile=quantile)
    modelParamValues = [a, k, tau, quantile]
    Utils.writeline(anl, "AR", res_surv, modelParamValues, Neoantigen.WEPS,
                        "AR", pof=pof_surv)
    anl.plot_survival(os.path.join(odir, "survival_plot_OS_q_"+str(quantile)+"_"+suffix+".pdf"), OS=True, pval=res_surv["pval_OS"], quantile=quantile)
    pof_surv.close()

    # AUC/ROC analysis
    response_values = set([el["response"] for el in mappingjs])
    perform_AUC = len(response_values) > 1
    print(set([el["response"] for el in mappingjs]))

    if perform_AUC:
        patient_class = {}
        if "CR" in response_values or "PR" in response_values:
            for el in mappingjs:
                patient_class[el["name"]] = 2 if el["response"] in ["CR", "PR"] else 1
        else:
            for el in mappingjs:
                patient_class[el["name"]] = el["response"]
        res_AUC = anl.classify_AUC(patient_class, beta=1, tau=tau, plot=True, ofile=os.path.join(odir, "ROC"+suffix+".pdf"))
        pof_AUC = open(os.path.join(odir, "pvalues_AUC"+suffix+".txt"), 'w')
        Utils.writelineHeader_AUC(pof_AUC, modelParamNames=modelParamNames)
        Utils.writeline_AUC(anl, "AR", res_AUC, modelParamValues, Neoantigen.WEPS,
                            "AR", pof=pof_AUC)
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

    outdir_trees = os.path.join(odir, 'Tree_plots')
    if not os.path.exists(outdir_trees):
        os.mkdir(outdir_trees)
    for patname in anl.patients:
        patient = anl.patients[patname]
        plottree = PlotTreeAll(patient, drivers=None, outdir=outdir_trees, tid=0, save_html=True, show=False)

