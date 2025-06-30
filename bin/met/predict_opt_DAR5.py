''''
Run as (examples):

hdir="/Users/mluksza/Dropbox/Research/10_PDAC/Modeling/Data_multi"
cdir="/Users/mluksza/Workspace/CIO/bin"
odir = $hdir/Results_03_12

python ${cdir}/predict.py -d $hdir -PDAC -tp_pref1 Primary -tp_pref2 Met -model AR -cleps 0.03 -o $odir
python ${cdir}/predict.py -d $hdir -PDAC -tp_pref1 Primary -tp_pref2 Met -model AR -cleps 0.03 -include_nested -o $odir
python ${cdir}/predict.py -d $hdir -PDAC -tp_pref1 Primary -tp_pref2 Met -model AR -cleps 0.0 -o $odir

python ${cdir}/predict.py -d $hdir -PDAC -tp_pref1 Primary -tp_pref2 Met -model DAR -cleps 0.03 -o $odir
python ${cdir}/predict.py -d $hdir -PDAC -tp_pref1 Primary -tp_pref2 Met -model DAR -cleps 0.03 -include_nested -o $odir
python ${cdir}/predict.py -d $hdir -PDAC -tp_pref1 Primary -tp_pref2 Met -model DAR -cleps 0.0 -include_A 0 -include_R 0 -o $odir

'''

import argparse
import json
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy.stats import spearmanr

from cio.fitness.DGCloneFitness import DGCloneFitness
from cio.fitness.ImmuneCloneFitness import ImmuneCloneFitness
from cio.fitness.neo_quality.DAR5NeoantigenQuality import DAR5NeoantigenQuality
from cio.patient.Sample import Sample
from cio.tree.SampleTree import SampleTree
from cio.tree.mutation.Neoantigen import Neoantigen
from cio.util.Analysis import Analysis
from cio.util.Utils import Utils

def plot_regr(c_x, c_y, xlab="X axis", ylab="Y axis", ofile="", labels=None, fig=None, ax=None):
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import pearsonr, linregress

    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    fig.subplots_adjust(right=0.95, left=0.15)

    c_pearson, c_pval = pearsonr(c_x, c_y)
    c_spearman, c_spval = spearmanr(c_x, c_y)
    if labels is None:
        ax.plot(c_x, c_y, 'k.')
    else:
        ulabs = set(labels)
        for ulab in ulabs:
            dat = [x for x in zip(c_x, c_y, labels) if x[2] == ulab]
            c_x1 = [x[0] for x in dat]
            c_y1 = [x[1] for x in dat]
            ax.plot(c_x1, c_y1, '.')
    slope_fit, intercept_fit, r_value_fit, p_value_fit, std_err_fit = linregress(np.array([c_x, c_y]).T)

    x_ax_arr = np.array([min([min(c_x) + min(c_x) / 50]), max(c_x) + max(c_x) / 50])
    y_ax_arr = np.array([max(c_y) + max(c_y) / 50])
    ax.plot(x_ax_arr, np.array([intercept_fit, intercept_fit]) + slope_fit * x_ax_arr, 'r',
            label='Regression fit: y = %.2ex+%.2e (R^2 = %.2f)\nPearson coeff: %.2f (pval = %.2e),\nSpearman coeff: %.2f (pval = %.2e)' % (
            slope_fit, intercept_fit, r_value_fit ** 2, c_pearson, c_pval, c_spearman, c_spval))

    rang = (max(c_x) - min(c_x)) / 5
    ux = list(set(c_x))
    ux.sort()
    yb = [[x, [y for x0, y in zip(c_x, c_y) if x - rang <= x0 <= x + rang]] for x in ux]
    yb = [el for el in yb if len(el[1]) > 0]
    xb = [el[0] for el in yb]
    yb = [np.mean(el[1]) for el in yb]
    ax.step(xb, yb)
    ax.legend(fontsize=10, frameon=False, loc=2, bbox_to_anchor=(0, 1.20))

    # plt.ylim(y_ax_arr)
    # plt.xlim(x_ax_arr)
    ax.set_xlabel(xlab, fontsize=14)
    ax.set_ylabel(ylab, fontsize=14)
    if ofile != "":
        plt.savefig(ofile)
    plt.close()


def plot_average_fitness(tumor_stats, odir1):
    cols = ["Immune_fitness_cost1", "Immune_fitness_cost2", "New_clones_immune_fitness_cost",
            "Private_volume", "IG"]
    for col in cols:
        ax = sns.swarmplot(x='Cohort', y=col, data=tumor_stats, dodge=True)
        sns.boxplot(x='Cohort', y=col, data=tumor_stats,
                    boxprops={'facecolor': 'None'})
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2], frameon=False)
        x1 = tumor_stats[tumor_stats.Cohort == 'STS'][col]
        x2 = tumor_stats[tumor_stats.Cohort != 'STS'][col]
        try:
            pval = scipy.stats.mannwhitneyu(x1, x2, alternative='two-sided').pvalue
        except:
            pval = 0.0
        ax.set_title("pvalue=" + str(pval))
        plt.savefig(os.path.join(odir1, col + ".pdf"))
        plt.close()

    for col in cols:
        plot_regr(tumor_stats.OS, tumor_stats[col], xlab="OS", ylab=col,
                  ofile=os.path.join(odir1, "OS_" + col + ".pdf"))
        plot_regr(tumor_stats.PFS, tumor_stats[col], xlab="PFS", ylab=col,
                  ofile=os.path.join(odir1, "PFS_" + col + ".pdf"))

def get_tree_mutation_frequency(tree, mid):
    '''
    Compute mutated gene frequency in the sample, by averaging over the trees.
    :param mid: str
            mutation identifier, <chr>_<pos>_<ref>_<alt>
    :param beta: float
            weight parameter
    :return: float
        the total tumor frequency at which the mutation is present in the tumor
    '''

    yy = [node.Y for node in tree.nodes.values() if tree.mutation_node_index[mid][node.id]]
    if len(yy) == 0:
        return 0
    else:
        return sum(yy)


def report_y(tpoint1, tpoint2, patient, odir1):

    trees1 = tpoint1.trees(num=0)
    trees2 = tpoint2.trees(num=0)
    nids = list(trees1[0].nodes.keys())

    ytab = [[np.mean([tree.nodes[nid].Y for tree in trees1]),
             np.mean([tree.nodes[nid].Y for tree in trees2]),
             np.mean([tree.nodes[nid].get_fitness() for tree in trees2])] for nid in nids]
    ytab = pd.DataFrame(ytab)
    tp1, tp2 = tpoint1.name, tpoint2.name
    ytab.columns = [tp1, tp2, "fitness"]
    ytab["Private"] = [y2 if y1 < 0.03 else 0 for y1, y2 in zip(ytab[tp1], ytab[tp2])]
    z = sum(ytab["Private"])
    if z > 0:
        ytab.Private = list(np.array(list(ytab.Private)) / z)
    ytab.to_csv(os.path.join(odir1, patient.name + "_y_" + tp1 + "_" + tp2 + ".txt"), sep="\t", index=False)


def compute_pair_stats(patient, tp1, tp2):
    tpoint1 = patient.timePoints[tp1]
    tpoint2 = patient.timePoints[tp2]

    # compute average immune fitness of tumors
    ave_immune_fit1 = tpoint1.fitness()  # immune fitness of the Primary tumor
    ave_immune_fit2 = tpoint2.fitness()  # immune fitness of the metastatic tumor
    ave_immune_new_fit2 = tpoint2.fitness(private=True)  # fitness of new clone in the met
    n_new_clones = tpoint2.number_of_new_clones()

    # iterate over all time-point pairs

    # TMB related numbers - 19.07.2021
    n_nsyn = tpoint1.average_over_samples(Sample.TMB_nsyn, filter_on_Taf=False)
    n_syn = tpoint1.average_over_samples(Sample.TMB_syn, filter_on_Taf=False)

    n_MHC = tpoint1.average_over_samples(Sample.TMB_nsyn, filter_on_Taf=False)

    n_nsyn1 = tpoint1.average_over_samples(Sample.TMB_nsyn, filter_on_Taf=True)
    n_nsyn2 = tpoint2.average_over_samples(Sample.TMB_nsyn, filter_on_Taf=True)

    n_syn1 = tpoint1.average_over_samples(Sample.TMB_syn, filter_on_Taf=True)
    n_syn2 = tpoint2.average_over_samples(Sample.TMB_syn, filter_on_Taf=True)

    n_MHC1 = tpoint1.average_over_samples(Sample.TMB_MHC, filter_on_Taf=True)
    n_MHC2 = tpoint2.average_over_samples(Sample.TMB_MHC, filter_on_Taf=True)

    n_nsyn_eff1 = tpoint1.average_over_samples(Sample.effective_TMB_nsyn, beta=1)
    n_nsyn_eff2 = tpoint2.average_over_samples(Sample.effective_TMB_nsyn, beta=1)
    n_syn_eff1 = tpoint1.average_over_samples(Sample.effective_TMB_syn, beta=1)
    n_syn_eff2 = tpoint2.average_over_samples(Sample.effective_TMB_syn, beta=1)

    n_MHC_eff1 = tpoint1.average_over_samples(Sample.effective_TMB_MHC, beta=1)
    n_MHC_eff2 = tpoint2.average_over_samples(Sample.effective_TMB_MHC, beta=1)

    # tmb = tpoint1.average_over_samples(lambda sample: len([mut for mut in sample.mutations.values()]))
    tmb = tpoint1.average_over_samples(Sample.TMB, filter_on_Taf=False)
    tmb1 = tpoint1.average_over_samples(Sample.TMB, filter_on_Taf=True)
    tmb2 = tpoint2.average_over_samples(Sample.TMB, filter_on_Taf=True)
    tnb = tpoint1.average_over_samples(Sample.presentation_score, kd0=500, strict=True, filter_on_Taf=False)
    tnb1 = tpoint1.average_over_samples(Sample.presentation_score, kd0=500, strict=True, filter_on_Taf=True)
    tnb2 = tpoint2.average_over_samples(Sample.presentation_score, kd0=500, strict=True, filter_on_Taf=True)

    tmb_eff1 = tpoint1.average_over_samples(Sample.effective_TMB, beta=1)
    tmb_eff2 = tpoint2.average_over_samples(Sample.effective_TMB, beta=1)
    tnb_eff1 = tpoint1.average_over_samples(Sample.effective_presentation_score, kd0=500, strict=True, beta=1)
    tnb_eff2 = tpoint2.average_over_samples(Sample.effective_presentation_score, kd0=500, strict=True, beta=1)

    mode = 'one'
    g_stat1 = tpoint1.average_over_samples(Sample.neo2mut_ratio, mode=mode, nsyn=True, kd0=500.,
                                           strict=True, filter_on_Taf=True)
    g_stat2 = tpoint2.average_over_samples(Sample.neo2mut_ratio, mode=mode, nsyn=True, kd0=500.,
                                           strict=True, filter_on_Taf=True)
    g_stat_eff1 = tpoint1.average_over_samples(Sample.effective_neo2mut_ratio, mode=mode, nsyn=True, kd0=500.,
                                               strict=True)
    g_stat_eff2 = tpoint2.average_over_samples(Sample.effective_neo2mut_ratio, mode=mode, nsyn=True, kd0=500.,
                                               strict=True)

    mode = 'all'
    gall_stat1 = tpoint1.average_over_samples(Sample.neo2mut_ratio, mode=mode, nsyn=True, kd0=500.,
                                              strict=True, filter_on_Taf=True)
    gall_stat2 = tpoint2.average_over_samples(Sample.neo2mut_ratio, mode=mode, nsyn=True, kd0=500.,
                                              strict=True, filter_on_Taf=True)
    gall_stat_eff1 = tpoint1.average_over_samples(Sample.effective_neo2mut_ratio, mode=mode, nsyn=True, kd0=500.,
                                                  strict=True)
    gall_stat_eff2 = tpoint2.average_over_samples(Sample.effective_neo2mut_ratio, mode=mode, nsyn=True, kd0=500.,
                                                  strict=True)

    ###
    heterogen1 = tpoint1.average_over_samples(Sample.entropy, beta=1)
    heterogen2 = tpoint2.average_over_samples(Sample.entropy, beta=1)
    try:
        heterogen_sh1 = tpoint1.average_over_samples(Sample.entropy, beta=1, shared_only=True)
    except:
        heterogen_sh1 = -1
    try:
        heterogen_sh2 = tpoint2.average_over_samples(Sample.entropy, beta=1, shared_only=True)
    except:
        heterogen_sh2 = -1

    private_volume = tpoint2.average_over_sample_trees(SampleTree.get_private_volume, beta=1)
    shared_volume = tpoint2.average_over_sample_trees(SampleTree.get_shared_volume, beta=1)

    # "trunkality" measure of the tumors - \sum_\alpha x_alpha(p)*y_alpha(met)
    # computed over shared clones only, excluding new clones from the met
    trunkality_shared1, trunkality_shared2 = patient.trunkality(tp1, tp2, beta=1., use_shared=True, eps=0.03)
    # computed over all clones
    trunkality_total1, trunkality_total2 = patient.trunkality(tp1, tp2, beta=1., use_shared=False)

    kl1, kl2, wd = patient.observed_distance(tp1, tp2, eps=0.0, beta=1., normalize=False)
    dstats = {
        "TMB": tmb,
        "TMB1": tmb1,
        "TMB2": tmb2,
        "TMB_eff1": tmb_eff1,
        "TMB_eff2": tmb_eff2,
        "TMB_MHC": n_MHC,
        "TMB_MHC1": n_MHC1,
        "TMB_MHC2": n_MHC2,
        "TMB_MHC_eff1": n_MHC_eff1,
        "TMB_MHC_eff2": n_MHC_eff2,
        "TMB_nsyn": n_nsyn,
        "TMB_nsyn1": n_nsyn1,
        "TMB_nsyn2": n_nsyn2,
        "TMB_syn": n_syn,
        "TMB_syn1": n_syn1,
        "TMB_syn2": n_syn2,  # separated in mutation classes, absolute
        "TMB_nsyn_eff1": n_nsyn_eff1,
        "TMB_nsyn_eff2": n_nsyn_eff2,
        "TMB_syn_eff1": n_syn_eff1,
        "TMB_syn_eff2": n_syn_eff2,  # and averaged over clones
        "TNB": tnb,
        "TNB1": tnb1,
        "TNB2": tnb2,
        "TNB_eff1": tnb_eff1,
        "TNB_eff2": tnb_eff2,
#        "g_stat1": g_stat1,
#        "g_stat2": g_stat2,
#        "g_stat_eff1": g_stat_eff1,
#        "g_stat_eff2": g_stat_eff2,
#        "gall_stat1": gall_stat1,
#        "gall_stat2": gall_stat2,
#        "gall_stat_eff1": gall_stat_eff1,
#        "gall_stat_eff2": gall_stat_eff2,
        "Immune_fitness_cost1": -ave_immune_fit1,
        "Immune_fitness_cost2": -ave_immune_fit2,
        "New_clones_immune_fitness_cost": -ave_immune_new_fit2,
        "Delta_immune_fitness": ave_immune_fit2 - ave_immune_fit1,
        "Trunkality": trunkality_total2,
        "Trunkality_shared_clones": trunkality_shared2,
        "New_clones": n_new_clones,
        "Private_volume": private_volume,
        "Shared_volume": shared_volume,
        "Entropy1": heterogen1,
        "Entropy2": heterogen2,
        "Entropy_shared1": heterogen_sh1,
        "Entropy_shared2": heterogen_sh2,
        "Observed_distance_KL1": kl1,
        "Observed_distance_KL2": kl2,
        "Observed_distance_WD": wd
    }
    return dstats


def compute_prediction_stats(patient, tp1, tp2, cleps, include_nested, odir1, sample_size_corr=1, args=None):
    '''

    :param patient: Patient
    :param tp1: str
    :param tp2: str
    :param cleps: float
    :param include_nested: bool
    :param odir1: str
    :param sample_size_corr: float
    :param args: dict
    :return:
    '''


    tpoint1 = patient.timePoints[tp1]
    tpoint2 = patient.timePoints[tp2]

#    trained_sigma = trained_params["sigma"]
#    trained_tau = trained_params["tau"]
#    trained_params = {"tau": trained_tau, "weights": {"immune": trained_sigma, "dg": 1-trained_sigma}}

    # EFFECTIVE population size for computing BIC and model comparison
    neff2 = tpoint2.average_over_samples(Sample.effective_sample_size, eps=0.0)
    norm_neff2 = neff2/sample_size_corr

    # exluding mutations not detected on the recurrent tumor sample
    neff_alt2 = tpoint2.average_over_samples(Sample.effective_sample_size, eps=0.0, filter_on_Taf=True)
    norm_neff_alt2 = neff_alt2 / sample_size_corr
    wnames = ["dg", "immune"]
    # Optimized predictions
    inparams = {"tau": None, "weights": {"dg": None, "immune": None}}
    nparams = [2, 1, 1, 0]

    if args is not None:
        if args.just_DG:
            inparams = {"tau": None, "weights": {"dg": 1, "immune": 0}}
            nparams = [1, 1, 1, 0]
        elif args.no_DG:
            inparams = {"tau": None, "weights": {"dg": 0, "immune": 1}}
            nparams = [1, 1, 1, 0]

    optparams, optdi = patient.optimized_prediction_distance(tp1, tp2, cleps, inparams, beta=1.,
                                                             include_nested=include_nested)
    #partial models
    inparams1 = {"tau": None, "weights": {"dg": 1., "immune": 0.}}
    optparams1, optdi1 = patient.optimized_prediction_distance(tp1, tp2, cleps, inparams1, beta=1.,
                                                               include_nested=include_nested)
    inparams2 = {"tau": None, "weights": {"dg": 0., "immune": 1.}}
    optparams2, optdi2 = patient.optimized_prediction_distance(tp1, tp2, cleps, inparams2, beta=1.,
                                                               include_nested=include_nested)

    zerodi = patient.distribution_distance(tp1, tp2, cleps, {'dg': 0, 'immune': 0}, 0.0,
                                           reference_time_point=2, beta=1., include_nested=include_nested,
                                           measure="KL")

    rec_difs = [optdi, optdi1, optdi2,  zerodi] #traineddi, opttraineddi,

    names = ["Optimized", "Optimized_immune", "Optimized_dg",  "Neutral"] #"Trained", "Trained_optimized",

    neffs = [neff2, norm_neff2, neff_alt2, norm_neff_alt2]
    neff_names = ["neff2", "norm_neff2", "neff_alt2", "norm_neff_alt2"]
    dcorr = {}
    for (neff_name, neff) in zip(neff_names, neffs):
        for (name, di, n) in zip(names, rec_difs, nparams):
            corrML = -neff*di - n *np.log(neff) / 2
            dcorr["ML_"+name+"_"+neff_name] = corrML

    # evaluate prediction
    tpoint1.compute_fitness(params=optparams, include_tau=True)
    ave_full_fit1 = patient.timePoints[tp1].fitness()  # fitness of the Primary tumor
    tpoint2.compute_fitness(params=optparams, include_tau=True)
    ave_full_fit2 = patient.timePoints[tp2].fitness()  # fitness of the Recurrent tumor

    patient.predict(tp1, tp2, cleps, optparams, beta=1, include_nested=include_nested)
    cdat = patient.clone_frequency_predictions(tp1, tp2, beta=1.0, cleps=0.03)
    cdat.to_csv(os.path.join(odir1, "Clone_optimized_predictions",
                             "Clone_data_tree_"+ patient.name + "_" + tp1 + "_" + tp2 +".txt"), sep="\t", index=False)

    patient.write_mutation_predictions(tp1, tp2, os.path.join(odir1, "Clone_optimized_predictions"), beta=1.0)

    #patient.write_mutation_fitness(tp1, tp2, os.path.join(odir1, "Mutation_optimized_fitness"), beta=1.0)
    patient.write_neoantigen_fitness(tp1, tp2, os.path.join(odir1, "Mutation_optimized_fitness"), beta=1.0)

    tpoint1.compute_fitness(params=optparams1, include_tau=True)
    tpoint2.compute_fitness(params=optparams1, include_tau=True)
    patient.predict(tp1, tp2, cleps, optparams1, beta=1, include_nested=include_nested)
    cdat = patient.clone_frequency_predictions(tp1, tp2, beta=1.0, cleps=0.03)
    cdat.to_csv(os.path.join(odir1, "Clone_dg_optimized_predictions",
                             "Clone_data_tree_" + patient.name + "_" + tp1 + "_" + tp2 + ".txt"), sep="\t", index=False)
    patient.write_mutation_predictions(tp1, tp2, os.path.join(odir1, "Clone_dg_optimized_predictions"), beta=1.0)

    tpoint1.compute_fitness(params=optparams2, include_tau=True)
    tpoint2.compute_fitness(params=optparams2, include_tau=True)
    patient.predict(tp1, tp2, cleps, optparams2, beta=1, include_nested=include_nested)
    cdat = patient.clone_frequency_predictions(tp1, tp2, beta=1.0, cleps=0.03)
    cdat.to_csv(os.path.join(odir1, "Clone_immune_optimized_predictions",
                             "Clone_data_tree_" + patient.name + "_" + tp1 + "_" + tp2 + ".txt"), sep="\t", index=False)
    patient.write_mutation_predictions(tp1, tp2, os.path.join(odir1, "Clone_immune_optimized_predictions"), beta=1.0)


    tpoint1.compute_fitness(params=optparams, include_tau=False)
    ave_full_fit1_no_tau = patient.timePoints[tp1].fitness()  # fitness of the Primary tumor
    tpoint2.compute_fitness(params=optparams, include_tau=False)
    ave_full_fit2_no_tau = patient.timePoints[tp2].fitness()  # fitness of the Recurrent tumor

    tau = optparams["tau"]
    tau1 = optparams1["tau"]
    tau2 = optparams2["tau"]

    predstats = {
        "Full_fitness1": ave_full_fit1,
        "Full_fitness2": ave_full_fit2,
        "Full_fitness_no_tau1": ave_full_fit1_no_tau,
        "Full_fitness_no_tau2": ave_full_fit2_no_tau,
        "Delta_full_fitness": ave_full_fit2 - ave_full_fit1,
        "KL0": zerodi,
        "N_eff2": neff2,
        "N_eff_alt2":neff_alt2,
        "N_norm_eff2": norm_neff2,
        "N_norm_eff_alt2": norm_neff_alt2,
        "Tau": tau, "Tau1": tau1, "Tau2": tau2,
        "MinKL": optdi,
        "MinKL_dg_only": optdi1,
        "MinKL_immune_only": optdi2
    }

    for col in dcorr:
        predstats[col] = dcorr[col]

    for wn in wnames:
        predstats[wn + "_1"] = optparams1["weights"][wn]
        predstats[wn + "_2"] = optparams2["weights"][wn]
        predstats[wn+"_opt"] = optparams["weights"][wn]
    return predstats

    #    # EFFECTIVE population size for computing BIC and model comparison
#    neff2 = tpoint2.average_over_samples(Sample.effective_sample_size, eps=0.0)
#
#    # exluding mutations not detected on the recurrent tumor sample
#    neff_alt2 = tpoint2.average_over_samples(Sample.effective_sample_size, eps=0.0, filter_on_Taf=True)
#
#    sample_ids1 = tpoint1.sampleIDs
#    sample_ids2 = tpoint2.sampleIDs
#
#    }


if __name__ == "__main__":

    anl = Analysis()

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-dir", help="directory")
    parser.add_argument("-mapping", default=None, help="mapping.json file path")
    parser.add_argument("-config", default=None, help="config.json file path")
    parser.add_argument("-ptab", help="path to score table")
    parser.add_argument("-pval_thr", default=0.01, type=float, help="threshold for including parameter combinations")
    parser.add_argument("-tp_pref1", default="name of the first time point, eg. Primary | Pre")
    parser.add_argument("-tp_pref2", default="name of the second time point, eg. Metastasis | Met | Rec ")
    parser.add_argument("-test", action='store_true', help="")
    parser.add_argument("-include_nested", action='store_true', help="")
    parser.add_argument("-with_synch", action='store_true', help="met dtaa specific argument, to include synchronous samples")
    parser.add_argument("-ep_dist_model_name", default='all_tcr_bl62_quad_model', help="Zach's model")
    parser.add_argument("-cleps", type=float, default=0.03, help="frequency threshold on new clones")
    parser.add_argument("-o", "--odir", default=None, help="output directory")

    #testing partial models
    parser.add_argument("-no_w", action='store_true', help="just parameters table")
    parser.add_argument("-just_w", action='store_true', help="just parameters table")
    parser.add_argument("-just_DG", action='store_true', help="just driver gene component")
    parser.add_argument("-no_DG", action='store_true', help="no driver gene")
    parser.add_argument("-just_R", action='store_true')
    parser.add_argument("-include_R", action='store_true', help="whether to include R")

    parser.add_argument("-TNB", action='store_true', help="whether to use neoantigen load model")
    parser.add_argument("-aggr_fun", default="max")
    # how to include optimal parameters
    parser.add_argument("-just_akw", action='store_true', help="use only a, k and wd parameters from training")

    # 1. Initialize
    args = parser.parse_args()
    cleps = args.cleps

    if args.odir is None:
        odir = os.path.join(args.dir, "Results_DAR5")
    else:
        odir = args.odir
    if not os.path.exists(odir):
        os.mkdir(odir)

    odir = os.path.join(odir, "cleps_" + str(cleps) + "_nested_" + str(args.include_nested))
    if not os.path.exists(odir):
        os.mkdir(odir)

    if args.config is None:
        configpath = os.path.join(args.dir, "config.json")
    else:
        configpath = args.config
    with open(configpath, 'r') as fp:
        config = json.load(fp)

    anl.npos = "1-9"
    if args.mapping is None:
        if args.test:
            mappingfile = os.path.join(args.dir, "mapping_test.json")
        else:
            mappingfile = os.path.join(args.dir, "mapping.json")
    else:
        mappingfile = args.mapping

    with open(mappingfile, 'r') as fp:
        mapping = json.load(fp)

    if not args.with_synch:
        mapping = [x for x in mapping if x["type"] != "Synchronous"]


    anl.ntrees = 5
    alndir = os.path.join(args.dir, config["aln_dir"])
    iedbfasta = None
    if "iedb_file" in config:
        iedbfasta = os.path.join(args.dir, config["iedb_file"])
    else:
        iedbfasta = os.path.join(alndir, "enemy.fasta")

    anl.initialize_config(config, mapping, dir=args.dir, kd_thr=500)

    model = "DAR5"
    model_name = model + "_dgw" + '_'+args.ep_dist_model_name
    print(model_name)
    odir1 = os.path.join(odir, model_name)

    if not os.path.exists(odir1):
        try:
            os.mkdir(odir1)
        except:
            pass

    AUC = False
    if args.ptab is not None:
        ptab = pd.read_csv(args.ptab, sep="\t")
        AUC = False
        if "score_OS" in list(ptab.columns):
            AUC=False
        elif "AUC" in list(ptab.columns):
            AUC = True

    # 2. Initialize NQ paramaters
    Neoantigen.WEPS = 0.0
    # optional - computes optimal parameters based on the pvalues.txt table from survival analysis
    Utils.set_ptab_params(ptab, args.pval_thr,
                          no_w=args.no_w, just_w=args.just_w,
                          no_DG=args.no_DG,
                          just_DG=args.just_DG, AUC=AUC)
    sigma_I = Utils.params["optimal_params"]["sigma"]
    tau = Utils.params["optimal_params"]["tau"]
    a = Utils.params["optimal_params"]["a"]
    k = 1
    wd = Utils.params["optimal_params"]["wd"]

    # write optimal parameters
    with open(os.path.join(odir1, "parameters_DAR.txt"), 'w') as of:
        of.write("just_w\t"+str(args.just_w)+"\n")
        of.write("no_w\t"+str(args.no_w)+"\n")
        of.write("just_dg\t"+str(args.just_DG)+"\n")
        of.write("no_dg\t"+str(args.no_DG)+"\n")
        for col in Utils.params["optimal_params"]:
            of.write(col+"\t"+str(Utils.params["optimal_params"][col])+"\n")
        if "optimal_score" in Utils.params:
            of.write("optimal_score"+"\t"+str(Utils.params["optimal_score"])+"\n")
        print(Utils.params["optimal_score"])



    # 2. set neoantigen quality model
    kwargs = {"include_R": args.include_R,
              "just_R": args.just_R,
              "TNB": args.TNB}
    Qmodel = DAR5NeoantigenQuality(alndir=alndir,
                                   iedbfasta=iedbfasta,
                                   ep_dist_model_name=args.ep_dist_model_name)
    if args.TNB:
        fitnessModelComp1 = ImmuneCloneFitness(aggrfun=sum)
    else:
        if args.aggr_fun == "max":
            fitnessModelComp1 = ImmuneCloneFitness(aggrfun=max)
        elif args.aggr_fun == "beta_max":
            fitnessModelComp1 = ImmuneCloneFitness(aggrfun=lambda x, beta: np.exp(beta*x-max(beta*x)))

    Qmodel.M = 1
    print("setting neoantigen qualities")
    anl.set_neantigen_quality_model(Qmodel)
    print("done.")

    fitnessModelComp2 = DGCloneFitness(genes=['TP53', 'KRAS', 'CDKN2A', 'SMAD4'])

    aggrname = str(max).split(" ")[-1].replace(">", "")
    print("computing neoantigen qualities")
    if args.just_akw:
        anl.compute_neoantigen_sample_qualities(a=a, k=k, wd=wd, kdthr=500, include_R=args.include_R, just_R=args.just_R)
    else:
        anl.compute_neoantigen_sample_qualities(weighted_params=Utils.params, kdthr=500, **kwargs)
    print("done.")

    # 3. Set fitness components, optimize component weights

    recompute_components = True
    anl.reset_fitness_model_components()
    anl.set_fitness_model_component(fitnessModelComp1, "immune", 1.)
    anl.set_fitness_model_component(fitnessModelComp2, "dg", 0.0)

    # compute fitness of all clones
    anl.compute_node_fitness(recompute_components=recompute_components)
    #anl.standardize_node_fitness()
    for patient in anl.patients.values():
        patient.write_trees(os.path.join(odir1, "Tree_stats"))

    anl.write_neoantigen_fitness(os.path.join(odir1, "NQ"),  exclusive=True, suffix="", longitudinal=True)
    print("Initializing tree pairs")
    anl.init_tree_pairs(args.tp_pref1, args.tp_pref2, cleps, beta=1, include_nested=args.include_nested)
    print("Done")

    dimmune_clone_pred = defaultdict(lambda: (0.,0.))
    sample_size_corr = {} #correct effective sample size by the number of recurrent samples
    immune_params = {"tau": 1, "weights": {"dg": 0., "immune": 1.}}

    for patient in anl.patients.values():
        tpairs = patient.get_time_pairs_by_prefix(args.tp_pref1, args.tp_pref2)
        sample_size_corr[patient.name] = len(tpairs)
        for tp1, tp2 in tpairs:
            tpoint1 = patient.timePoints[tp1]
            tpoint1.compute_fitness(params=immune_params, include_tau=True)
            patient.write_mutation_fitness(tp1, tp2, os.path.join(odir1, "Mutation_immune_fitness"), beta=1.0)
            patient.write_neoantigen_fitness(tp1, tp2, os.path.join(odir1, "Mutation_immune_fitness"), beta=1.0)
            patient.predict(tp1, tp2, cleps, immune_params, beta=1, include_nested=args.include_nested)
            accres = patient.clone_frequency_predictions(tp1, tp2, beta=1.0)
            patient.write_mutation_predictions(tp1, tp2, os.path.join(odir1, "Clone_immune_predictions"), beta=1.0)
            dimmune_clone_pred[patient.name+"|"+tp1+"|"+tp2] = accres
    recompute_components = False

    # standardize fitness components (there is still tau parameter to regulate the fitness amplitude)
    #anl.standardize_node_fitness()

    params = anl.Qmodel.get_parameters()
    ptab = pd.DataFrame([[p, params[p]] for p in params.keys()])
    ptab.to_csv(os.path.join(odir1, "parameters.txt"), sep="\t", index=False)

    # information gain on predictions
    tab = []
    trunkality_shared2 = 0.0
    tables = defaultdict(list)

    paired_stats = []
    prim_stats = []

    first_elem = True
    odiry = os.path.join(odir1, "Report_y")
    if not os.path.exists(odiry):
        os.mkdir(odiry)

    if not os.path.exists(os.path.join(odir1, "Clone_immune_predictions")):
        os.mkdir(os.path.join(odir1, "Clone_immune_predictions"))

    for patient in anl.patients.values():
        # get time point pairs, primary-recurrent tumors
        tpairs = patient.get_time_pairs_by_prefix(args.tp_pref1, args.tp_pref2)
        print(patient.name + " n tpairs=" + str(len(tpairs)))

        # Sets the mutation - node index, for faster access to node mutation content
        patient.set_mutation_node_index()

        #iterate over all pairs
        for tp1, tp2 in tpairs:
            tpoint1 = patient.timePoints[tp1]
            tpoint2 = patient.timePoints[tp2]

            sample_ids1 = tpoint1.sampleIDs
            sample_ids2 = tpoint2.sampleIDs

            sites1 = "_".join(tpoint1.tissues)
            sites2 = "_".join(tpoint2.tissues)

            report_y(tpoint1, tpoint2, patient, odiry)

            inparams1 = {"tau": 1, "weights": {"dg": 1., "immune": 0.}}
            inparams2 = {"tau": 1, "weights": {"dg": 0., "immune": 1.}}

            tpoint1.compute_fitness(params=inparams2, include_tau=False)
            tpoint2.compute_fitness(params=inparams2, include_tau=False)
            #mark which clones are new
            patient.mark_clone_sharing(tp1, tp2, args.cleps)

            general_stats = compute_pair_stats(patient, tp1, tp2)

            trained_params = Utils.params["optimal_params"]

            pred_stats = compute_prediction_stats(patient, tp1, tp2, cleps, args.include_nested,
                                                  odir1,
                                                  sample_size_corr=sample_size_corr[patient.name],args=args
                                                    )

            #acc_cols = [("Accuracy_immune", "acc"),
            #            ("Weighted_accuracy_immune", "wacc"),
            #            ("Accuracy_immune_X", "acc_x"),
            #            ("Weighted_accuracy_immune_X", "wacc_x"),
            #            ("Weighted_rec_accuracy_immune", "wacc_rec"),
            #            ("Weighted_rec_accuracy_immune_X", "wacc_rec_x")]

            #for i, (acc_col, acc_key) in enumerate(acc_cols):
            #    pred_stats[acc_col] = dimmune_clone_pred[patient.name+"|"+tp1+"|"+tp2][acc_key]

            if first_elem:
                paired_columns = ["Patient", "Sample_ID1", "Sample_ID2", "Cohort", "OS", "PFS", "TP1", "TP2", "Sites_prim", "Sites_met"]
                gcols = [col for col in general_stats]
                predcols = [col for col in pred_stats]
                paired_columns += gcols
                paired_columns += predcols
                first_elem = False

            row_data = [patient.name, sample_ids1, sample_ids2,
                        patient.cohort, patient.OS, patient.PFS,
                        tp1, tp2, sites1, sites2]
            row_data += [general_stats[gstat] for gstat in general_stats]
            row_data += [pred_stats[pstat] for pstat in pred_stats]
            paired_stats.append(row_data)

        #patient.write_trees(os.path.join(odir1, "Tree_stats"))

    paired_stats = pd.DataFrame(paired_stats)
    paired_stats.columns = paired_columns

    paired_stats["IG"] = paired_stats.KL0 - paired_stats.MinKL
    paired_stats["IGn"] = paired_stats.IG / paired_stats.KL0

    paired_stats["IG_dg_only"] = paired_stats.KL0 - paired_stats.MinKL_dg_only
    paired_stats["IGn_dg_only"] = paired_stats.IG_dg_only / paired_stats.KL0

    paired_stats["IG_immune_only"] = paired_stats.KL0 - paired_stats.MinKL_immune_only
    paired_stats["IGn_immune_only"] = paired_stats.IG_immune_only / paired_stats.KL0

    cols = ["TMB", "TMB_nsyn", "TMB_syn",
            "TMB_eff", "TMB_nsyn_eff", "TMB_syn_eff",
            "TMB_MHC", "TMB_MHC_eff",
            "TNB", "TNB_eff", #"g_stat", "g_stat_eff","gall_stat", "gall_stat_eff",
            "Entropy",
            "Immune_fitness_cost"
            ]

    for col in cols:
        y1 = list(paired_stats[col + "1"])
        y2 = list(paired_stats[col + "2"])
        eps = (max(y1 + y2) - min(y1 + y2)) / 100

        paired_stats["Ratio_" + col] = (paired_stats[col + "2"] + eps) / (paired_stats[col + "1"] + eps)
        paired_stats["logRatio_" + col] = np.log((paired_stats[col + "2"] + eps) / (paired_stats[col + "1"] + eps))
        paired_stats["Delta_" + col] = paired_stats[col + "2"] - paired_stats[col + "1"]

    paired_stats.Cohort = [c.replace("V", "") for c in paired_stats.Cohort]
    paired_stats.to_csv(os.path.join(odir1, "fitness_stats.txt"), sep="\t", index=False)
    plot_average_fitness(paired_stats, odir1)
