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

from cfit.patient.Sample import Sample
from cfit.tree.SampleTree import SampleTree
from cfit.tree.mutation.Neoantigen import Neoantigen
from cfit.util.Analysis import Analysis
from cfit.util.Utils import Utils

if __name__ == "__main__":
    anl = Analysis()

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-dir", help="directory")
    parser.add_argument("-mapping", default=None, help="mapping.json file path")
    parser.add_argument("-config", default=None, help="config.json file path")
    parser.add_argument("-with_synch", action='store_true', help="met dtaa specific argument, to include synchronous samples")

    # 1. Initialize
    args = parser.parse_args()

#    if args.odir is None:
#        odir = os.path.join(args.dir, "Results_DAR5")
#    else:
#        odir = args.odir
#    if not os.path.exists(odir):
#        os.mkdir(odir)

#    odir = os.path.join(odir, "cleps_" + str(cleps) + "_nested_" + str(args.include_nested))
#    if not os.path.exists(odir):
#        os.mkdir(odir)

    if args.config is None:
        configpath = os.path.join(args.dir, "config.json")
    else:
        configpath = args.config
    with open(configpath, 'r') as fp:
        config = json.load(fp)

    if args.mapping is None:
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

    anl.initialize_config(config, mapping, args.dir, kd_thr=500)

    purs = []
    for pname in anl.patients:
        pat = anl.patients[pname]
        cohort = pat.cohort.replace("V", "")
        for tp in pat.timePoints:
            tpoint = pat.timePoints[tp]
            for sample in tpoint.samples.values():
                purity = sample.average_tree_function(tree_fun=lambda tree: tree.purity)
                purs.append([pname, cohort, tp, sample.tissue, sample.name, purity])
    purs = pd.DataFrame(purs)
    purs.columns = ["Patient", "Cohort", "TP", "Tissue", "Sample", "Purity"]

    ax = sns.swarmplot(x='Cohort', y="Purity", data=purs, dodge=True, order=["STS", "LTS"], palette={"STS": "r", "LTS": 'b'})
    sns.boxplot(x='Cohort', y="Purity", data=purs,
                boxprops={'facecolor': 'None'}, order=["STS", "LTS"])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False)
    x1 = purs[purs.Cohort == 'STS']["Purity"]
    x2 = purs[purs.Cohort != 'STS']["Purity"]
    try:
        pval = scipy.stats.mannwhitneyu(x1, x2, alternative='two-sided').pvalue
    except:
        pval = 0.0
    ax.set_title("p=" + str(round(pval,4)))
    plt.savefig(os.path.join(args.dir, "purity.pdf"))
    plt.close()

    ax = sns.swarmplot(x='Tissue', y="Purity", data=purs, dodge=True, hue="Cohort", palette={"STS": 'r', "LTS": 'b'})
    sns.boxplot(x='Tissue', y="Purity", data=purs, boxprops={'facecolor': 'None'})
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(args.dir, "purity_by_tissue.pdf"))
    plt.close()




    plt.hist(purs.Purity)
    plt.savefig(os.path.join(args.dir, "purity_hist.pdf"))
    plt.close()

    purs.to_csv(os.path.join(args.dir, "purity.txt"), sep="\t", index=False)