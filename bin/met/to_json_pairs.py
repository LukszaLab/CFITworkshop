import argparse
import json
import os

from cfit.util.Analysis import Analysis
from cfit.util.Log import Log

if __name__ == "__main__":

    '''
    Run as python compute_ntau_AR.py -d <data_directory> -p -PDAC

    
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d", "--dir", help="data directory")
    parser.add_argument("-config", default=None, help="config.json file path")
    parser.add_argument("-mapping", default=None, help="config.json file path")
    parser.add_argument("-netMHC", default="34", help="netMHC version, 34 or 40")
    parser.add_argument('-ns', default='9', help="peptide lengths, e.g '9', '8,9', '8-14")

    args = parser.parse_args()

    if "-" in args.ns:
        ns = [int(x) for x in  args.ns.split("-")[:2]]
        ns = list(range(ns[0], ns[1]+1))
    else:
        ns = [int(x) for x in (args.ns).split(",")]

    odir = os.path.join("/Users/mluksza/Workspace/NeoantigenEditing/data/Patient_data_pairs")
    if not os.path.exists(odir):
        os.mkdir(odir)

    anl = Analysis()
    anl.set_MHC_version(args.netMHC)

    if args.mapping is None:
        mappingfile = os.path.join(args.dir, "mapping.json")
    else:
        mappingfile = args.mapping
    with open(mappingfile) as f:
        mappingjs = json.load(f)

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
    anl.initialize_config(configjs, mappingjs, args.dir, kd_thr=500, ns=ns)

    anl.clonal = 1  # whether to take the trees into account in the analysis, 1 - yes
    anl.ntrees = 5  # number of top scoring trees
    anl.model = "AR"

    for patient in anl.patients.values():
        if patient.type != "Metachronous":
            continue

        tpairs = patient.get_time_pairs_by_prefix("Primary", "Met")
        for (tp1, tp2) in tpairs:
            patient.mark_clone_sharing(tp1, tp2, eps=0.03, dir=1)
            tpoint1 = patient.timePoints[tp1]
            tpoint2 = patient.timePoints[tp2]
            fname = patient.name+"_"+tp1+"_"+tp2
            podir = os.path.join(odir, fname)
            podir = os.path.join(odir, fname)
            #if not os.path.exists(podir):
            #    os.mkdir(podir)

            pairjs = {}
            for (tp, tpoint) in [("primary", tpoint1), ("recurrent", tpoint1)]:
                pairjs[tp] = []
                print(len(tpoint.samples))
                for sname in tpoint.samples: #there is only one sample
                    sample = tpoint.samples[sname]
                    js = sample.toJSON()
                    js["patient"] = patient.name
                    js["cohort"] = (patient.cohort).replace("V", "")
                    js["OS"] = patient.OS
                    js["PFS"] = patient.PFS
                    js["status"] = patient.dead
                    js['HLA_genes'] = patient.HLAS
                    js["mutations"] = [mut.toJSON() for mut in patient.mutations.values()]
                    js["neoantigens"] = [neo.toJSON() for neo in patient.neoantigens.values()]
                    pairjs[tp] = js
            with open(os.path.join(odir, fname+".json"), 'w') as of:
                json.dump(pairjs, of, indent=True)

#        js = patient.toJSON()
#
#
