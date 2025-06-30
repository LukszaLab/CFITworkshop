import argparse
import json
import os

from cfit.util.Analysis import Analysis
from cfit.util.Log import Log

if __name__ == "__main__":

    '''
    Export to JSON

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

    odir = os.path.join("../NeoantigenEditing/data/Patient_data")
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
    print('args', args, args.dir)
    anl.initialize_config(configjs, mappingjs, dir=args.dir, kd_thr=500, ns=ns)
    print(anl.samples.values())

def karina():

    anl.clonal = 1  # whether to take the trees into account in the analysis, 1 - yes
    anl.ntrees = 5  # number of top scoring trees
    anl.model = "AR"

#    js = anl.toJSON()
#    with open(os.path.join(odir, "data.json"), 'w') as of:
#        json.dump(js, of, indent=True)

    for patient in anl.patients.values():
        if patient.type != "Metachronous":
            continue
        patdir = os.path.join(odir, patient.name)
        primdir = os.path.join(patdir, "Primary")
        recdir = os.path.join(patdir, "Recurrent")
        if not os.path.exists(patdir):
            os.mkdir(patdir)
            os.mkdir(primdir)
            os.mkdir(recdir)

        tpoint1 = patient.timePoints["Primary"]
        sodir = primdir
        for sname in tpoint1.samples:
            fname = sname+".json"
            sample = tpoint1.samples[sname]
            js = sample.toJSON()
            js["patient"] = patient.name
            js["cohort"] = (patient.cohort).replace("V", "")
            js["OS"] = patient.OS
            js["PFS"] = patient.PFS
            js["status"] = patient.dead
            js['HLA_genes'] = patient.HLAS
            js["mutations"] = [mut.toJSON() for mut in patient.mutations.values()]
            js["neoantigens"] = [neo.toJSON() for neo in patient.neoantigens.values()]
            with open(os.path.join(sodir, fname), 'w') as of:
                json.dump(js, of, indent=True)

        tpairs = patient.get_time_pairs_by_prefix("Primary", "Met")
        for (tp1, tp2) in tpairs:
            tpoint1 = patient.timePoints[tp1]
            tpoint2 = patient.timePoints[tp2]

            tree_pairs = patient.get_tree_pairs(tp1, tp2, beta=1)
            for tree1, tree2, w in tree_pairs:
                new_clones = patient.identify_new_clones(tree1, tree2, eps=0.03)
                tree1.soft_remove_nodes(new_clones, include_nested=False, overwrite=False)
                tree2.soft_remove_nodes(new_clones, include_nested=False, overwrite=False)
            patient.mark_clone_sharing(tp1, tp2, eps=0.03, dir=1)
            sodir = recdir
            for sname in tpoint2.samples:
                for psname in tpoint1.samples:
                    fname = "paired_primary_tumor_"+sname+".json"
                    sample = tpoint1.samples[psname]
                    js = sample.toJSON()
                    js["patient"] = patient.name
                    js["cohort"] = (patient.cohort).replace("V", "")
                    js["OS"] = patient.OS
                    js["PFS"] = patient.PFS
                    js["status"] = patient.dead
                    js['HLA_genes'] = patient.HLAS
                    js["mutations"] = [mut.toJSON() for mut in patient.mutations.values()]
                    js["neoantigens"] = [neo.toJSON() for neo in patient.neoantigens.values()]
                    with open(os.path.join(sodir, fname), 'w') as of:
                        json.dump(js, of, indent=True)

                sample = tpoint2.samples[sname]
                fname = sname+".json"
                js = sample.toJSON()
                js["patient"] = patient.name
                js["cohort"] = (patient.cohort).replace("V", "")
                js["OS"] = patient.OS
                js["PFS"] = patient.PFS
                js["status"] = patient.dead
                js['HLA_genes'] = patient.HLAS
                js["mutations"] = [mut.toJSON() for mut in patient.mutations.values()]
                js["neoantigens"] = [neo.toJSON() for neo in patient.neoantigens.values()]
                with open(os.path.join(sodir, fname), 'w') as of:
                    json.dump(js, of, indent=True)


#        js = patient.toJSON()
#        with open(os.path.join(odir, patient.name+".json"), 'w') as of:
#            json.dump(js, of, indent=True)
