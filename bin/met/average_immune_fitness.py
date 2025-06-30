import json
import os

from cfit.util.Analysis import Analysis
from cfit.util.Log import Log
from cfit.util.Utils import Utils
#from cfit.fitness.neo_quality.CARNeoantigenQuality import
if __name__ == "__main__":

    '''
    Run as python compute_ntau_AR.py -d <data_directory> -p -PDAC

    
    '''

    parser = Utils.make_cohort_parser()
    args = parser.parse_args()

    if "-" in args.ns:
        ns = [int(x) for x in  args.ns.split("-")[:2]]
        ns = list(range(ns[0], ns[1]+1))
    else:
        ns = [int(x) for x in (args.ns).split(",")]

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
    anl.initialize_config(configjs, mappingjs, args=args)

    anl.clonal = 1  # whether to take the trees into account in the analysis, 1 - yes
    anl.ntrees = 5  # number of top scoring trees
    anl.model = "AR"

    js = anl.toJSON()
    with open(os.path.join(args.dir, "data.json"), 'w') as of:
        json.dump(js, of, indent=True)

    for patient in anl.patients.values():
        if patient.type != "Metachronous":
            continue
        patdir = os.path.join(args.odir, patient.name)
        primdir = os.path.join(patdir, "Primary")
        recdir = os.path.join(patdir, "Recurrent")
        if not os.path.exists(patdir):
            os.mkdir(patdir)
            os.mkdir(primdir)
            os.mkdir(recdir)
        for tp in patient.timePoints:
            tpoint = patient.timePoints[tp]
            sodir = primdir if "Primary" in tp else recdir
            for sname in tpoint.samples:
                fname = sname+".json"
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
                with open(os.path.join(sodir, fname), 'w') as of:
                    json.dump(js, of, indent=True)


#        js = patient.toJSON()
#        with open(os.path.join(odir, patient.name+".json"), 'w') as of:
#            json.dump(js, of, indent=True)
