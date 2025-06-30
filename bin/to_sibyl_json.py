import argparse
import json
import os

from cfit.util.Analysis import Analysis
from cfit.util.Log import Log
from cfit.util.Utils import Utils

if __name__ == "__main__":

    '''
    Export to JSON

    '''
    parser = Utils.make_cohort_parser()
#    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#    parser.add_argument("-d", "--dir", help="data directory")
#    parser.add_argument("-config", default=None, help="config.json file path")
#    parser.add_argument("-mapping", default=None, help="config.json file path")
#    parser.add_argument("-netMHC", default="34", help="netMHC version, 34 or 40")
#    parser.add_argument('-ns', default='9', help="peptide lengths, e.g '9', '8,9', '8-14")

    args = parser.parse_args()

    if "-" in args.ns:
        ns = [int(x) for x in  args.ns.split("-")[:2]]
        ns = list(range(ns[0], ns[1]+1))
    else:
        ns = [int(x) for x in (args.ns).split(",")]

    odir = os.path.join("./sibyl_data")
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

    # for patient in anl.patients.values():
    #     print(patient.to_sibyl_JSON())
    # Step 1: Create an array to store JSON representations of patients
    patients_json_array = []

    # Step 1: Iterate over all patients and convert them to JSON
    for patient in anl.patients.values():
        #timePoints = patient.timePoints
        patient.tr
        #print(patient.name)
        # patient.print_attributes()
        patient_json = patient.get_patient_to_sibyl()
        tree_nodes_json = patient.get_tree_nodes_to_sibyl()
        sample_nodes_json = patient.get_sample_tree_nodes_to_sibyl()
        mutations_json = patient.get_mutations_to_sibyl()
        neoantigens_json = patient.get_neoantigens_to_sibyl()
     

        patient_dir = os.path.join(odir, patient.name)
   
        if not os.path.exists(patient_dir):
            os.makedirs(patient_dir)
            
         
        # Constructing the filename by appending '.json' to the patient ID
        patient_filename = 'patient_'+patient.name + '.json'
        tree_node_filename = 'tree_nodes_' + patient.name + '.json'
        sample_tree_node_filename = 'sample_tree_nodes_'+patient.name+'.json'
        mutations_filename = 'mutations_'+patient.name+'.json'
        neoantigens_filename = 'neoantigens_'+patient.name+'.json'
        
        
        patient_file_path = os.path.join(patient_dir, patient_filename)
        tree_node_file_path = os.path.join(patient_dir, tree_node_filename)
        sample_tree_node_file_path = os.path.join(patient_dir, sample_tree_node_filename)
        mutations_file_path = os.path.join(patient_dir, mutations_filename)
        neoantigens_file_path = os.path.join(patient_dir, neoantigens_filename)
        
            
        with open(patient_file_path, 'w') as json_file:
            json.dump(patient_json, json_file, indent=2)
        
        with open(tree_node_file_path, 'w') as json_file:
            json.dump(tree_nodes_json, json_file, indent=2)
                    
        with open(sample_tree_node_file_path, 'w') as json_file:
            json.dump(sample_nodes_json, json_file, indent=2)

        with open(mutations_file_path, 'w') as json_file:
            json.dump(mutations_json, json_file, indent=2)

        with open(neoantigens_file_path, 'w') as json_file:
            json.dump(neoantigens_json, json_file, indent=2)
            


        
       