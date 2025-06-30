import glob
import os
from Bio import SeqIO
import pandas as pd

print(glob.glob("Neoantigens_all/34/*"))
for file in glob.glob("Neoantigens_all/34/*"):
    fname = os.path.basename(file)
    dat = pd.read_csv(file, sep="\t")
    print(dat.shape)
    idx = [len(pep)==9 for pep in dat.MT_Peptide]
    dat = dat[idx]
    dat.to_csv(os.path.join("Neoantigens","34", fname), sep="\t",index=False)