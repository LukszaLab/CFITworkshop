# Cancer FITness Modeling 

Package implementing the neoantigen fitness model and evolutionary predictions.

Installation:

```
cd CFIT
pip install . --user
```

To run n(tau) calculation with AxR model: 
```
hdir=test/Data_orig
python bin/compute_ntau_AR.py -d $hdir -PDAC -c $hdir/config.json -netMHC 3.4 -ns 9 -kd_thr 500
```
and CAR model using score landscape from precomputed pvalues.txt file:
```
hdir=test/Data_orig
python bin/compute_ntau_CAR.py -d $hdir -config $hdir/config_2022_01.json -netMHC 3.4 -ns 9 -kd_thr 500
#or
python bin/compute_ntau_CAR.py -d $hdir -config $hdir/config_2022_01.json -ptab $hdir/Results/Results_CAR_2022_01/pvalues.txt -netMHC 3.4 -ns 9 -kd_thr 500

```
