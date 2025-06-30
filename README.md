# Cancer FITness Modeling 

Package implementing the neoantigen fitness model and evolutionary predictions.

Installation:

```
cd CFIT
pip install . --user
```
Download test data from https://www.dropbox.com/s/9d85d83ctgew1r1/test.zip?dl=0
and unzip into the main folder (CFIT).

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

JSON output
```
hdir=test/Data_orig
python bin/to_json.py -d $hdir
```
For later:


Model training (takes longer and can be parallelized)
1. AxR model
```
hdir=test/Data_orig
iedb=2016
for patdir in ${hdir}/VCF/*
do
    patient="$(basename $patdir)"
    echo $patient
    python ${cdir}/bin/align_neoantigens_per_patient.py -patient $patient -dir $hdir -iedb $iedb
done
```

```
hdir=test/Data_orig
python bin/score_landscape_AR.py -d $hdir -PDAC -config $hdir/config.json -netMHC 3.4 -kd_thr 500

```
2. CAR model (w logC + (1-w)log A)xR
```
hdir=test/Data_orig_fs
python bin/score_landscape_CAR.py -d $hdir -PDAC -config $hdir/config_2022_01.json -netMHC 3.4 -kd_thr 500

```
