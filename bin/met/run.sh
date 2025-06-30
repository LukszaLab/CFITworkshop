#hdir=$HOME/Dropbox/Research/10_PDAC/Modeling/Data_met_fs/
#python score_landscape_AR.py --dir $hdir -PDAC
#python score_landscape_TMB.py --dir $hdir

#python to_json.py --dir $hdir

#python to_json_pairs.py --dir $hdir

hdir=/Users/veramazeeva/MtSinai/Data_met_fs_corrected
cdir=$HOME/MtSinai/CFIT
python3 $cdir/bin/predict_opt.py -d $hdir -tp_pref1 Prim -tp_pref2 Met -kd_thr 500

