#!/bin/bash

while read f; do
    python make_muon_calib_plots.py $f
    echo ""
done < <(ls ../muon_calibration/*_Z_*)  

