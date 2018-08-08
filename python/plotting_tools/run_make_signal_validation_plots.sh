#!/bin/bash

while read f; do
    python make_signal_validation_plots.py $f
    echo ""
done < <(ls ../validation/*b0*)  

