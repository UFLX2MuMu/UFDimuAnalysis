#!/bin/bash

while read f; do
    echo "FILE $f"
    python make_validation_plots.py $f
    echo ""
done < <(ls ../validation/*b0*)  

