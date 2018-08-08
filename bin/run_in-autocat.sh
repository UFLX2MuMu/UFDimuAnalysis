#!/bin/bash

TREEDIR="/home/puno/h2mumu/bdt/studies/h2mumu/trees/"
#find ${TREEDIR} -maxdepth 1 -name "*.xml"

while read f; do
    echo "Categorizing $f ..."
    echo ""
    ./categorize --categories=$f --binning=-4 --sig_xlumi=0
    echo ""
done < <(find ${TREEDIR} -maxdepth 1 -name "*.xml")

