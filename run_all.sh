#!/bin/bash
# Script to run all analysis scripts.

function clean {
    rm -f BCI.pdb R.py *.pdf *.Rout file*.pdb network* .Rhistory
}

for dir in 7pv*
do
    clean
    cd $dir
    R CMD BATCH --no-save --no-restore consensusCNA.R
    cd ..
done