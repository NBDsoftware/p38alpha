# Community Network Analysis for P38-alpha (MAPK14) with non-canonical ligands #

This repository contains folders and files analyzed with Grant's lab. [Bio3D R package](http://thegrantlab.org/bio3d/) using replicate unbiased molecular dynamics trajectories (obtained using the BioExcel Building Blocks (BioBB's) [BioBB's API documentation](http://biobb-md.readthedocs.io/en/latest/))   

Trajectories have been converted from the gromacs xtc format to dcd format using catdcd.

Once Bio3D is fully installed with the dependencies needed for full Community Network
Analysis the scripts can be easily invoked via command line with:  

```
R CMD BATCH --no-save --no-restore consensusCNA.R
```

Or they can all be run at once using the bash script in this repo:

```
./run_all.sh
```
