#!/bin/bash

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
#setup dunesw v09_89_01d01 -q e26:prof 
setup dunesw v09_91_03d00 -q e26:prof
setup gallery v1_22_06 -q e26:prof
# if changing versions, check what versions work with each other by checking dependencies with ups depend <packagename> <version> -q <qualifier>
