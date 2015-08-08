#!/bin/bash

# This file defines where WEST and NAMD can be found
# Modify to taste

#attempt to ensure correct paths
export PATH=/usr/local/anaconda/bin:$PATH

# Inform WEST where to find Python and our other scripts where to find WEST
export WEST_PYTHON=$(which python)
#export WEST_ROOT=$(readlink -f $PWD/../../..)
export WEST_ROOT=$HOME/westpa

# Explicitly name our simulation root directory
if [[ -z "$WEST_SIM_ROOT" ]]; then
    export WEST_SIM_ROOT="$PWD"
fi
# Set environment variable for NAMD for convenience
#export NAMD=$(which namd2)
#set up paths for tablemc
export tablemc=~/tabulation-virus/tablemc-virus-10-27-14e/tablemc
#export morphanal=~/tabulation-virus/morph-anal-3-31-15/morph-anal
#after 500th iteration, fix bug involving last fragment
export morphanal=~/tabulation-virus/morph-anal-4-8-15/morph-anal
export files=$WEST_SIM_ROOT/tablemc-config/
#export on_compute_node=no
#Try copying the tables to the tmp space, to speed up I/O.
if [ ! -z "$PBS_NODEFILE" ]; then
	export tables=/tmp/tables
#	export on_compute_node=yes
	if [ ! -d $tables ]; then
		#pcmd mkdir -p $tables
		#pcmd cp -v ~/tabulation-virus/hepb-new-go5/tables/*-0.5-5-15-10.dat $tables
		#pbsdsh -u -v mkdir -p $tables
		#pbsdsh -u -v cp -v ~/tabulation-virus/hepb-new-go5/tables/*-0.5-5-15-10.dat $tables
		mkdir -p $tables
                cp -v ~/tabulation-virus/hepb-new-go5/tables/*-0.5-5-15-10.dat $tables
	fi
	echo "tables copied at `date`" >> special-log
	#scratch drive folder to work in
	export SCRDIR=/scr/$PBS_JOBID

	#if the scratch drive doesn't exist (it shouldn't) make it.
	if [[ ! -e $SCRDIR ]]; then
        	mkdir $SCRDIR
	fi
	echo scratch drive ${SCRDIR}

else
	export tables=~/tabulation-virus/hepb-new-go5/tables/
fi

# Set simulation name
export SIM_NAME=$(basename $WEST_SIM_ROOT)
echo "simulation $SIM_NAME root is $WEST_SIM_ROOT"
