#!/bin/bash
# node.sh
#   Runs WESTPA on a single node
#   Written by Karl Debiec on 13-09-30, last updated by Karl Debiec on 14-05-06
################################################# MODULES AND SETTINGS #################################################
source /etc/profile.d/sam-modules.sh
module purge
module load env/intel-11.1-openmpi-1.5
#module load westpa
cd $PBS_O_WORKDIR  #needed in order to find env.sh
source env.sh || exit 1
cd     $WEST_SIM_ROOT
 
######################################################### MAIN #########################################################
#Try copying the tables to the tmp space, to speed up I/O.
if [ ! -z "$PBS_NODEFILE" ]; then
        export tables=/tmp/tables
        #if [ ! -d $tables ]; then
                pcmd mkdir -p $tables
                pcmd cp -v ~/tabulation-virus/hepb-new-go5/tables/*-0.5-5-15-10.dat $tables
        #fi
        echo "tables copied at `date`" >> special-log
else
        export tables=~/tabulation-virus/hepb-new-go5/tables/
fi

$WEST_ROOT/bin/w_run "$@"
