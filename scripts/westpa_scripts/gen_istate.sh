#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -p $(dirname $WEST_ISTATE_DATA_REF)
#cd $(dirname $WEST_ISTATE_DATA_REF)

#Generate all new state (lattice with random orientations) and put in $WEST_ISTATE_DATA_REF
ln -vs $files/hepb-capsid-wt-dimer-*-ca-only-some-sites.pdb .
sed -e "s@:new;@$WEST_ISTATE_DATA_REF@g" $files/generate.txt > generate.txt
$tablemc generate.txt

