#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT
mkdir -p temp
cd temp

tablemc=~/tabulation-virus/tablemc-9-9-14/tablemc
files=$WEST_SIM_ROOT/tablemc-config/
script=native-contacts.txt

ln -vs $files/contact-map* .
ln -vs $files/hepb-capsid-wt-dimer-*-ca-only-some-sites.pdb .
sed -e "s@seg.rest@$WEST_STRUCT_DATA_REF@g" $WEST_SIM_ROOT/tablemc-config/$script > $script
if [ -n "$SEG_DEBUG" ] ; then
        ls -l
fi

# Calculate progress coordinate and output coordinates
#115 is the number of native contacts in a pentamer
$tablemc $script > native.out     
grep -i "Total number of native contacts:" native.out | awk '{print $NF/115}' >> $WEST_PCOORD_RETURN

