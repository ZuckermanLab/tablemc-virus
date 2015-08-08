#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT
#This really should be replaced by an actual measurement of the fraction of native contacts for each initial state.
echo "60 0" >> $WEST_PCOORD_RETURN
