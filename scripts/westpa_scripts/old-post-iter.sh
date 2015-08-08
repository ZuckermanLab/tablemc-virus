#!/bin/bash

cp -v west.h5 $PBS_O_WORKDIR/
cp -v west-$PBS_JOBID.log .
srcdir=`printf $WEST_SIM_ROOT/traj_segs/%06d $WEST_CURRENT_ITER`
destdir=`printf $PBS_O_WORKDIR/traj_segs/%06d $WEST_CURRENT_ITER`
cp -vr $srcdir/* $destdir/
logfiles=`printf $WEST_SIM_ROOT/seg_logs/%06d-\* $WEST_CURRENT_ITER`
cp -v $logfiles $PBS_O_WORKDIR/seg_logs/
