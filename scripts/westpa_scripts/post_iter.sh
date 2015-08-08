#!/bin/bash

cp -v west.h5 $PBS_O_WORKDIR/
cp -v west-$PBS_JOBID.log $PBS_O_WORKDIR/
srcdir=`printf $WEST_SIM_ROOT/traj_segs/%06d $WEST_CURRENT_ITER`
destdir=`printf $PBS_O_WORKDIR/traj_segs/%06d $WEST_CURRENT_ITER`
#cp -vr $srcdir/* $destdir/
#copy back only the files we absolutely need.
cd $srcdir
mkdir $destdir
for seg in `ls`
do          
	mkdir $destdir/$seg
	#trajectory, sim. output, restart file
	cp -v $srcdir/$seg/seg.* $destdir/$seg/
	#analysis report and pdb file
	cp -v $srcdir/$seg/report-current* $destdir/$seg/
done
cd $WEST_SIM_ROOT
logfiles=`printf $WEST_SIM_ROOT/seg_logs/%06d-\* $WEST_CURRENT_ITER`
cp -v $logfiles $PBS_O_WORKDIR/seg_logs/
