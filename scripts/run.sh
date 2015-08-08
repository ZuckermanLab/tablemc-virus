#!/bin/bash

source env.sh

if [ ! -z "$PBS_NODEFILE" ]; then
	logfile=west-$PBS_JOBID.log
else
	logfile=west.log
fi
rm -f $logfile
$WEST_ROOT/bin/w_run  "$@" &> $logfile
