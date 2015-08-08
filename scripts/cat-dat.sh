#!/bin/bash

if [ "$#" -lt "2" ]; then
	echo "usage: cat-dat tracefile outputfile [freq]"
	exit
fi
#This generates a "dat" file (centers/orientations) by concatenating individual segments.
trace=$1 #output by w_trace
output=$2
if [ "$#" -ge "3" ]; then 
	freq=$3
else
	freq=1
fi
datoutput=`echo "$output.dat"`
dcdoutput=`echo "$output.dcd"`
rm $datoutput
awk -v output=$datoutput '{if (NF==0) next; if (substr($0,1,1)=="#") next; iter=$1; seg=$2; if (iter<=0) next; print iter,seg; system(sprintf("cat traj_segs/%06d/%06d/seg.dat >> %s\n",iter,seg,output))}' $trace

#This expands the trajectory into atomic coordinates.
nframes=`awk 'BEGIN {n=0} ($2==1) {n++} END {print n}' $datoutput`
#get the box size from one of the output files, and the number of fragments from a restart file
boxsize=`grep -i "box" traj_segs/000001/000000/seg.out | awk '{print $(NF-1)}'`
echo $nframes $boxsize
#Copy FRAG declarations from master script and add the directory
awk '/FRAG/ {printf("%s %s tablemc-config/%s\n",$1,$2,$3)}' tablemc-config/setup.txt > expand.txt
echo "INSERT ab 150" >> expand.txt #this depends on actual system
echo "INSERT cd 150" >> expand.txt
echo "EXPAND $datoutput $nframes $freq $boxsize $dcdoutput" >> expand.txt
echo "END" >> expand.txt
~/tabulation-virus/tablemc-virus-10-27-14e/tablemc expand.txt
