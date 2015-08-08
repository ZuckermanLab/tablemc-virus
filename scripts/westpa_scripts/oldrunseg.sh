#!/bin/bash

if [ -n "$SEG_DEBUG" ] ; then
    set -x
    env | sort
fi

cd $WEST_SIM_ROOT

mkdir -pv $WEST_CURRENT_SEG_DATA_REF || exit 1
cd        $WEST_CURRENT_SEG_DATA_REF || exit 1

echo "Starting iteration $WEST_CURRENT_ITER segment $WEST_CURRENT_SEG_ID at `date`"
free -m
ps -aux | sort -n -k4 -r | head
#clean out directory to begin with
rm -v *
#these need to be done no matter what kind of segment we are doing
script=seg.txt 
ln -vs $files/contact-map* .  
ln -vs $tables/*-0.5-5-15-10.dat .
ln -vs $files/hepb-capsid-wt-dimer-*-ca-only-all-sites.pdb .
ln -vs $files/$script .
if [ -n "$SEG_DEBUG" ] ; then
	ls -l
fi

#analysis parameters
rcutoff=80.0    
acutoff=30.0


# Set up the run
case $WEST_CURRENT_SEG_INITPOINT_TYPE in
    SEG_INITPOINT_CONTINUES)
        # A continuation from a prior segment
        INPUT_COOR=parent.rest
	ln -vs $WEST_PARENT_DATA_REF/seg.rest $INPUT_COOR
 	$tablemc $script > seg.out || exit 1
        #obtain information for parent from previous analysis, saves 90 seconds
        ln -vs $WEST_PARENT_DATA_REF/report-current report-parent
    ;;

    SEG_INITPOINT_NEWTRAJ)
        # Initiation of a new trajectory; $WEST_PARENT_DATA_REF contains the reference to the
        # appropriate basis state or generated initial state
	ln -vs $files/setup.* .
	ln -vs $files/capsid.rest .
	./setup.sh
        $tablemc $script > seg.out || exit 1
        #need to calculate reaction coordinates for parent also, since they are not previously available.
        awk '(NF==3) {step=$3} (NF>3) {print step,$0}' parent.rest > parent.dat
        $morphanal 300 1258.39 $rcutoff $acutoff capsid.rest parent.dat report-parent report-parent.pdb
    ;;

    *)
        echo "unknown init point type $WEST_CURRENT_SEG_INITPOINT_TYPE"
        exit 2
    ;;
esac

echo "Done with Monte Carlo, $WEST_CURRENT_ITER segment $WEST_CURRENT_SEG_ID at `date`"
free -m
ps -aux | sort -n -k4 -r | head
# Calculate progress coordinate and output coordinates
#Ernesto puts the progress coordinate of the previous segment first, then the current segment.
wait
#create pseudo trajectory file for the morph anal program and analyze current position
awk '(NF==3) {step=$3} (NF>3) {print step,$0}' seg.rest > lastframe.dat
ln -vs $files/capsid.rest .
#parameters: number of fragments, box size, distance cutoff, angle cutoff, template, trajectory, report file, pdb file
#the time data will come out on main output, not in the tee file
time $morphanal 300 1258.39 $rcutoff $acutoff capsid.rest lastframe.dat report-current report-current.pdb | tee morph-anal.out
nw=`grep -ci "Warning: failed to converge" morph-anal.out`
if [ "$nw" -gt "0" ]; then
	echo "RMSD analysis convergence failure, iteration $WEST_CURRENT_ITER segment $WEST_CURRENT_SEG_ID" >> $WEST_SIM_ROOT/convergence-failures
fi
#pick out reaction coordinates from reports
cat /dev/null > $WEST_PCOORD_RETURN
navgfrag=1 #Number of fragments whose nearest distances are averaged.
for x in `echo "parent current"`        
do
        #report is a list of clusters from largest to smallest. $3 is the size, $5 the min. distance to edges in the largest cluster
        clustersize=`awk '/cluster/ {print $4}' report-$x | head -n1`
        dist=`awk '($1~/fragment/) && ($4!=1)' report-$x | sort -n -k5 | head -n$navgfrag | awk -v n=$navgfrag 'BEGIN {x=0} {x+=$5} END {print x/n}'`
	#If the capsid is complete, then the distance is the largest distance inside the fragment
	if [ "$clustersize" -ge "120" ]; then
		dist=`awk '($1~/fragment/) && ($4==1)' report-$x | sort -n -k5 | tail -n1 | awk '{print $5}'`
	fi
        echo "$clustersize $dist" >> $WEST_PCOORD_RETURN
done

echo "Done with analysis, $WEST_CURRENT_ITER segment $WEST_CURRENT_SEG_ID at `date`"
free -m
ps -aux | sort -n -k4 -r | head
cat $WEST_PCOORD_RETURN


# Clean up
#rm -f $INPUT_COOR $INPUT_VEL *.conf *.dcd *.prm *.psf *.xsc [0-9]*
