#!/bin/bash

$tablemc setup.txt
head -n1 capsid.rest > parent.rest
#fragments 1-20 and 41-50 correspond to the three "pentamers" we wish to exclude
awk '(NF>3) && ((($1>=21) && ($1<=40)) || ($1>=51))' capsid.rest > temp
awk '(NF>3)' lattice.rest >> temp
awk '($2=="ab")' temp > temp2
awk '($2=="cd")' temp >> temp2
awk 'BEGIN {n=1} {$1=n; print $0; n++}' temp2 >> parent.rest
