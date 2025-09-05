#!/bin/bash

dumpsDir="../../All/dumps/"
target="DBS"
nummols=100
timestep="10" ## time step of snapshots, ps
eqtime="5" ## equilibration time, ns

mkdir -p ../results/evaluate
../bin/computeAngle $dumpsDir $target $nummols $timestep $eqtime  
