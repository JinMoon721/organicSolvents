#!/bin/bash

dumpsDir="../../All/dumps/"
target="DBS"
nummols=10
numatoms=20
terminal=
boxX="29.84" ## box size A
boxY="29.84" ## box size A
boxZ="59.68" ## box size A
timestep="0.05" ## time step of snapshots, ps
eqtime="2.8652" ## equilibration time, ns

mkdir -p ../results/evaluate
../bin/evaluate $target $density $field $cutoffin $cutoffout $boxX $boxY $boxZ $timestep $eqtime 
