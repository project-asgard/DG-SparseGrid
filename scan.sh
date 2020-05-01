#!/bin/bash
MATLAB=/Applications/MATLAB_R2019a.app/bin/matlab
for lev in {3..3}
do
for deg in {1..4}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(continuity1,'lev',$lev,'deg',$deg,'implicit',false);" &
done
done
