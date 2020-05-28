#!/bin/bash
MATLAB=/Applications/MATLAB_R2019a.app/bin/matlab
for lev in {3..3}
do
for deg in {1..4}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method', 'BE','dt', 1, 'num_steps', 50,'grid_type', 'SG','lev',$lev,'deg',$deg, 'time_independent_A', true);" &
done
done
