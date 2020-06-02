#!/bin/bash
MATLAB=/usr/local/bin/matlab
for lev in {5..5}
do
for deg in {3..7}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method', 'BE','dt', 0.05, 'num_steps', 600,'grid_type', 'SG','lev',$lev,'deg',$deg);" &
done
done
