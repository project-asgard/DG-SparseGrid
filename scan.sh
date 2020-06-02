#!/bin/bash
MATLAB=/usr/local/bin/matlab
for lev in {3..5}
do
for deg in {3..7}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method', 'BE','dt', 1, 'num_steps', 50,'grid_type', 'FG','lev',$lev,'deg',$deg);" &
done
done
