#!/bin/bash
MATLAB=/usr/local/bin/matlab
for lev in {3..3}
do
for deg in {3..8}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'lev',$lev,'deg',$deg,'implicit',true, 'dt', 0.1, 'num_steps', 8, 'grid_type', 'FG', 'time_independent_A', true);" &
done
done
