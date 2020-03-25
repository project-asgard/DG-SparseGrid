#!/bin/bash
MATLAB=/usr/local/bin/matlab
for lev in {3..3}
do
for deg in {3..7}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'lev',$lev,'deg',$deg,'implicit',true, 'dt', 0.05, 'num_steps', 16, 'grid_type', 'FG', 'time_independent_A', true);" &
done
done
