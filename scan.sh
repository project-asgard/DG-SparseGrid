#!/bin/bash
MATLAB=/usr/local/bin/matlab
for lev in {5..5}
do
for deg in {3..7}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'lev',$lev,'deg',$deg,'implicit',true, 'dt', 0.1, 'num_steps', 8, 'grid_type', 'SG', 'adapt', true, 'adapt_initial_condition', true, 'adapt_threshold', 1e-4);" &
done
done
