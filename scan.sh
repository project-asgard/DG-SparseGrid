#!/bin/bash
MATLAB=/usr/local/bin/matlab
for x in {3..7}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./')); asgard(fokkerplanck2_complete, 'lev', 4, 'deg', $x, 'dt', 0.2, 'num_steps', 250, 'grid_type', 'FG', 'timestep_method', 'BE')" &
done
