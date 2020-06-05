#!/bin/bash
MATLAB=/usr/local/bin/matlab
<<<<<<< HEAD
for lev in {5..5}
do
for deg in {6..7}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method', 'BE','dt', 1, 'num_steps', 50,'grid_type', 'FG','lev',$lev,'deg',$deg);" &
done
=======
for x in {3..7}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./')); asgard(fokkerplanck2_complete, 'lev', 4, 'deg', $x, 'dt', 0.2, 'num_steps', 250, 'grid_type', 'FG', 'timestep_method', 'BE')" &
>>>>>>> bdbfe7c86722d7397ea9677bc55ad3b8b1a354bb
done
