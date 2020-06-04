#!/bin/bash
MATLAB=/usr/local/bin/matlab
<<<<<<< HEAD
for x in {3..7}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./')); asgard(fokkerplanck2_complete, 'lev', 4, 'deg', $x, 'dt', 0.2, 'num_steps', 250, 'grid_type', 'FG', 'timestep_method', 'BE')" &

=======
for lev in {5..5}
do
for deg in {6..7}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method', 'BE','dt', 0.125, 'num_steps', 400,'grid_type', 'SG','lev',$lev,'deg',$deg);" &
done
>>>>>>> 187aabc00617d0e839210c168e43b2041b0802e6
done
