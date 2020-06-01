#!/bin/bash
MATLAB=/usr/local/bin/matlab
for lev in {3..5}
do
for deg in {3..7}
do
<<<<<<< HEAD
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method', 'BE','dt', 1, 'num_steps', 50,'grid_type', 'FG','lev',$lev,'deg',$deg);" &
=======
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method', 'BE','dt', 1, 'num_steps', 50,'grid_type', 'SG','lev',$lev,'deg',$deg, 'time_independent_A', true);" &
>>>>>>> 3a76ff22a4b57fa6a3ef839d6362d88db6f5caa9
done
done
