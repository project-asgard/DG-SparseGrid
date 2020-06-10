#!/bin/bash
MATLAB=/Applications/MATLAB_R2020a.app/bin/matlab
for lev in {4..5}
do
for deg in {4..8}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method','BE','lev',$lev,'deg',$deg,'dt',1,'num_steps',50,'grid_type','SG','time_independent_build_A',true,'quiet',false,'save_output',true,'output_filename_id','-l$lev-d$deg')
;" &
done
done
