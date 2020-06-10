#!/bin/bash
MATLAB=matlab
for lev in {5..6}
do
for deg in {4..9}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method','BE','lev',$lev,'deg',$deg,'dt',1,'num_steps',50,'grid_type','FG','time_independent_build_A',true,'quiet',false,'save_output',true,'output_filename_id','-l$lev-d$deg-FG','plot_freq',100,'save_freq',100);" &
done
done
