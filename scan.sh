#!/bin/bash
MATLAB=matlab

grid_type=FG
dt=1
num_steps=50
PDE=fokkerplanck2_complete
timestep_method=BE

for lev in {5..6}
do
for deg in {5..9}
do

ID=-$PDE-l$lev-d$deg-$grid_type-$dt-$num_steps-$timestep_method

$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard($PDE,'timestep_method','$timestep_method','lev',$lev,'deg',$deg,'dt',$dt,'num_steps',$num_steps,'grid_type','$grid_type','time_independent_build_A',true,'quiet',false,'save_output',true,'output_filename_id','$ID','plot_freq',100,'save_freq',100);" > log-asgard$ID.txt 2>&1 &

done
done
