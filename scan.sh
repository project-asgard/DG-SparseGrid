#!/bin/bash
MATLAB=/usr/local/bin/matlab
for deg in {5..8}
do
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(fokkerplanck2_complete,'timestep_method', 'BE','dt', 0.2, 'num_steps', 300, 'lev', 5, 'deg', $deg, 'save_output', true, 'output_filename_id', 'f2d_SG_lev5_deg$deg', 'grid_type', 'SG')" &
done
