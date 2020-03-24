#!/bin/bash
MATLAB=/Applications/MATLAB_R2019a.app/bin/matlab
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(continuity1,'lev',3,'deg',1,'implicit',false);" &
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(continuity1,'lev',3,'deg',2,'implicit',false);" &
$MATLAB -noFigureWindows -batch "addpath(genpath('./'));asgard(continuity1,'lev',3,'deg',3,'implicit',false);" 
