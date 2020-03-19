# Matlab version of ASGarD

## Quickstart

1.  Clone repository.
```
https://code.ornl.gov/lmm/DG-SparseGrid.git
```
2. Add the `DG-SparseGrid` folder (and subfolders) to your matlab path. In Matlab do 
```
cd DG-SparseGrid
addpath(genpath('./'))
```
3. Run a PDE, e.g., 
```
asgard(diffusion2)
```

## Run options
Run as 
```
[err,f_wSpace,f_rSpace] = asgard(PDE,name,value)
```
where the `name`,`value` pair options are defined in `runtime_defaults.m` ...
```
'lev', 3 | [3,4,5];                 % Scalar or num_dimensions long vector
'deg' , 2;                          % Globally defined basis function order
'num_steps', 5;                     % Number of time steps
'quiet', true | **false**;          % Display output and plots
'implicit', true | **false**;       % Implicit or explicit time advance
'grid_type', **'SG'** | 'FG';       % Sparse-grid (SG) or Full-grid (FG)
'CFL' = 0.01;                       % CFL number
'adapt' = true | **false**;         % Enable adaptivity
'use_oldhash' = true | **false**;   % Use original hash table or not
```
e.g., 
```
asgard(continuity1)
asgard(continuity1,'lev',4,'deg',3,'implicit',false)
asgard(continuity1,'implicit',true)
asgard(continuity1,'implicit',true,'CFL',0.1)
```

## Running the tests
Only the PDE tests
```
cd DG-SparseGrid
addpath(genpath('./'))
close all; clear all
runtests
```
All the tests (PDE, two_scale, kronmult)
```
cd DG-SparseGrid
addpath(genpath('./'))
close all; clear all
runtests('IncludeSubfolders',true)
```

## Producing the C++ gold testing data
```
cd DG-SparseGrid
addpath(genpath('./'))
cd gold
close all; clear all
make_all_gold
```

## Running with -nodesktop
This will enable running with no GUI for say running within a `screen` session.
```
matlab -nodesktop -nosplash -noFigureWindows
cd DG-SparseGrid
addpath(genpath('./'))
asgard(continuity1,'lev',4,'deg',3,'implicit',false)
```


