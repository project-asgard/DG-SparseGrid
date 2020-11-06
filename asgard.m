%% MATLAB (reference) version of the ASGarD solver

function [err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard(pde_handle,varargin)

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Check for valid PDE function handle
if ~exist('pde_handle','var')
    if ~isa(pde_handle,'functionhandle')
        error('Invalid PDE function handle, exiting.');
    end
end

%% Create options
opts = OPTS(varargin);

%% Create PDE
pde = pde_handle(opts);

%% Launch ASGarD with this PDE
[err,fval,fval_realspace,nodes,err_realspace,outputs] = asgard_run_pde(opts,pde);

end

