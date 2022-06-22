%% MATLAB (reference) version of the ASGarD solver

function [err,fval,fval_realspace,nodes,err_realspace,outputs,opts] = asgard(pde_handle,varargin)

%format short e
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
output_filename = create_output_filename(opts);
if exist(output_filename,'file') && opts.save_output
    err=-1;
    fval=-1;
    fval_realspace=-1;
    nodes=-1;
    err_realspace=-1;
    outputs=-1;
    disp('Output already exists ... skipping');
else
    [err,fval,fval_realspace,nodes,err_realspace,outputs,opts] = asgard_run_pde(opts,pde);
end

end

