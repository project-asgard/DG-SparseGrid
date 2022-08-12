%% MATLAB (reference) version of the ASGarD pde system solver

function [ err, opts, pde_system ] = asgard_pde_system( pde_system_handle, varargin )

    format short e
    folder = fileparts(which(mfilename));
    addpath(genpath(folder));

    %% Check for valid PDE function handle
    if ~exist('pde_system_handle','var')
        if ~isa(pde_system_handle,'functionhandle')
            error('Invalid PDE function handle, exiting.');
        end
    end

    %% Create options
    opts = OPTS(varargin);

    %% Create PDE system
    pde_system = pde_system_handle(opts);

    %% Launch ASGarD with this PDE
    output_filename = create_output_filename(opts);
    if exist(output_filename,'file') && opts.save_output
        disp('Output already exists ... skipping');
    else
        [ err, opts ] = asgard_run_pde_system( opts, pde_system );
    end

end