%% 
% Load the PDE
if ~exist('pde','var') || isempty(pde)
    pde = Vlasov8;
end

nTerms = numel(pde.terms);
nDims = numel(pde.dimensions);

%%
% Simulation end time
if ~exist('TEND','var') || isempty(TEND)   
    TEND = 0.001;
end
opts.TEND = TEND;

%%
% Number of levels
if exist('lev','var')
    for d=1:nDims
        pde.dimensions{d}.lev = lev;
    end
end

%%
% Polynomial degree
% Deg = 2 Means Linear Element
if exist('deg','var')

    for d=1:nDims
        pde.dimensions{d}.deg = deg;
    end
end

%%
% Enable / disable print statements
if ~exist('quiet','var') || isempty(quiet)
    quiet = 0;
end
opts.quiet = quiet;

%%
% Use or not the compression reference version
if ~exist('compression','var') || isempty(compression)
    compression = 4;
end
opts.compression = compression;

%%
% Sparse or Full grid
if ~exist('gridType','var') || isempty(gridType)
    gridType = 'SG'; % Use 'SG' or 'FG'
else
    if strcmp(gridType,'SG') || strcmp(gridType,'FG')
    else
        error("gridType must be set to 'SG' or 'FG'");
    end
end
opts.gridType = gridType;

%%
% Implicit or explicit time advance
if ~exist('implicit','var') || isempty(implicit)
    opts.implicit = 0;
else
    opts.implicit = implicit;
end

%%
% Use connectivity or not
if ~exist('useConnectivity','var') || isempty(useConnectivity)
    opts.useConnectivity = 0;
else
    opts.useConnectivity = useConnectivity;
end

%%
% CFL number
if ~exist('CFL','var') || isempty(CFL)
    if ~isfield(pde,'CFL')
        pde.CFL = 0.1;
    end
else
    pde.CFL = CFL;
end
