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
runTimeOpts.TEND = TEND;

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
runTimeOpts.quiet = quiet;

%%
% Use or not the compression reference version
if ~exist('compression','var') || isempty(compression)
    compression = 4;
end
runTimeOpts.compression = compression;

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
runTimeOpts.gridType = gridType;

%%
% Implicit or explicit time advance
if ~exist('implicit','var') || isempty(implicit)
    runTimeOpts.implicit = 0;
else
    runTimeOpts.implicit = implicit;
end

%%
% Use connectivity or not
if ~exist('useConnectivity','var') || isempty(useConnectivity)
    runTimeOpts.useConnectivity = 0;
else
    runTimeOpts.useConnectivity = useConnectivity;
end

%%
% CFL number
if ~exist('CFL','var') || isempty(CFL)
    pde.CFL = 0.1;
else
    pde.CFL = CFL;
end
