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

%%
% Use or not the compression reference version
if ~exist('compression','var') || isempty(compression)
    compression = 4;
end

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

%%
% Implicit or explicit time advance
if ~exist('implicit','var') || isempty(implicit)
    pde.implicit = 0;
else
    pde.implicit = implicit;
end

%%
% If implicit select we must choose a particular tensor prod compression
if pde.implicit
    compression = 1;
end
