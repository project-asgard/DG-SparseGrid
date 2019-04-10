function default_dim = checkDimension(dim);

% Set defaults for dimension

default_dim.name = 'x';
default_dim.BCL = 2; % neumann
default_dim.BCL_fList = [];
default_dim.BCR = 2; % neumann
default_dim.BCR_fList = [];
default_dim.domainMin = 0;
default_dim.domainMax = 1;
default_dim.lev = 3;
default_dim.deg = 2;
default_dim.FMWT = []; % Gets filled in later
default_dim.init_cond_fn = @(x,p) x.*0;

% Overwrite with present entries from specified dimension

if isfield(dim,'name')    
    default_dim.name = dim.name;
end

if isfield(dim,'BCL')    
    default_dim.BCL = dim.BCL;
end

if isfield(dim,'BCL_fList')
    default_dim.BCL_fList = dim.BCL_fList;
end

if isfield(dim,'BCR')    
    default_dim.BCR = dim.BCR;
end

if isfield(dim,'BCR_fList')
    default_dim.BCR_fList = dim.BCR_fList;
end

if isfield(dim,'domainMin')
    default_dim.domainMin = dim.domainMin;
end

if isfield(dim,'domainMax')
    default_dim.domainMax = dim.domainMax;
end

if isfield(dim,'lev')
    default_dim.lev = dim.lev;
end

if isfield(dim,'deg')
    default_dim.deg = dim.deg;
end

if isfield(dim,'FMWT')
    default_dim.FMWT = dim.FMWT;
end

if isfield(dim,'init_cond_fn')
    default_dim.init_cond_fn = dim.init_cond_fn;
end

end