function default_dim = checkDimension(nDims,dim)

% Set defaults for dimension

default_dim.name = 'x';
default_dim.BCL = 'N'; % neumann
default_dim.BCR = 'N'; % neumann

for d=1:nDims % BC variation in all dimensions
    default_dim.BCL_fList{d} = @(x,p,t) x.*0;
    default_dim.BCR_fList{d} = @(x,p,t) x.*0;
end
default_dim.BCL_fList{nDims+1} = @(t,p) 1;  % time variation of BCS
default_dim.BCR_fList{nDims+1} = @(t,p) 1;

default_dim.domainMin = 0;
default_dim.domainMax = 1;
default_dim.lev = 3;
default_dim.FMWT = []; % Gets filled in later
default_dim.init_cond_fn = @(x,p) x.*0;

% Check to make sure all fields exist.
% If not, use default.

fn = fieldnames(default_dim);
for k=1:numel(fn)
    if isfield(dim,fn{k})
        default_dim.(fn{k}) = dim.(fn{k});
    end
end

%%
% Check for erroneous fields

fn = fieldnames(dim);
for k=1:numel(fn)
    if ~isfield(default_dim,fn{k})
        error(strcat('Unrecognized field in dim', fn{k} ));
    end
end

end