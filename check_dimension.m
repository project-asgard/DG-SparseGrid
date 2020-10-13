function default_dim = check_dimension(opts,num_dimensions,dim)

% Set defaults for dimension

default_dim.name = 'x';
default_dim.domainMin = 0;
default_dim.domainMax = 1;
default_dim.lev = opts.lev;
default_dim.init_cond_fn = @(x,p) x.*0;
default_dim.jacobian = @(x,p,t) x.*0 + 1;

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
        error(strcat('Unrecognized field in dim: ', fn{k}));
    end
end

end
