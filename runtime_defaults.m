input_parser = inputParser;

input_parser.KeepUnmatched = true;

default_start_time = 0;
default_lev = 3;
default_deg = 2;
default_num_steps = 5;
default_quiet = false;
default_implicit = false;
default_grid_type = 'SG';
valid_grid_types = {'SG','FG'};
check_grid_type = @(x) any(validatestring(x,valid_grid_types));
default_CFL = 0.01;
default_dt = 1.0;
default_adapt = false;
default_use_oldhash = false;
default_use_oldcoeffmat = false;
default_timestep_method = 'RK3';
valid_timestep_methods = {'BE','CN','BDFm','RK45m','RK3'};
check_timestep_method = @(x) any(validatestring(x,valid_timestep_methods));
default_time_independent_A = false;
default_many_solution_capable = false;
default_max_lev = 8;
default_adapt_threshold = 1e-3;
default_refinement_method = 1;
default_adapt_initial_condition = false;

addRequired(input_parser, 'pde', @isstruct);
addParameter(input_parser,'lev',default_lev, @isnumeric);
addParameter(input_parser,'deg',default_deg, @isnumeric);
addOptional(input_parser,'start_time',default_start_time, @isnumeric);
addOptional(input_parser,'num_steps',default_num_steps, @isnumeric);
addOptional(input_parser,'quiet',default_quiet,@islogical);
addOptional(input_parser,'implicit',default_implicit, @islogical);
addOptional(input_parser,'timestep_method',default_timestep_method, check_timestep_method);
addOptional(input_parser,'grid_type',default_grid_type, check_grid_type);
addOptional(input_parser,'CFL',default_CFL, @isnumeric);
addOptional(input_parser,'dt',default_dt, @isnumeric);
addOptional(input_parser,'adapt',default_adapt, @islogical);
addOptional(input_parser,'use_oldhash',default_use_oldhash, @islogical);
addOptional(input_parser,'use_oldcoeffmat',default_use_oldcoeffmat, @islogical);
addOptional(input_parser,'time_independent_A',default_time_independent_A,@islogical);
addOptional(input_parser,'many_solution_capable',default_many_solution_capable,@islogical);
addOptional(input_parser,'max_lev',default_max_lev, @isnumeric);
addOptional(input_parser,'adapt_threshold',default_adapt_threshold, @isnumeric);
addOptional(input_parser,'refinement_method',default_refinement_method, @isnumeric);
addOptional(input_parser,'adapt_initial_condition',default_adapt_initial_condition,@islogical);

if numel(varargin) == 0 && ~exist('pde','var')
    
    num_parameters = numel(input_parser.Parameters);
    
    disp('ASGarD - Adaptive Sparse Grid Discrization');
    disp('');
    disp('Run with ...');
    disp('');
    disp("asgard(pde_name,'opt_name',opt_val)");
    disp('');
    disp('e.g.,');
    disp('');
    disp("asgard(continuity1,'lev',4,'deg',3,'implicit',true,'CFL',0.1,'adapt',true)");
    disp('');
    disp('Available options are ...');
    
    for p=1:num_parameters
        disp(input_parser.Parameters(p));
    end
    
    error('Run with no arguments, exiting.');
    
end

parse(input_parser,pde,varargin{:})

num_terms       = numel(pde.terms);
num_dimensions  = numel(pde.dimensions);
num_steps       = input_parser.Results.num_steps;
t               = input_parser.Results.start_time;

if numel(input_parser.Results.lev) == num_dimensions
    %%
    % Specify lev_vec at runtime to have dimension dependent lev
    for d=1:num_dimensions
        pde.dimensions{d}.lev = input_parser.Results.lev(d);
    end
else
    %%
    % Specify a single lev which all dimensions get
    assert(numel(input_parser.Results.lev)==1);
    for d=1:num_dimensions
        pde.dimensions{d}.lev = input_parser.Results.lev;
    end
end

pde.deg = input_parser.Results.deg;

% CFL priority
% low        : default_CFL
% med        : pde.CFL
% high       : command line CFL
% extra high : set dt at command line

CFL_set = true;
opts.dt_set_at_runtime = true;
for c=input_parser.UsingDefaults
    if strcmp(c{1},'CFL')
        CFL_set = false;
    end
    if strcmp(c{1},'dt')
        opts.dt_set_at_runtime = false;
    end
end

opts.CFL = default_CFL;

if isfield(pde,'CFL')
    opts.CFL = pde.CFL;
end
if CFL_set
    opts.CFL = input_parser.Results.CFL;
end

opts.dt = input_parser.Results.dt;
opts.quiet = input_parser.Results.quiet;
opts.grid_type = input_parser.Results.grid_type;
opts.timestep_method = input_parser.Results.timestep_method;
opts.implicit = false;
if sum(strcmp(opts.timestep_method,{'BE','CN','BDFm'}))
    opts.implicit = true;
end

opts.adapt = input_parser.Results.adapt;
opts.use_oldhash = input_parser.Results.use_oldhash;
opts.use_oldcoeffmat = input_parser.Results.use_oldcoeffmat;
opts.time_independent_A = input_parser.Results.time_independent_A;
opts.many_solution_capable = input_parser.Results.many_solution_capable;
opts.max_lev = input_parser.Results.max_lev;
opts.adapt_threshold = input_parser.Results.adapt_threshold;
opts.refinement_method = input_parser.Results.refinement_method;
opts.adapt_initial_condition = input_parser.Results.adapt_initial_condition;

if opts.adapt
    opts.use_oldhash = false;
end

if ~isempty(fieldnames(input_parser.Unmatched))
   disp('Extra inputs:')
   disp(input_parser.Unmatched)
   error('Unrecognised input.')
end
