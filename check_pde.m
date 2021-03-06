function default_pde = check_pde(pde,opts)

%% Other workflow options that should perhpas not be in the PDE?

default_pde.set_dt = []; % Function which accepts the pde (after being updated with CMD args).
default_pde.solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
default_pde.applySpecifiedE = 0; % Controls the "workflow" ... something we still don't know how to do generally. 
default_pde.checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.
default_pde.CFL = 0.01;
default_pde.dimensions = {};
default_pde.terms = {};
default_pde.params = {};
default_pde.sources = {};
default_pde.termsLHS = {};
[~, default_pde.transform_blocks] = OperatorTwoScale_wavelet2(opts.deg, opts.max_lev);

if opts.many_solution_capable
    default_pde.solutions = {};
    default_pde.initial_conditions = {};
else
    default_pde.analytic_solutions_1D = {};
end

% Check to make sure all fields exist.
% If not, use default.

fn = fieldnames(default_pde);
for k=1:numel(fn)
    if isfield(pde,fn{k})
        default_pde.(fn{k}) = pde.(fn{k});
    end
end

%%
% Check if there are erroneous field names

fn = fieldnames(pde);
for k=1:numel(fn)
    if ~isfield(default_pde,fn{k})
        error(strcat('Unrecognized term in PDE: ', fn{k} ));
    end
end

end
