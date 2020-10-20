classdef PDE
    properties        
        % these all need to be removed or moved to opts
        solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally.
        applySpecifiedE = 0; % Controls the "workflow" ... something we still don't know how to do generally.
        checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.
        
        set_dt = @(pde,CFL) 1; % Function which accepts the pde (after being updated with CMD args).
        dimensions = {};
        terms = {};
        params = {};
        sources = {};
        termsLHS = {};
        transform_blocks = {}; % will be updated in asgard.m
        analytic_solutions_1D = {}; % to be removed when transition all to many_solution_capable
        
        % MSC = Many Solution Capable fields
        % (all PDEs will transition to this soon)
        solutions = {};
        initial_conditions = {};        
    end
    
    methods
        
        function pde = PDE(opts,dimensions,terms,LHS_terms,sources,params,set_dt,analytic_solutions_1D,initial_conditions_MSC,solutions_MSC)
            
            pde.dimensions = dimensions;
            pde.terms = terms;
            pde.termsLHS = LHS_terms;
            pde.sources = sources;
            pde.analytic_solutions_1D = analytic_solutions_1D;
            pde.params = params;
            pde.set_dt = set_dt;
            [~,pde.transform_blocks] = OperatorTwoScale_wavelet2(opts.deg, opts.max_lev);
            if nargin > 8
                pde.solutions = solutions_MSC;
                pde.initial_conditions = initial_conditions_MSC;
            end
            
        end 
    end  
end
