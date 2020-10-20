classdef PDE
    properties
        
        % these all need to be removed or moved to opts
        solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally.
        applySpecifiedE = 0; % Controls the "workflow" ... something we still don't know how to do generally.
        checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.
        %         CFL = 0.01;
        set_dt = @(pde,CFL) 1; % Function which accepts the pde (after being updated with CMD args).
        dimensions = {};
        terms = {};
        params = {};
        sources = {};
        termsLHS = {};
        transform_blocks = {}; % will be updated in asgard.m
        solutions = {};
        analytic_solutions_1D = {};
        
    end
    
    methods
        
        function pde = PDE(opts,dimensions,terms,LHS_terms,sources,params,set_dt,analytic_solutions)
            
            pde.dimensions = dimensions;
            pde.terms = terms;
            pde.termsLHS = LHS_terms;
            pde.sources = sources;
            pde.analytic_solutions_1D = analytic_solutions; 
            pde.params = params;
            pde.set_dt = set_dt;
            [~,pde.transform_blocks] = OperatorTwoScale_wavelet2(opts.deg, opts.max_lev);       
            
        end
        
    end
    
    
end
