classdef PDE
    properties
        
        % these all need to be removed or moved to opts
        solvePoisson = 0; % Controls the "workflow" ... something we still don't know how to do generally.
        applySpecifiedE = 0; % Controls the "workflow" ... something we still don't know how to do generally.
        checkAnalytic = 1; % Will only work if an analytic solution is provided within the PDE.
        CFL = 0.01;
        
        set_dt = []; % Function which accepts the pde (after being updated with CMD args).
        dimensions = {};
        terms = {};
        params = {};
        sources = {};
        termsLHS = {};
        transform_blocks = {}; % will be updated in asgard.m
        solutions = {};
        initial_conditions = {};
        analytic_solutions_1D = {};
                       
    end
    
    methods
        
        function pde = PDE(dimensions)
            
            if nargin == 1
                pde.dimensions = dimensions;
                
                [~,pde.transform_blocks] = OperatorTwoScale_wavelet2(opts.deg, opts.max_lev); 
                
            end
            
        end
    end
    
end