classdef EQUATION < handle
    
    properties
        unknown;
        terms;
        type = 'evolution';
        LHS_func;
        LHS_term;
    end
    
    methods
        
        function equation = EQUATION( unknown, terms, type, LHS_func )
            
            assert(strcmp(type,'evolution') || strcmp(type,'closure'),...
                   "equation type must be 'evolution' or 'closure' " )
            
            equation.unknown = unknown;
            equation.terms   = terms;
            equation.type    = type;
            
            if( isempty(LHS_func) )
                
                num_dims = numel( unknown.dimensions );
                
                equation.LHS_func = cell( 1, num_dims + 1 );
                
                for i = 1 : num_dims
                    
                    equation.LHS_func{i} = @(x,p,t,d) 0*x+1;
                    
                end
                
                equation.LHS_func{num_dims+1} = @(t,p) 0*t+1;
                
            else
                
                equation.LHS_func = LHS_func;
                
            end
            
            sd_terms = cell( 1, num_dims ); % Does not include temporal component;
            
            for i = 1 : num_dims
                
                sd_terms{i} = SD_TERM( {MASS( equation.LHS_func{i} )} );
                
            end
            
            md_term = MD_TERM( num_dims, sd_terms );
            
            equation.LHS_term = TERM( unknown, {unknown}, {md_term} );
            
        end
        
    end
end

