classdef MICRO_MACRO_1X1V
    
    properties
        dim_x;
        dim_v;
        lev_x;
        lev_v;
        deg;
        BCs_X = 'PERIODIC';
        BCs_V = 'PERIODIC';
    end
    
    methods
        
        function obj = MICRO_MACRO_1X1V( opts, dimensions, BCs_X, BCs_V )
            
            assert( numel(dimensions) == 2, 'Only Two Dimensions (1X and 1V)' )
            
            obj.dim_x = dimensions{1};
            obj.lev_x = dimensions{1}.lev;
            obj.dim_v = dimensions{2};
            obj.lev_v = dimensions{2}.lev;
            obj.deg   = opts.deg;
            
            if exist( 'BCs_X', 'var' )
                
                obj.BCs_X = BCs_X;
                
            end
            
            if exist( 'BCs_V', 'var' )
                
                obj.BCs_V = BCs_V;
                
            end
            
        end
        
        function rhs_Maxwellian = evaluate_rhs_Maxwellian( obj, opts, Q, t )
            
            N_x   = 2^obj.lev_x;
            N_v   = 2^obj.lev_v;
            dof_x = obj.deg * N_x;
            dof_v = obj.deg * N_v;
            dof   = dof_x * dof_v;
            
            rhs_Maxwellian = zeros( dof, 1 );
            
        end
        
        function rhs_vDotGradMaxwellian = evaluate_rhs_vDotGradMaxwellian( obj, opts, Q, t )
            
            N_x   = 2^obj.lev_x;
            N_v   = 2^obj.lev_v;
            dof_x = obj.deg * N_x;
            dof_v = obj.deg * N_v;
            dof   = dof_x * dof_v;
            
            rhs_vDotGradMaxwellian = zeros( dof, 1 );
            
        end
        
    end
end

