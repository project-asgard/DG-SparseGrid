classdef PARTIAL_SD_TERM
    
    properties
        
        type
        dat
        
        % RHS g_func, matrix, and vector 
        % (vector presently unused but should be where we but the surface integrals) 
        g
        mat
        
        % LHS g_func, matrix, and vector
        LHS_g
        LHS_mat
        
        % for debugging purposes only
        mat_unrotated
        LHS_mat_unrotated
    end
    
    methods
        
        function pt = PARTIAL_SD_TERM(type_,g_,dat_)           
            if nargin<1
                type_ = 'mass';
            end
            if nargin<2
                g_ = @(x,p,t,dat) x.*0+1;
            end
            if nargin<3
                dat_ = [];
            end
            
            pt.type = type_;
            pt.g = g_;
            pt.dat = dat_;
        end
        
    end
    
end
