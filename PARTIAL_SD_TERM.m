classdef PARTIAL_SD_TERM
    
    properties
        
        type
        dat
        
        % RHS g_func, matrix, and vector 
        g
        mat
        
        % LHS g_func, matrix, and vector
        LHS_mass_g
        LHS_mass_mat
        
        % dV / jacobian / volume_element
        dV
        
        % for debugging purposes only
        mat_unrotated
        LHS_mass_mat_unrotated
    end
    
    methods
        
        function pterm = PARTIAL_SD_TERM(type_,g_,LHS_mass_g_,dat_,dV_)           
            if nargin<1
                type_ = 'mass';
            end
            if nargin<2
                g_ = @(x,p,t,dat) x.*0+1;
            end
            if nargin<3
                LHS_mass_g_ = @(x,p,t,dat) x.*0+1;
            end
            if nargin<4
                dat_ = [];
            end
            if nargin<5
                dV_ = @(x,p,t,d) x.*0+1;
            end
            
            pterm.type = type_;
            pterm.g = g_;
            pterm.LHS_mass_g = LHS_mass_g_;
            pterm.dat = dat_;
            pterm.dV = dV_;
        end
        
    end
    
end
