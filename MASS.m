classdef MASS < PARTIAL_SD_TERM
    
    properties
    end
    
    methods
        
        function m = MASS(g_,LHS_mass_g_,dat_,dV_)
            if nargin<1
                g_ = @(x,p,t,dat) x.*0+1;
            end
            if nargin<2
                LHS_mass_g_ = @(x,p,t,dat) x.*0+1;
            end
            if nargin<3
                dat_ = [];
            end
            if nargin<4
                dV_ = @(x,p,t,dat) x.*0+1;
            end
            
            m@PARTIAL_SD_TERM('mass',g_,LHS_mass_g_,dat_,dV_);
            
        end
    end
    
end
