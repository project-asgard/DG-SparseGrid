classdef MASS < PARTIAL_SD_TERM
    
    properties
    end
    
    methods
        
        function m = MASS(g_,LHS_mass_g_,dat_,dV_)
            
            g = @(x,p,t,dat) x.*0+1;
            if exist('g_','var')
                if ~isempty(g_)
                    g = g_;
                end
            end
            
            LHS_mass_g = @(x,p,t,dat) x.*0+1;
            if exist('LHS_mass_g_','var')
                if ~isempty(LHS_mass_g_)
                    LHS_mass_g = LHS_mass_g_;
                end
            end
            
            if nargin<3
                dat_ = [];
            end
            if nargin<4
                dV_ = @(x,p,t,dat) x.*0+1;
            end
            
            m@PARTIAL_SD_TERM('mass',g,LHS_mass_g,dat_,dV_);
            
        end
    end
    
end
