classdef MASS < PARTIAL_SD_TERM
    
    properties
    end
    
    methods
        
        function m = MASS(g_,dat_)
            if nargin<1
                g_ = @(x,p,t,dat) x.*0+1;
            end
            if nargin<2
                dat_ = [];
            end
            
            m@PARTIAL_SD_TERM('mass',g_,dat_);
            
            m.g = g_;
            m.dat = dat_;
        end
    end
    
end
