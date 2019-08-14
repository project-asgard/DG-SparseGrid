classdef mass < partial_term_1D
    
    properties
    end
    
    methods
        
        function m = mass(g_,dat_)
            if nargin<1
                g_ = @(x,p,t,dat) x.*0+1;
            end
            if nargin<2
                dat_ = [];
            end
            
            m@partial_term_1D('mass',g_,dat_);
            
            m.g = g_;
            m.dat = dat_;
        end
    end
    
end