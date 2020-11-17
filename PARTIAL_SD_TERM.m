classdef PARTIAL_SD_TERM
    
    properties
        type
        g
        dat
        mat
        mat_unrotated
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
