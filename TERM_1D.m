classdef TERM_1D
    
    properties
        pterms
        time_dependent {is_valid(time_dependent)} = false;
        mat
        mat_unrotated
    end
      
    methods
        function t = TERM_1D(pterms_,time_dependent_)
            if nargin<1
                pterms_ = {mass()};
            end
            if nargin<2
                time_dependent_ = false;
            end
            
            t.pterms = pterms_;
            t.time_dependent = time_dependent_;
        end
    end
    
end

function is_valid(td)
assert(islogical(td));
end
