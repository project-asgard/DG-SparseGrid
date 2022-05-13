classdef GLOBAL_SOLUTION_VECTOR < handle
    
    properties
        fvec;
        fvec_size;
        num_unknowns;
        unknown_index_lo;
        unknown_index_hi;
    end
    
    methods
        
        function gsv = GLOBAL_SOLUTION_VECTOR( unknowns )
            
            gsv.num_unknowns = numel(unknowns);
            
            gsv.fvec_size = 0;
            for i = 1 : gsv.num_unknowns
                gsv.fvec_size = gsv.fvec_size + unknowns{i}.size();
            end
            
            gsv.fvec = zeros( gsv.fvec_size, 1 );
            
            [ gsv.unknown_index_lo, gsv.unknown_index_hi ]...
                = gsv.GetUnknownBounds( unknowns );
            
        end
        
        function [ lo, hi ] = GetUnknownBounds( obj, unknowns )
            
            lo = uint64(zeros(obj.num_unknowns,1));
            hi = uint64(zeros(obj.num_unknowns,1));
            
            os = 0;
            for i = 1 : obj.num_unknowns
                lo(i) = os + 1;
                hi(i) = os + unknowns{i}.size();
                os = hi(i);
            end
            
        end
        
        function out = copy( obj )
            
            out = obj.fvec;
            
        end
        
        function out = zeros( obj )
            
            out = zeros(size(obj.fvec));
            
        end
        
    end
end

