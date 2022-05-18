classdef GLOBAL_SOLUTION_VECTOR < handle
    
    properties
        fvec;
        fvec_size;
        unknowns;
        num_unknowns;
        lbounds;
        ubounds;
    end
    
    methods
        
        function gsv = GLOBAL_SOLUTION_VECTOR( unknowns )
            
            gsv.unknowns = unknowns;
            
            gsv.num_unknowns = numel(gsv.unknowns);
            
            gsv.fvec_size = 0;
            for i = 1 : gsv.num_unknowns
                gsv.fvec_size = gsv.fvec_size + gsv.unknowns{i}.size();
            end
            
            gsv.fvec = zeros( gsv.fvec_size, 1 );
            
            gsv.set_bounds();
            
            for i = 1 : numel(gsv.unknowns)
                gsv.unknowns{i}.set_bounds( gsv.lbounds(i), gsv.ubounds(i) );
            end
            
        end
        
        function set_bounds( obj )
            
            obj.lbounds = uint64(zeros(obj.num_unknowns,1));
            obj.ubounds = uint64(zeros(obj.num_unknowns,1));
            
            os = 0;
            for i = 1 : obj.num_unknowns
                obj.lbounds(i) = os + 1;
                obj.ubounds(i) = os + obj.unknowns{i}.size();
                os = obj.ubounds(i);
            end
            
        end
        
        function [ sz ] = size_evolution_unknowns( obj )
            
            sz = 0;
            for i = 1 : numel(obj.unknowns)
                if(strcmp( obj.unknowns{i}.type, 'evolution' ) )
                    sz = sz + obj.unknowns{i}.size();
                end
            end
            
        end
        
        function [ sz ] = size_closure_unknowns( obj )
            
            sz = 0;
            for i = 1 : numel(obj.unknowns)
                if(strcmp( obj.unknowns{i}.type, 'closure' ) )
                    sz = sz + obj.unknowns{i}.size();
                end
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

