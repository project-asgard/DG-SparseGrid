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
            
            gsv.fvec_size = gsv.size();
            
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
        
        function [ sz ] = size( obj )
            
            sz = 0;
            for i = 1 : numel(obj.unknowns)
                sz = sz + obj.unknowns{i}.size();
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
        
        function out = copy_evolution_unknowns( obj )
            
            out = zeros(obj.size_evolution_unknowns(),1);
            
            os = 0;
            for i = 1 : numel(obj.unknowns)
                
                unknown = obj.unknowns{i};
                
                if( strcmp( unknown.type, 'evolution') )
                    
                    unknown_size = unknown.size();
                    lo = unknown.lo_global;
                    hi = unknown.hi_global;
                    out(os+1:os+unknown_size) = obj.fvec(lo:hi);
                    os = os + unknown_size;
                    
                end
                
            end
            
        end
        
        function insert_evolution_unknowns( obj, u_evol )
            
            os = 0;
            for i = 1 : numel(obj.unknowns)
                
                unknown = obj.unknowns{i};
                
                if( strcmp( unknown.type, 'evolution' ) )
                    
                    unknown_size = unknown.size();
                    lo = unknown.lo_global;
                    hi = unknown.hi_global;
                    obj.fvec(lo:hi) = u_evol(os+1:os+unknown_size);
                    os = os + unknown_size;
                    
                end
                
            end
            
        end
        
        function out = zeros( obj )
            
            out = zeros(size(obj.fvec));
            
        end
        
        function out = zeros_evolution( obj )
            
            out = zeros(obj.size_evolution_unknowns(),1);
            
        end
        
        function [ Qvec ] = get_input_unknowns( obj, Q, x )
            
            if( nargin < 3 )
                use_input_vector = false;
            else
                use_input_vector = true;
                assert( numel( obj.fvec ) == numel( x ) )
            end
            
            num_Q = numel(Q);
            
            Qvec = cell(num_Q,1);
            
            for i = 1 : num_Q
                
                lo = Q{i}.lo_global;
                hi = Q{i}.hi_global;
                
                if( use_input_vector )
                    Qvec{i} = x(lo:hi);
                else
                    Qvec{i} = obj.fvec(lo:hi);
                end
                
            end
            
        end
        
        function put_output_unknowns( obj, Q, Qvec )
            
            num_Q = numel(Q);
            
            for i = 1 : num_Q
                
                lo = Q{i}.lo_global;
                hi = Q{i}.hi_global;
                
                obj.fvec(lo:hi) = Qvec{i};
                
            end
            
        end
        
    end
end

