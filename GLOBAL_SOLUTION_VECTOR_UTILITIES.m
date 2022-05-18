classdef GLOBAL_SOLUTION_VECTOR_UTILITIES
    
    methods
        
        function global_solution_vector_utilities = GLOBAL_SOLUTION_VECTOR_UTILITIES()
            
        end
        
        function [ gsv_evol ] = checkout_evolution( obj, gsv, u_evol )
            
            gsv_evol = GLOBAL_SOLUTION_VECTOR( gsv.unknowns );
            
            index_lo = gsv.unknown_index_lo;
            index_hi = gsv.unknown_index_hi;
            
            os = 0;
            for i = 1 : numel(gsv.unknowns)
                
                if( strcmp(gsv.unknowns{i}.type,'evolution') )
                    
                    unknown_size = gsv.unknowns{i}.size();
                    
                    gsv_evol.fvec(index_lo(i):index_hi(i))...
                        = u_evol(os+1:os+unknown_size);
                    
                    os = os + unknown_size;
                    
                end
                
            end
            
        end
        
        function [ u_evol ] = extract_evolution( obj, gsv_evol )
            
            u_evol = zeros(gsv_evol.size_evolution(),1);
            
            index_lo = gsv_evol.unknown_index_lo;
            index_hi = gsv_evol.unknown_index_hi;
            
            os = 0;
            for i = 1 : numel(gsv_evol.unknowns)
                
                if( strcmp(gsv_evol.unknowns{i}.type,'evolution') )
                    
                    unknown_size = gsv_evol.unknowns{i}.size();
                    
                    u_evol(os+1:os+unknown_size)...
                        = gsv_evol.fvec(index_lo(i):index_hi(i));
                    
                    os = os + unknown_size;
                    
                end
                
            end
            
            gsv_evol.delete();
            
        end
        
    end
end

