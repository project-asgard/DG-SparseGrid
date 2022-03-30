classdef PDE_TERMS
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        linear_terms = {};
        driver;
    end
    
    methods
        
        function [ solutions ] = pde_driver( obj, solutions )
            
            for i = 1 : numel(obj.linear_terms)
                
                data = obj.linear_terms{i};
                solutions{data.output_id(1)}.F_fval(:,data.output_id(2))...
                  = solutions{data.output_id(1)}.F_fval(:,data.output_id(2))...
                  + apply_A( data.term, solutions{data.input_id(1)}.fval(:,data.input_id(2)) );
              
            end
            
        end
        
        function obj = untitled5(inputArg1,inputArg2)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function pde_terms = PDE_TERMS( data )
            
            pde_terms.linear_terms = data;
            
        end
        
    end
end

