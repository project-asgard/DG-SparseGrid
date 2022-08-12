classdef VELOCITY_MOMENT_TERM < PARTIAL_SD_TERM
    
    properties
    end

    methods

        function obj = VELOCITY_MOMENT_TERM( g_, dV_ )
            
            g = @(x,p,t,dat) x.*0+1;
            if exist( 'g_', 'var' )
                if ~isempty(g_)
                    g = g_;
                end
            end
            
            dV = @(x,p,t,d) x.*0+1;
            if exist( 'dV_', 'var' )
                if ~isempty(dV_)
                    dV = dV_;
                end
            end
            
            obj@PARTIAL_SD_TERM( 'moment', g, '', '', dV )
            
        end

    end
end