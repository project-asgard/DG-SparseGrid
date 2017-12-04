% EXAMPLE_INTERFACE Example MATLAB/Octave class wrapper to an underlying C++ class
classdef example_interface < cppinterface

    methods
        %% Constructor
        function this = example_interface()

            % initialise the cppinterface parent class by passing the 
            % mexfunction to the superclass constructor
            this = this@cppinterface(@MMeshMaker);
        end

        function refine (this,a)
            this.cppcall ('refine',a);
        end
        
        function print (this)
            this.cppcall ('print');
        end
        
         function init (this,a,b,c,d)
            this.cppcall ('init',a,b,c,d);
         end
         
         function x = findv (this, a )
             x = this.cppcall ('findv',a );
         end
         
          function x = findk (this, a )
             x = this.cppcall ('findk',a );
         end
    end

end
