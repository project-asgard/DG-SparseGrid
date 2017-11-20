% EXAMPLE_INTERFACE Example MATLAB/Octave class wrapper to an underlying C++ class
classdef Mesh_interface < cppinterface

    methods
        %% Constructor
        function this = Mesh_interface()

            % initialise the cppinterface parent class by passing the 
            % mexfunction to the superclass constructor
            this = this@cppinterface(@MMeshMaker);
        end

        function refine (this,a,b,c)
            this.cppcall ('refine',a,b,c);
        end
        
        function print (this)
            this.cppcall ('print');
        end
        
         function init (this,a,b,c,d)
            this.cppcall ('init',a,b,c,d);
         end
         
         function x = findvs (this, a )
             x = this.cppcall ('findvs',a );
         end
         
          function x = findk (this, a )
             x = this.cppcall ('findk',a );
          end
         
          function x = leafv (this )
             x = this.cppcall ('leafv');
          end
           function x = allv (this )
             x = this.cppcall ('allv');
          end 
         function x = size (this )
             x = this.cppcall ('size' );
         end
    end

end
