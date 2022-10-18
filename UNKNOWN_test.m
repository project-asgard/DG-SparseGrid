function [ tests ] = UNKNOWN_test( )
    tests = functiontests( localfunctions );
end

function test_UNKNOWN_1( testCase )

    args = {'lev',[1,4],'deg',3,'grid_type','FG'};
    opts = OPTS(args);
    
    dim_x = DIMENSION(0,1,opts.lev(1));
    dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
    
    dim_v = DIMENSION(-6,+6,opts.lev(2));
    dim_v.moment_dV = @(v,p,t,dat) 0*v+1;

    dimensions_xv = {dim_x,dim_v};

    f = UNKNOWN( opts, dimensions_xv, {}, {} );
    
end