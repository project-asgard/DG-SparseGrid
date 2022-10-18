function tests = asgard_pde_system_test()

    fh = localfunctions;
    tests = functiontests(fh);

end

function asgard_diffusion_system1_explicit_test(testCase)
    addpath(genpath(pwd));
    disp( 'UNKNOWN_test' )
    UNKNOWN_test();
    disp( 'Testing diffusion_system1 (SSPRK3)' );
    [ err ] = asgard_pde_system...
                ( @diffusion_system1, 'lev', 4, 'deg', 3, 'CFL', 0.1, 'quiet', true, 'grid_type', 'FG',...
                  'fast_FG_matrix_assembly', false, 'timestep_method', 'SSPRK3', 'num_steps', 10 );
    verifyLessThan(testCase, err, 1e-3);
    disp( 'Testing euler_system1 (LinearAdvection, SSPRK3)' );
    [ err ] = asgard_pde_system...
                ( @euler_system1, 'lev', 4, 'deg', 3, 'CFL', 0.1, 'quiet', true, 'grid_type', 'FG',...
                  'fast_FG_matrix_assembly', false, 'timestep_method', 'SSPRK3', 'num_steps', 10,...
                  'case', 1 );
    verifyLessThan(testCase, err, 1e-3);
end