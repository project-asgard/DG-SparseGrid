function spitzer_test()

% modify PDE
m_D_cgs = 3.3443*10^-24; %Deuterium mass in g
m_He_cgs = 6.7*10^-24; %helium 4 mass in g
m_B_cgs = 1.82*10^-23; %Boron 11 mass in g
m_vals = 0.001.*[m_D_cgs,m_He_cgs,m_B_cgs];
z_vals = [1,2,5];
% run PDE

for i = 1:length(m_vals)
    args = {'timestep_method','BE','case',3,'dt',1e-5,'lev',4, 'deg', 3, 'num_steps', 5,'quiet',true,'save_output',true, 'update_params_each_timestep', true};
    opts = OPTS(args);
    pde = mirror2_velocity_div(opts);
    pde.params.b.m = m_vals(i);
    pde.params.b.Z = z_vals(i);
   [err,fval,fval_realspace,nodes,err_realspace] = asgard_run_pde(opts,pde);
   load('asgard-out-l4-d3-SG-dt1e-05-adapt-n.mat');
    plot_mirror_output(nodes,outputs,pde,opts);
    hold on
end
hold off
end