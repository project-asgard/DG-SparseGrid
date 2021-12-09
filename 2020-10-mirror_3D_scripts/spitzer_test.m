function spitzer_test()

% modify PDE
m_D_cgs = 3.3443*10^-24; %Deuterium mass in g
m_He_cgs = 6.7*10^-24; %helium 4 mass in g
m_B_cgs = 1.82*10^-23; %Boron 11 mass in g
m_vals = 0.001.*[m_D_cgs,m_He_cgs,m_B_cgs];
z_vals = [1,2,5];
vals.E = 9.17e-05;
vals.E_D = 73.84;
vals.Z = 1;
vals.m = m_vals(1);
vals.dt = 3.9e-6;
args = {'timestep_method','BE','case',3,'lev',4, 'deg', 4, 'num_steps', 5,'quiet',true,'save_output',true, 'cmd_args', vals,'update_params_each_timestep', true};
opts = OPTS(args);
pde = mirror2_velocity_div(opts);
ratio_max = 0.1;
N = 15;
ratio = linspace(0.02,ratio_max,N);

for i = 1:length(ratio)
    vals.E = ratio(i)*pde.params.E_dreicer_si;
    vals.dt = 3.9e-6*vals.E/vals.E_D;
    vals.num_steps = 12*3.9e-6/vals.dt;
    args{16} = vals;
    opts = OPTS(args);
    pde = mirror2_velocity_div(opts);
   [~,~,~,nodes,~,outputs(i)] = asgard_run_pde(opts,pde);
   alpha_nr(i,vals.Z) = outputs(i).alpha_t{end};
    %plot_mirror_output(nodes,outputs,pde,opts);
    %legend( sprintf('Z = %g', vals.Z) );
end
semilogy(ratio,alpha_nr(:,1),'LineWidth',5,'DisplayName','ASGarD (nr, Z=1)','color','red');

end