%% Generate gold data for C++ testing of time advance component

% time advance
data_dir = strcat("generated-inputs", "/", "time_advance", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

% diffusion1 l4d4
out_format = strcat(data_dir,'diffusion1_sg_l4_d4_t%d.dat');
run_time_advance(diffusion1,out_format,...
    'CFL',0.01,'lev',4,'deg',4,'grid_type','SG','use_oldhash',true,'timestep_method','BE');
% advection1 l4d4
out_format = strcat(data_dir,'advection1_sg_l4_d4_t%d.dat');
run_time_advance(advection1,out_format,...
    'CFL',0.01,'lev',4,'deg',4,'grid_type','SG','use_oldhash',true,'timestep_method','BE');

% continuity1 sg l2d2
out_format = strcat(data_dir, 'continuity1_sg_l2_d2_t%d.dat');
run_time_advance(continuity1,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');
% continuity1 fg l2d2
out_format = strcat(data_dir, 'continuity1_fg_l2_d2_t%d.dat');
run_time_advance(continuity1,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','FG','use_oldhash',true,'timestep_method','RK3');
% continuity1 sg l4d3
out_format = strcat(data_dir, 'continuity1_sg_l4_d3_t%d.dat');
run_time_advance(continuity1,out_format,...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');

% continuity2 sg l2d2
out_format = strcat(data_dir, 'continuity2_sg_l2_d2_t%d.dat');
run_time_advance(continuity2,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');
% continuity2 fg l2d2
out_format = strcat(data_dir, 'continuity2_fg_l2_d2_t%d.dat');
run_time_advance(continuity2,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','FG','use_oldhash',true,'timestep_method','RK3');
% continuity2 sg l4d3
out_format = strcat(data_dir, 'continuity2_sg_l4_d3_t%d.dat');
run_time_advance(continuity2,out_format,...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');

% continuity3 sg l2d2
out_format = strcat(data_dir,'continuity3_sg_l2_d2_t%d.dat');
run_time_advance(continuity3,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');
% continuity3 sg l4d3
out_format = strcat(data_dir,'continuity3_sg_l4_d3_t%d.dat');
run_time_advance(continuity3,out_format,...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');

% continuity6 sg l2d2
out_format = strcat(data_dir, 'continuity6_sg_l2_d2_t%d.dat');
run_time_advance(continuity6,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');
% continuity6 sg l2d3
out_format = strcat(data_dir, 'continuity6_sg_l2_d3_t%d.dat');
run_time_advance(continuity6,out_format,...
    'CFL',0.01,'lev',2,'deg',3,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');

% fokkerplanck1_4p2 sg l2d2
out_format = strcat(data_dir,'fokkerplanck1_4p2_sg_l2_d2_t%d.dat');
run_time_advance(fokkerplanck1_4p2,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');
% fokkerplanck1_4p2 fg l2d2
out_format = strcat(data_dir,'fokkerplanck1_4p2_fg_l2_d2_t%d.dat');
run_time_advance(fokkerplanck1_4p2,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','FG','use_oldhash',true,'timestep_method','RK3');

% fokkerplanck1_4p1a sg l2d2
out_format = strcat(data_dir,'fokkerplanck1_4p1a_sg_l2_d2_t%d.dat');
run_time_advance(fokkerplanck1_4p1a,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');

% fokkerplanck1_4p3 sg l2d2
out_format = strcat(data_dir,'fokkerplanck1_4p3_sg_l2_d2_t%d.dat');
run_time_advance(fokkerplanck1_4p3,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','use_oldhash',true,'timestep_method','RK3');


%% implicit

% fokkerplanck2_complete sg l3d3

out_format = strcat(data_dir,'fokkerplanck2_complete_implicit_sg_l3_d3_t%d.dat');
% any other timestep method will give diverging results
run_time_advance(fokkerplanck2_complete,out_format,...
    'CFL',0.01,'lev',3,'deg',3,'grid_type','SG','use_oldhash',true,...
    'timestep_method', 'BE');

out_format = strcat(data_dir,'fokkerplanck2_complete_implicit_sg_l4_d3_t%d.dat');
run_time_advance(fokkerplanck2_complete,out_format,...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','use_oldhash',true,...
    'timestep_method', 'BE');

out_format = strcat(data_dir,'fokkerplanck2_complete_implicit_sg_l5_d3_t%d.dat');
run_time_advance(fokkerplanck2_complete,out_format,...
    'CFL',0.01,'lev',5,'deg',3,'grid_type','SG','use_oldhash',true,...
    'timestep_method', 'BE');

% continuity1 sg l2d2
out_format = strcat(data_dir,'continuity1_implicit_l2_d2_t%d.dat');
run_time_advance(continuity1,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','use_oldhash',true,'timestep_method','BE');
% continuity1 sg l4d3
out_format = strcat(data_dir,'continuity1_implicit_l4_d3_t%d.dat');
run_time_advance(continuity1,out_format,...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','use_oldhash',true,'timestep_method','BE');

% continuity2 sg l2d2
out_format = strcat(data_dir,'continuity2_implicit_l2_d2_t%d.dat');
run_time_advance(continuity2,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','use_oldhash',true,'timestep_method','BE');
% continuity2 sg l4d3
out_format = strcat(data_dir,'continuity2_implicit_l4_d3_t%d.dat');
run_time_advance(continuity2,out_format,...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','use_oldhash',true,'timestep_method','BE');

function run_time_advance(pde,out_format,varargin)

  runtime_defaults
    
  pde = check_pde( pde, opts );
  dt = pde.set_dt(pde, pde.CFL );

  lev_vec = zeros(numel(pde.dimensions),1);
  for d = 1 : numel(pde.dimensions)
    lev_vec(d) = pde.dimensions{d}.lev;
  end
  [HASH,HASHInv] = hash_table_nD(lev_vec, opts.grid_type);

  for i=1:length(pde.dimensions)
      pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
  end

  t = 0;
  TD = 0;
  pde = get_coeff_mats(pde,t,TD,opts.use_oldcoeffmat);

  A_data = global_matrix(pde,opts,HASHInv);

  Vmax = 0;
  Emax = 0;
  out = initial_condition_vector(pde,opts,HASHInv,0);

  for i=0:4
      time = i*dt;
      out = time_advance(pde,opts,A_data,out,time,dt,pde.deg,HASHInv,Vmax,Emax);
      write_octave_like_output(sprintf(out_format,i), out);
  end
end

function generate_data(pde, output_prefix, varargin)

runtime_defaults

opts.use_oldcoeffmat = 0;
opts.use_oldhash = 1;
opts.compression = 4;
opts.useConnectivity = 0;

pde = check_pde(pde, opts);

dt = pde.set_dt(pde, opts.CFL); 

% time advance
data_dir = strcat("generated-inputs/time_advance/", output_prefix, "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/', char(data_dir)]);

level = pde.dimensions{1}.lev;
degree = pde.deg;
grid_type= opts.grid_type;

% order of sprintf args: output_prefix, implicit, grid_type.to_lower(), level, degree, time
out_format = strcat(data_dir, '%s_%s_%s_l%d_d%d_t%d.dat');

if opts.implicit == 0
  implicit_string = 'e';
else
  implicit_string = 'i';
end

for i=1:length(pde.dimensions)
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

for i=1:numel(pde.dimensions)
  lev_vec( i ) = pde.dimensions{i}.lev;
end

%lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD(lev_vec,grid_type);

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD,opts.use_oldcoeffmat);

A_data = global_matrix(pde,opts,HASHInv);

Vmax = 0;
Emax = 0;
out = initial_condition_vector(pde,opts,HASHInv,0);
deg = pde.deg;
for i=0:4
    time = i*dt;
    out = time_advance(pde,opts,A_data,out,time,dt,deg,HASHInv,Vmax,Emax);
    write_octave_like_output(...
    sprintf(out_format, output_prefix, implicit_string, lower(grid_type), level, degree, i),...
    out);
end

end
