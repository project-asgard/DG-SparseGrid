%% Generate gold data for C++ testing of time advance component

% fokkerplanck2_complete l3d3
out_format = strcat(data_dir, 'fokkerplanck2_complete_sg_l%i_d%i_t%d.dat');
fokkerplanck2_complete_cfl = 1e-10;
generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete',...
                  'CFL', fokkerplanck2_complete_cfl, 'use_oldhash', true, 'lev', 2, 'deg', 3,...
               'timestep_method', 'RK3' );

generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete',...
                  'CFL', fokkerplanck2_complete_cfl, 'use_oldhash', true, 'lev', 3, 'deg', 4,...
               'timestep_method', 'RK3' );

generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete',...
                  'CFL', fokkerplanck2_complete_cfl, 'use_oldhash', true, 'lev', 4, 'deg', 3,...
               'timestep_method', 'RK3' );

generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete',...
                  'CFL', fokkerplanck2_complete_cfl, 'use_oldhash', true, 'lev', 2, 'deg', 3,...
               'timestep_method', 'BE' );

generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete',...
                  'CFL', fokkerplanck2_complete_cfl, 'use_oldhash', true, 'lev', 3, 'deg', 4,...
               'timestep_method', 'BE' );

generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete',...
                  'CFL', fokkerplanck2_complete_cfl, 'use_oldhash', true, 'lev', 4, 'deg', 3,...
               'timestep_method', 'BE' );

%Temporarily cut off data generation below this line - needs updating
return;

%diffusion1 implicit
diffusion_cfl = 0.003;
generate_data( diffusion1, 'diffusion1', ...
               'lev', 2, 'deg', 2, 'grid_type', 'SG', 'timestep_method', 'BE',...
               'CFL', diffusion_cfl );
generate_data( diffusion1, 'diffusion1', ...
               'lev', 3, 'deg', 3, 'grid_type', 'SG', 'timestep_method', 'BE',...
               'CFL', diffusion_cfl );
generate_data( diffusion1, 'diffusion1', ...
               'lev', 4, 'deg', 4, 'grid_type', 'SG', 'timestep_method', 'BE',...
               'CFL', diffusion_cfl );

%diffusion1 explicit
generate_data( diffusion1, 'diffusion1', 'lev', 2, 'deg', 2, 'grid_type', 'SG',...
               'CFL', diffusion_cfl );
generate_data( diffusion1, 'diffusion1', 'lev', 3, 'deg', 3, 'grid_type', 'SG',...
               'CFL', diffusion_cfl );
generate_data( diffusion1, 'diffusion1', 'lev', 4, 'deg', 4, 'grid_type', 'SG',...
               'CFL', diffusion_cfl );

%diffusion2 implicit
generate_data( diffusion2, 'diffusion2', ...
               'lev', 2, 'deg', 2, 'grid_type', 'SG', 'timestep_method', 'BE',...
               'CFL', diffusion_cfl );

generate_data( diffusion2, 'diffusion2', ...
               'lev', 3, 'deg', 3, 'grid_type', 'SG', 'timestep_method', 'BE',...
               'CFL', diffusion_cfl );

generate_data( diffusion2, 'diffusion2', ...
               'lev', 4, 'deg', 4, 'grid_type', 'SG', 'timestep_method', 'BE',...
               'CFL', diffusion_cfl );

%diffusion2 explicit
generate_data( diffusion2, 'diffusion2', 'lev', 2, 'deg', 2, 'grid_type', 'SG',...
               'CFL', diffusion_cfl );
generate_data( diffusion2, 'diffusion2', 'lev', 3, 'deg', 3, 'grid_type', 'SG',...
               'CFL', diffusion_cfl );
generate_data( diffusion2, 'diffusion2', 'lev', 4, 'deg', 4, 'grid_type', 'SG',...
               'CFL', diffusion_cfl );


opts.use_oldcoeffmat = 0;
opts.use_oldhash = 1;

% time advance
data_dir = strcat("generated-inputs", "/", "time_advance", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

% diffusion1 l4d4
out_format = strcat( data_dir, 'diffusion1_sg_l4_d4_t%d.dat' );
run_time_advance(diffusion1, out_format, ...
'CFL', 0.01, 'lev', 4, 'deg', 4, 'grid_type', 'SG', 'implicit', false, 'use_oldhash', true);

% advection1 l4d4
out_format = strcat( data_dir, 'advection1_sg_l4_d4_t%d.dat');
run_time_advance( advection1, out_format, ...
'CFL', 0.01, 'lev', 4, 'deg', 4, 'grid_type', 'SG', 'implicit', false, 'use_oldhash', true);

% continuity1

%sg l2d2
out_format = strcat(data_dir, 'continuity1_sg_l2_d2_t%d.dat');
pde = continuity1;
level = 2;
degree = 2;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,i), out);
end

%fg l2d2
out_format = strcat(data_dir, 'continuity1_fg_l2_d2_t%d.dat');
pde = continuity1;

level = 2;
degree = 2;
grid_type='FG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD(lev_vec,grid_type);

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD,opts.use_oldcoeffmat);

A_data = global_matrix(pde,opts,HASHInv);

Vmax = 0;
Emax = 0;
out = initial_condition_vector(pde,opts,HASHInv,0);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);
deg = pde.deg;
for i=0:4
    time = i*dt;
    out = time_advance(pde,opts,A_data,out,time,dt,deg,HASHInv,Vmax,Emax);
    write_octave_like_output(sprintf(out_format,i), out);
end

%sg l4d3
out_format = strcat(data_dir, 'continuity1_sg_l4_d3_t%d.dat');
pde = continuity1;

level = 4;
degree = 3;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD(lev_vec,grid_type);

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD,opts.use_oldcoeffmat);

A_data = global_matrix(pde,opts,HASHInv);

Vmax = 0;
Emax = 0;
out = initial_condition_vector(pde,opts,HASHInv,0);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);
deg = pde.deg;
for i=0:4
    time = i*dt;
    out = time_advance(pde,opts,A_data,out,time,dt,deg,HASHInv,Vmax,Emax);
    write_octave_like_output(sprintf(out_format,i), out);
end

% continuity2

%sg l2d2
out_format = strcat(data_dir, 'continuity2_sg_l2_d2_t%d.dat');
pde = continuity2;

level = 2;
degree = 2;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,i), out);
end

%fg l2d2
out_format = strcat(data_dir, 'continuity2_fg_l2_d2_t%d.dat');
pde = continuity2;

level = 2;
degree = 2;
grid_type='FG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD(lev_vec,grid_type);

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD,opts.use_oldcoeffmat);

A_data = global_matrix(pde,opts,HASHInv);

Vmax = 0;
Emax = 0;
out = initial_condition_vector(pde,opts,HASHInv,0);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);
deg = pde.deg;
for i=0:4
    time = i*dt;
    out = time_advance(pde,opts,A_data,out,time,dt,deg,HASHInv,Vmax,Emax);
    write_octave_like_output(sprintf(out_format,i), out);
end

%sg l4d3
out_format = strcat(data_dir, 'continuity2_sg_l4_d3_t%d.dat');
pde = continuity2;
level = 4;
degree = 3;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,i), out);
end

% continuity3

%sg l2d2
out_format = strcat(data_dir, 'continuity3_sg_l2_d2_t%d.dat');
pde = continuity3;

level = 2;
degree = 2;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,i), out);
end

%sg l4d3
out_format = strcat(data_dir, 'continuity3_sg_l4_d3_t%d.dat');
pde = continuity3;
level = 4;
degree = 3;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,i), out);
end

% continuity6

%sg l2d2
out_format = strcat(data_dir, 'continuity6_sg_l2_d2_t%d.dat');
pde = continuity6;

level = 2;
degree = 2;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD,opts.use_oldcoeffmat);
pde.CFL=0.1;

lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD(lev_vec,grid_type);

A_data = global_matrix(pde,opts,HASHInv);

Vmax = 0;
Emax = 0;
out = initial_condition_vector(pde,opts,HASHInv,0);
deg = pde.deg;
for i=0:4
    time = i*dt;
    out = time_advance(pde,opts,A_data,out,time,dt,deg,HASHInv,Vmax,Emax);
    write_octave_like_output(sprintf(out_format,i), out);
end

%sg l2d3
out_format = strcat(data_dir, 'continuity6_sg_l2_d3_t%d.dat');
pde = continuity6;

level = 2;
degree = 3;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD,opts.use_oldcoeffmat);

lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD(lev_vec,grid_type);

A_data = global_matrix(pde,opts,HASHInv);

Vmax = 0;
Emax = 0;
out = initial_condition_vector(pde,opts,HASHInv,0);
deg = pde.deg;
for i=0:4
    time = i*dt;
    out = time_advance(pde,opts,A_data,out,time,dt,deg,HASHInv,Vmax,Emax);
    write_octave_like_output(sprintf(out_format,i), out);
end

% fokkerplanck1_4p2

%sg l2d2
out_format = strcat(data_dir, 'fokkerplanck1_4p2_sg_l%i_d%i_t%d.dat');
pde = fokkerplanck1_4p2;
level = 2;
degree = 2;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,level,degree,i), out);
end

%sg l2d2
out_format = strcat(data_dir, 'fokkerplanck1_4p2_fg_l%i_d%i_t%d.dat');
pde = fokkerplanck1_4p2;
level = 2;
degree = 2;
grid_type='FG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;

lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,level,degree,i), out);
end


% fokkerplanck1_4p1a

% sg l2d2
lev = 2;
deg = 2;
for i=1:5
    [err,fval,fval_realspace] = asgard (fokkerplanck1_4p1a, ...
        'lev',lev, 'deg',deg, 'CFL', 0.01, 'use_oldhash', true, 'num_steps', i, 'quiet', true);
    out_format = strcat(data_dir, 'fokkerplanck1_4p1a_sg_l%i_d%i_t%d.dat');
    write_octave_like_output(sprintf(out_format,lev,deg,i-1), full(fval));
end

% fokkerplanck1_4p3

% sg l2d2
lev = 2;
deg = 2;
for i=1:5
    [err,fval,fval_realspace] = asgard (fokkerplanck1_4p3, ...
        'lev',lev, 'deg',deg, 'CFL', 0.01, 'use_oldhash', true, 'num_steps', i, 'quiet', true);
    out_format = strcat(data_dir, 'fokkerplanck1_4p3_sg_l%i_d%i_t%d.dat');
    write_octave_like_output(sprintf(out_format,lev,deg,i-1), full(fval));
end

% fokkerplanck2_complete

% sg l3d3
lev = 3;
deg = 3;
for i=1:5
    [err,fval,fval_realspace] = asgard (fokkerplanck2_complete, ...
        'lev',lev, 'deg',deg, 'CFL', 1e-10, 'use_oldhash', true, 'num_steps', i, 'quiet', true);
    out_format = strcat(data_dir, 'fokkerplanck2_complete_sg_l%i_d%i_t%d.dat');
    write_octave_like_output(sprintf(out_format,lev,deg,i-1), full(fval));
end


%% implicit

% continuity1

%sg l2d2
out_format = strcat(data_dir, 'continuity1_implicit_l2_d2_t%d.dat');
pde = continuity1;
level = 2;
degree = 2;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 1;
opts.implicit_method="BE";
lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,i), out);
end

%sg l4d3
out_format = strcat(data_dir, 'continuity1_implicit_l4_d3_t%d.dat');
pde = continuity1;
level = 4;
degree = 3;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 1;
opts.implicit_method="BE";
lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,i), out);
end

% continuity2

%sg l2d2
out_format = strcat(data_dir, 'continuity2_implicit_l2_d2_t%d.dat');
pde = continuity2;
level = 2;
degree = 2;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 1;
opts.implicit_method="BE";

lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,i), out);
end

%sg l4d3
out_format = strcat(data_dir, 'continuity2_implicit_l4_d3_t%d.dat');
pde = continuity2;
level = 4;
degree = 3;
grid_type='SG';

for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
CFL = 0.01;
dt = pde.set_dt(pde,CFL);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 1;
opts.implicit_method="BE";

lev_vec = zeros(numel(pde.dimensions),1)+level;
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
    write_octave_like_output(sprintf(out_format,i), out);
end

function run_time_advance( pde, out_format, varargin )

  runtime_defaults
  pde = check_pde( pde, opts );
  dt = pde.set_dt(pde, pde.CFL );

  opts.compression = 4;
  opts.useConnectivity = 0;
  opts.implicit_method="BE";

  lev_vec = zeros(numel(pde.dimensions),1);
  for d = 1 : numel(pde.dimensions)
    lev_vec(d) = pde.dimensions{d}.lev
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

if( strcmp(opts.timestep_method, 'RK3') )
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
