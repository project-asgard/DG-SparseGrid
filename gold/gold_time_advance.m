%% Generate gold data for C++ testing of time advance component

opts.use_oldcoeffmat = 0;
opts.use_oldhash = 1;

% time advance
data_dir = strcat("generated-inputs", "/", "time_advance", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

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
CFL = .1;
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
CFL=.1;
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
CFL = 0.1;
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
CFL = .1;
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
CFL = 0.1;
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
CFL = .1;
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
CFL = .1;
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
CFL = 0.1;
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
CFL = .1;
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
CFL = .1;
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
CFL = .1;
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
CFL = .1;
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
CFL = .1;
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
CFL = .1;
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
CFL = .1;
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

