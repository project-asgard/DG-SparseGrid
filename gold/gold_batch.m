%% Generate gold data for C++ testing of batch component

% batch
batch_dir = strcat(pwd, "/", "generated-inputs", "/", "batch", "/");
mkdir (batch_dir);

% continuity2 - sg
out_format = strcat(batch_dir, 'continuity2_sg_l2_d2_t%d.dat');
pde = continuity2;
level = 2;
degree = 2;
gridType='SG';

for i=1:length(pde.dimensions)
  pde.dimensions{i}.lev = level;
  pde.deg = degree;
  pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end
pde = check_pde(pde);
pde = check_terms(pde);

num_dimensions = length(pde.dimensions);
lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD(lev_vec,gridType);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;
opts.use_oldhash = 1;

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD);

A_data = global_matrix (pde,opts,HASHInv);

Vmax = 0;
Emax = 0;

f = ones(size(HASHInv,2) * degree^num_dimensions);

out = apply_A(pde,opts,A_data,f,degree,Vmax,Emax) ;
write_octave_like_output(sprintf(out_format,1), out);

% continuity2 - fg
out_format = strcat(batch_dir, 'continuity2_fg_l3_d4_t%d.dat');
pde = continuity2;
level = 3;
degree = 4;
gridType='FG';

for i=1:length(pde.dimensions)
  pde.dimensions{i}.lev = level;
  pde.deg = degree;
  pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end
pde = check_pde(pde);
pde = check_terms(pde);

num_dimensions = length(pde.dimensions);
lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD (lev_vec, gridType);

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;
opts.use_oldhash = 1;

A_data = global_matrix (pde,opts,HASHInv);

Vmax = 0;
Emax = 0;

f = ones(size(HASHInv,2) * degree^num_dimensions);

out = apply_A(pde,opts,A_data,f,degree,Vmax,Emax);
write_octave_like_output(sprintf(out_format,1), out);

% continuity3 - sg
out_format = strcat(batch_dir, 'continuity3_sg_l3_d4_t%d.dat');
pde = continuity3;
level = 3;
degree = 4;
gridType='SG';

for i=1:length(pde.dimensions)
  pde.dimensions{i}.lev = level;
  pde.deg = degree;
  pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end
pde = check_pde(pde);
pde = check_terms(pde);

num_dimensions = length(pde.dimensions);
lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD (lev_vec, gridType);

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;
opts.use_oldhash = 1;

A_data = global_matrix (pde,opts,HASHInv);

Vmax = 0;
Emax = 0;

f = ones(size(HASHInv,2) * degree^num_dimensions);

out = apply_A(pde,opts,A_data,f,degree,Vmax,Emax);
write_octave_like_output(sprintf(out_format,1), out);

% continuity6 - sg 
% TODO : why is this the old version / using old coeff_mat
% also when the rest of the tests above do not?

out_format = strcat(batch_dir, 'continuity6_sg_l2_d3_t%d.dat');
pde = continuity6;  
level = 2;
degree = 3;
gridType='SG';

for i=1:length(pde.dimensions)
  pde.dimensions{i}.lev = level;
  pde.deg = degree;
  pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
end

pde = check_pde(pde);
pde = check_terms(pde);

num_dimensions = length(pde.dimensions);
lev_vec = zeros(numel(pde.dimensions),1)+level;
[HASH,HASHInv] = hash_table_nD(lev_vec,gridType);

opts.compression = 4;
opts.useConnectivity = 0;
opts.implicit = 0;
opts.use_oldhash = 1;
opts.use_oldcoeffmat = 0;

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD,opts.use_oldcoeffmat);
pde.CFL=0.1;

A_data = global_matrix (pde,opts,HASHInv);

Vmax = 0;
Emax = 0;

f = ones(size(HASHInv,2) * degree^num_dimensions);

out = apply_A(pde,opts,A_data,f,degree,Vmax,Emax);
write_octave_like_output(sprintf(out_format,1), out);