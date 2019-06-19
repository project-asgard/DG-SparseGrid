%% Generate gold data for C++ testing of transformations component

% transformations testing 
transformations_dir = strcat(pwd, "/", "generated-inputs", "/", "transformations", "/");
mkdir (transformations_dir);

% multiwavelet file generation
out_base = strcat(transformations_dir, "multiwavelet_1_");
[h0,h1,g0,g1,scale_co,phi_co]=MultiwaveletGen(1);
save(strcat(out_base, "h0.dat"), 'h0');
save(strcat(out_base, "h1.dat"), 'h1');
save(strcat(out_base, "g0.dat"), 'g0');
save(strcat(out_base, "g1.dat"), 'g1');
save(strcat(out_base, "scale_co.dat"), 'scale_co');
save(strcat(out_base, "phi_co.dat"), 'phi_co');


out_base = strcat(transformations_dir, "multiwavelet_3_");
[h0,h1,g0,g1,scale_co,phi_co]=MultiwaveletGen(3);
save(strcat(out_base, "h0.dat"), 'h0');
save(strcat(out_base, "h1.dat"), 'h1');
save(strcat(out_base, "g0.dat"), 'g0');
save(strcat(out_base, "g1.dat"), 'g1');
save(strcat(out_base, "scale_co.dat"), 'scale_co');
save(strcat(out_base, "phi_co.dat"), 'phi_co');

% combine dimensions file generation
out_base = strcat(transformations_dir, "combine_dim_");

filename = strcat(out_base, "dim2_deg2_lev3_sg.dat");
lev = 3;
dim = 2;
deg = 2;

[fwd, back] = HashTable(lev, dim, 'SG', 1);
fd = {[1:16], [17:32]};
ft = 2.0;

ans = combine_dimensions_D(fd, ft, back, deg);
save(filename, "ans");

filename = strcat(out_base, "dim3_deg3_lev2_fg.dat");
lev = 2;
dim = 3;
deg = 3;

[fwd, back] = HashTable(lev, dim, 'FG', 1);
fd = {[1:12], [13:24], [25:36]};
ft = 2.5;

ans = combine_dimensions_D(fd, ft, back, deg);
save(filename, "ans");


% operator two scale file generation
out_base = strcat(transformations_dir, "operator_two_scale_");


filename = strcat(out_base, "2_2.dat");
degree = 2;
level = 2;
vect = OperatorTwoScale(degree, 2^level);
save(filename, 'vect');

filename = strcat(out_base, "2_3.dat");
degree = 2;
level = 3;
vect = OperatorTwoScale(degree, 2^level);
save(filename, 'vect');

filename = strcat(out_base, "4_3.dat");
degree = 4;
level = 3;
vect = OperatorTwoScale(degree, 2^level);
save(filename, 'vect');


filename = strcat(out_base, "5_5.dat");
degree = 5;
level = 5;
vect = OperatorTwoScale(degree, 2^level);
save(filename, 'vect');

filename = strcat(out_base, "2_6.dat");
degree = 2;
level = 6;
vect = OperatorTwoScale(degree, 2^level);
save(filename, 'vect');

% forward MWT test generation
out_base = strcat(transformations_dir, "forward_transform_");

filename = strcat(out_base, "2_2_neg1_pos1_double.dat");
degree = 2;
level = 2;
l_min = -1.0;
l_max = 1.0;
double_it = @(n,x) (n*2);
vect = forwardMWT(level, degree, l_min, l_max, double_it, 0);
save(filename, 'vect');


filename = strcat(out_base, "3_4_neg2_pos2_doubleplus.dat");
degree = 3;
level = 4;
l_min = -2.0;
l_max = 2.0;
double_plus = @(n,x) (n + n*2);
vect = forwardMWT(level, degree, l_min, l_max, double_plus, 0);
save(filename, 'vect');