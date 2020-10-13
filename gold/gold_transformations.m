%% Generate gold data for C++ testing of transformations component

% transformations testing 
data_dir = strcat("generated-inputs", "/", "transformations", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

degree = 1;
generate_multiwavelet_data( data_dir, degree );
degree = 2;
generate_multiwavelet_data( data_dir, degree );
degree = 3;
generate_multiwavelet_data( data_dir, degree );
degree = 4;
generate_multiwavelet_data( data_dir, degree );

% combine dimensions file generation
out_base = strcat(data_dir, "combine_dim_");
max_lev = 8;
filename = strcat(out_base, "dim2_deg2_lev3_sg.dat");


lev = 3;
dim = 2;
deg = 2;
lev_vec = zeros(dim,1)+lev;
grid_type = 'SG';
[elements, elements_idx]    = hash_table_sparse_nD (lev_vec, max_lev, grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;

fd = {[1:16], [17:32]};
ft = 2.0;
use_oldhash = false;
ans = combine_dimensions_D(deg, fd, ft, hash_table, use_oldhash);
write_octave_like_output(filename, full(ans));

filename = strcat(out_base, "dim3_deg3_lev2_fg.dat");
lev = 2;
dim = 3;
deg = 3;

lev_vec = zeros(dim,1)+lev;
grid_type = 'FG';
[elements, elements_idx]    = hash_table_sparse_nD (lev_vec, max_lev, grid_type);
hash_table.elements         = elements;
hash_table.elements_idx     = elements_idx;

fd = {[1:12], [13:24], [25:36]};
ft = 2.5;
use_oldhash = false;
ans = combine_dimensions_D(deg, fd, ft, hash_table, use_oldhash);
write_octave_like_output(filename, full(ans));


% operator two scale file generation
out_base = strcat(data_dir, "operator_two_scale_");


filename = strcat(out_base, "2_2.dat");
degree = 2;
level = 2;
vect = OperatorTwoScale(degree, level, 'wavelet');
write_octave_like_output(filename, full(vect));

filename = strcat(out_base, "2_3.dat");
degree = 2;
level = 3;
vect = OperatorTwoScale(degree, level, 'wavelet');
write_octave_like_output(filename, full(vect));

filename = strcat(out_base, "4_3.dat");
degree = 4;
level = 3;
vect = OperatorTwoScale(degree, level, 'wavelet');
write_octave_like_output(filename, full(vect));


filename = strcat(out_base, "5_5.dat");
degree = 5;
level = 5;
vect = OperatorTwoScale(degree, level, 'wavelet');
write_octave_like_output(filename, full(vect));

filename = strcat(out_base, "2_6.dat");
degree = 2;
level = 6;
vect = OperatorTwoScale(degree, level, 'wavelet');
write_octave_like_output(filename, full(vect));

% forward MWT test generation
out_base = strcat(data_dir, "forward_transform_");
max_lev = 8;

filename = strcat(out_base, "2_2_neg1_pos1_double.dat");
degree = 2;
level = 2;
l_min = -1.0;
l_max = 1.0;
double_it = @(x,p,t) (x*2);
[~, transform_blocks] = OperatorTwoScale_wavelet2(degree, max_lev);
vect = forward_wavelet_transform(degree, level, l_min, l_max, double_it,...
                                 [], transform_blocks, 0);
write_octave_like_output(filename, full(vect));


filename = strcat(out_base, "3_4_neg2_pos2_doubleplus.dat");
degree = 3;
level = 4;
l_min = -2.0;
l_max = 2.0;
double_plus = @(x,p,t) (x + x*2);
[~, transform_blocks] = OperatorTwoScale_wavelet2(degree, max_lev);
vect = forward_wavelet_transform(degree, level, l_min, l_max, double_plus,...
                                 [], transform_blocks, 0);
write_octave_like_output(filename, full(vect));

% multiwavelet file generation
function generate_multiwavelet_data( output_prefix, degree )

out_base = sprintf('%smultiwavelet_%d_', output_prefix, degree);
[h0,h1,g0,g1,scale_co,phi_co]=MultiwaveletGen(degree);
write_octave_like_output(strcat(out_base, "h0.dat"), h0);
write_octave_like_output(strcat(out_base, "h1.dat"), h1);
write_octave_like_output(strcat(out_base, "g0.dat"), g0);
write_octave_like_output(strcat(out_base, "g1.dat"), g1);
write_octave_like_output(strcat(out_base, "scale_co.dat"), scale_co);
write_octave_like_output(strcat(out_base, "phi_co.dat"), phi_co);

end





