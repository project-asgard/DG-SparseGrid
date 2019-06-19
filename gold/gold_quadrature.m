%% Write C++ testing data for the quadrature component

quad_dir = strcat(pwd, "/", "generated-inputs", "/", "quadrature", "/");
mkdir (quad_dir);

% testing legendre poly/deriv

out_format = strcat(quad_dir, "legendre_");

x = [-1.0];
degree = 2;
[deriv, poly] = lin_dlegendre2(x, degree);
write_octave_like_output(strcat(out_format, 'deriv_neg1_2.dat'), deriv);
write_octave_like_output(strcat(out_format, 'poly_neg1_2.dat'), poly);

x = linspace(-2.5, 3.0, 11);
degree = 5;
[deriv, poly] = lin_dlegendre2(x, 5);
write_octave_like_output(strcat(out_format, 'deriv_linspace_5.dat'), deriv);
write_octave_like_output(strcat(out_format, 'poly_linspace_5.dat'), poly);

% testing legendre_weights

out_format = strcat(quad_dir, "lgwt_");
[roots, weights] = lgwt(10, -1, 1);
write_octave_like_output(strcat(out_format, 'roots_10_neg1_1.dat'), roots);
write_octave_like_output(strcat(out_format, 'weights_10_neg1_1.dat'), weights);

out_format = strcat(quad_dir, "lgwt_");
[roots, weights] = lgwt(2^5, -5, 2);
write_octave_like_output(strcat(out_format, 'roots_32_neg5_2.dat'), roots);
write_octave_like_output(strcat(out_format, 'weights_32_neg5_2.dat'), weights);

clear