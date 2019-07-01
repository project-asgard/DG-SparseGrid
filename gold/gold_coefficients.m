% coefficient testing
coeff_dir = strcat(pwd, "/", "generated-inputs", "/", "coefficients", "/");
[stat,msg] = mkdir (coeff_dir);

% continuity1 term
pde = check_pde(continuity1_old);
pde.dimensions{1}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{1}.lev,'wavelet');
time = 0;
mat = coeff_matrix_old(pde.deg, time, pde.dimensions{1}, pde.terms{1}{1});
write_octave_like_output(strcat(coeff_dir,'continuity1_coefficients.dat'), full(mat));

% continuity2 terms
pde = check_pde(continuity2_old);
level = 4;
degree = 3;
for i=1:length(pde.dimensions)
  pde.dimensions{i}.lev = level;
  pde.deg = degree;
  pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, 'continuity2_coefficients_l4_d3_%d_%d.dat');
%doesn't matter, the term is time independent...
time = 1.0;
for t=1:length(pde.terms)
  for d=1:length(pde.dimensions)
    coeff_mat = coeff_matrix_old(pde.deg,time,pde.dimensions{d},pde.terms{t}{d});
    write_octave_like_output(sprintf(out_format,t,d), full(coeff_mat));
  end
end

% continuity3 terms
pde = check_pde(continuity3_old);
level = 3;
degree = 5;
for i=1:length(pde.dimensions)
  pde.dimensions{i}.lev = level;
  pde.deg = degree;
  pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, 'continuity3_coefficients_l3_d5_%d_%d.dat');
%doesn't matter, the term is time independent...
time = 1.0;
for t=1:length(pde.terms)
  for d=1:length(pde.dimensions)
    coeff_mat = coeff_matrix_old(pde.deg,time,pde.dimensions{d},pde.terms{t}{d});
    write_octave_like_output(sprintf(out_format,t,d), full(coeff_mat));
  end
end

% continuity6 terms
pde = check_pde(continuity6_old);
level = 2;
degree = 4;
for i=1:length(pde.dimensions)
  pde.dimensions{i}.lev = level;
  pde.deg = degree;
  pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, 'continuity6_coefficients_l2_d4_%d_%d.dat');
%doesn't matter, the term is time independent...
time = 1.0;
for t=1:length(pde.terms)
  for d=1:length(pde.dimensions)
    coeff_mat = coeff_matrix_old(pde.deg,time,pde.dimensions{d},pde.terms{t}{d});
    write_octave_like_output(sprintf(out_format,t,d), full(coeff_mat));
  end
end
