% coefficient testing
coeff_dir = strcat("generated-inputs", "/", "coefficients", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(coeff_dir)]);

% continuity1 term
pde = check_pde(continuity1);
pde.dimensions{1}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{1}.lev,'wavelet');
time = 0;
sd_term = pde.terms{1}.terms_1D{1};
sd_term_out = coeff_matrix(1,pde.deg, time, pde.dimensions{1}, sd_term, pde.params);
write_octave_like_output(strcat(coeff_dir,'continuity1_coefficients.dat'), full(sd_term_out.mat));

% continuity2 terms
pde = check_pde(continuity2);
level = 4;
degree = 3;
for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, 'continuity2_coefficients_l%i_d%i_%d_%d.dat');
%doesn't matter, the term is time independent...
time = 1.0;
for t=1:length(pde.terms)
    for d=1:length(pde.dimensions)
        sd_term = pde.terms{t}.terms_1D{d};
        sd_term_out = coeff_matrix(2,pde.deg,time,pde.dimensions{d},sd_term,pde.params);
        coeff_mat = sd_term_out.mat;
        write_octave_like_output(sprintf(out_format,level,degree,t,d), full(coeff_mat));
    end
end

% continuity3 terms
pde = check_pde(continuity3);
level = 4;
degree = 4;
for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, 'continuity3_coefficients_l%i_d%i_%d_%d.dat');
out_format0 = strcat(coeff_dir, 'continuity3_coefficients_norotate_l%i_d%i_%d_%d.dat');
%doesn't matter, the term is time independent...
time = 1.0;
for t=1:length(pde.terms)
    for d=1:length(pde.dimensions)
        sd_term = pde.terms{t}.terms_1D{d};
        sd_term_out = coeff_matrix(3,pde.deg,time,pde.dimensions{d},sd_term,pde.params);
        mat = sd_term_out.mat;
        mat0 = sd_term_out.mat_unrotated;
        write_octave_like_output(sprintf(out_format,level,degree,t,d), full(mat));
        write_octave_like_output(sprintf(out_format0,level,degree,t,d), full(mat0));
    end
end

% continuity6 terms
pde = check_pde(continuity6);
level = 2;
degree = 4;
for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, 'continuity6_coefficients_l%i_d%i_%d_%d.dat');
%doesn't matter, the term is time independent...
time = 1.0;
for t=1:length(pde.terms)
    for d=1:length(pde.dimensions)
        sd_term = pde.terms{t}.terms_1D{d};
        sd_term_out = coeff_matrix(6,pde.deg,time,pde.dimensions{d},sd_term,pde.params);
        mat = sd_term_out.mat;
        write_octave_like_output(sprintf(out_format,level,degree,t,d), full(mat));
    end
end

% fokkerplanck1_4p2 terms
pde = check_pde(fokkerplanck1_4p2);
level = 3;
degree = 4;
for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, 'fokkerplanck1_4p2_coefficients_l%i_d%i_%d_%d.dat');
out_format0 = strcat(coeff_dir, 'fokkerplanck1_4p2_coefficients_norotate_l%i_d%i_%d_%d.dat');
%doesn't matter, the term is time independent...
time = 1.0;
for t=1:length(pde.terms)
    for d=1:length(pde.dimensions)
        sd_term = pde.terms{t}.terms_1D{d};
        sd_term_out = coeff_matrix(numel(pde.dimensions),pde.deg,time,pde.dimensions{d},sd_term,pde.params);
        mat = sd_term_out.mat;
        mat0 = sd_term_out.mat_unrotated;
        write_octave_like_output(sprintf(out_format,level,degree,t,d), full(mat));
        write_octave_like_output(sprintf(out_format0,level,degree,t,d), full(mat0));
    end
end

% fokkerplanck2_complete terms
pde = check_pde(fokkerplanck2_complete);
level = 3;
degree = 4;
for i=1:length(pde.dimensions)
    pde.dimensions{i}.lev = level;
    pde.deg = degree;
    pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, 'fokkerplanck2_complete_coefficients_l%i_d%i_%d_%d.dat');
out_format0 = strcat(coeff_dir, 'fokkerplanck2_complete_coefficients_norotate_l%i_d%i_%d_%d.dat');
%doesn't matter, the term is time independent...
time = 1.0;
for t=1:length(pde.terms)
    for d=1:length(pde.dimensions)
        sd_term = pde.terms{t}.terms_1D{d};
        sd_term_out = coeff_matrix(numel(pde.dimensions),pde.deg,time,pde.dimensions{d},sd_term,pde.params);
        mat = sd_term_out.mat;
        mat0 = sd_term_out.mat_unrotated;
        write_octave_like_output(sprintf(out_format,level,degree,t,d), full(mat));
        write_octave_like_output(sprintf(out_format0,level,degree,t,d), full(mat0));
    end
end
