% coefficient testing

%diffusion1
generate_data( diffusion1, 'diffusion1', 'lev', 2, 'deg', 2 );
generate_data( diffusion1, 'diffusion1', 'lev', 4, 'deg', 4 );
generate_data( diffusion1, 'diffusion1', 'lev', 5, 'deg', 5 );

%diffusion2
generate_data( diffusion2, 'diffusion2', 'lev', 2, 'deg', 2 );
generate_data( diffusion2, 'diffusion2', 'lev', 4, 'deg', 4 );
generate_data( diffusion2, 'diffusion2', 'lev', 5, 'deg', 5 );

%continuity1
generate_data( continuity1, 'continuity1', 'lev', 2, 'deg', 2 );

coeff_dir = strcat("generated-inputs", "/", "coefficients", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(coeff_dir)]);

%fokkerplanck2_complete
generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete', 'lev', 2, 'deg', 2 );
generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete', 'lev', 4, 'deg', 4 );
generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete', 'lev', 5, 'deg', 5 );

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

function generate_data(pde, output_prefix, varargin)

runtime_defaults

opts.use_oldcoeffmat = 0;
opts.use_oldhash = 1;
opts.compression = 4;
opts.useConnectivity = 0;

pde = check_pde(pde, opts);

CFL = 0.01;
dt = pde.set_dt(pde,CFL);

% coefficient testing
coeff_dir = strcat("generated-inputs/coefficients/", output_prefix, "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(coeff_dir)]);

% diffusion1 terms
for i=1:length(pde.dimensions)
  pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, 'coefficients_l%i_d%i_%d_%d.dat');
%doesn't matter, the term is time independent...
time = 1.0;
level = pde.dimensions{1}.lev;
degree = pde.deg;
for t=1:length(pde.terms)
    for d=1:length(pde.dimensions)
        sd_term = pde.terms{t}.terms_1D{d};
        sd_term_out = ...
        coeff_matrix(numel(pde.dimensions),pde.deg,time,pde.dimensions{d},sd_term,pde.params);
        coeff_mat = sd_term_out.mat;
        write_octave_like_output(sprintf(out_format,level,degree,t,d), full(coeff_mat));
    end
end
end
