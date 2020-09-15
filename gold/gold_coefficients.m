generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete', 'lev', 3, 'deg', 3 );
generate_data( fokkerplanck2_complete, 'fokkerplanck2_complete', 'lev', 4, 'deg', 4 );
generate_data( continuity1, 'continuity1', 'lev', 2, 'deg', 2 );
generate_data( continuity2, 'continuity2', 'lev', 4, 'deg', 3 );
generate_data( continuity3, 'continuity3', 'lev', 4, 'deg', 4 );
generate_data( continuity6, 'continuity6', 'lev', 2, 'deg', 4 );
generate_data( diffusion1, 'diffusion1', 'lev', 5, 'deg', 6 );
generate_data( diffusion2, 'diffusion2', 'lev', 3, 'deg', 5 );
generate_data( fokkerplanck1_4p1a, 'fokkerplanck1_4p1a', 'lev', 4, 'deg', 3 );
generate_data( fokkerplanck1_4p2, 'fokkerplanck1_4p2', 'lev', 5, 'deg', 2 );
generate_data( fokkerplanck1_4p3, 'fokkerplanck1_4p3', 'lev', 2, 'deg', 5 );
generate_data( fokkerplanck1_4p4, 'fokkerplanck1_4p4', 'lev', 5, 'deg', 3 );
generate_data( fokkerplanck1_4p5, 'fokkerplanck1_4p5', 'lev', 3, 'deg', 5 );

function generate_data(pde, pde_name, varargin)
% coefficient testing
coeff_dir = "generated-inputs/coefficients/";
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(coeff_dir)]);

runtime_defaults

opts.use_oldcoeffmat = 0;
opts.use_oldhash = 0;

pde = check_pde(pde, opts);

CFL = 0.01;
dt = pde.set_dt(pde,CFL);

for i=1:length(pde.dimensions)
  pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev,'wavelet');
end

out_format = strcat(coeff_dir, pde_name, '_coefficients_l%i_d%i_%d_%d.dat');
out_format_unrot = strcat(coeff_dir, pde_name, '_coefficients_norotate_l%i_d%i_%d_%d.dat');

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
        unrotated = sd_term_out.mat_unrotated;
        write_octave_like_output(sprintf(out_format,level,degree,t,d), full(coeff_mat));
        write_octave_like_output(sprintf(out_format_unrot,level,degree,t,d), full(unrotated));
    end
end
end
