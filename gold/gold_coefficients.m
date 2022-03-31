levels = [4, 4];
generate_data( @fokkerplanck2_complete, levels, 'fokkerplanck2_complete', 'case', 4, 'lev', 4, 'deg', 4 );

% uniform level
% note fp2 must be calculated with modified x domain [1,20] to avoid singularity problem
levels = [3, 3];
generate_data( @fokkerplanck2_complete, levels, 'fokkerplanck2_complete', 'case', 4, 'deg', 3 );
levels = [4, 4];
generate_data( @fokkerplanck2_complete, levels, 'fokkerplanck2_complete', 'case', 4, 'lev', 4, 'deg', 4 );
levels = [2];
generate_data( @continuity1, levels, 'continuity1', 'deg', 2 );
levels = [4, 4];
generate_data( @continuity2, levels, 'continuity2', 'deg', 3 );
levels = [4, 4, 4];
generate_data( @continuity3, levels, 'continuity3', 'deg', 4 );
levels = [2, 2, 2, 2, 2, 2];
generate_data( @continuity6, levels, 'continuity6', 'deg', 4 );
levels = [5];
generate_data( @diffusion1, levels, 'diffusion1', 'deg', 6 );
levels = [3, 3];
generate_data( @diffusion2, levels, 'diffusion2', 'deg', 5 );
levels = [4];
generate_data( @fokkerplanck1_pitch_E, levels, 'fokkerplanck1_4p1a', 'case', 1, 'deg', 3 );
generate_data( @fokkerplanck1_pitch_E, levels, 'fokkerplanck1_pitch_E_case2', 'case', 2, 'deg', 3 );
levels = [5];
generate_data( @fokkerplanck1_pitch_C, levels, 'fokkerplanck1_4p2', 'deg', 2 );
levels = [2];
generate_data( @fokkerplanck1_pitch_R, levels, 'fokkerplanck1_4p3', 'deg', 5 );
levels = [5];
generate_data( @fokkerplanck1_pitch_CE, levels, 'fokkerplanck1_4p4', 'deg', 3 );
levels = [3];
generate_data( @fokkerplanck1_pitch_CER, levels, 'fokkerplanck1_4p5', 'deg', 5 );

% non-uniform level
% note fp2 must be calculated with modified x domain [1,20] to avoid singularity problem
levels = [2, 3];
generate_data( @fokkerplanck2_complete, levels, 'fokkerplanck2_complete', 'case', 4, 'deg', 3 );
levels = [4, 2];
generate_data( @fokkerplanck2_complete, levels, 'fokkerplanck2_complete', 'case', 4, 'lev', 4, 'deg', 4 );
levels = [4, 5];
generate_data( @continuity2, levels, 'continuity2', 'deg', 3 );
levels = [2, 3, 2];
generate_data( @continuity3, levels, 'continuity3', 'deg', 4 );
levels = [2, 3, 3, 3, 2, 4];
generate_data( @continuity6, levels, 'continuity6', 'deg', 4 );
levels = [2, 3];
generate_data( @diffusion2, levels, 'diffusion2', 'deg', 5 );

function generate_data(pde_handle, lev_vec, pde_name, varargin)
coeff_dir = "generated-inputs/coefficients/";
root = get_root_folder();
[~,~] = mkdir ([root,'/gold/',char(coeff_dir)]);

opts = OPTS(varargin);  
opts.quiet = 1;
opts.use_oldcoeffmat = 0;
opts.use_oldhash = 0;
opts.max_lev_coeffs = 0;
pde = pde_handle( opts );

for i=1:length(pde.dimensions)
  pde.dimensions{i}.lev = lev_vec(i);
end

[~, transform_blocks] = OperatorTwoScale_wavelet2(opts.deg,...
                                                  opts.max_lev);

out_format = strcat(coeff_dir, pde_name, '_coefficients_l%sd%i_%d_%d.dat');
out_format_unrot = strcat(coeff_dir, pde_name, '_coefficients_norotate_l%sd%i_%d_%d.dat');

lev_string = "";
for d=1:length(pde.dimensions)
    lev_string = lev_string + int2str(lev_vec(d)) + "_";
end

pde = compute_dimension_mass_mat(opts, pde);
pde = get_coeff_mats(pde, opts, 1.0, 0);

degree = opts.deg;
time = 1.0;
for t=1:length(pde.terms)
    for d=1:length(pde.dimensions)
        sd_term = pde.terms{t}.terms_1D{d};

        coeff_mat = sd_term.mat;
        unrotated = sd_term.mat_unrotated;

        write_octave_like_output(sprintf(out_format,lev_string,degree,t,d), full(coeff_mat));
        write_octave_like_output(sprintf(out_format_unrot,lev_string,degree,t,d), full(unrotated));
    end
end

end
