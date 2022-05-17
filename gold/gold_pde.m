% pde testing
data_dir = strcat("generated-inputs", "/", "pde", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

t = 5;
domain = [ 0.1, 0.2, 0.3, 0.4, 0.5 ];

out_format = strcat( data_dir, "diffusion_1_");
run_pde(@diffusion1, out_format, domain, t, 'lev', 3, 'deg', 2, 'CFL', 1);

out_format = strcat( data_dir, "diffusion_2_");
run_pde(@diffusion2, out_format, domain, t, 'lev', 3, 'deg', 2, 'CFL', 1);

out_format = strcat( data_dir, "advection_1_");
run_pde(@advection1, out_format, domain, t, 'lev', 3, 'deg', 2, 'CFL', 1);

% note fp2 must be calculated with modified x domain [1, 20] to avoid
% singularity problem
out_format = strcat( data_dir, "fokkerplanck2_complete_" );
run_pde(@fokkerplanck2_complete, out_format, domain, t, 'case', 4, 'lev', 5, 'deg', 4);

out_format = strcat( data_dir, "continuity_1_" );
run_pde(@continuity1, out_format, domain, t, 'lev', 2, 'deg', 2);

out_format = strcat( data_dir, "continuity2_" );
run_pde(@continuity2, out_format, domain, t, 'lev', 5, 'deg', 4);

out_format = strcat( data_dir, "continuity_3_" );
run_pde(@continuity3, out_format, domain, t, 'lev', 5, 'deg', 4);

out_format = strcat( data_dir, "continuity_6_" );
run_pde(@continuity6, out_format, domain, t);

out_format = strcat( data_dir, "fokkerplanck2_complete_" );
gfunc_test(@fokkerplanck2_complete, out_format, domain, t, 'case', 4);

out_format = strcat( data_dir, "diffusion_2_" );
gfunc_test(@diffusion2, out_format, domain, t);

out_format = strcat( data_dir, "fokkerplanck2_complete_" );
lhsmass_test(@fokkerplanck2_complete, out_format, domain, t, 'case', 4, 'lev', [4,4], 'deg', 4);


function run_pde(pde_handle, out_format, x, t, varargin)

  opts = OPTS(varargin);  
  opts.quiet = 1;
  opts.CFL = 1;
  pde = pde_handle( opts );
  dt = pde.set_dt(pde,opts.CFL);
  
  write_octave_like_output(strcat(out_format, 'dt.dat'), dt);
  
  for d=1:length(pde.dimensions)
    init_cond = pde.initial_conditions{1}{d};
    y_init = init_cond(x, pde.params, 0);
    write_octave_like_output(strcat(out_format, sprintf('initial_dim%d.dat', d-1)), y_init);

    if ~isempty(pde.solutions)
        exact_sol = pde.solutions{1}{d};
        y_exact = exact_sol(x, pde.params, t);
        write_octave_like_output(strcat(out_format, sprintf('exact_dim%d.dat', d-1)), y_exact);
    end
  end
  
  if ~isempty(pde.solutions)
    y_exact_time = pde.solutions{1}{length(pde.dimensions)+1}(t);
    write_octave_like_output(strcat(out_format, 'exact_time.dat'), y_exact_time);
  end
  
  for s=1:length(pde.sources)
    for d=1:length(pde.dimensions)
      y_source = pde.sources{s}{d}(x);
      write_octave_like_output(strcat(out_format, sprintf('source%d_dim%d.dat',s-1,d-1)), ...
                               y_source);
      y_source_t = pde.sources{s}{length(pde.sources{s})}(t);
      write_octave_like_output(strcat(out_format, sprintf('source%d_time.dat',s-1)), y_source_t);
    end
  end
end

function gfunc_test(pde_handle, out_format, x, t, varargin)
  opts = OPTS(varargin);
  opts.quiet = 1;
  opts.CFL = 1;
  pde = pde_handle( opts );
  
  pde = get_coeff_mats(pde, opts, t, 0);
  
  gfuncs = [];
  dvfuncs = [];
  
  for d=1:length(pde.dimensions)
    for t=1:length(pde.terms)
        md_term = pde.terms{t};
        sd_term = md_term.terms_1D{d};
        for p=1:length(sd_term.pterms)
            gfunc = sd_term.pterms{p}.g(x, pde.params, t, sd_term.pterms{p}.dat);
            dvfunc = sd_term.pterms{p}.dV(x, pde.params, t, sd_term.pterms{p}.dat);
            
            gfuncs = [gfuncs; gfunc];
            dvfuncs = [dvfuncs; dvfunc];
        end
    end
  end
  
  write_octave_like_output(strcat(out_format, 'gfuncs.dat'), gfuncs);
  write_octave_like_output(strcat(out_format, 'dvfuncs.dat'), dvfuncs);
end

function lhsmass_test(pde_handle, out_format, x, t, varargin)
  opts = OPTS(varargin);
  opts.quiet = 1;
  opts.CFL = 1;
  pde = pde_handle( opts );

  pde = get_coeff_mats(pde, opts, t, 0);

  pterm_mats = [];
  lhs_mass = [];
  for d=1:length(pde.dimensions)
    for t=1:length(pde.terms)
        md_term = pde.terms{t};
        sd_term = md_term.terms_1D{d};
        for p=1:length(sd_term.pterms)
            mass = sd_term.pterms{p}.LHS_mass_mat;

            lhs_mass = [lhs_mass; mass];

            mat = sd_term.pterms{p}.mat;
            pterm_mats = [pterm_mats; mat];
        end
    end
  end

  write_octave_like_output(strcat(out_format, 'lhsmass.dat'), lhs_mass);
  write_octave_like_output(strcat(out_format, 'pterm_mats.dat'), pterm_mats);
end
