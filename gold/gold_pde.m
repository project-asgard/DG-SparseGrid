% pde testing

data_dir = strcat("generated-inputs", "/", "pde", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

out_format = strcat( data_dir, "diffusion_1_");
run_pde(diffusion1, out_format, 4.2, 0, 'lev', 3, 'deg', 2, 'CFL', 1);

out_format = strcat( data_dir, "advection_1_");
run_pde(advection1, out_format, 4.2, 0, 'lev', 3, 'deg', 2, 'CFL', 1);

out_format = strcat( data_dir, "fokkerplanck2_complete_" )
domain = [ 0.1, 0.2, 0.3, 0.4, 0.5 ];
t = 0;
%Captain
run_pde( fokkerplanck2_complete, out_format, domain, t, 'lev', 5, 'deg', 4, 'CFL', 1 )

% continuity 1
out_format = strcat(data_dir, "continuity_1_");
pde = continuity1;
pde = check_pde(pde);
pde.CFL = 1.0;
dt = pde.set_dt(pde,pde.CFL);
x = 1.1;
y_init = pde.dimensions{1}.init_cond_fn(x);
y_source0_x = pde.sources{1}{1}(x);
y_source0_t = pde.sources{1}{2}(x);
y_source1_x = pde.sources{2}{1}(x);
y_source1_t = pde.sources{2}{2}(x);
y_exact_x = pde.analytic_solutions_1D{1}(x);
y_exact_time = pde.analytic_solutions_1D{2}(x);

write_octave_like_output(strcat(out_format, 'dt.dat'), dt);
write_octave_like_output(strcat(out_format, 'initial_dim0.dat'), y_init);
write_octave_like_output(strcat(out_format, 'source0_dim0.dat'), y_source0_x);
write_octave_like_output(strcat(out_format, 'source0_time.dat'), y_source0_t);
write_octave_like_output(strcat(out_format, 'source1_dim0.dat'), y_source1_x);
write_octave_like_output(strcat(out_format, 'source1_time.dat'), y_source1_t);
write_octave_like_output(strcat(out_format, 'exact_dim0.dat'), y_exact_x);
write_octave_like_output(strcat(out_format, 'exact_time.dat'), y_exact_time);

% continuity 2
out_format = strcat(data_dir, "continuity_2_");
pde = continuity2;
pde = check_pde(pde);
x = 2.2;
for d=1:length(pde.dimensions)
  y_init = pde.dimensions{d}.init_cond_fn(x);
  write_octave_like_output(strcat(out_format, sprintf('initial_dim%d.dat', d-1)), y_init);
  y_exact = pde.analytic_solutions_1D{d}(x);
  write_octave_like_output(strcat(out_format, sprintf('exact_dim%d.dat', d-1)), y_exact);
end
y_exact_time = pde.analytic_solutions_1D{length(pde.analytic_solutions_1D)}(x);
write_octave_like_output(strcat(out_format, 'exact_time.dat'), y_exact_time);

for s=1:length(pde.sources)
  for d=1:length(pde.dimensions)
    y_source = pde.sources{s}{d}(x);
    write_octave_like_output(strcat(out_format, sprintf('source%d_dim%d.dat',s-1,d-1)), y_source);
  y_source_t = pde.sources{s}{length(pde.sources{s})}(x);
  write_octave_like_output(strcat(out_format, sprintf('source%d_time.dat',s-1)), y_source_t);
  end
end
pde.CFL=1;
dt = pde.set_dt(pde,pde.CFL);
write_octave_like_output(strcat(out_format, 'dt.dat'), dt);

% continuity 3
out_format = strcat(data_dir, "continuity_3_");
pde = continuity3;
pde = check_pde(pde);
pde.dimensions{1}.lev = 2;
pde.dimensions{2}.lev = 2;
pde.dimensions{3}.lev = 2;
x = 3.3;
for d=1:length(pde.dimensions)
  y_init = pde.dimensions{d}.init_cond_fn(x);
  write_octave_like_output(strcat(out_format, sprintf('initial_dim%d.dat', d-1)), y_init);
  y_exact = pde.analytic_solutions_1D{d}(x);
  write_octave_like_output(strcat(out_format, sprintf('exact_dim%d.dat', d-1)), y_exact);
end
y_exact_time = pde.analytic_solutions_1D{length(pde.analytic_solutions_1D)}(x);
write_octave_like_output(strcat(out_format, 'exact_time.dat'), y_exact_time);

for s=1:length(pde.sources)
  for d=1:length(pde.dimensions)
    y_source = pde.sources{s}{d}(x);
    write_octave_like_output(strcat(out_format, sprintf('source%d_dim%d.dat',s-1,d-1)), y_source);
  y_source_t = pde.sources{s}{length(pde.sources{s})}(x);
  write_octave_like_output(strcat(out_format, sprintf('source%d_time.dat',s-1)), y_source_t);
  end
end

pde.CFL=1;
dt = pde.set_dt(pde,pde.CFL);
write_octave_like_output(strcat(out_format, 'dt.dat'), dt);

% continuity 6
out_format = strcat(data_dir, "continuity_6_");
pde = check_pde(continuity6);
x = 6.6;
for d=1:length(pde.dimensions)
  y_init = pde.dimensions{d}.init_cond_fn(x);
  write_octave_like_output(strcat(out_format, sprintf('initial_dim%d.dat', d-1)), y_init);
  y_exact = pde.analytic_solutions_1D{d}(x);
  write_octave_like_output(strcat(out_format, sprintf('exact_dim%d.dat', d-1)), y_exact);
end
y_exact_time = pde.analytic_solutions_1D{length(pde.analytic_solutions_1D)}(x);
write_octave_like_output(strcat(out_format, 'exact_time.dat'), y_exact_time);

for s=1:length(pde.sources)
  for d=1:length(pde.dimensions)
    y_source = pde.sources{s}{d}(x);
    write_octave_like_output(strcat(out_format, sprintf('source%d_dim%d.dat',s-1,d-1)), y_source);
  y_source_t = pde.sources{s}{length(pde.sources{s})}(x);
  write_octave_like_output(strcat(out_format, sprintf('source%d_time.dat',s-1)), y_source_t);
  end
end
pde.CFL=1;
dt = pde.set_dt(pde,pde.CFL);
write_octave_like_output(strcat(out_format, 'dt.dat'), dt);

% fokkerplanck2_complete
out_format = strcat(data_dir, "fokkerplanck2_complete_");
pde = check_pde(fokkerplanck2_complete);
x = 0.5;
for d=1:length(pde.dimensions)
  y_init = pde.dimensions{d}.init_cond_fn(x);
  write_octave_like_output(strcat(out_format, sprintf('initial_dim%d.dat', d-1)), y_init);
end

pde.CFL=1;
dt = pde.set_dt(pde,pde.CFL);
write_octave_like_output(strcat(out_format, 'dt.dat'), dt);

function run_pde(pde, out_format, x, t, varargin)

  runtime_defaults;
  pde = check_pde( pde, opts );
  pde.CFL=1;
  dt = pde.set_dt(pde,pde.CFL);
  write_octave_like_output(strcat(out_format, 'dt.dat'), dt);
  for d=1:length(pde.dimensions)
    y_init = pde.dimensions{d}.init_cond_fn(x, 0, t);
    write_octave_like_output(strcat(out_format, sprintf('initial_dim%d.dat', d-1)), y_init);
    y_exact = pde.analytic_solutions_1D{d}(x, 0, 0);
    write_octave_like_output(strcat(out_format, sprintf('exact_dim%d.dat', d-1)), y_exact);
  end
  y_exact_time = pde.analytic_solutions_1D{length(pde.analytic_solutions_1D)}(x);
  write_octave_like_output(strcat(out_format, 'exact_time.dat'), y_exact_time);
  for s=1:length(pde.sources)
    for d=1:length(pde.dimensions)
      y_source = pde.sources{s}{d}(x);
      write_octave_like_output(strcat(out_format, sprintf('source%d_dim%d.dat',s-1,d-1)), ...
                               y_source);
      y_source_t = pde.sources{s}{length(pde.sources{s})}(x);
      write_octave_like_output(strcat(out_format, sprintf('source%d_time.dat',s-1)), y_source_t);
    end
  end
end
