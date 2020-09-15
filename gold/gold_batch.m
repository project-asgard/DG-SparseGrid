%% Generate gold data for C++ testing of batch component

% batch
batch_dir = strcat("generated-inputs", "/", "batch", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(batch_dir)]);

% continuity1 - sg
out_format = strcat(batch_dir, 'continuity1_sg_l2_d2_t%d.dat');
run_batch(continuity1,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');

% continuity1 - sg
out_format = strcat(batch_dir, 'continuity1_sg_l3_d4_t%d.dat');
run_batch(continuity1,out_format,...
    'CFL',0.01,'lev',3,'deg',4,'grid_type','SG','timestep_method','RK3');

% continuity2 - sg
out_format = strcat(batch_dir, 'continuity2_sg_l2_d2_t%d.dat');
run_batch(continuity2,out_format,...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');

% continuity2 - fg
out_format = strcat(batch_dir, 'continuity2_fg_l3_d4_t%d.dat');
run_batch(continuity2,out_format,...
    'CFL',0.01,'lev',3,'deg',4,'grid_type','FG','timestep_method','RK3');

% continuity3 - sg
out_format = strcat(batch_dir, 'continuity3_sg_l3_d4_t%d.dat');
run_batch(continuity3,out_format,...
    'CFL',0.01,'lev',3,'deg',4,'grid_type','SG','timestep_method','RK3');

% continuity6 - sg 
out_format = strcat(batch_dir, 'continuity6_sg_l2_d3_t%d.dat');
run_batch(continuity6,out_format,...
    'CFL',0.01,'lev',2,'deg',3,'grid_type','SG','timestep_method','RK3');

function run_batch(pde,out_format,varargin)

  runtime_defaults
    
  pde = check_pde( pde, opts );
  
  [elements, elements_idx]    = hash_table_sparse_nD (pde.lev_vec, pde.max_lev, opts.grid_type);
  hash_table.elements         = elements;
  hash_table.elements_idx     = elements_idx;

  for i=1:length(pde.dimensions)
      pde.dimensions{i}.FMWT = OperatorTwoScale(pde.deg,pde.dimensions{i}.lev);
  end

  t = 0;
  TD = 0;
  pde = get_coeff_mats(pde,t,TD,0);

  A_data = global_matrix(pde,opts,hash_table);

  Vmax = 0;
  Emax = 0;

  out = initial_condition_vector(pde,opts,hash_table,0);
  out = ones(numel(out),1);
  
  out = apply_A(pde,opts,A_data,out,pde.deg,Vmax,Emax) ;
  write_octave_like_output(sprintf(out_format,1), out);
end


