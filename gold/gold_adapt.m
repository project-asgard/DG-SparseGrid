% adaptivity functions
data_dir = "generated-inputs/adapt/";
out_format = strcat(data_dir, 'continuity1_l4_d3_');
run_adapt(@continuity1,@fval_gen_1, out_format, ...
    'lev',4,'deg',3,'grid_type','SG');
out_format = strcat(data_dir, 'continuity2_l5_d2_');
run_adapt(@continuity2,@fval_gen_zeroes, out_format, ...
    'lev',5,'deg',2,'grid_type','SG');
out_format = strcat(data_dir, 'continuity3_l4_d4_');
run_adapt(@continuity3,@fval_gen_2, out_format, ...
    'lev',4,'deg',4,'grid_type','SG');


% initial condition adapt
out_string = strcat(data_dir, 'diffusion2_l2_d3_initial.dat');
run_initial(@diffusion2,out_string, ...
     'lev',[2,2],'deg',3,'grid_type','SG');
out_string = strcat(data_dir, 'diffusion1_l3_d4_initial.dat');
run_initial(@diffusion1,out_string, ...
     'lev',3,'deg',4,'grid_type','SG');  
 
function [x] = fval_gen_1(num_entries, elem_size)
    x = zeros(num_entries, 1);
    num_elems = num_entries/elem_size;
    for i=1:num_elems-1
        entry = (i*elem_size)+1;
        if mod(i, 2) == 0
            x(entry) = 0.5e-5; % should be coarsened
        end
        if mod(i, 3) == 0
            x(entry) = 1e-3; % should be refined
        end
    end
    x(num_entries-elem_size+1:num_entries) = 1.0;
end

function [x] = fval_gen_2(num_entries, elem_size)
    x = zeros(num_entries, 1);
    num_elems = num_entries/elem_size;
    for i=1:num_elems-1
        entry = (i*elem_size)+1;
        if mod(i, 3) == 0
            x(entry) = 1e-5; % should be coarsened
            x(entry + elem_size) = 1e-5;
        end
        if mod(i, 5) == 0
            x(entry) = 1; % should be refined
            x(entry + elem_size) = 0.5;
        end
    end
    x = x(1:num_entries);
    x(num_entries-elem_size+1:num_entries) = 1.0;
end

function [x] = fval_gen_zeroes(num_entries, elem_size)
    x = zeros(num_entries, 1);
    x(num_entries-elem_size+1:num_entries) = 1.0;
end

function run_adapt(pde_handle, fval_gen, out_format, varargin)

  opts = OPTS(varargin);  
  opts.quiet = 1;
  opts.adapt_threshold = 1e-4;
  opts.max_lev_coeffs = 0;
  
  pde = pde_handle( opts );
  pde = compute_dimension_mass_mat(opts, pde);
  pde = get_coeff_mats(pde,opts,0,0,0);
  lev_vec = get_lev_vec(pde);
  num_dims = length(lev_vec);

  [elements, elements_idx]    = hash_table_sparse_nD (lev_vec, ...
                                          opts.max_lev, opts.grid_type);
  hash_table.elements         = elements;
  hash_table.elements_idx     = elements_idx;
  
  hash_table_orig = hash_table;
  
  elem_size = opts.deg^num_dims;
  fval_size = numel(hash_table.elements_idx)*elem_size;
  fval_orig = fval_gen(fval_size, elem_size);
  fval = fval_orig;
  write_octave_like_output(strcat(out_format,'orig.dat'), fval_orig);  
  
  % coarsen
  coarsen = 1;
  refine = 0;
  [~,fval,hash_table,~,~,~,~,~,~,~,~,~] ...
    = adapt(pde,opts,[],fval,hash_table_orig,[],[],[],[],[],coarsen,refine);
  write_octave_like_output(strcat(out_format,'coarse.dat'), fval);  
  [~, ~, cell_val] = find(hash_table.elements.pos_p1);
  num_entries = size(cell_val,1) / num_dims; 
  table = zeros(num_entries, num_dims * 2);
  for j=1:num_entries
     id = hash_table.elements_idx(j);
     table(j, :) = horzcat(full(hash_table.elements.lev_p1(id,:)), ...
         full(hash_table.elements.pos_p1(id,:))) - 1;
  end
  write_octave_like_output(strcat(out_format, 'coarse_table.dat'),table);
        
  
  % refine
  coarsen = 0;
  refine = 1;
  [~,fval,hash_table,~,~,~,~,~,~,~,~,~] ...
    = adapt(pde,opts,[],fval_orig,hash_table_orig,[],[],[],[],[],coarsen,refine);
  write_octave_like_output(strcat(out_format, 'refine.dat'), fval);  
  [~, ~, cell_val] = find(hash_table.elements.pos_p1);
  num_entries = size(cell_val,1) / num_dims; 
  table = zeros(num_entries, num_dims * 2);
  for j=1:num_entries
     id = hash_table.elements_idx(j);
     table(j, :) = horzcat(full(hash_table.elements.lev_p1(id,:)), ...
         full(hash_table.elements.pos_p1(id,:))) - 1;
  end
  write_octave_like_output(strcat(out_format,'refine_table.dat'),table);
end


function run_initial(pde_handle, out_string, varargin)
  opts = OPTS(varargin);  
  opts.quiet = 1;
  opts.adapt_threshold = 1e-4;
  opts.max_lev_coeffs = 1;

  pde = pde_handle( opts );
  t = 0;
  TD = 0;
  pde = compute_dimension_mass_mat(opts, pde);
  pde = get_coeff_mats(pde,opts,t,TD,opts.use_oldcoeffmat);

  lev_vec = get_lev_vec(pde);

  [elements, elements_idx]    = hash_table_sparse_nD (lev_vec, ...
                                          opts.max_lev, opts.grid_type);
  hash_table.elements         = elements;
  hash_table.elements_idx     = elements_idx;
  
  fval = initial_condition_vector(pde, opts, hash_table, t);

  keep_adapting_initial_condition = true;
  while keep_adapting_initial_condition
      num_pre_adapt = numel(fval);
      % first refine
      [pde,fval_tmp,hash_table,~,~,~,~,~,~,~,~,~] ...
                = adapt(pde,opts,[],fval,hash_table,[],[],[],[], ...
                [],0,1);
      if num_pre_adapt == numel(fval_tmp)
          keep_adapting_initial_condition = false;
      end
      clear fval_tmp;
      % reproject initial condition onto refined basis
      fval = initial_condition_vector(pde, opts, hash_table, t);
  end
  % coarsen
  [pde,~,hash_table,~,~,~,~,~,~,~,~,~] ...
    = adapt(pde,opts,[],fval,hash_table,[],[],[],[], ...
      [],1,0);
  % reproject onto coarsend basis
  fval = initial_condition_vector(pde, opts, hash_table, t);
  write_octave_like_output(out_string,fval);
end
