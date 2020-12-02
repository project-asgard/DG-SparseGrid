% continuity1 sg l4d3
data_dir = "generated-inputs/adapt/";
out_format = strcat(data_dir, 'continuity1_l4_d3_');
run_adapt(@continuity1,@fval_gen_1, out_format, ...
    'lev',4,'deg',3,'grid_type','SG');
out_format = strcat(data_dir, 'continuity2_l5_d2_');
run_adapt(@continuity1,@fval_gen_zeroes, out_format, ...
    'lev',5,'deg',2,'grid_type','SG');


function [x] = fval_gen_1(num_entries, elem_size)
    x = zeros(num_entries, 1);
    num_elems = num_entries/elem_size;
    for i=1:num_elems-1
        entry = i*elem_size;
        if mod(i, 2) == 0
            x(entry:entry+elem_size-1) = 1e-7;
        end
        if mod(i, 3) == 0
            x(entry:entry+elem_size-1) = (1/elem_size)*2e-5;
        end
    end
    x(num_entries-elem_size+1:num_entries) = 1.0;
end

function [x] = fval_gen_zeroes(num_entries, elem_size)
    x = zeros(num_entries, 1);
    x(num_entries-elem_size+1:num_entries) = 1.0;
end

function run_adapt(pde_handle, fval_gen, out_format, varargin)

  opts = OPTS(varargin);  
  opts.quiet = 1;
  opts.adapt_threshold = 1e-6;
  opts.max_lev_coeffs = 0;
  
  pde = pde_handle( opts );
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
