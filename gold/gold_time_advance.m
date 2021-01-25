
%% Generate gold data for C++ testing of time advance component

% time advance
data_dir = strcat("generated-inputs", "/", "time_advance", "/");
root = get_root_folder();
[stat,msg] = mkdir ([root,'/gold/',char(data_dir)]);

%% adaptive time advance testing

out_format = strcat(data_dir,'fokkerplanck1_4p1a_ad_sg_l4_d4_t%d.dat');
run_time_advance(@fokkerplanck1_pitch_E,out_format, ...
    'CFL',0.01,'lev',4,'deg',4,'grid_type','SG','timestep_method','RK3', 'adapt', true, ...
    'adapt_initial_condition', true, 'adapt_threshold', 1e-4);

out_format = strcat(data_dir,'continuity2_ad_sg_l3_d4_t%d.dat');
run_time_advance(@continuity2,out_format, ...
    'CFL',0.01,'lev',3,'deg',4,'grid_type','SG','timestep_method','RK3', 'adapt', true, ...
    'adapt_initial_condition', true, 'adapt_threshold', 1e-3);

out_format = strcat(data_dir,'diffusion2_ad_sg_l3_d4_t%d.dat');
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',3,'deg',4,'grid_type','SG','timestep_method','RK3', 'adapt', true, ...
    'adapt_initial_condition', true, 'adapt_threshold', 0.5e-1);

out_format = strcat(data_dir,'diffusion2_ad_implicit_sg_l3_d4_t%d.dat');
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',3,'deg',4,'grid_type','SG','timestep_method','BE', 'adapt', true, ...
    'adapt_initial_condition', true, 'adapt_threshold', 0.5e-1);

%% non-adaptive

%% uniform level, explicit

% diffusion2 l2d2 RK3
out_format = strcat(data_dir,'diffusion2_sg_l2_d2_t%d.dat');
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');
% diffusion2 l3d3 RK3
out_format = strcat(data_dir,'diffusion2_sg_l3_d3_t%d.dat');
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',3,'deg',3,'grid_type','SG','timestep_method','RK3');
% diffusion2 l4d4 RK3
out_format = strcat(data_dir,'diffusion2_sg_l4_d4_t%d.dat');
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',4,'deg',4,'grid_type','SG','timestep_method','RK3');

% diffusion1 l2d2 RK3
out_format = strcat(data_dir,'diffusion1_sg_l2_d2_t%d.dat');
run_time_advance(@diffusion1,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');
% diffusion1 l3d3 RK3
out_format = strcat(data_dir,'diffusion1_sg_l3_d3_t%d.dat');
run_time_advance(@diffusion1,out_format, ...
    'CFL',0.01,'lev',3,'deg',3,'grid_type','SG','timestep_method','RK3');
% diffusion1 l4d4 RK3
out_format = strcat(data_dir,'diffusion1_sg_l4_d4_t%d.dat');
run_time_advance(@diffusion1,out_format, ...
    'CFL',0.01,'lev',4,'deg',4,'grid_type','SG','timestep_method','RK3');

% continuity1 sg l2d2
out_format = strcat(data_dir, 'continuity1_sg_l2_d2_t%d.dat');
run_time_advance(@continuity1,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');
% continuity1 fg l2d2
out_format = strcat(data_dir, 'continuity1_fg_l2_d2_t%d.dat');
run_time_advance(@continuity1,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','FG','timestep_method','RK3');
% continuity1 sg l4d3
out_format = strcat(data_dir, 'continuity1_sg_l4_d3_t%d.dat');
run_time_advance(@continuity1,out_format, ...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','timestep_method','RK3');

% continuity2 sg l2d2
out_format = strcat(data_dir, 'continuity2_sg_l2_d2_t%d.dat');
run_time_advance(@continuity2,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');
% continuity2 fg l2d2
out_format = strcat(data_dir, 'continuity2_fg_l2_d2_t%d.dat');
run_time_advance(@continuity2,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','FG','timestep_method','RK3');
% continuity2 sg l4d3
out_format = strcat(data_dir, 'continuity2_sg_l4_d3_t%d.dat');
run_time_advance(@continuity2,out_format, ...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','timestep_method','RK3');

% continuity3 sg l2d2
out_format = strcat(data_dir,'continuity3_sg_l2_d2_t%d.dat');
run_time_advance(@continuity3,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');
% continuity3 sg l4d3
out_format = strcat(data_dir,'continuity3_sg_l4_d3_t%d.dat');
run_time_advance(@continuity3,out_format, ...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','timestep_method','RK3');

% continuity6 sg l2d2
out_format = strcat(data_dir, 'continuity6_sg_l2_d2_t%d.dat');
run_time_advance(@continuity6,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');
% continuity6 sg l2d3
out_format = strcat(data_dir, 'continuity6_sg_l2_d3_t%d.dat');
run_time_advance(@continuity6,out_format, ...
    'CFL',0.01,'lev',2,'deg',3,'grid_type','SG','timestep_method','RK3');

% fokkerplanck1_4p2 sg l2d2
out_format = strcat(data_dir,'fokkerplanck1_4p2_sg_l2_d2_t%d.dat');
run_time_advance(@fokkerplanck1_pitch_C,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');

% fokkerplanck1_4p2 fg l2d2
out_format = strcat(data_dir,'fokkerplanck1_4p2_fg_l2_d2_t%d.dat');
run_time_advance(@fokkerplanck1_pitch_C,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','FG','timestep_method','RK3');

% fokkerplanck1_4p1a sg l2d2
out_format = strcat(data_dir,'fokkerplanck1_4p1a_sg_l2_d2_t%d.dat');
run_time_advance(@fokkerplanck1_pitch_E,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');

%fokkerplanck1_4p3 sg l2d2
out_format = strcat(data_dir,'fokkerplanck1_4p3_sg_l2_d2_t%d.dat');
run_time_advance(@fokkerplanck1_pitch_R,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','RK3');


%% uniform level, implicit
% note fp2 must be calculated with modified x domain [1, 20] to avoid
% singularity problem
% fokkerplanck2_complete sg l3d3
out_format = strcat(data_dir,'fokkerplanck2_complete_implicit_sg_l3_d3_t%d.dat');
% any explicit timestep method will give diverging results
run_time_advance(@fokkerplanck2_complete,out_format, ...
    'CFL',0.01,'lev',3,'deg',3,'grid_type','SG','case', 4,...
    'timestep_method', 'BE');

out_format = strcat(data_dir,'fokkerplanck2_complete_implicit_sg_l4_d3_t%d.dat');
run_time_advance(@fokkerplanck2_complete,out_format, ...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','case', 4,...
    'timestep_method', 'BE');

out_format = strcat(data_dir,'fokkerplanck2_complete_implicit_sg_l5_d3_t%d.dat');
run_time_advance(@fokkerplanck2_complete,out_format, ...
    'CFL',0.01,'lev',5,'deg',3,'grid_type','SG','case', 4,...
    'timestep_method', 'BE');

% continuity1 sg l2d2
out_format = strcat(data_dir,'continuity1_implicit_l2_d2_t%d.dat');
run_time_advance(@continuity1,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','BE');
% continuity1 sg l4d3
out_format = strcat(data_dir,'continuity1_implicit_l4_d3_t%d.dat');
run_time_advance(@continuity1,out_format, ...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','timestep_method','BE');

% continuity2 sg l2d2
out_format = strcat(data_dir,'continuity2_implicit_l2_d2_t%d.dat');
run_time_advance(@continuity2,out_format, ...
    'CFL',0.01,'lev',2,'deg',2,'grid_type','SG','timestep_method','BE');
% continuity2 sg l4d3
out_format = strcat(data_dir,'continuity2_implicit_l4_d3_t%d.dat');
run_time_advance(@continuity2,out_format, ...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG','timestep_method','BE');

% diffusion1 l4d4 BE
out_format = strcat(data_dir,'diffusion1_implicit_sg_l4_d4_t%d.dat');
run_time_advance(@diffusion1,out_format, ...
    'CFL',0.01,'lev',4,'deg',4,'grid_type','SG','timestep_method','BE');

% diffusion 2 implicit
out_format = strcat(data_dir,'diffusion2_implicit_sg_l3_d3_t%d.dat');
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',3,'deg',3,'grid_type','SG',...
    'timestep_method', 'BE');

out_format = strcat(data_dir,'diffusion2_implicit_sg_l4_d3_t%d.dat');
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',4,'deg',3,'grid_type','SG',...
    'timestep_method', 'BE');

out_format = strcat(data_dir,'diffusion2_implicit_sg_l5_d3_t%d.dat');
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',5,'deg',3,'grid_type','SG',...
    'timestep_method', 'BE');

% advection1 l4d4 implicit BE
out_format = strcat(data_dir,'advection1_implicit_sg_l4_d4_t%d.dat');
run_time_advance(@advection1,out_format, ...
    'CFL',0.01,'lev',4,'deg',4,'grid_type','SG','timestep_method','BE');


%% non-uniform level, explicit

% diffusion2 d2
out_format = strcat(data_dir,'diffusion2');
levels = [4, 5];
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',levels,'deg',2,'grid_type','SG','timestep_method','RK3');

% continuity2 fg d3
out_format = strcat(data_dir, 'continuity2');
levels = [3, 4];
run_time_advance(@continuity2,out_format, ...
    'CFL',0.01,'lev',levels,'deg',3,'grid_type','FG','timestep_method','RK3');


% continuity3 sg d4
out_format = strcat(data_dir,'continuity3');
levels = [3, 4, 2];
run_time_advance(@continuity3,out_format, ...
    'CFL',0.01,'lev',levels,'deg',4,'grid_type','SG','timestep_method','RK3');

% continuity6 sg d2
out_format = strcat(data_dir, 'continuity6');
levels = [2, 3, 2, 3, 3, 2];
run_time_advance(@continuity6,out_format, ...
    'CFL',0.01,'lev',levels,'deg',2,'grid_type','SG','timestep_method','RK3');

%% non-uniform level, implicit

% diffusion2 sg d2
out_format = strcat(data_dir,'diffusion2_implicit');
levels = [4, 5];
run_time_advance(@diffusion2,out_format, ...
    'CFL',0.01,'lev',levels,'deg',2,'grid_type','SG','timestep_method','BE');

% continuity2 fg d3
out_format = strcat(data_dir, 'continuity2_implicit');
levels = [3, 4];
run_time_advance(@continuity2,out_format, ...
    'CFL',0.01,'lev',levels,'deg',3,'grid_type','FG','timestep_method','BE');

% fokkerplanck2_complete sg d3
% note fp2 must be calculated with modified x domain [1, 20] to avoid
% singularity problem
out_format = strcat(data_dir,'fokkerplanck2_complete_implicit');
% any explicit timestep method will give diverging results
levels = [2, 3];
run_time_advance(@fokkerplanck2_complete,out_format, ...
    'CFL',0.01,'lev',levels,'deg',3,'grid_type','SG','case', 4,...
    'timestep_method', 'BE');

function run_time_advance(pde_handle, out_format, varargin)

  opts = OPTS(varargin);  
  opts.quiet = 1;
  pde = pde_handle( opts );
  dt = pde.set_dt( pde, opts.CFL );
  if length(pde.dimensions) > 3
      opts.max_lev = 4;
  end
  opts.max_lev_coeffs = 1;
  lev_vec = get_lev_vec(pde);
  
  % non-uniform case
  if ~(length(unique(lev_vec)) == 1)
      lev_string = "";
      for i=1:length(pde.dimensions)
         lev_string = lev_string + int2str(lev_vec(i)) + "_";
      end
      out_format = strcat(out_format, '_', lower(opts.grid_type),...
          '_l', lev_string, 'd', int2str(opts.deg), '_t%d.dat');
  end
  
  [elements, elements_idx]    = hash_table_sparse_nD (lev_vec, ...
                                          opts.max_lev, opts.grid_type);
  hash_table.elements         = elements;
  hash_table.elements_idx     = elements_idx;

  [~, pde.transform_blocks] = OperatorTwoScale_wavelet2(opts.deg,opts.max_lev);

  t = 0;
  TD = 0;
  pde = get_coeff_mats(pde,opts,t,TD,opts.use_oldcoeffmat);

  A_data = global_matrix(pde,opts,hash_table);

  fval = initial_condition_vector(pde,opts,hash_table,0);

  
  %% Check to see if initial resolution meets requested accuracy
if opts.adapt
    if opts.adapt_initial_condition        
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
        [pde,~,hash_table,A_data,~,~,~,~,~,~] ...
            = adapt(pde,opts,[],fval,hash_table,[],[],[],[], ...
            [],1,0);
        % reproject onto coarsend basis
        fval = initial_condition_vector(pde, opts, hash_table, t);
    end
end
  
  
  for i=0:4
    t = i*dt;
    % Coarsen Grid
    if opts.adapt
        [pde,fval,hash_table,A_data,~,~,~,~,~,~] ...
            = adapt(pde,opts,[],fval,hash_table,[],[],[],[],[],1,0);
    end
    
    needs_adapting = true;
    while needs_adapting

        
        % construct the TD coeff_mats.
        TD = 1;
        pde = get_coeff_mats(pde,opts,t,TD);

        
        fval_unstepped = fval;
        fval = time_advance(pde,opts,A_data,fval,t,dt,opts.deg,hash_table,[],[]);
        
        %%
        % Refine Grid - determine which elements to add, but reset fval to
        % fval_previous with those new elements added and time_advance
        % again
        if opts.adapt
            
            num_elements_0 = numel(fval);
            
            [pde,~,hash_table,A_data,~,~,~,~,~,~,fval_unstepped_adapted] ...
                = adapt(pde,opts,[],fval,hash_table,[],[],[],[], ...
                [],0,1,fval_unstepped);
            
            num_elements_adapted = numel(fval_unstepped_adapted);
            
            if num_elements_0 == num_elements_adapted
                needs_adapting = false;
                fval = time_advance(pde,opts,A_data,fval_unstepped_adapted,t,dt,opts.deg,hash_table,[],[]);
            else
                fval = fval_unstepped_adapted;
            end
        else
            needs_adapting = false;
        end
        
    end
    write_octave_like_output(sprintf(out_format,i), fval);
  end
end
