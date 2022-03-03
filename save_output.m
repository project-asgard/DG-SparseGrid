function [outputs] = save_output(outputs,L,pde,opts,num_dims,fval,fval_realspace,f_realspace_analytic_nD,nodes,nodes_nodups,nodes_count,t,dt,t1,root_directory,hash_table)

outputs.fval_t{L+1} = fval;
outputs.nodes_t{L+1} = nodes_nodups;
outputs.time_array(L+1) = t+dt;
outputs.wall_clock_time(L+1) = t1;
outputs.dt = dt;
outputs.hash_table = hash_table;

f_realspace_nD = singleD_to_multiD(num_dims,fval_realspace,nodes);
if strcmp(opts.output_grid,'fixed') || strcmp(opts.output_grid,'elements')
    f_realspace_nD = ...
        remove_duplicates(num_dims,f_realspace_nD,nodes_nodups,nodes_count);
end

if num_dims <= 3
    outputs.f_realspace_nD_t{L+1} = f_realspace_nD;
    outputs.f_realspace_analytic_nD_t{L+1} = f_realspace_analytic_nD;
else
    error('Save output for num_dimensions >3 not yet implemented');
end

store_alpha_only = false;
if store_alpha_only
    outputs_all = outputs;
    clear outputs;
    outputs.time_array = outputs_all.time_array;
    if isfield(outputs_all,'alpha_t')
        outputs.alpha_t = outputs_all.alpha_t;
    end
    outputs.fval_t = {};
    outputs.fval_t{1} = outputs_all.fval_t{end};
    clear pde;
    pde = 1;
end

if opts.save_output && (mod(L,opts.save_freq)==0 || L==opts.num_steps)
    [status, msg, msgID] = mkdir([root_directory,'/output']);
    
    outputs.output_file_name = create_output_filename(opts);
    disp(['    Saving output ',char(outputs.output_file_name),' ...']);
    save(outputs.output_file_name,'pde','opts','outputs','-nocompression');
    disp(['    DONE']);
end

end