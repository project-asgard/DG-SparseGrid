function [outputs] = save_output(outputs,L,pde,opts,num_dims,fval,fval_realspace,f_realspace_analytic_nD,nodes,nodes_nodups,nodes_count,t,dt,toc,root_directory,hash_table)

    outputs.fval_t{L+1} = fval;
    outputs.nodes_t{L+1} = nodes_nodups;
    outputs.time_array(L+1) = t+dt;
    outputs.wall_clock_time(L+1) = toc;
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
    
    if opts.save_output && (mod(L,opts.save_freq)==0 || L==opts.num_steps)
        [status, msg, msgID] = mkdir([root_directory,'/output']);

        outputs.output_file_name = create_output_filename(opts);
        save(outputs.output_file_name,'pde','opts','outputs');
        
    end

end