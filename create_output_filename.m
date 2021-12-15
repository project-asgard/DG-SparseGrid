function [filename] = create_output_filename(opts)

grid_type_str = opts.grid_type;
if opts.adapt
    adapt_threshold_str = num2str(opts.adapt_threshold,'%1.1e');
    grid_type_str = 'ASG';
end

filename = [ ...
    '-',grid_type_str, ...   
    '-l',replace(num2str(opts.lev),' ','') ...
    '-d',num2str(opts.deg), ...
    '-dt',num2str(opts.dt,'%1.1e')];

if opts.adapt
    filename = [filename,'-at',adapt_threshold_str];
end

if ~isempty(opts.output_filename_id)
    filename = ['-',opts.output_filename_id,filename];
end

root_folder = get_root_folder();

filename = append(root_folder,"/output/asgard-out",filename,".mat");

end