function [root] = get_root_folder()

% Get the DG-SparseGrid ROOT level folder

dir = pwd();
top_level = 'DG-SparseGrid';
root = [];
if contains(dir,top_level)
    loc = strfind(dir,top_level);
    root = extractBefore(dir,loc+strlength(top_level));
else
    error('ERROR: must be run within the DG-SparseGrid name repo');
end

assert(~isempty(root));

end