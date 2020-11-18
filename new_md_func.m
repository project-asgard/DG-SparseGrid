function md_func = new_md_func(num_dims,functions)

% Return a multi-dimensional (md) function list (md_func)
%
% e.g.,
% fun_x = @(x,p,t) x.^2+1;
% fun_y = @(y,p,t) y.*0+1;
% fun_t = @(t,p) exp(-t);
%
% md_term = md_func(2,{fun_x,fun_y,fun_t}); % provide all arguments
% md_term = md_func(2,{[],[],fun_t}); % rely on defaults for space
% md_term = md_func(2,{fun_x,fun_y,[]}); rely on defaults for time
% md_term = md_func(2,{fun_x,[]}); % do not specify time at all and rely on default

% Fill a default with 1
for d=1:num_dims
    md_func{d} = @(x,p,t) x.*0+1;
end
md_func{d+1} = @(t,p) t.*0+1;

if nargin>1
    
    num_funcs = numel(functions);
    assert(num_funcs>=num_dims)
    
    % Overwrite with inputted functions if provided
    for d=1:num_funcs
        if ~isempty(functions{d})
            md_func{d} = functions{d};
        end
    end
    
    if num_funcs == num_dims+1
        if ~isempty(functions{num_dims+1})
            md_func{num_dims+1} = functions{num_dims+1};
        end
    end
    
end

end