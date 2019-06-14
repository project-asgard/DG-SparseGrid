function [keystr] = get_hash_key_str(key)

hash_format =  'i%04.4d';

% Inputs
%
%   key : for 1D (2x1) 
%       (lev )
%       (cell)
%
%       : for 2D (2x2)
%       (lev_dim1,  lev_dim2 )
%       (cell_dim1, cell_dim2)

    keystr = sprintf(hash_format,key(:));

end