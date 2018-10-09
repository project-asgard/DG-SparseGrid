function [keystr] = getHashKeyStr(key)

hash_format =  'i%04.4d_';

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