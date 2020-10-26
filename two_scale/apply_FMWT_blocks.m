function Y = apply_FMWT_blocks(level, blocks, X, transpose_side)
%
% Y = apply_FMWT( level, blocks, X, [transpose_side] )
% 
% Y = FMWT * X (using dense blocks rather than matrix storage)
%
% transpose_side =  'LN'  left multiply, no-transpose
% transpose_side =  'RN'  right multiply, no-transpose
% transpose_side =  'LT'  left multiply, transpose 
% transpose_side =  'RT'  right multiply, transpose
%

assert(level >= 0);
if(level == 0)
    Y = X;
    return
end

assert(length(blocks) > 0);
assert(mod(length(blocks),2) == 0);
max_level = length(blocks) / 2;
assert(level <= max_level, "provided dense blocks do not support level arg");

% we will check proper block dims later
degree = size(blocks{1}, 1);
n = degree * 2^level;

rows_X = size(X,1);
cols_X = size(X,2);

trans_side_default = 'LN';
TS = trans_side_default;
if (nargin >= 4),
    TS = transpose_side;
end;

if (strcmp(TS,'RN')),
        multiply_left = 0; do_transpose = 0;
elseif (strcmp(TS,'LT')),
        multiply_left = 1; do_transpose = 1;
elseif (strcmp(TS,'RT')),
        multiply_left = 0; do_transpose = 1;
elseif (strcmp(TS,'LN')),
        multiply_left = 1; do_transpose = 0;
else
   error('transpose_side string not recognized; see apply_FMWT_blocks.m');
end;

multiply_right = ~multiply_left;

% -------------------------
% determine shape of output
% -------------------------
isok = (multiply_left  && (rows_X == n)) || ...
       (multiply_right && (cols_X == n));
if (~isok),
     error('apply_FMWT: multiply_left=%d, rows_X=%d,cols_X=%d,n=%d', ...
           multiply_left, rows_X, cols_X, n);
end;

if (multiply_left),
    rows_Y = n; cols_Y = cols_X;
else
    rows_Y = rows_X; cols_Y = n;
end;
Y = zeros(rows_Y, cols_Y);

% -----------------------------------
% do math - start with coarsest block
% -----------------------------------
first_block = blocks{(max_level - level) * 2 + 1};
assert(size(first_block, 1) == degree);
assert(size(first_block, 2) == n);
if(multiply_left)
    if(do_transpose)
        Y(1:n, 1:cols_X) = first_block' * X(1:degree, 1:cols_X);
    else
        Y(1:degree, 1:cols_X) = first_block * X(1:n, 1:cols_X);
    end
else % right multiplying
    if(do_transpose)
        Y(1:rows_X, 1:degree) = X(1:rows_X, 1:n) * first_block';
    else
        Y(1:rows_X, 1:n) = X(1:rows_X, 1:degree) * first_block;
    end
end

% ------------------
% rest of the blocks
% ------------------
block_offset = (max_level - level) * 2 + 1; 
degree_start = degree + 1;

for i=0:(level-1)
   num_cells = 2^i;
   cell_size = n/num_cells;
   block = blocks{block_offset + i * 2 + 1};

   for j=0:(num_cells-1)
       cell_start = j * cell_size + 1;
       cell_end = cell_start + cell_size - 1; 
       degree_end = degree_start + degree - 1;
       
       if(multiply_left)
           if(do_transpose)
             Y(cell_start:cell_end, 1:cols_X) = ...
             Y(cell_start:cell_end, 1:cols_X) + ...
             block' * X(degree_start:degree_end, 1:cols_X);     
           else
             Y(degree_start:degree_end, 1:cols_X) = ...
             Y(degree_start:degree_end, 1:cols_X) + ...
             block * X(cell_start:cell_end, 1:cols_X); 
           end
       else
           if(do_transpose)
             Y(1:rows_X, degree_start:degree_end) = ...
             Y(1:rows_X, degree_start:degree_end) + ...
             X(1:rows_X, cell_start:cell_end) * block';
           else
             Y(1:rows_X, cell_start:cell_end) = ...
             Y(1:rows_X, cell_start:cell_end) + ...
             X(1:rows_X, degree_start:degree_end) * block;
           end
       end
       
       degree_start = degree_end + 1;
   end
end


