function fval = initial_condition_vector(f,HASHInv,pde)

% fxList = {fx};
% fvList = {fv};
fList = {f};
ftList = {1};

fval = combine_dimensions_2(fList,ftList,HASHInv,pde);

end