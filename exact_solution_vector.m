function fval = exact_solution_vector(HASHInv,pde,time)

% Returns the wavelet transformed exact solution

% LevX = pde.dimensions{1}.lev;
% LevV = pde.dimensions{2}.lev;
% Deg  = pde.dimensions{1}.deg;

nDims = numel(pde.dimensions);

%%
% Loop over the number of dimensions (+time) to construct the analytic solution.

for d=1:nDims
    fList{d} = forwardMWT(pde.dimensions{d}.lev,pde.dimensions{d}.deg,...
        pde.dimensions{d}.domainMin,pde.dimensions{d}.domainMax,...
        pde.analytic_solutions_1D{d},pde.params);
end
% AS_d{nDims+1} = 

% fx1 = AS_d{1};
% fv1 = AS_d{2};
% ft1 = AS_d{3};

% fx1 = forwardMWT(LevX,Deg,pde.dimensions{1}.domainMin,pde.dimensions{1}.domainMax,pde.ExactFx,pde.params);
% fv1 = forwardMWT(LevV,Deg,pde.dimensions{2}.domainMin,pde.dimensions{2}.domainMax,pde.ExactFv,pde.params);
% ft1 = pde.ExactFt(time);

% fxList = {fx1};
% fvList = {fv1};
% ftList = {ft1};

% fval = combine_dimensions_2(fxList,fvList,ftList,HASHInv,pde);

% fList = {AS_d{1},AS_d{2}};
ft = pde.analytic_solutions_1D{nDims+1}(time);
fval = combine_dimensions_D(fList,ft,HASHInv,pde);

end