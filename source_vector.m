function fval = source_vector(HASHInv,pde,time)

% Returns the wavelet transformed source

% LevX = pde.dimensions{1}.lev;
% LevV = pde.dimensions{2}.lev;
% Deg = pde.dimensions{1}.deg;
% 
% Lmin = pde.dimensions{1}.domainMin;
% Lmax = pde.dimensions{1}.domainMax;
% Vmin = pde.dimensions{2}.domainMin;
% Vmax = pde.dimensions{2}.domainMax;

nDims = numel(pde.dimensions);
nSources = numel(pde.sources);

%%
% Loop over the number of sources, each of which has nDims + time elements.

fval = 0;
for s=1:nSources
    for d=1:nDims
        fList{d} = forwardMWT(pde.dimensions{d}.lev,pde.dimensions{d}.deg,...
            pde.dimensions{d}.domainMin,pde.dimensions{d}.domainMax,...
            pde.sources{s}{d},pde.params);
    end
    fs_d{s}{nDims+1} = pde.sources{s}{nDims+1}(time);
    
    ft = pde.sources{s}{nDims+1}(time);
    fval = fval + combine_dimensions_D(fList,ft,HASHInv,pde);
end

% fx1 = forwardMWT(LevX,Deg,Lmin,Lmax,pde.source1x,pde.params);
% fv1 = forwardMWT(LevV,Deg,Vmin,Vmax,pde.source1v,pde.params);
% ft1 = pde.source1t(time);
% 
% fx2 = forwardMWT(LevX,Deg,Lmin,Lmax,pde.source2x,pde.params);
% fv2 = forwardMWT(LevV,Deg,Vmin,Vmax,pde.source2v,pde.params);
% ft2 = pde.source2t(time);
% 
% fx3 = forwardMWT(LevX,Deg,Lmin,Lmax,pde.source3x,pde.params);
% fv3 = forwardMWT(LevV,Deg,Vmin,Vmax,pde.source3v,pde.params);
% ft3 = pde.source3t(time);

% fx1 = fs_d{1}{1};
% fv1 = fs_d{1}{2};
% ft1 = fs_d{1}{3};
% 
% fx2 = fs_d{2}{1};
% fv2 = fs_d{2}{2};
% ft2 = fs_d{2}{3};
% 
% fx3 = fs_d{3}{1};
% fv3 = fs_d{3}{2};
% ft3 = fs_d{3}{3};
% 
% fxList = {fx1,fx2,fx3};
% fvList = {fv1,fv2,fv3};
% ftList = {ft1,ft2,ft3};
% 
% fvalO = combine_dimensions_2(fxList,fvList,ftList,HASHInv,pde);

end