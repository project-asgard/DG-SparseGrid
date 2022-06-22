function [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data)
%Create matrix to take (sparse) wavelet space to full-grid
%realspace with of same poly degree (B) and to constants only (B0)

assert(numel(pde.dimensions) == 2,'Works for 2 dimensions only');

[perm,iperm,pvec] = sg_to_fg_mapping_2d(pde,opts,A_data);
I = speye(numel(iperm));
I_FG = sparse(size(pvec,1),numel(iperm)); I_FG(pvec,:) = I(perm(pvec),:);
FMWT_x = OperatorTwoScale_wavelet2(opts.deg,pde.dimensions{1}.lev);
FMWT_v = OperatorTwoScale_wavelet2(opts.deg,pde.dimensions{2}.lev);
FMWT_2D = kron(FMWT_x,FMWT_v);

%Get constant coeffs data
B0 = kronrealspace2DtoDG(pde,opts,1)*FMWT_2D'*I_FG;
%Conver to cell average
B0 = B0*sqrt(2^(pde.dimensions{1}.lev-1)*2^(pde.dimensions{2}.lev-1));


B = kronrealspace2DtoDG(pde,opts,opts.deg)*FMWT_2D'*I_FG;

end

