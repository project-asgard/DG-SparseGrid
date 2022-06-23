function [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data)
%Create matrix to take (sparse) wavelet space to full-grid
%realspace with of same poly degree (B) and to constants only (B0)

persistent FG2DG FG2DG0

assert(numel(pde.dimensions) == 2,'Works for 2 dimensions only');

[perm,iperm,pvec] = sg_to_fg_mapping_2d(pde,opts,A_data);
I = speye(numel(iperm));
I_FG = sparse(size(pvec,1),numel(iperm)); I_FG(pvec,:) = I(perm(pvec),:);

if isempty(FG2DG)
    FMWT_x = OperatorTwoScale_wavelet2(opts.deg,pde.dimensions{1}.lev);
    FMWT_v = OperatorTwoScale_wavelet2(opts.deg,pde.dimensions{2}.lev);
    FMWT_2D = kron(FMWT_x,FMWT_v);

    lev_vec = [pde.dimensions{1}.lev,pde.dimensions{2}.lev];
    FG2DG0 = kronrealspace2DtoDG(lev_vec,opts.deg,1)*FMWT_2D';
    FG2DG  = kronrealspace2DtoDG(lev_vec,opts.deg,opts.deg)*FMWT_2D';
end

%Get constant coeffs data
B0 = FG2DG0*I_FG*sqrt(2^(pde.dimensions{1}.lev-1)*2^(pde.dimensions{2}.lev-1));
%B0_2 = kronrealspace2DtoDG(pde,opts,1)*FMWT_2D'*I_FG;
%Convert to cell average
%B0 = B0_2*sqrt(2^(pde.dimensions{1}.lev-1)*2^(pde.dimensions{2}.lev-1));

B = FG2DG*I_FG;
%B = kronrealspace2DtoDG(pde,opts,opts.deg)*FMWT_2D'*I_FG;

end

