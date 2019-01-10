function [Mat,Mat1,Mat2] = MatrixDiff_Momentum(Lev,Deg,LInt,LEnd,FunCoef,q_bcL ,q_bcR ,f_bcL ,f_bcR )
%function Mat to compute the Diff operator
% d/dx[(1-x^2)df/dx]
% Denote : q = (1-x^2)df/dx
% (d/dx q,v) =>               [  0  ]  [Q]
% (q,w) - ((1-x^2)df/dx,w) => [ I   ]  [F]
% Mat1 denotes the matrix for:
% (q,v) = Mat1*f
% ()
% This is the non-divergence form
%-----------------------------------------------------
if ~exist('FluxVal','var') || isempty(FluxVal)
    FluxVal = 1;
end

% q_bcL = 0;q_bcR = 1;f_bcL = 1;f_bcR = 0;
Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,@(x)1,   @(x)0,f_bcL,f_bcR);% equation for q

Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd, 1,FunCoef,@(x)0,q_bcL,q_bcR); % equation for f
Mat = Mat2*Mat1;

end
