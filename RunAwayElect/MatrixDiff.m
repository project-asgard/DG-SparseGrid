function [Mat,Mat1,Mat2] = MatrixDiff(Lev,Deg,LInt,LEnd,FunCoef)
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

% - derivative of FunCoef

% FunCoef2 = @(x)( 2*x);

% FunCoef2 = @(x)(-2*x);
FunCoef2 = @(x)-dFunCoef(x);


% % % FullPitchAngleDyn
% % FunCoef2 = @(x)( 2*x);
% % % % use alternating flux
% % Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd, 1,@(x)1,@(x)0,0,0);% equation for q
% % Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,FunCoef,FunCoef2,1,1); % equation for f
% % % 
% % % % central flux
% % % % Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,0,@(x)1,@(x)0,0,2);
% % % % Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd,0,FunCoef,FunCoef2,2,0);
% % % 
% % % 
% % Mat = Mat1*Mat2;

% % % use alternating flux
% % Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd, 1,@(x)1,@(x)0,0,2);% equation for q
% % Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,FunCoef,FunCoef2,2,0); % equation for f
% % 
% % % central flux
% % % Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,0,@(x)1,@(x)0,0,2);
% % % Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd,0,FunCoef,FunCoef2,2,0);
% % 
% % 
% % Mat = Mat1*Mat2;

q_bcL = 0;q_bcR = 0;f_bcL = 1;f_bcR = 1;
Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,1,@(x)1,   @(x)0,f_bcL,f_bcR);% equation for q

Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,FunCoef,@(x)0,q_bcL,q_bcR); % equation for f
Mat = Mat2*Mat1;

end

function y = dFunCoef(x)
psi = @(x)(1./x.^2.*(erf(x)-2*x/sqrt(pi).*exp(-x.^2)));
dpsi = @(x)(2*exp(-x.^2)/sqrt(pi)-(erf(x)-2*x.*exp(-x.^2)/sqrt(pi))./x.^3);
if abs(x)<1e-5
    y = 0;
else
    y = -(1/2)*(erf(x)-2*x.*exp(-x.^2)/sqrt(pi))./x.^2+2*x.*exp(-x.^2)/sqrt(pi);
end
end