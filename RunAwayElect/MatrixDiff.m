function Mat = MatrixDiff(Lev,Deg,LInt,LEnd,FunCoef)
%function Mat to compute the Diff operator
% d/dx[(1-x^2)df/dx]
% Denote : q = (1-x^2)df/dx
% (d/dx q,v) =>                         [ 0  ]  [Q]
% (q,w) - ((1-x^2)df/dx,w) => [ I   ]  [F]
% This is the non-divergence form
%-----------------------------------------------------
if ~exist('FluxVal','var') || isempty(FluxVal)
    FluxVal = 1;
end

% - derivative of FunCoef

FunCoef2 = @(x)( 2*x);

% FunCoef2 = matlabFunction( diff(FunCoef(x)) );

% use alternating flux
Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,1);
Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,FunCoef,FunCoef2);

% % central flux
% Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,0);
% Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd,0,FunCoef,FunCoef2);


Mat = Mat1*Mat2;
% figure;spy(Mat)

end