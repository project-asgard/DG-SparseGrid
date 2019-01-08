function Mat = MatrixDiff(Lev,Deg,LInt,LEnd,FunCoef,FluxVal)
%function Mat to compute the Diff operator
% d/dx[(1-x^2)df/dx]
% Denote : q = (1-x^2)df/dx
% (d/dx q,v) =>                         [ 0  ]  [Q]
% (q,w) - ((1-x^2)df/dx,w) => [ I   ]  [F]
%-----------------------------------------------------
if ~exist('FluxVal','var') || isempty(FluxVal)
    FluxVal = 1;
end
L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

% Mat = sparse(2*DoF,2*DoF);

% quad_num = 10;

Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,@(x)(1));
Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd,FunCoef);

II = speye(DoF);

% Mat = [II Mat1;Mat2 sparse(DoF,DoF)];

Mat = Mat2*Mat1;
figure;spy(Mat)

end