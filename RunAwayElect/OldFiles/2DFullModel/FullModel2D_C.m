% This is a test for FullModel 2D problem
% Gamma^C
% Term 1
% q = Ca * p^2 * df/dp
% Q = (Gp x I ) * F
% Term 2
% r = Cb * (1-x^2) * df/dx
% R = (I x Gx ) * F
% Term 3
% 1/p^2*d/dp [q] + 1/p^2 d/dp [Cf*p^2*f] + 1/p^4 d/dx[r]


Gam = 1;
f_bcL = 1; f_bcR = 0;
q_bcL = 0; q_bcR = 1;

Mat_Mass_p = MatrixMass(Lev,Deg,pInt,pEnd,@(x)(x.^2));

% Term 1
% q = Ca * p^2 * df/dp
% Q = (Gp x I ) * F
f_bcL = 1; f_bcR = 0;
q_bcL = 0; q_bcR = 1;

% PDE_FP2;
[Mat_Term1_p] = MatrixDiff_Momentum(Lev,Deg,pInt,pEnd,@(x)(x.^2.*Ca(x)),q_bcL,q_bcR,f_bcL,f_bcR);

Mat_Term1 = kron(Mat_Term1_p,speye(DoFs,DoFs));

% Term 2
% r = 1/p^2 * d/dp*[p^2*Cf*f]
% R = (Gp * I ) * F
% Here we should use BC for x variable
Mat_Term2_p = MatrixGrad(Lev,Deg,pInt,pEnd,1,@(x)x.^2.*Cf(x),@(x)0,f_bcL,f_bcR); 
Mat_Term2 = kron(Mat_Term2_p,speye(DoFs,DoFs));

% Term 3
% r = Cb * (1-x^2) * df/dx
% R = (I x Gx ) * F
% Here we should use BC for x variable
f_bcL = 1; f_bcR = 1;
q_bcL = 0; q_bcR = 0;

[Mat_Term3_x] = MatrixDiff_Momentum(Lev,Deg,LInt,LEnd,@(x)(1-x.^2),q_bcL,q_bcR,f_bcL,f_bcR);

Mat_Term3_p = MatrixMass(Lev,Deg,pInt,pEnd,@(x)(Cb(x)./x.^2));
Mat_Term3 = kron(Mat_Term3_p,Mat_Term3_x);


Mat_C = Mat_Term1 + Mat_Term2;% + Mat_Term3;