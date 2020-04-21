
% Gamma^C term

Mat_Mass_p = MatrixMass(Lev,Deg,pInt,pEnd,@(x)(x.^2));

% Term 1
% q = Ca * p^2 * df/dp
% Q = (Gp x I ) * F
[Mat_Term1_p] = MatrixDiff_Momentum(Lev,Deg,pInt,pEnd,@(x)(x.^2.*Ca(x)),qp_bcL,qp_bcR,fp_bcL,fp_bcR);
Mat_Term1_x = speye(DoFs,DoFs);
Mat_Term1 = kron(Mat_Term1_x,Mat_Term1_p);

% Term 2
% r = 1/p^2 * d/dp*[p^2*Cf*f]
% R = (Gp * I ) * F
% Here we should use BC for x variable
Mat_Term2_p = MatrixGrad(Lev,Deg,pInt,pEnd,1,@(x)x.^2.*Cf(x),@(x)0,fp_bcL,fp_bcR); 
Mat_Term2_x = speye(DoFs,DoFs);
Mat_Term2 = kron(Mat_Term2_x,Mat_Term2_p);
% Term 3
% r = Cb * (1-x^2) * df/dx
% R = (I x Gx ) * F
% Here we should use BC for x variable
[Mat_Term3_x] = MatrixDiff_Momentum(Lev,Deg,LInt,LEnd,@(x)(1-x.^2),qx_bcL,qx_bcR,fx_bcL,fx_bcR);

Mat_Term3_p = MatrixMass(Lev,Deg,pInt,pEnd,@(x)(Cb(x)./x.^2));
Mat_Term3 = kron(Mat_Term3_x,Mat_Term3_p);


Mat_C = Mat_Term1 + Mat_Term2 + Mat_Term3;