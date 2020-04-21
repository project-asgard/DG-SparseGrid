% Gamma^R term


Mat_Mass_p = MatrixMass(Lev,Deg,pInt,pEnd,@(x)(x.^2));

% Term 1
% x*d/dp(p^2*E*f)
Mat_Term1_p = MatrixGrad(Lev,Deg,pInt,pEnd,0,@(x)x.^3.*gamma(x)/tau,@(x)0,fp_bcL,fp_bcR); 
Mat_Term1_x = MatrixMass(Lev,Deg,LInt,LEnd,@(x)(1-x.^2));
Mat_Term1 =  kron(Mat_Term1_x,Mat_Term1_p);

% Term 2
% r = 1/p^2 * d/dp*[p^2*Cf*f]
% R = (Gp * I ) * F
% Here we should use BC for x variable
Mat_Term2_p = MatrixMass(Lev,Deg,pInt,pEnd,@(x)(-x.^2./(tau*gamma(x))));
Mat_Term2_x = MatrixGrad(Lev,Deg,LInt,LEnd,0,@(x)x.*(1-x.^2),@(x)0,fx_bcL,fx_bcR); 
Mat_Term2 =  kron(Mat_Term2_x,Mat_Term2_p);

Mat_R = (Mat_Term1 + Mat_Term2);