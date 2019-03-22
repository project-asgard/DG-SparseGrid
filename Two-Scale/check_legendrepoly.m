% simple script to check legendrepoly
%
% the first few Legendre polynomials are
% P_0(x) = 1
% P_1(x) = x
% P_2(x) = (1/2) * (3 * x^2 - 1)
% P_3(x) = (1/2) * (5 * x^3 - 3*x)
% P_4(x) = (1/8) * (35 * x^4 - 30 * x^2 + 3)
% P_5(x) = (1/8) * (63 * x^5 - 70*x^3 + 15*x)
% P_6(x) = (1/16) * (231 * x^6 - 315 * x^4 + 105 * x^2 - 5)
%
n = 1024;
x = linspace(-1,1,n);
x = reshape(x,n,1);

x2 = x.*x;
x3 = x2 .* x;
x4 = x2 .* x2;
x5 = x2.*x3;
x6 = x3.*x3;

P0 = ones(n,1);
P1 = x;
P2 = (1/2) * (3 * x2 - 1);
P3 = (1/2) * (5 * x3 - 3*x);
P4 = (1/8) * (35 * x4 - 30*x2 + 3);
P5 = (1/8) * (63 * x5 - 70*x3 + 15*x);
P6 = (1/16) * (231 * x6 - 315 * x4 + 105 * x2 - 5);

Fok = [P0,P1,P2,P3,P4,P5,P6];

ndeg = 6;
F = legendrepoly( ndeg, x );

err = norm(F-Fok,1);
relerr = err/norm(Fok,1);

tol = 1e-7;
isok = (relerr < tol) && (err < tol);
if (~isok),
    disp(sprintf('err=%g, relerr=%g', ...
                  err, relerr ));
    figure();
    subplot(2,1,1); plot( x, Fok); title('analytic expression');
    subplot(2,1,2); plot( x, F);   title('legendrepoly');
else
    disp(sprintf('ALL OK, err=%g, relerr=%g', err, relerr));
end;
