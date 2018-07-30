C=[2.0000e+00   4.0000e+00   2.9530e-01   1.8886e+00
2.0000e+00   5.0000e+00   1.7203e-01   1.5001e+00
2.0000e+00   6.0000e+00   6.2190e-02   8.2806e-01
2.0000e+00   7.0000e+00   2.9634e-02   4.7486e-01
2.0000e+00   8.0000e+00   1.0839e-02   2.5014e-01
2.0000e+00   9.0000e+00   3.8700e-03   1.2554e-01];     %cfl=0.0068


A=[2.0000e+00   4.0000e+00   1.6057e-01   1.2770e+00
2.0000e+00   5.0000e+00   5.9727e-02   9.1575e-01
2.0000e+00   6.0000e+00   3.0425e-02   4.0578e-01
2.0000e+00   7.0000e+00   9.5843e-03   1.8581e-01
2.0000e+00   8.0000e+00   3.0692e-03   7.5016e-02
2.0000e+00   9.0000e+00   7.4479e-04   2.7476e-02];     %cfl=0.0045


U=[2.0000e+00   4.0000e+00   1.5184e-01   4.5193e-01
2.0000e+00   5.0000e+00   5.8034e-02   2.2244e-01
2.0000e+00   6.0000e+00   1.1166e-02   1.0084e-01
2.0000e+00   7.0000e+00   2.8075e-03   3.4250e-02
2.0000e+00   8.0000e+00   1.4295e-03   1.1689e-02
2.0000e+00   9.0000e+00   1.0272e-03   4.2937e-03];     %cfl=0.00625

for i =2:6
    c(i,1)=log2(C(i-1,4)/C(i,4));
    a(i,1)=log2(A(i-1,4)/A(i,4));
    u(i,1)=log2(U(i-1,4)/U(i,4));
end

n=A(:,2);
e1=C(:,4);
e2=A(:,4);
e3=U(:,4);
plot(n,e1,n,e2,n,e3,n,2.^(-n),n,2.^(-1.5*n))
legend('central flux with CFL=0.0068','alternating flux with CFL=0.0045','up-winding flux with CFL=0.00625')
