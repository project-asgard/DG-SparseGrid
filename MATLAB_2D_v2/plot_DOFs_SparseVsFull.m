
% N=10000;
D=1:0.1:6;

fact=10*8/(1024^5)

figure
N1=100;
semilogy(D,N1.^D*fact,'r-o','LineWidth',2,'MarkerSize',8);
hold on
plot(D,N1.*(log10(N1)).^(D-1)*fact,'b-+','LineWidth',2,'MarkerSize',8)
N2=1000;
plot(D,N2.^D*fact,'r-d','LineWidth',2,'MarkerSize',8);
hold on
plot(D,N2.*(log10(N2)).^(D-1)*fact,'b-*','LineWidth',2,'MarkerSize',8)
N3=10000;
plot(D,N3.^D*fact,'r-^','LineWidth',2,'MarkerSize',8);
hold on
plot(D,N3.*(log10(N3)).^(D-1)*fact,'b--','LineWidth',2,'MarkerSize',8)
set(gca,'fontsize',20)
legend(['Full Grid ',num2str(N1)],['Sparse Grid ',num2str(N1)],...
    ['Full Grid ',num2str(N2)],['Sparse Grid ',num2str(N2)],...
    ['Full Grid ',num2str(N3)],['Sparse Grid ',num2str(N3)]...
    )


