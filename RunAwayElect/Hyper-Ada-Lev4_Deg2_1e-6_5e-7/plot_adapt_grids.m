Lstart = -1;
Lend = 1;
Lev = 4;

sigma = 0.1;
f0 = @(x)( exp(-x.^2/sigma^2) );
exactf = @(x,t)(...
    (1-phi(x,t).^2)./(1-x.^2).*f0(phi(x,t)) ...
    );
CFL = 0.01;
h = (Lend-Lstart)/2^Lev;
dt = CFL*h^((Deg)/3);

% num_In = 2^4;%2^(10-Lev);
% plot_start = 2^(10-Lev)/2;

figure(1);
% hold on
% figure(2);
% hold on
for I = 1:ceil(1200/40)
    i = 40*I;
    load(['IHash_t',num2str(i),'.mat']);
    figure(1);subplot(1,2,1)
    plotgrid(IHash,Lstart,Lend,I,'b-x');hold on;
    
    Id = [];
    for j = 1:size(IHash,2)
        l1 = IHash{j};
        key = [l1(1),l1(2)];
        lr = l1(3);%Hash.(sprintf('i%g_',key));
        
        Id(Deg*(j-1)+[1:Deg]) = Deg*(lr-1)+[1:Deg];
    end
    
    MMeval = Meval(plot_start:num_In:end,Id);
    sol_num = MMeval*fval_MW;
       
    time = i*dt;
    sol_exact = exactf(xx_node,time);
    L2error = norm(sol_num(1:end-1)-sol_exact(1:end-1));
    
    figure(1);subplot(1,2,2)
    plot(xx_node,sol_num,'r-o',xx_node,sol_exact,'b--','LineWidth',2)
    title(['L2error is ',num2str(L2error),' Number of Grid is ', num2str(size(IHash,2))])
    hold off
    pause
end

    
    
    