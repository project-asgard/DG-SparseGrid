Lstart = -1;
Lend = 1;
Lev = 3;
figure(1);
% hold on
% figure(2);
% hold on
[Hash,IHash,FineIndex] = HashTable1D(Lev);
option.l = 'r-o';option.c = 'r';
plotgrid(IHash,Lstart,Lend,0,option);hold on;

for I = 1:ceil(1320/40)
    i = 40*I;
    load(['Damp_IHash_t',num2str(i),'.mat']);
    figure(1);%subplot(1,2,1)
    if mod(I,2)==0
        option.l = 'r-o';option.c = 'r';
        plotgrid(IHash,Lstart,Lend,I,option);hold on;
    else
        option.l = 'k-o';option.c = 'k';
        plotgrid(IHash,Lstart,Lend,I,option);hold on;
    end
        
    
%     Id = [];
%     for j = 1:size(IHash,2)
%         l1 = IHash{j};
%         key = [l1(1),l1(2)];
%         lr = l1(3);%Hash.(sprintf('i%g_',key));
%         
%         Id(Deg*(j-1)+[1:Deg]) = Deg*(lr-1)+[1:Deg];
%     end
%     
%     MMeval = Meval(plot_start:num_In:end,Id);
%     sol_num = MMeval*fval_MW;
%        
%     time = i*dt;
%     sol_exact = exactf(xx_node,time);
%     L2error = norm(sol_num(1:end-1)-sol_exact(1:end-1));
%     
%     figure(1);subplot(1,2,2)
%     plot(xx_node,sol_num,'r-o',xx_node,sol_exact,'b--','LineWidth',2)
%     title(['L2error is ',num2str(L2error),' Number of Grid is ', num2str(size(IHash,2))])
%     hold off
%     pause
end