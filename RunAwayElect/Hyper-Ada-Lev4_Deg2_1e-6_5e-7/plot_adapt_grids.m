figure;hold on
for i = 20:40:1200
    load(['IHash_t',num2str(i),'.mat']);
    plotgrid(IHash,Lstart,Lend,1/30*i,'b-x');
end