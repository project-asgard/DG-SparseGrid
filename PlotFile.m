load('Data_Plot.mat')
figure;
subplot(1,2,1)
mesh(x_2D_plot,y_2D_plot,FucVal);


[f_nd_t] = converttoRealSpace(2,lev,deg,'SG',[0,0],[1,1],F0,6);
subplot(1,2,2)
mesh(reshape(f_nd_t,64,64))