numfiles = 3;
colors = ['m', 'k'];
for i=1:2
myfilename = sprintf('f1d_x-deg_%d_lev_5', i+5);
load(myfilename, 'x', 'f1d');
semilogy(x,f1d,colors(i), 'LineWidth', 2);
hold on
end
hold off