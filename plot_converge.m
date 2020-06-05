colors = ['b', 'r', 'g', 'm', 'k'];
for i=1:5
myfilename = sprintf('f1d_x_FG-deg_%d_lev_4', i+2);
load(myfilename, 'x', 'f1d');
semilogy(x,f1d,colors(i), 'LineWidth', 2);
hold on
end
hold off
