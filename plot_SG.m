<<<<<<< HEAD
colors = ['b', 'r', 'g', 'm', 'k'];
for i=1:5
myfilename = sprintf('f1d_x-deg_%d_lev_4', i+2);
=======
numfiles = 3;
colors = ['m', 'k'];
for i=1:2
myfilename = sprintf('f1d_x-deg_%d_lev_5', i+5);
>>>>>>> 187aabc00617d0e839210c168e43b2041b0802e6
load(myfilename, 'x', 'f1d');
semilogy(x,f1d,colors(i), 'LineWidth', 2);
hold on
end
hold off