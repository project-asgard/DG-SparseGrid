<<<<<<< HEAD
numfiles = 3;
colors = ['m', 'k'];
for i=1:2
myfilename = sprintf('f1d_x-deg_%d_lev_5', i+5);
=======
colors = ['b', 'r', 'g', 'm', 'k'];
for i=1:5
myfilename = sprintf('f1d_x-deg_%d_lev_4', i+2);
>>>>>>> bdbfe7c86722d7397ea9677bc55ad3b8b1a354bb
load(myfilename, 'x', 'f1d');
semilogy(x,f1d,colors(i), 'LineWidth', 2);
hold on
end
<<<<<<< HEAD
hold off
=======
hold off
>>>>>>> bdbfe7c86722d7397ea9677bc55ad3b8b1a354bb
