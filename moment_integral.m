function MomentValue = moment_integral(pde, fval_realspace, gfunc, t)

xmin = pde.dimensions{1,1}.domainMin;
xmax = pde.dimensions{1,1}.domainMax;
Lev = pde.lev_vec;
h = (xmax - xmin)/2^(Lev(1));
deg = pde.deg;
num_dimensions = length(pde.dimensions);
p = pde.params;

[quad_xx, quad_ww] = lgwt(deg, -1, 1);

quad_ww = 2^(-Lev(1))/2*quad_ww;

ww = repmat(quad_ww, 2^Lev(1), 1); 

if num_dimensions >= 2
    for i = 2:num_dimensions
         domainMin = pde.dimensions{1,i}.domainMin;
         domainMax = pde.dimensions{1,i}.domainMax;
         ww = kron(ww,ww)*(domainMax - domainMin);
    end
end

ww = ww.*(xmax - xmin);

[x, w] = lgwt(deg, 0, h);

points = [];
num = 1;

for j = 1:length(Lev)
    num = num*2^Lev(j)*deg;
    for i = 0:2^Lev(j)-1
        points = [points; xmin + x + i*h];
    end
end    
points2 = repmat(points, [2^Lev(1)*deg/2 1]);

%points2 = reshape(reshape(points2,num,num)',num*num,1);
mag = length(fval_realspace);
if num_dimensions == 1
    MomentValue = sum(ww.*fval_realspace.*gfunc(points,p,t));
else
    MomentValue = sum(ww.*fval_realspace.*gfunc(points2,p,t));
end
end
