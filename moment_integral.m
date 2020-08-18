function MomentValue = moment_integral(pde, fval_realspace, gfunc)

xmin = pde.dimensions{1,1}.domainMin;
xmax = pde.dimensions{1,1}.domainMax;
Lev = pde.lev_vec;
h = (xmax - xmin)/2^(Lev(1));
deg = pde.deg;
num_dimensions = length(pde.dimensions);

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

%[x, w] = lgwt(deg, 0, h);

%points = [];
%num = 2^Lev*deg;


%for i = 0:2^Lev(i)-1
%    points = [points; xmin + x + i*h];
%end

%points2 = repmat(points, [num 1]);

%points2 = reshape(reshape(points2,num,num)',num*num,1);
mag = length(fval_realspace);
MomentValue = sum(ww.*fval_realspace(2:mag-1).*gfunc);

end
