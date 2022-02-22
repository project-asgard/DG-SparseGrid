function [F] = buildNonSeparableSource(r,th,k,f)
%Builds (f(x,v),phi_j)_{\W}

num_r = numel(r)-1;
num_th = numel(th)-1;

jac_r = (r(2)-r(1))/2;
jac_th = (th(2)-th(1))/2;

[quad_ref, w_ref]  = lgwt(10,-1,1);
quad_ref = quad_ref';

[leg_vals,~,~,~] = buildLegendre(10,k);
%test_ref = repmat(w_ref',k+1,1).*leg_vals; %Weights included with the test functions

weight = w_ref*w_ref';

%Create local function values and store in matrix
test_ref = zeros((k+1)^2,100);
for i=1:k+1
    for j=1:k+1
        %Order is [1,th,th^2,x,x th, x th^2,...];
        test_ref((i-1)*(k+1)+j,:) = ...
            reshape( ((leg_vals(j,:)/sqrt(jac_th))'*(leg_vals(i,:)/sqrt(jac_r))).*weight*jac_r*jac_th,[],1)';
    end
end

F = zeros((k+1)^2*num_r*num_th,1);
count = 1;
for i=1:num_r
    quad_r = quad_ref*(r(i+1)-r(i))/2 + (r(i+1)+r(i))/2;
    for j=1:num_th
        quad_th = quad_ref*(th(j+1)-th(j))/2 + (th(j+1)+th(j))/2;
        [R,TH] = meshgrid(quad_r,quad_th);
        f_val = reshape(f(R,TH),[],1);
        F(count:count+(k+1)^2-1) = test_ref*f_val;
        count = count+(k+1)^2;
    end
end

end

