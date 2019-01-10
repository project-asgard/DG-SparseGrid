function [x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot)

if ~exist('num_plot','var') || isempty(num_plot)
    num_plot = Deg;
end

% Plotting 
% num_plot = Deg;
[quad_x,quad_w]=lgwt(num_plot,-1,1);

% quad_x = [-1,1]';
h = (LEnd - LInt)/2^Lev;
p_val = legendre(quad_x,Deg);
for L=0:2^Lev-1
    %---------------------------------------------
    % Generate the coefficients for DG bases
    %---------------------------------------------
    
    Iu = [Deg*L+1:Deg*(L+1)];
    
    Iv = [num_plot*L+1:num_plot*(L+1)];
%         xi=h*(quad_x/2+1/2+L);
    
    x0 = LInt+L*h;
    x1 = x0+h;
    xi = quad_x*(x1-x0)/2+(x1+x0)/2;
    
    
    Meval(Iv,Iu)=sqrt(1/h)*p_val;
    x_node(Iv,1)=xi;
    
end