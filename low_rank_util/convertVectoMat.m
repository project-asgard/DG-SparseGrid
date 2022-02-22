function [S] = convertVectoMat(x,v,k,u)
%Converts vector in dg-basis to a matrix such that
%u(x,v) = PHI(x)*S*PSI(v)'
%where PHI(x) = [phi_1(x),...,phi_{n_x}(x)] are the x-dofs and
%      PSI(v) = [psi_1(v),...,phi_{n_v}(v)] are the v-dofs.

%Ex: k=1, then each (k+1)^2=4 chunk of u:[u1;u2;u3;u4] needs to be rewitten
%as [u1 u2;u3 u4]. 

num_x = numel(x)-1;
num_v = numel(v)-1;

S = zeros((k+1)*num_x,(k+1)*num_v);

for i=1:num_x
   for j=1:num_v
       %Get (k+1)^2 u chunk we care about
       u_temp = u((k+1)^2*((i-1)*num_v+(j-1))+1:(k+1)^2*((i-1)*num_v+j));
       %Convert to (k+1)x(k+1) matrix and transpose it
       s_temp = reshape(u_temp,k+1,k+1)';
       %Store it
       S((i-1)*(k+1)+1:i*(k+1),(j-1)*(k+1)+1:j*(k+1)) = s_temp;
   end
end
       
end

