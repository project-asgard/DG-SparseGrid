function [u] = convertMattoVec(x,v,k,S)
%Does the opposite of convertMattoVec.  Look for the documentation there.
num_x = numel(x)-1;
num_v = numel(v)-1;

u = zeros((k+1)^2*num_x*num_v,1);

for i=1:num_x
   for j=1:num_v
       %Get (k+1)\times(k+1) S chunk we care about
       s_temp = S((i-1)*(k+1)+1:i*(k+1),(j-1)*(k+1)+1:j*(k+1));
       %Transpose it and convert to vector
       u_temp = reshape(s_temp',[],1);
       %Store it
       u((k+1)^2*((i-1)*num_v+(j-1))+1:(k+1)^2*((i-1)*num_v+j)) = u_temp;
   end
end
       
end

