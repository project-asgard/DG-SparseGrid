function cell = AllCell(Lev)
%===================================================
% Compute all the cell information from Lev
% Lev is 1xDim array
% cell is *xDim array, and each row corresponding to
%       varying Lev(i), i=1~Dim
%===================================================
dim = length(Lev);

% nMax = zeros(1,dim);
% for ii =1:dim
%     nMax(ii) = max(0,2^max(0,Lev(ii)-1)-1);
% end
% nz = prod(nMax+1);
%     
% 
% cell = zeros(nz,dim);

% for ii=1:dim
%     % Both of the following two methods can work
%     % Method 1
%     % tmp=repmat([0:nMax(ii)]',nz/(nMax(ii)+1),1);
%     % cell(:,ii)=tmp(:);
%     % Method 2
% 
%     cell(:,ii)=repmat([0:nMax(ii)]',nz/(nMax(ii)+1),1);
% 
% end


for i = 1:dim
    nMax = max(0,2^max(0,Lev(i)-1)-1);
    if nMax == 0
        value{i} = 0;
    else
    value{i}=[0:nMax];
    end
end
cell = allcomb(value{:});





end
