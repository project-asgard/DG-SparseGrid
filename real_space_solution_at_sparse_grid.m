function fval_r = real_space_solution_at_sparse_grid(pde,fval)

N = numel(pde.elementsIDX);
nDims = numel(pde.dimensions);
deg = pde.dimensions{1}.deg;

%%
% Get center coordinates for all elements

% allCoords = zeros(N,nDims);
% for n=1:N
%     idx = pde.elementsIDX(n);
%     allCoords(n,:) = getMyRealSpaceCoord(pde,idx);
% end

%%
% For now just use a simple, ordered set of coordinates (the above would
% give the centers of all the elements)

allCoords(:,1) = [-1:0.01:+1];

nPts = numel(allCoords(:,1));

fval_r_D = zeros(nPts,nDims);
fval_r = zeros(nPts,1);

%%
% Loop over elements, and evaluate which other elements each contributes to
% and evaluate and add its contribution.

for n=1:N
    
    %%
    % Get this element coordinate vector
    
    idx = pde.elementsIDX(n);
    
    [myCoord,myCoordL,myCoordR] = getMyRealSpaceCoord(pde,idx);
    myLevVec = pde.elements.lev_p1(idx,:)-1;
    myPosVec = pde.elements.pos_p1(idx,:)-1;
    
    %%
    % Get this element coefficients
    
    elementDOF = deg^nDims;
    i1 = (n-1)*elementDOF+1;
    i2 = n*elementDOF;
    myCoeffs = fval(i1:i2);
    
    for d=1:nDims
        
        %%
        % Get this element spatial range for each dim
        
        myLev = myLevVec(d)
        myLev_n = 2^myLev;
        myLev_min = pde.dimensions{d}.domainMin;
        myLev_max = pde.dimensions{d}.domainMax;
        myLev_h = (myLev_max-myLev_min)/myLev_n;    
        
        %%
        % Find the idx's of all the elements this element contributes to
        
        xMin = myCoordL(d); % left side of element
        xMax = myCoordR(d); % right side of element
        
        allCoord = allCoords(:,d);
        
        iiList = find(allCoord >= xMin & allCoord <= xMax);
        
        %%
        % Get value of this element at the center of all those element it
        % contributes to (mapped to the [-1,1] range
        
        allCoordNorm = (allCoord - xMin)/(xMax-xMin)*2-1;
        
        p_val = lin_legendre(allCoordNorm(iiList),deg)*sqrt(1/myLev_h);
        
        %%
        % Now multiply the basis functions by their coefficients for all
        % the points we need them at
        
        fval_r_D(iiList,d) = fval_r_D(iiList,d) + p_val * myCoeffs;
        
    end
    
    %%
    % Now combine all the dimensions via multiplication (as this is point
    % wise)
    
    fval_r = fval_r_D(:,1);
    for d=2:nDims
        fval_r = fval_r .* fval_r_D(:,d);
    end
    
    figure(99);
    hold on
    x = allCoords(:,1);
    [X,I] = sort(x);
    plot(X,fval_r(I),'-o')
    
end

%%
% Test plot

x = allCoords(:,1);
[X,I] = sort(x);
plot(X,fval_r(I),'o')

end