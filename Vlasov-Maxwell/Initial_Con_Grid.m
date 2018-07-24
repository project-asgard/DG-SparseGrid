function InCond = Initial_Con_Grid(HASH,inverseHash,InCond,pde)
%==========================================================
% This code is to generate initial conditions 
% which depend on Grids
% Input::
%       HASH and inverseHash
%       InCond
% Output::
%       finit
% Note::
%   Initial condition for f is on the sparse grids
%   Initial condition for E/B is on the full grids for 
%               scaling functions
%==========================================================
dof = size(inverseHash,1);
f0 = zeros(dof,1);

Dim = pde.DimX+pde.DimV;

DimX = pde.DimX;
DimV = pde.DimV;


if DimX == 1
    I{2} = 1;
elseif DimX ==2
    I{3} = 1;
end
if DimV == 1
    I{5} = 1;
elseif DimV ==2
    I{6} = 1;
end

fsep = InCond.fsep;
% Compute the coefficients in f on the sparse grids
for i = 1:dof
    ll = inverseHash{i};
    Ind_1D = ll(Dim*2+1:end);
    Index = Deg^(Dim)*(i-1)+1:Deg^(Dim)*i;
    for dim = 1:DimX
        loc = Ind_1D(dim);
        I{dim} = fsep(Deg*(loc-1):Deg*loc,dim);
    end
    for dim = 1:DimV
        loc = Ind_1D(dim+DimX);
        I{DimX+dim}=fsep(Deg*(loc-1):Deg*loc,dim+DimX);
    end
    f0(Index) = kron6D(I);
end

if Solver == 'VM' % compute the coefficients for E and B
    EDim = pde.E_Dim;
    BDim = pde.B_Dim;

    Esep = InCond.Esep;
    Bsep = InCond.Bsep;
    dof = 2^(maxLev);
    E0 = zeros(dof*3,1);
    B0 = zeros(dof*3,1);
    
    % need a smart way to full grids from 
    % (i1,i2,i3)-->global Index and
    % global -->(i1,i2,i3)
    
    dof = 2^(maxLev)^DimX;
    
    I = {};
    I{1} = 1;I{2} = 1;I{3} = 1;
    I{4} = 1;I{5} = 1;I{6} = 1;
    I{7} = 1;I{8} = 1;I{9} = 1;
    
    for dim = 1:DimX
        I{dim} = Esep(1:2^maxLev,dim);
    end
    for dim = 1:EDim
        E0((dim-1)*dof+[1:dof]) = kron3D(I{dim+[1:3]});
    end
    
    I = {};
    I{1} = 1;I{2} = 1;I{3} = 1;
    I{4} = 1;I{5} = 1;I{6} = 1;
    I{7} = 1;I{8} = 1;I{9} = 1;
    for bdim = 1:BDim
        for dim = 1:DimX
            I{dim} = Esep((bdim-1)*2^maxLev+[1:2^maxLev],dim);
        end
    
        B0((bdim-1)*dof+[1:dof]) = kron3D(I);
    end
end

if solver == 'VM'
    InCond = struct('f0',f0,'E0',E0,'B0',B0);
else
    InCond = struct('f0',f0);
end