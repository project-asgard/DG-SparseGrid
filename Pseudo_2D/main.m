function main
%----------------------------------------------------------
% This code implements the generation of global matrix for
% 2D Poisson equation
% Maxmium Level for mesh: Np=5;
% Polynomial Degree: k=3; (quadratic element)
%----------------------------------------------------------
% The following is the most inefficient implementation
% load 2D Mesh: H Index-->(nx,ix,px,ny,iy,py)
% Property for Mesh: nx+ny<=Np
%                    ix=[0:2^max(nx-1,0)-1],iy=[0:2^max(ny-1,0)-1]
%                    px=[0:k],py=[0:k]
% load 1D Mesh: M (n,i,p)--> Index for x or y component
% load Stiffness matrix S with size (DOF1DxDOF1D)
% load Mass matrix I (identity) with size (DOF1DxDOF1D)
% dof denotes the size for the 2D solver

Np=5;k=3;Dim=2;

load('Data.mat','H','S');
I=speye(2^Np*k);

dof=size(H,1);
A_s=sparse(dof,dof);

for IndexI=1:dof
    % 
    % for x-component
    nx=H(IndexI,2);
    ix=H(IndexI,3);
    px=H(IndexI,4);
    % for y-component
    ny=H(IndexI,5);
    iy=H(IndexI,6);
    py=H(IndexI,7);
    for IndexJ=1:dof
            % for x-component
            mx=H(IndexJ,2);
            jx=H(IndexJ,3);
            kx=H(IndexJ,4);
            % for y-component
            my=H(IndexJ,5);
            jy=H(IndexJ,6);
            ky=H(IndexJ,7);
            
            % pull from S and I to generate global matrix
            %  S*I+I*S
            Ix=Map_1D(nx,ix,px,k);Jx=Map_1D(mx,jx,kx,k);
            Iy=Map_1D(ny,iy,py,k);Jy=Map_1D(my,jy,ky,k);

            val=S(Ix,Jx)*I(Iy,Jy)+I(Ix,Jx)*S(Iy,Jy);
            A_s=A_s+sparse(IndexI,IndexJ,val,dof,dof);
            
        
    end
end

spy(A_s)


function z=Map_1D(level,cell,poly,k)
    if level==0
        z=poly;
    else
        z=k*2^(level-1)+k*cell+poly;
    end
end

end
