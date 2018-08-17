
nx = 128;
ny = 128;

Nx = 16;
Ny = 16;

ratio_x=floor(nx/Nx);
ratio_y=floor(ny/Ny);

index=[1:(nx+1)*(ny+1)];
mapping=reshape(index,nx+1,ny+1);
Map=[];


count_y=1;
count=1;
Block = zeros(Nx,Ny,2);
Block_ind = zeros(Nx,Ny);

for I=1:Nx
    
    count_x=1;
    if I ~= Nx
        y_InEnd=min(count_y+ratio_y-1,1+ny);
    else
        y_InEnd=min(count_y+ratio_y,1+ny);
    end
    
    for J=1:Ny
        Block(I,J,1)=count;
        
        if J ~= Ny
            x_InEnd=min(count_x+ratio_x-1,1+nx);
            
        elseif J == Ny
            x_InEnd=min(count_x+ratio_x,1+nx);
%             y_InEnd=min(count_y+ratio_y-1,1+ny);
        end
           
        tmp=mapping(count_x:x_InEnd,count_y:y_InEnd);
        tmp=tmp(:);
        Map=[Map;tmp];
        
        count=count+size(tmp,1);
        count_x=x_InEnd+1;
        Block(I,J,2)=count-1;
        
        Block_ind(I,J) = Ny * (I-1) + J;
        
    end
    
    count_y=y_InEnd+1;

end
