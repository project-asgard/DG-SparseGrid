function f = TimeAdvance(A,f,dt,slow)
%-------------------------------------------------
% Time Advance Method
% Input: Matrix:: A
%        Vector:: f
%        Time Step:: dt
% Output: Vector:: f
%-------------------------------------------------

f = RungeKutta3(A,f,dt,slow);

end

function fval = RungeKutta3(A,f,dt,slow)
%----------------------------------
% 3-rd Order Runge Kutta Method
%----------------------------------
ftmp = ApplyA(A,f,slow);
f_1 = f +dt *ftmp;
ftmp = ApplyA(A,f_1,slow);
f_2 = 3/4*f+1/4*f_1+1/4*dt*(ftmp);
ftmp = ApplyA(A,f_2,slow);
fval = 1/3*f+2/3*f_2+2/3*dt*(ftmp);

end

function ftmp = ApplyA(A,f,slow)

%-----------------------------------
% Multiply Matrix A by Vector f
%-----------------------------------
dof = size(f,1);
ftmp=sparse(dof,1);


if slow ~= 1
    
    for i=1:size(A,2)
        
        
        tmpA=A{i}.A1;
        tmpB=A{i}.A2;
        IndexI=A{i}.IndexI;
        IndexJ=A{i}.IndexJ;
        
        ftmp(IndexI)=ftmp(IndexI)+kron_mult2(tmpA,tmpB,f(IndexJ));
        
    end
    
else
    
    dof = size(A{1}.element_global_row_index);
    
    dofCnt = 1;
    conCnt = 1;
    
    for i=1:dof
        
        nConnected = A{1}.element_n_connected(i);
        
        for j=1:nConnected
            
            Index_I1 = A{1}.element_local_1_index(dofCnt);
            Index_I2 = A{1}.element_local_2_index(dofCnt);
            Index_J1 = A{1}.connected_local_1_index(conCnt);
            Index_J2 = A{1}.connected_local_2_index(conCnt);
            
            IndexI = A{1}.element_global_row_index(dofCnt);
            IndexJ = A{1}.connected_global_col_index(conCnt);
            
            % Apply term1 v.d_dx (vMassV . GradX) 
            
            tmpA = A{1}.vMassV(Index_I1,Index_J1);
            tmpB = A{1}.GradX(Index_I2,Index_J2);
            
            ftmp(IndexI)=ftmp(IndexI)+tmpA*tmpB*f(IndexJ);
            
            % Apply term 2 E.d_dv (EMassX . GradV)
            
            tmpA = A{1}.GradV(Index_I1,Index_J1);
            tmpB = A{1}.EMassX(Index_I2,Index_J2);
            
            ftmp(IndexI)=ftmp(IndexI)+tmpA*tmpB*f(IndexJ);
            
            conCnt = conCnt+1;
            
        end
        
        dofCnt = dofCnt + 1;
        
    end
    
end

end
