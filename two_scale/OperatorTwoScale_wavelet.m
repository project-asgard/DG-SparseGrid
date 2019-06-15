function FMWT_COMP = OperatorTwoScale_wavelet(maxDeg,maxLev)
% FMWT_COMP = OperatorTwoScale_wavelet(maxDeg,maxLev)
%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: maxDeg
%        Level: Np
% Output: Convert Matrix: FMWT_COMP
%**********************************

% % Load G0 and H0 from file
idebug = 0;
maxCells = 2^maxLev;

fileName = ['two_scale/two_scale_rel_',num2str(maxDeg),'.mat'];

if exist(fileName,'file') == 2
    load(fileName);
else
    disp('Generating two-scale file');
    [H0,G0,scale_co,phi_co]=MultiwaveletGen(maxDeg);
    save(fileName,'H0','G0','scale_co','phi_co');
end

H0(find(abs(H0)<1e-5))=0; % Why are we doing this?
G0(find(abs(G0)<1e-5))=0;

H1 = zeros(maxDeg);
G1 = zeros(maxDeg);

for j_x = 1:maxDeg
    for j_y = 1:maxDeg
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2)  )*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(maxDeg+j_x+j_y-2))*G0(j_x,j_y);
    end
end

use_sparsify = 1;
FMWT = zeros(maxDeg*maxCells);
FMWT2 = zeros(maxDeg*maxCells);


% Unroll the matlab for easier porting
porting = 1;

if porting
    for j=1:maxCells/2
        
        %FMWT( maxDeg*(j-1)+1 : maxDeg*j, 2*maxDeg*(j-1)+1 : 2*maxDeg*j )=[H0 H1];
        
        rs = maxDeg*(j-1)+1;
        cs = 2*maxDeg*(j-1)+1;
        
        for j_x = 1:maxDeg
            for j_y = 1:maxDeg
                
                FMWT2( rs+j_x-1, cs+j_y-1 ) = H0(j_x,j_y);
                FMWT2( rs+j_x-1, cs+maxDeg+j_y-1 ) = H1(j_x,j_y);
            end
        end
        
        %FMWT( maxDeg*(j+maxCells/2-1)+1 : maxDeg*(j+maxCells/2), 2*maxDeg*(j-1)+1 : 2*maxDeg*j) = [G0 G1];
        
        rs = maxDeg*(j+maxCells/2-1)+1;
        cs = 2*maxDeg*(j-1)+1;
        
        for j_x = 1:maxDeg
            for j_y = 1:maxDeg
                
                FMWT2( rs+j_x-1, cs+j_y-1 ) = G0(j_x,j_y);
                FMWT2( rs+j_x-1, cs+maxDeg+j_y-1 ) = G1(j_x,j_y);
            end
        end
    end
    
end

for j=1:maxCells/2
    % The reverse order from Lin
    FMWT(maxDeg*(j-1)+1:maxDeg*j,2*maxDeg*(j-1)+1:2*maxDeg*j)=[H0 H1];
    FMWT(maxDeg*(j+maxCells/2-1)+1:maxDeg*(j+maxCells/2),2*maxDeg*(j-1)+1:2*maxDeg*j) = [G0 G1];
end

if (idebug >= 1),
   figure();
   spy( FMWT );
   title(sprintf('line 86:N=%d,maxDeg=%d,maxCells=%d',...
		 size(FMWT,1),maxDeg,maxCells));
end;



if porting; assert(isequal(FMWT,FMWT2)); end

if (use_sparsify),
  FMWT_COMP = speye(maxDeg*maxCells,maxDeg*maxCells);
  FMWT_COMP2 = speye(maxDeg*maxCells,maxDeg*maxCells);
else
  FMWT_COMP = eye(maxDeg*maxCells);
  FMWT_COMP2 = eye(maxDeg*maxCells);
end;

n = floor( log2(maxCells) );
% n = maxCells;

for j=1:n
    cFMWT = FMWT;
    cFMWT2 = FMWT;
    % Reverse the index in matrix from Lin
    if j>1
      if (use_sparsify),
        nzmax = (maxDeg*maxCells)*maxDeg;
        cFMWT = sparse([],[],[], maxDeg*maxCells, maxDeg*maxCells,nzmax);
        cFMWT2 = sparse([],[],[],maxDeg*maxCells, maxDeg*maxCells,nzmax);
      else
        cFMWT = zeros(maxDeg*maxCells);
        cFMWT2 = zeros(maxDeg*maxCells);
      end;
        
        cn = 2^(n-j+1)*maxDeg;
        cnr=maxCells*maxDeg-cn;
        
        if porting
            
            % cFMWT(cn+1:maxDeg*maxCells,cn+1:maxDeg*maxCells)=eye(maxCells*maxDeg-cn);
            
            rs = cn+1;
            cs = cn+1;

            use_for_loop = 0;
            if (use_for_loop),
              for ii=0:maxDeg*maxCells - (cn+1)
                for jj=0:maxDeg*maxCells - (cn+1)
                    if (ii==jj)
                        cFMWT2(rs+ii,cs+jj) = 1;
                    end
                end
              end
            else
              % ---------------------------
              % short-cut loop to assign 1 
              % ---------------------------
              for ii=0:maxDeg*maxCells - (cn+1)
                jj = ii;
                cFMWT2(rs+ii,cs+jj) = 1;
              end
            end;
            
            % cFMWT(1:cn/2,1:cn)=FMWT(1:cn/2,1:cn);
            
            for ii=1:cn/2
                for jj=1:cn
                    cFMWT2(ii,jj) = FMWT(ii,jj);
                end
            end
            
            % cFMWT(cn/2+1:cn,1:cn)=FMWT(maxDeg*maxCells/2+1:maxDeg*maxCells/2+cn/2,1:cn);
            
            rs = maxDeg*maxCells/2+1;
            for ii=0:cn/2-1
                for jj=1:cn
                    cFMWT2(cn/2+1+ii,jj) = FMWT(rs+ii,jj);
                end
            end
            
        end
        
	if (use_sparsify),
          cFMWT(cn+1:maxDeg*maxCells,cn+1:maxDeg*maxCells)=speye(maxCells*maxDeg-cn,maxCells*maxDeg-cn);
	else
          cFMWT(cn+1:maxDeg*maxCells,cn+1:maxDeg*maxCells)=eye(maxCells*maxDeg-cn);
	end;
        cFMWT(1:cn/2,1:cn)=FMWT(1:cn/2,1:cn);
        cFMWT(cn/2+1:cn,1:cn)=FMWT(maxDeg*maxCells/2+1:maxDeg*maxCells/2+cn/2,1:cn);

	if (idebug >= 1),
	  figure();
          subplot(2,1,1);
          spy(cFMWT);
	  title(sprintf('line 178: j=%d, cFMWT, N=(%d,%d)',...
			 j, size(cFMWT,1),size(cFMWT,2) ));
          subplot(2,1,2);
          spy(FMWT_COMP);
	  title(sprintf('FMWT COMP,N=(%d,%d)', ...
			size(FMWT_COMP,1),size(FMWT_COMP,2)  ));
	end;
        
        if porting; assert(isequal(cFMWT,cFMWT2)); end
        
    end

    if (use_sparsify),
      cFMWT = sparsify_matrix(cFMWT);
      FMWT_COMP = sparsify_matrix(FMWT_COMP); 
    end;

    FMWT_COMP = cFMWT*FMWT_COMP;
    
    if porting
        
      if (use_sparsify),
        cFMWT2 = sparsify_matrix(cFMWT2);
        FMWT_COMP2 = sparsify_matrix(FMWT_COMP2);
      end;

        FMWT_COMP2 = cFMWT2*FMWT_COMP2;
        
        assert(isequal(cFMWT,cFMWT2));
        assert(isequal(FMWT_COMP,FMWT_COMP2));
        
    end
    
end

