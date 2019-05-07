function [f_nd_t] = converttoRealSpace(pde,Dim,Lev_solution,Deg,gridType,Lmin,Lmax,fval_t,Lev_viz)
%
% [f_nd_t] = converttoRealSpace(Dim,Lev_solution,Deg,gridType,Lmin,Lmax,fval,Lev_viz)
% convert from wavelet coefficient to real space  on grid at level Lev_viz
% Note grid at level Lev_viz may be finer than solution grid at level Lev_solution
% -----------------------------------------------------------------------
if ~exist('gridType','var') || isempty(gridType)
    gridType = 'SG';
else
    if strcmp(gridType,'SG') || strcmp(gridType,'FG')
    else
        error("gridType must be set to 'SG' or 'FG'"); 
    end
end

Lev = Lev_viz;

idebug = 0;
use_find = 0;

if (idebug >= 1),
    disp(sprintf('converttoRealSpace:Dim=%d,Deg=%d,Lev_solution=%d, Lev_viz=%d', ...
		  Dim,Deg,   Lev_solution, Lev_viz));
end;


[nW,nT] = size(fval_t);
if (idebug >= 1),
   disp(sprintf('converttoRealSpace:nW=%d,nT=%d,size(fval_t)', ...
                                    nW,   nT));
   size(fval_t)

   for k=1:numel(Lmin),
           disp(sprintf('Lmin(%d)=%g,Lmax(%d)=%g', ...
                              k,Lmin(k),   k,Lmax(k) ));
   end;
end;

LevX = Lev;
LevV = Lev;


[HASH,HASHInv] = HashTable(pde,Lev_solution,Dim,gridType);
if (idebug >= 1),
    disp(sprintf('numel(HASH)=%d, numel(HASHInv)=%d', ...
        numel(HASH),    numel(HASHInv) ));
    disp(sprintf('Dim=%g,gridType=%s,Lev_solution', ...
                  Dim,   gridType));
    Lev_solution
end;

n = 2^Lev*(Deg-1);

for idim=1:Dim,
% ----------------------------------------------
% Note that fxv_loc and fv_loc are small matrices
% may consider fully dense storage
% sparsity ratio about 18.75 percent
% ----------------------------------------------
   xv{idim} = linspace(Lmin(idim),Lmax(idim),n );
   fxv_loc{idim} = EvalWavPoint4(Lmin(idim),Lmax(idim),Lev,Deg,xv{idim});
end;

if (idebug >= 1),
    for idim=1:Dim,
       disp(sprintf('idim=%d,size(fxv_loc{idim})=(%d,%d)', ...
                     idim,   size(fxv_loc{idim},1), size(fxv_loc{idim},2) ));
    end;
end;

%%
% Plot basis functions in 1D

% figure(11);
% bFns1D = fxv_loc{1};
% [nB,nP] = size(bFns1D);
% iii=1;
% for l=1:numel(bFns1D(:,1))
%     ax=subplot(nB,1,iii);
%     plot(bFns1D(iii,:));
%     iii=iii+1;
% end

steps = 1:nT;


nS = numel(steps);
nsteps = numel(steps);
if (idebug >= 1),
     disp(sprintf('nS=%d,nsteps=%d', nS,nsteps));
end;

%n_t = zeros(nX,nS);

%% nD Reverse MWT

% Evaluation starts from here
% Method :: evaluate every basis at the given points (xx) and (vv)
% fx_loc = phi_x (xx) with size DoF1DxNum_xx
% fv_loc = phi_y (yy) with size DoF1DxNum_vv
% Loop over every point with (xx_i,yy_j) and combine the value
%   for fx_loc(:,i) and fv_loc(:,j)

% if we have computed fval, we can just start with following
% Input:: Lstart, Lend, Lev, Deg, xx, yy


%% Main loop

    
    time_main = tic();
    
    f_nd_t = zeros( [n^Dim, nsteps]);
    
    nHash = numel(HASHInv);
    for i=1:nHash,
        ll = HASHInv{i};
        levels = ll(1:Dim);
        icells = ll( Dim + (1:Dim) );
        
           

        for k=1:Dim,
            index_Ik{k} = LevCell2index( levels(k), icells(k));
            index_IkDeg{k} = (index_Ik{k}-1)*Deg + (1:Deg);
        end;
        
        %% 
        % Check LevCell2index output against approach in other chunks of code. 
        
        nDims = Dim;
        for d=1:nDims
            this_idx1D = ll(nDims*2+d);          
            if (this_idx1D ~= index_Ik{d}),
              disp(sprintf('d=%d,this_idx1d=%d,index_Ik{d}', ...
                            d,   this_idx1d,   index_Ik{d}));
            end;
            assert(this_idx1D==index_Ik{d});
        end

        Index = (i-1)*(Deg^Dim) + (1:(Deg^Dim));
       
        % --------------------------------------------------
        % Note K1 and K2 should have contiguous index values
        % --------------------------------------------------
        sizes = zeros(1,Dim);
        for k=1:Dim,
          if (use_find),
            Kindex{k} = find( sum(fxv_loc{k}(index_Ik{k},:) ~= 0,1));
            K_start = min(Kindex{k});
            K_end = max(Kindex{k});
          else
            [K_start,K_end] = calc_overlap_index( levels(k),icells(k),Lev,Deg);
            Kindex{k} = K_start:K_end;
          end;
          Kstart(k) = K_start;
          Kend(k) = K_end;
        end;

        sizes = n * ones(1,Dim);
        Jindex = index_nd( Dim, Kstart,Kend, sizes );
        Jindex = reshape( Jindex, numel(Jindex),1);
        if (idebug >= 1),
                disp(sprintf('converttoRealSpace:numel(Jindex)=%g', ...
                                                 numel(Jindex)));
        end;

        
        % --------------------------------------------
        % note convert to full storage to use efficient
        % DGEMM operations
        %
        % Note in C = kron(A,B), index in B vary faster
        %
        % C( [ib,ia], [jb,ja]) = B(ib,jb) * A(ia,ja)
        % --------------------------------------------
        nkron = Dim;
        for k=1:Dim,
           kk = Dim-k+1;
           Acell{kk} = transpose( full( fxv_loc{k}(index_IkDeg{k},Kindex{k}) ) );
           if (idebug >= 1),
              disp(sprintf('index_Ik{%d}',k));
              index_Ik{k}
              disp(sprintf('Kindex{%d}=%d:%d',...
                   k,min(Kindex{k}), max(Kindex{k})));
           end;
        end;
        Xmat = fval_t(  Index, 1:nsteps);

        if (idebug >= 1),
                disp(sprintf('converttoRealSpace:size(Xmat) '));
                size(Xmat)
        end;
        
        % -----------------------------------------------
        % Ymat = Bmat * Xmat * transpose(Amat)
        % Ymat(K2,K1) = transpose(fx_loc(index_I2, K2)) *
        %                  Xmat(Index_I2,Index_I1) *
        %                    ( fv_loc(index_I1,K1) )
        % -----------------------------------------------
        
        Ymat = kron_multd(nkron,Acell,Xmat);
        if (idebug >= 1),
          disp(sprintf('nkron=%d',nkron));
          for k=1:nkron,
             disp(sprintf('k=%d,size(Acell{k})=(%d,%d)', ...
                   k,        size(Acell{k},1), size(Acell{k},2)));
          end;
        disp(sprintf('after kron_multd,nsteps=%d',...
                      nsteps));
        disp(sprintf('size(Ymat)'));
        size(Ymat)
        disp(sprintf('size(f_nd_t)'));
        size(f_nd_t)
        end;
        
        f_nd_t(Jindex, 1:nsteps) = f_nd_t(Jindex, 1:nsteps) + ...
            reshape( Ymat, [numel(Jindex), nsteps]);
        
    end;
    
    elapsed_time = toc( time_main );
    disp(sprintf('elapsed time for %d time steps is %g sec', ...
        nsteps,   elapsed_time ));
    
end


