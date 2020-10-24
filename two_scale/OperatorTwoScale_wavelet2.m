function [FMWT, blocks] = OperatorTwoScale_wavelet2(deg,nLev)
% FMWT_COMP = OperatorTwoScale_wavelet2(maxDeg,nLev)
%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: maxDeg
%        Level: Np
% Output: Convert Matrix: FMWT_COMP
%         Dense blocks also output for testing
%**********************************

% % Load G0 and H0 from file
idebug = 0;
assert(nLev >= 0);

asgard_root = get_root_folder();

fileName = [asgard_root,'/two_scale/two_scale_rel_',num2str(deg),'.mat'];

if exist(fileName,'file') == 2
    load(fileName);
else
    disp('Generating two-scale file');
    [H0,H1,G0,G1,scale_co,phi_co]=MultiwaveletGen(deg);
    save(fileName,'H0','G0','scale_co','phi_co');
end

if(nLev == 0)
   FMWT = eye(deg, deg);
   blocks = {};
   return
end

% ---------------------------------------------
% some entries should  theoretically be zero
% but may be very small due to numerical roundoff
% not critical to computations
% ---------------------------------------------
tol = 10^4 * eps;
H0(find(abs(H0) < tol))=0; 
G0(find(abs(G0) < tol))=0;

H1 = zeros(deg,deg);
G1 = zeros(deg,deg);

for j_x = 1:deg
    for j_y = 1:deg
        % -----------------------------------------
        % note  (-1)^k  is  (is_even(k)) ? 1 : -1,  
	% where is_even(k) = (mod(k,2) == 0);
        % no need to call std::pow()
        % -----------------------------------------
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2)  )*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(deg+j_x+j_y-2))*G0(j_x,j_y);
    end
end

% nLev = round( log2( maxLev) );
% isok = (2^nLev == maxLev);
% if (~isok),
%    error('Operator_TwoScale_wavelet2: maxLev=%d is not a power of 2', ...
% 	                              maxLev );
%    return;
% end;
n =  deg * 2^(nLev);

% -----------------------------------------------
% need index space  -1..(nLev-1) so add ioff=2
% so  blocks{ ioff + (-1)} is accessing blocks{1}
% -----------------------------------------------
ioff = 2;

use_portable = 1;
if (use_portable),
   Gmat = zeros(deg, n);
   Hmat = zeros(deg, n);
   Hmat(1:deg,1:deg) = eye(deg,deg);
else
Hmat =  [ eye(deg,deg)];
end;


for j=(nLev-1):-1:0,
   ncells = 2^(j);
   isize = n/ncells;
   isizeh = isize/2;

   if (use_portable),
     Htmp = zeros(deg, isizeh);
     Htmp(1:deg, 1:isizeh) = Hmat(1:deg,1:isizeh);


     Gmat(1:deg, 1:isizeh) = G0(1:deg,1:deg) * Htmp(1:deg, 1:isizeh );
     Gmat(1:deg, isizeh + (1:isizeh)) = G1(1:deg,1:deg) * Htmp(1:deg, 1:isizeh);


     Hmat(1:deg, 1:isizeh) = H0(1:deg,1:deg) * Htmp(1:deg,1:isizeh);
     Hmat(1:deg, isizeh + (1:isizeh)) = H1(1:deg,1:deg) * Htmp(1:deg,1:isizeh);
   else
     Gmat = [G0 * Hmat, G1 * Hmat];
     Hmat = [H0 * Hmat, H1 * Hmat];
   end;
   blocks{(j+1) * 2} = Gmat(1:deg,1:isize);
   if (j == 0)
       h_cols = n;
   else
       h_cols = isize;
   end;
       blocks{j*2 + 1} = Hmat(1:deg,1:h_cols);
end;

if (idebug >= 1),
   for j=1:(nLev*2),
      Fmat = blocks{j};
      nrow = size(Fmat,1);
      ncol = size(Fmat,2);
      disp(sprintf( 'j=%d, size(blocks{j})=(%d,%d)', ...
		     j,    nrow,ncol  ));
   end;
end;

% ---------------------------------------------------
% form sparse matrix from small dense blocks stored in blocks{-1..(nLev-1)}
% ---------------------------------------------------
nzmax = deg * (nLev)*n;
FMWT = sparse( [], [], [], n,n, nzmax );
irow = 1;
j = -1;
FMWT( 1:deg, 1:n) = blocks{ioff+(-1)};

irow = 1 + deg;
for j=0:(nLev-1),
    ncells = 2^j;
    isize = n/ncells;
    Fmat = blocks{ioff + j*2};
    for icell=1:ncells,
       j1 = (icell-1)*isize  + 1;
       j2 = j1 + isize - 1;
       i1 = (icell-1)*deg + irow;
       i2 = i1 + deg-1;
       FMWT( i1:i2, j1:j2) = Fmat(1:deg, 1:isize);
    end;
    irow = irow + deg * ncells;
end;

end








