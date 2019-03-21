function isok = check_FMWT( kdeg, Lev, FMWT )
% isok = check_FMWT( kdeg, Lev, FMWT )
% --------------------------------------------------
% script to check identical subblocks in FMWT matrix
% --------------------------------------------------
n = kdeg * 2^Lev;

      gmaxerr = 0;
      ip = 2*kdeg + 1;
      for lev=1:(Lev-1),
          ncells = 2^lev;
          isize = n/ncells;

          icell = 1;
          ipend = ip + kdeg-1;
          col1 = 1 + (icell-1)*isize;
          col2 = col1 + isize-1;
          Fmat = FMWT(ip:ipend,col1:col2);

          maxerr = 0;
          for icell=1:ncells,
             ipend = ip + kdeg-1;
                  col1 = 1 + (icell-1)*isize;
                  col2 = col1 + isize-1;
                  abserr = norm( FMWT(ip:ipend,col1:col2) - Fmat,1);
                  maxerr = max( maxerr,abserr );
             ip = ipend + 1;
          end;
          disp(sprintf('lev=%d, maxerr=%g', ...
                        lev,    maxerr));
        gmaxerr = max( gmaxerr, maxerr);
      end;

      tol = 1e-10;
      isok = (gmaxerr < tol);

end
