function v=lin_dlegendre(x,k)

% recursive dlegendre routine
    function dp = dlg_dp(x,n)
        if n == 0
            dp = x.* 0;
        else
            pmn_1 = legendre(n-1,x);
            pn_1 = pmn_1(1,:)';
            dp_1 = dlg_dp(x,n-1);
            dp = n * pn_1 + x.*dp_1;
        end
    end

use_matlab_legendre = true;
if use_matlab_legendre
    
    for kk=1:k
        n=kk-1;
        v(:,kk) = dlg_dp(x,n)*sqrt((n+1)*2-1);
    end
    
else
    
    % Legendre Polynomials with degree k on [-1,1]
    
    if k>7
        error('ERROR: max supported deg=9');
    end
    
    v(:,1)=x-x;
    v(:,2)=x-x+1*sqrt(2*2-1);
    v(:,3)=1/2*(3*2*x)*sqrt(3*2-1);
    v(:,4)=1/2*(5*3*x.^2-3)*sqrt(4*2-1);
    v(:,5)=1/8*(35*4*x.^3-30*2*x)*sqrt(5*2-1);
    v(:,6)=1/8*(63*5*x.^4-70*3*x.^2+15)*sqrt(6*2-1);
    v(:,7)=1/16*(231*6*x.^5-315*4*x.^3+105*2*x)*sqrt(7*2-1);
    v=v(:,1:k);
    
end

end
