function p = fokkerplanck_parameters(opts)

%%
% Define a few relevant functions

nuEE = 1;
vT = 1;
switch opts.case_
    case 1
        delta = 0.042;
        Z = 1;
        E = 0.0025;
        tau = 10^5;
    case 2
        delta = 0.042;
        Z = 1;
        E = 0.25;
        tau = 10^5;
    case 3 % 6p1
        delta = 0.042;
        Z = 1;
        E = 0.0025;
        tau = 10^5;
    case 4
        delta = 0.3;
        Z = 5;
        E = 0.4;
        tau = 10^5;
    case 5 % for the runaway production rate electric field scan
        delta = 0.044;
        Z = 1;
        E = opts.cmd_args.E;
        tau = 1.2e5;    
end

phi = @phi_;
psi = @psi_;
f0_z = @f0_z_;
f0_p = @f0_p_;
soln_z = @soln_z_;
soln_p = @soln_p_;

gamma = @(p)sqrt(1+(delta*p).^2);
vx = @(p)1/vT*(p./gamma(p));

Ca = @(p)nuEE*vT^2*(psi(vx(p))./vx(p));

Cb = @(p)1/2*nuEE*vT^2*1./vx(p).*(Z+phi(vx(p))-psi(vx(p))+delta^4*vx(p).^2/2);

Cf = @(p)2*nuEE*vT*psi(vx(p));

    function ret = phi_(x)
        ret = erf(x);
    end

    function ret = psi_(x,t)     
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);
        ret = 1./(2*x.^2) .* (phi(x) - x.*dphi_dx);   
        ix = find(abs(x)<1e-5); % catch singularity at boundary
        ret(ix) = 0;
    end

    function ret = f0_z_(x)
        test = opts.case_;
        ret = zeros(size(x));
        switch test
            case 1
                ret = x.*0+1;
            case 2
                ret = x.*0+1;
            case 3
                h = [3,0.5,1,0.7,3,0,3];
                
                for l=1:numel(h)
                    
                    L = l-1;
                    P_m = legendre(L,x); % Use matlab rather than Lin's legendre.
                    P = P_m(1,:)';
                    
                    ret = ret + h(l) * P;
                    
                end
            case 4
                ret = x.*0 + 1;
            case 5
                ret = x.*0 + 1;
        end
    end

    function ret = f0_p_(x)
        test = opts.case_;       
        ret = zeros(size(x));
        switch test
            
            case 1
                for i=1:numel(x)
                    if x(i) <= 5
                        ret(i) = 3/(2*5^3);
                    else
                        ret(i) = 0;
                    end
                end
                
            case 2 
                a = 2;
                ret = 2/(sqrt(pi)*a^3) * exp(-x.^2/a^2);
            case 3
                ret = 2/(3*sqrt(pi)) * exp(-x.^2);
            case 4
                N = 1000;
                h = 20/N;
                Q = 0;
                Fun = @(p)exp(-2/delta^2*sqrt(1+delta^2*p.^2));
                for i = 1:N
                    x0 = (i-1)*h;
                    x1 = i*h;
                    [xi,w] = lgwt(20,x0,x1);
                    Q = Q+sum(w.*Fun(xi).*xi.^2);
                 end
                ret = exp(-2/delta^2*sqrt(1+delta^2*x.^2))/(2*Q);
            case 5 
                a = 1;
                ret = 2/(sqrt(pi)*a^3) * exp(-x.^2/a^2);
        end
    end

    function ret = soln_z_(x,~)
        ret = x.*0+1;
    end
    function ret = soln_p_(x,t)
        ret = 2/sqrt(pi) * exp(-x.^2);
    end

% throw all these into a parameters structure for use elsewhere

save_filename = 'parameters.mat';
save(save_filename);
p = load(save_filename);
delete(save_filename);

end
