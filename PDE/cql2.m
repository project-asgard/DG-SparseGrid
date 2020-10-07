function pde = cql2

% Velocity dynamics for eq 2.4.27 from Killen
%
% df/dt = 1/v^2 d/dv * (A*f + B*df/dv + C*df/dth) 
%          + 1/(v^2*sin(th)) * (D*f + E*df/dv + F*df/dth) 
%          + K*f + J
%
% Run with
%
% asgard(cql2,'timestep_method','BE','dt',1e-6);

    function m = maxwellian(v,vth)
        m = sqrt(2/pi) * (1/vth)^2 .* v^2 * exp(-v.^2/(2*vth^2));
    end

%% Dimensions
%

dim_v.domainMin = 0;
dim_v.domainMax = 1e6;
dim_v.init_cond_fn = @(v,p,t) maxwellian(v,vth)

dim_th.domainMin = 0;
dim_th.domainMax = pi;
div_th.init_cond_fn = @(th,p,t) th.*0+1;

end