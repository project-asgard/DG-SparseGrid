function pde=VlasovMaxwell1
%=================================================
% Numerical Example for Vlasov Maxwell Equation
% Straming Weibel Instability Model
% 1X2V system with f(x2,v1,v2,t)
% E=[E1(x2,t),E2(x2,t),0] and B=[0,0,B3(x2,t)]
%=================================================

% parameters
k_0 = 0.2; 
delta = 0.5; 
beta = 0.01;
xi_1 = 0.3;
xi_2 = 0.3;
b = 0.001;

DimX = 1; DimV = 2;
Lmax = 2*pi/k_0; 
Vmax = 1.2;

E_Dim = 2; B_Dim = 1;

%-------------------------
% Initial condition for f
%-------------------------
    function f=Fx_0(x)
        % Initial condition for x variable
        
        nz = length(x);
        %f = zeros(nz,DimX);
        f = ones(nz,DimX);
    end

    function f=Fv_0(v)
        % Initial condition for v variable
        nz = length(v);
        f = zeros(nz,DimV);
        f(:,1) = exp(   delta *exp(-(v-xi_1).^2/beta)+...
                     (1-delta)*exp(-(v-xi_2).^2/beta)   );
        f(:,2) = exp(-v.^2/beta)/(pi*beta);
    end
    function f=Fxv_0(x,v)
        % Initial condition for x and v variable
        nz = length(x);
        f = zeros(nz,Dim);
        f=Fv_0(v).*Fx_0(x);
    end

%-------------------------
% Initial condition for E
%-------------------------
    function f=E_0(x)
        % Initial condition for E=[E1(x2,t=0),E2(x2,t=0),0]
        nz = length(x);
        f = zeros(nz,E_Dim);      
    end

%-------------------------
% Initial condition for B
%-------------------------
    function f=B_0(x)
        % Initial condition for B=[0,0,B3(x2,t=0)]
        nz = length(x);
        f = zeros(nz,B_Dim);
        f = b*sin(k_0*x);      
    end

    
pde = struct('Fx_0',@Fx_0,'Fv_0',@Fv_0,'Fxv_0',@Fxv_0,...
             'E_0',@E_0,'B_0',@B_0,'Vmax',Vmax,'Lmax',Lmax,...
             'DimX',DimX,'DimV',DimV,'E_Dim',E_Dim,'B_Dim',B_Dim);
end