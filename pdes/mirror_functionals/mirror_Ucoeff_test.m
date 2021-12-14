function tests = mirror_Ucoeff_test()
tests = functiontests(localfunctions);
end

function test_Ucoeff(testCase)

uVal = 1;
rel_tol = 1e-6;
testVal = [];
func = @(x,y) x.*y.^2;

%testing U_coeff for l = 0 Legendre polynomial
lIndex = 0;
gold_M = 2.9348;
testVal = mirror_Ucoeff(func,uVal,lIndex);
rel_err = abs(testVal - gold_M)/abs(gold_M);
verifyLessThan(testCase, rel_err, rel_tol);

end

function test_Ucoeff_l1(testCase)

uVal = 3.5314*10^6;
temp = 100; %temperature in eV
n_o = 8e20; %equilibrium density in eV
e = 1.602*10^-19; %charge in Coulombs
m_e = 9.109*10^-31; %electron mass in kg
params = mirror_parameters();
params.a.m = m_e; %beam is electrons
params.b.m = m_e; %background is electrons
params.ln_delt = 15;
u_th = sqrt(2*temp*e/m_e); %thermal velocity in m/s
rel_tol = 1e-6;
func = @(u,z) n_o.*exp(-u.^2./u_th^2).*(z.*0 + 1)./(pi^(3/2)*u_th^3);

%testing V_coeff for l = 1 Legendre polynomial
lIndex = [0,1,2,3];
gold_M = func(uVal,0);
for i = 1:length(lIndex)
    testVals(i) = mirror_Ucoeff(func,uVal,lIndex(i));
end
testVal = sum(testVals(:));
rel_err = abs(testVal - gold_M)/abs(gold_M);
verifyLessThan(testCase,rel_err, rel_tol);



end