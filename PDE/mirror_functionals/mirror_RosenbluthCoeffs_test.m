function tests = mirror_RosenbluthCoeffs_test()

fh = localfunctions;
tests = functiontests(fh);

end

function test_RosenbluthCoeffs_Mawellian(testCase)

uVal = 3.5314*10^6; %test velocity in m/s
temp = 100; %temperature in eV
n_o = 8e20; %equilibrium density in eV
e = 1.602*10^-19; %charge in Coulombs
m_e = 9.109*10^-31; %electron mass in kg
u_th = sqrt(2*temp*e/m_e); %thermal velocity in m/s
rel_tol = 1e-6;
func = @(u,z) n_o.*exp(-u.^2./u_th^2).*(z.*0 + 1)./(pi^(3/2)*u_th^3);
testVal = [];

%testing Rosenbluth Coeffs for l = 0
lIndex = 0;
gold_M = 2.9798*10^-3;
testVal = mirror_RosenbluthCoeffs(func,uVal,lIndex);
rel_err = abs(testVal(1) - gold_M)/abs(gold_M);
verifyLessThan(testCase, rel_err, rel_tol);

end

function test_RosenbluthCoeffs_Maxwellian_l1(testCase)

uVal = 7.315*10^6; %test velocity in m/s
temp = 152; %temperature in eV
m_e = 9.109*10^-31; %electron mass in kg
u_th = sqrt(2*temp/m_e); %thermal velocity in m/s
rel_tol = 1e-6;
func = @(u) exp(-u.^2/u_th^2);
testVal = [];

%testing Rosenbluth Coeffs functional for l = 1
lIndex = 1;
gold_M = 0.8862;
testVal = mirror_RosenbluthCoeffs(func,uVal,lIndex);
rel_err = abs(testVal(1) - gold_M)/abs(gold_M);
verifyLessThan(testCase,rel_err, rel_tol);


end