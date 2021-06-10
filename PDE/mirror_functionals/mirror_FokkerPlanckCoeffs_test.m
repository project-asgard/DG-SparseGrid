function tests = mirror_FokkerPlanckCoeffs_test()

fh = localfunctions;
tests = functiontests(fh);

end

function test_FokkerPlanckCoeffs_Mawellian(testCase)

uVals = [3.5314*10^6, 7.3150*10^6, 1.1099*10^7, 1.4882*10^7, 1.8666*10^7, 2.2450*10^7, 2.6233*10^7, 3.0017*10^7, 3.3800*10^7, 3.7584*10^7]; %test velocity in m/s
temp = 100; %temperature in eV
n_o = 8e20; %equilibrium density in eV
e = 1.602*10^-19; %charge in Coulombs
m_e = 9.109*10^-31; %electron mass in kg
u_th = sqrt(2*temp*e/m_e); %thermal velocity in m/s
rel_tol = 1e-6;
func = @(u,z) n_o.*exp(-u.^2./u_th^2).*(z.*0 + 1)./(pi^(3/2)*u_th^3);
test_Avals = [];
testVals = [];

%testing Rosenbluth Coeffs for l = 0
lIndex = 0;
z = pi/4;
gold_M = 2.9798*10^-3;
for i = 1:length(uVals)
    uVal = uVals(i);
    testVals = mirror_FokkerPlanckCoeffs(func,uVal,lIndex,z);
    test_Avals(i) = testVals(1);
end
rel_err = abs(testVal(1) - gold_M)/abs(gold_M);
verifyLessThan(testCase, rel_err, rel_tol);

end

function test_FokkerPlanckCoeffs_Maxwellian_l1(testCase)

uVal = 7.315*10^6; %test velocity in m/s
temp = 152; %temperature in eV
m_e = 9.109*10^-31; %electron mass in kg
u_th = sqrt(2*temp/m_e); %thermal velocity in m/s
rel_tol = 1e-6;
func = @(u) exp(-u.^2/u_th^2);
testVal = [];

%testing Rosenbluth Coeffs functional for l = 1
lIndex = 1;
z = pi/4;
gold_M = 0.8862;
testVal = mirror_FokkerPlanckCoeffs(func,uVal,lIndex,z);
rel_err = abs(testVal(1) - gold_M)/abs(gold_M);
verifyLessThan(testCase,rel_err, rel_tol);


end