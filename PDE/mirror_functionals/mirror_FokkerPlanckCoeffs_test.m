function tests = mirror_FokkerPlanckCoeffs_test()

fh = localfunctions;
tests = functiontests(fh);

end

function test_FokkerPlanckCoeffs_Mawellian(testCase)

uVals = [3.5314*10^8, 7.3150*10^8, 1.1099*10^9, 1.4882*10^9, 1.8666*10^9, 2.2450*10^9, 2.6233*10^9, 3.0017*10^9, 3.3800*10^9, 3.7584*10^9]; %test velocity in m/s
%uVals = [3.5314*10^6];
temp = 1.6022e-10; %temperature in erg
k_e = 1.3807e-16; %Boltzmann constant in cgs 
n_o = 8e14; %equilibrium density in cm.^-3
e = 4.803*10^-10; %charge in Fr
m_e = 9.109*10^-28; %electron mass in g
c = 3*10^10; %cm/s

%setting up integrand matrix

params = mirror_parameters();
params.gamma = @(u) sqrt(1 + u.^2./c^2); %relativistic correction
params.a.Z = -1;
params.b.Z = -1;
params.a.m = m_e; %beam is electrons
params.b.m = m_e; %background is electrons
params.ln_delt = 15;
u_th = sqrt(2*temp/m_e); %thermal velocity in cm/s
rel_tol = 1e-6;
func = @(u,z) n_o.*exp(-u.^2./u_th^2).*(z.*0 + 1)./(pi^(3/2)*u_th^3);
test_Avals = [];
testVals = [];
numer_func = interpolate_numerical_distribution();

%testing Rosenbluth Coeffs for l = 0
lIndex = [0 1];
z = pi/4;
gold_M = 2.9798*10^-3;
for i = 1:length(uVals)
    uVal = uVals(i);
    for j = 1:length(lIndex)
        gamma_a = 4*pi*(params.a.Z)^2*e^4/((params.a.m)^2);
        numerVals = mirror_FokkerPlanckCoeffs(numer_func,uVal,lIndex(j),z,params);
        testVals = mirror_FokkerPlanckCoeffs(func,uVal,lIndex(j),z,params);
        %get FP coefficients from numerical function
        numer_Avals(i,j) = numerVals(1);
        numer_Bvals(i,j) = numerVals(2);
        numer_Cvals(i,j) = numerVals(3);
        numer_Dvals(i,j) = numerVals(4);
        numer_Evals(i,j) = numerVals(5);
        numer_Fvals(i,j) = numerVals(6);
        %get FP coefficients for analytic Maxwellian
        test_Avals(i,j) = testVals(1);
        test_Bvals(i,j) = testVals(2);
        test_Cvals(i,j) = testVals(3);
        test_Dvals(i,j) = testVals(4);
        test_Evals(i,j) = testVals(5);
        test_Fvals(i,j) = testVals(6);
    end
    total_numA(i) = sum(numer_Avals(i,:));
    total_numB(i) = sum(numer_Bvals(i,:));
    total_numC(i) = sum(numer_Cvals(i,:));
    total_numD(i) = sum(numer_Dvals(i,:));
    total_numE(i) = sum(numer_Evals(i,:));
    total_numF(i) = sum(numer_Fvals(i,:));
    total_A(i) = sum(test_Avals(i,:));
    total_B(i) = sum(test_Bvals(i,:));
    total_C(i) = sum(test_Cvals(i,:));
    total_D(i) = sum(test_Dvals(i,:));
    total_E(i) = sum(test_Evals(i,:));
    total_F(i) = sum(test_Fvals(i,:));
end
data = [uVals; total_A; total_B; total_C; total_D; total_E; total_F];
numer_data = [uVals; total_numA; total_numB; total_numC; total_numD; total_numE; total_numF];
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