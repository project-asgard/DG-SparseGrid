function tests = mirror_Vcoeff_test()
tests = functiontests(localfunctions);
end

function test_Vcoeff(testCase)

uVal = 1;
rel_tol = 1e-6;
testVal = [];
func = @(x,y) x.*y.^2;

%testing V_coeff for l = 0 Legendre polynomial
lIndex = 0;
gold_M = 2.9348;
testVal = mirror_Vcoeff(func,lIndex,uVal);
rel_err = abs(testVal - gold_M)/abs(gold_M);
verifyLessThan(testCase, rel_err, rel_tol);

end

function test_Vcoeff_l1(testCase)

uVal = 1;
rel_tol = 5e-6;
testVal = [];
func = @(x,y) x.*y.^2;

%testing V_coeff for l = 1 Legendre polynomial
lIndex = 1;
gold_M = -6.5735;
testVal = mirror_Vcoeff(func,uVal,lIndex);
rel_err = abs(testVal - gold_M)/abs(gold_M);
verifyLessThan(testCase,rel_err, rel_tol);



end