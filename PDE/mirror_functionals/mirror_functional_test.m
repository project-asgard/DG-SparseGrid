function tests = mirror_functional_test()
tests = functiontests(localfunctions);
end

function test_functional_Maxwellian(testCase)

minVal = 0;
rel_tol = 7e-5;
func = @(x) exp(-x.^2);
testVal = [];

%testing M functional for l = 0
lIndex = 0;
gold_M = 0.5;
testVal = mirror_functional(func,minVal,lIndex);
rel_err = abs(testVal(1) - gold_M)/abs(gold_M);
verifyLessThan(testCase, rel_err, rel_tol);

%testing M functional for l = 1
lIndex = 1;
gold_M = 0.8862;
testVal = mirror_functional(func,minVal,lIndex);
rel_err = abs(testVal(1) - gold_M)/abs(gold_M);
verifyLessThan(testCase,rel_err, rel_tol);


minVal = 0.01; %resetting minVal for R functional

%testing R functional for l = 0
lIndex = 0;
gold_R = 0.5;
testVal = mirror_functional(func,minVal,lIndex);
rel_err = abs(testVal(3) - gold_R)/abs(gold_R);
verifyLessThan(testCase, rel_err, rel_tol);

%testing R functional for l = 1
lIndex = 1;
gold_R = 0.44313;
testVal = mirror_functional(func,minVal,lIndex);
rel_err = abs(testVal(3) - gold_R)/abs(gold_R);
verifyLessThan(testCase,rel_err, rel_tol);



end

function test_functional_Simple(testCase)

func = @(x) x.*0 + 1;
maxVal = 2;
rel_tol = 7e-5;

%testing N functional for l = 0
lIndex = 0;
gold_N = 8/3;
testVal = mirror_functional(func,maxVal,lIndex);
rel_err = abs(testVal(2) - gold_N)/abs(gold_N);
verifyLessThan(testCase, rel_err, rel_tol);

%testing N functional for l = 1
lIndex1 = 1;
gold_N = 4;
testVal = mirror_functional(func,maxVal,lIndex1);
rel_err = abs(testVal(2) - gold_N)/abs(gold_N);
verifyLessThan(testCase, rel_err, rel_tol);


maxVal = 0.5; %setting maxVal for E functional 

%testing E functional for l = 0
lIndex = 0;
gold_E = 0.00625;
testVal = mirror_functional(func,maxVal,lIndex);
rel_err = abs(testVal(4) - gold_E)/abs(gold_E);
verifyLessThan(testCase, rel_err, rel_tol);

%testing E functional for l = 1
lIndex1 = 1;
gold_E = 0.002604;
testVal = mirror_functional(func,maxVal,lIndex1);
rel_err = abs(testVal(4) - gold_E)/abs(gold_E);
verifyLessThan(testCase, rel_err, rel_tol);

end