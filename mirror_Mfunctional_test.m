function tests = mirror_Mfunctional_test()
tests = functiontests(localfunctions);
end

function test_Mfunctional_Maxwellian(testCase)

func = @(x) exp(-x.^2);
minVal = 0;
lIndex = 0;
testVal = mirror_Mfunctional(func,minVal,lIndex);
diff = abs(testVal - 0.5);
verifyLessThan(testCase, diff, 1e-6);
lIndex1 = 1;
testVal1 = mirror_Mfunctional(func,minVal,lIndex1);
diff1 = abs(testVal1 - 0.8862);
verifyLessThan(testCase,diff1, 5e-5);

end