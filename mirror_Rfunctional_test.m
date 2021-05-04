function tests = mirror_Rfunctional_test()
tests = functiontests(localfunctions);
end

function test_Rfunctional_Maxwellian(testCase)

func = @(x) exp(-x.^2);
minVal = 0.01;
lIndex = 0;
testVal = mirror_Rfunctional(func,minVal,lIndex);
diff = abs(testVal - 0.50);
verifyLessThan(testCase, diff, 3e-5);
lIndex1 = 1;
testVal1 = mirror_Rfunctional(func,minVal,lIndex1);
diff1 = abs(testVal1 - 0.44313);
verifyLessThan(testCase,diff1, 3e-5);

end