function tests = mirror_Efunctional_test()
tests = functiontests(localfunctions);
end

function test_Efunctional_Simple(testCase)

func = @(x) x.*0 + 1;
maxVal = 0.5;
lIndex = 0;
testVal = mirror_Efunctional(func,maxVal,lIndex);
diff = abs(testVal - 0.00625);
verifyLessThan(testCase, diff, 1e-6);
lIndex1 = 1;
testVal1 = mirror_Efunctional(func,maxVal,lIndex1);
diff1 = abs(testVal1 - 0.002604);
verifyLessThan(testCase,diff1, 1e-6);

end