function tests = mirror_Nfunctional_test()
tests = functiontests(localfunctions);
end

function test_Nfunctional_Simple(testCase)

func = @(x) x.*0 + 1;
maxVal = 2;
lIndex = 0;
testVal = mirror_Nfunctional(func,maxVal,lIndex);
diff = abs(testVal - 8/3);
verifyLessThan(testCase, diff, 1e-6);
lIndex1 = 1;
testVal1 = mirror_Nfunctional(func,maxVal,lIndex1);
diff1 = abs(testVal1 - 4.00000);
verifyLessThan(testCase,diff1, 1e-6);

end