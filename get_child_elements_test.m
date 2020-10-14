function tests = get_child_elements_test
tests = functiontests(localfunctions);
end

function test_2D_1(testCase)

[lev,pos,cnt] = get_child_elements ([3,2], [1,0], 6);

gold_lev = [    [4,2]' ...
                [4,2]' ...
                [3,3]' ...
                [3,3]' ];
            
gold_pos = [    [2,0]' ...
                [3,0]' ...
                [1,0]' ...
                [1,1]' ];

assert( norm(gold_lev'-lev) == 0);

end

function test_2D_2(testCase)

[lev,pos,cnt] = get_child_elements ([3,0], [1,0], 6);

gold_lev = [    [4,0]' ...
                [4,0]' ...
                [3,1]' ];
            
gold_pos = [    [2,0]' ...
                [3,0]' ...
                [1,0]' ];

assert( norm(gold_lev'-lev) == 0);

end
