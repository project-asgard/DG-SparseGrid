# Runaway Implementation

### Changes in Runaway-electron problem

The file /PDE/Vlasov_RE.m contains the specific functions (listed at the end) 
and weak-formulation components (listed in the beginning) 
for the runaway-electron problem. 

### Comparing Function Outputs

Matlab:
```
run Print_Test.m
```

C++: (from FK6D branch 'feature/adams-thing')
```
cd FK6D/include
g++ -Wall -std=c++11 runaway-test.cpp
./a.out
```