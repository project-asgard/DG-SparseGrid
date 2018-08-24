# Runaway Implementation

### Changes in Runaway-electron problem

The file /PDE/Vlasov_RE.m contains the specific functions (listed at the end) 
and weak-formulation components (listed in the beginning) 
for the runaway-electron problem. 

### Comparing Function Outputs

Matlab: (from this folder)
```
run Print_Test.m
```

C++:
```
cd FK6D/include
g++ -Wall -std=c++11 runaway-test.cpp
./a.out
```