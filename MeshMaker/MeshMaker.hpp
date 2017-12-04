//
//  MeshMaker.hpp
//  MMeshMaker
//
//  Created by CJiao on 7/15/17.
//
//

#ifndef MeshMaker_hpp
#define MeshMaker_hpp

#include <stdio.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <vector>
#include <numeric>
#include <cmath>
#include <unordered_map>
using namespace std;

//hash function
namespace std {
    template<>
    struct hash<vector<int> >{
        size_t operator()(vector<int> const& vec) const {
            size_t seed = vec.size();
            for(auto& i : vec) {
                seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
            
        }
    };
}

// The class that we are interfacing to, this class would likely be
// defined elsewhere such as a library in ny real situation
class MeshMaker
{
public:
    
    MeshMaker () { }
    MeshMaker(int dim, int np, int d, int dm);
    ~MeshMaker () { }

    
    void init(int dim, int np, int d, int dm);
    
    vector<int> findk(int value);
    vector<int> findvs(vector<vector<int> > &keys);
    vector<vector<int> > findks(vector<int> values);
    void refine(vector<int> adds, vector<vector<int>> dims, vector<int> removes, int Nmax);
    void refine1(vector<int> adds, vector<int> removes);

    
    
    vector<vector<int> > GenerateKeys();
    vector<vector<int> > GenerateKeysOnlyLeaf();
    
    unordered_map<vector<int>, unsigned long> GenerateWholeMap();
    unordered_map<vector<int>, unsigned long> GenerateLeafMap();
    
    vector<vector<int> > get_LCD()const;
    int get_Dim()const;
    int get_numKeys()const;
    int get_bigmapsize();
    
    unordered_map<vector<int>, unsigned long> get_WholeMap()const;
    unordered_map<vector<int>, unsigned long> get_LeafMap()const;
    
    //compute level
    //input:Dim,Np
    //combinations of k numbers sums less than n
    vector<vector<int> > Compute_Level(int k,//Dim
                                       int n,//Np
                                       int isLeaf);
    //compute cell from level
    vector<vector<int> > Compute_LevelCell(vector<vector<int> > keys_L);
    
    //compute degrees from level&cell
    vector<vector<int> > Compute_LevelCellDegree(vector<vector<int> > keys_LC);
    
    void CombinationSumN0( int k,
                          int n,
                          vector<vector<int> > &res,
                          vector<int> &sol);
    void CombinationSumN1( int k,
                          int n,
                          vector<vector<int> > &result,
                          vector<int> &sol);
    
    //for one level, compute combinations of cells
    vector<vector<int> > combinationCell(const vector<vector<int> > keys_Level);
    
    //Descartes product of multi vectors
    void cart_product(
                      vector<vector<int> >& rvvi,  // final result
                      vector<int>&  rvi,   // current result
                      vector<vector<int> >::const_iterator me, // current input
                      vector<vector<int> >::const_iterator end); // final input
    
    void PrintKeys();
    void PrintWholeMap();
    void PrintLeafMap();
    bool comp(pair<vector<int>,int> a, pair<vector<int>,int> b) {
        return a.second < b.second;
    }
    
    friend ostream & operator<<(ostream& os, const MeshMaker& mm){
        for(const auto & row : mm.get_LCD()){
            for(const auto & s : row){
                os<<" "<<s;
            }
        }
        return os;
    }
    
    bool operator==(const MeshMaker& mm) const{
        return (mm.get_LCD() == keys_LCD);
    }
    vector<vector<int>> findParents(vector<int> key);
    
private:
    //init parameters
    int Dim;
    int Np;
    int D;
    int Dm;
    
    //keys for elements
    vector<vector<int> > keys_LCD;
    vector<vector<int> > keys_LC;
    vector<vector<int> > keys_L;
    
    //number of initial keys
    int num_Keys;
    unordered_map<vector<int>, unsigned long> WholeMap;
    unordered_map<vector<int>, unsigned long> LeafMap;
    
    int total;//for generating values of map
    
};

#endif /* MeshMaker_hpp */

