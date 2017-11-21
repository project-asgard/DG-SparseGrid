//
//  MeshMaker.cpp
//  MMeshMaker
//
//  Created by CJiao on 7/15/17.
//
//

#include "MeshMaker.hpp"

MeshMaker::MeshMaker(int dim, int np, int d, int dm){
    
    Dim = dim;
    Np = np;
    D = d;
    Dm = dm;
}

void MeshMaker::init(int dim, int np, int d, int dm){
    //Meaning of mGrid(Dim, Np, D, Dm, Is_Leaf_Only);
    //    MeshMaker mm(2, 2, 1, 4);
    Dim = dim;
    Np = np;
    D = d;
    Dm = dm;
    
    
    //A hashtable containnng all the elements
    WholeMap = GenerateWholeMap();
    //     PrintWholeMap();
    
    //A hashtable only containnng leaf elements(no child)
    LeafMap = GenerateLeafMap();
    //     PrintLeafMap();
}


vector<vector<int> > MeshMaker::GenerateKeys(){
    keys_L = Compute_Level(Dim, Np, 0);
    keys_LC = Compute_LevelCell(keys_L);
    keys_LCD = Compute_LevelCellDegree(keys_LC);
    return keys_LCD;
}
vector<vector<int> > MeshMaker::GenerateKeysOnlyLeaf(){
    keys_L = Compute_Level(Dim, Np, 1);
    keys_LC = Compute_LevelCell(keys_L);
    keys_LCD = Compute_LevelCellDegree(keys_LC);
    return keys_LCD;
    
}
unordered_map<vector<int>, unsigned long> MeshMaker::GenerateWholeMap(){
    vector<vector<int> > keys = GenerateKeys();
    unsigned long count = keys.size();
    for (unsigned long i = 0; i < count; ++i){
        total = i;
        WholeMap[keys[i]] = total;
    }
    return WholeMap;
}
unordered_map<vector<int>, unsigned long> MeshMaker::GenerateLeafMap(){
    vector<vector<int> > keys = GenerateKeysOnlyLeaf();
    unsigned long count = keys.size();
    for (unsigned long i = 0; i < count; ++i){
        LeafMap[keys[i]] = WholeMap[keys[i]];
    }
    return LeafMap;
}

vector<vector<int> > MeshMaker::get_LCD()const{
    return keys_LCD;
}
int MeshMaker::get_numKeys()const{
    return num_Keys;
}
int MeshMaker::get_Dim()const{
    return Dim;
}
unordered_map<vector<int>, unsigned long> MeshMaker::get_WholeMap()const{
    return WholeMap;
}
unordered_map<vector<int>, unsigned long> MeshMaker::get_LeafMap()const{
    return LeafMap;
}

vector<int> MeshMaker::findvs(vector<vector<int> > &keys){
    vector<int> values;
    for(const auto & key:keys){
        // Check if key  exists in the map
        if (WholeMap.count(key) > 0){
            unsigned long tmp = WholeMap[key];
            values.push_back(tmp);
        }
        else{
            values.push_back(-1);
//            cout<<"key not exist. "<<endl;
        }
        
    }
    return values;
}
/*
 int MeshMaker::findv(vector<int>  &key){
 
 // Check if key  exists in the map
 if (WholeMap.count(key) > 0){
 return( WholeMap[key] );
 }
 else{
 cout<<"key not exist. "<<endl;
 return  -1;
 }
 }
 */

//
void MeshMaker::refine(vector<int> adds, vector<vector<int>> dims, vector<int> removes, int Nmax){
    
    //first, remove nodes from 2 maps
    //set value = -2
    //make sure not remove all nodes added last time
    //Otherwise, have to move a node from wholemap to leafmap.a little difficult
    vector<vector<int> > keys_remove = findks(removes);
    

    
    for(auto &key : keys_remove){
        // Check if key  exists in the map
        if (LeafMap.count(key) <= 0)
        {
//            cout<<"leafmap  delete error. no key= "<<LeafMap[key]<<endl;
        }
        else{
            //LeafMap[key] = -2;
            LeafMap.erase(key);
//            cout<<"leafmap  value changed.  key= -2"<<endl;
        }
    }
    
    for(auto &key : keys_remove){
        // Check if key  exists in the map
        if (WholeMap.count(key) <= 0)
        {
//            cout<<"WholeMap delete error. no key= "<<WholeMap[key]<<endl;
        }
        else{
            //WholeMap[key] = -2;
            WholeMap.erase(key);
//            cout<<"WholeMap value changed.  key= -2"<<endl;
        }
    }
    
    
    //second, add node to 2maps
    //remove a node from leaf node, because it is not a leaf anymore
    vector<vector<int> > keys = findks(adds);

    for(int i=0; i<keys.size(); ++i){
        
        //Level

        
        vector<vector<int>> newkey_L;
        vector<int> tmpL;
        for (int j=0; j<get_Dim(); ++j){
            tmpL.push_back(keys[i][j]);
        }
        for (int j=0; j<get_Dim(); ++j){
            newkey_L.push_back(tmpL);
        }
        
        int count = 0;//record the number how many dimition not need to refine
        for (int j=0; j<get_Dim(); ++j){
            if(newkey_L[j][j] < Nmax && dims[i][j] == 1){
                count ++;
                newkey_L[j][j] = newkey_L[j][j] + 1;
            }
            
        }
        
        //Cell
        
        vector<vector<int>> newkey_C;
        vector<int> tmpC;
        for (int j=get_Dim(); j<2*get_Dim(); ++j){
            tmpC.push_back(keys[i][j]);
        }
        for (int j=0; j<2*get_Dim(); ++j){
            newkey_C.push_back(tmpC);
        }
        
        for (int j=0; j<get_Dim(); ++j){
            if(newkey_L[j][j] < Nmax && dims[i][j] == 1){
                newkey_C[2*j][j] = 2*newkey_C[2*j][j];
                if(keys[i][j] == 0){
                    newkey_C[2*j+1][j] = 2*newkey_C[2*j+1][j];
                }
                else{
                    newkey_C[2*j+1][j] = 2*newkey_C[2*j+1][j] + 1;
                }
            }
        }

        //Level+Cell
        vector<vector<int>> newkey_LC;
        for(int j=0; j<newkey_L.size(); ++j){
            vector<int> row_LC;
            row_LC.reserve(newkey_L[0].size() + newkey_C[0].size());
            row_LC.insert( row_LC.end(), newkey_L[j].begin(), newkey_L[j].end() );
            row_LC.insert( row_LC.end(), newkey_C[2*j].begin(), newkey_C[2*j].end() );
            newkey_LC.push_back(row_LC);

            vector<int> row_LC2;
            row_LC2.reserve(newkey_L[0].size() + newkey_C[0].size());
            row_LC2.insert( row_LC2.end(), newkey_L[j].begin(), newkey_L[j].end() );
            row_LC2.insert( row_LC2.end(), newkey_C[2*j+1].begin(), newkey_C[2*j+1].end() );
            newkey_LC.push_back(row_LC2);

        }

        // Level + Cell + Degree
        vector<vector<int> > newkey_LCD;
        newkey_LCD = Compute_LevelCellDegree(newkey_LC);
        

        for(auto &key : newkey_LCD){
            // Check if key  exists in the map
            if (WholeMap.count(key) <= 0)
            {
                WholeMap[key] = total+1;/*WholeMap.size();*/
                total = total + 1 ;
                LeafMap[key] = WholeMap[key];
                //                cout<<"LeafMap  new key= "<<WholeMap[key]<<endl;
                //                cout<<"WholeMap new key= "<<LeafMap[key]<<endl;
                
            }
            else{
                //                cout<<"WholeMap add error. key already exists.  key= "<<WholeMap[key]<<endl;
                //                cout<<"LeafMap  add error. key already exists.  key= "<<LeafMap[key]<<endl;
                
            }
            //search parents of the new node in leaf map
            //move them to the big map
            
/********************************************************
 delete elements that are not leaf any more from leaf map
            vector<vector<int>> parents = findParents(key);
            for(auto &row : parents){
                

                
                if (LeafMap.count(row) > 0)
                {
                    LeafMap.erase(row);
                    //            cout<<"leafmap  delete error. no key= "<<LeafMap[key]<<endl;
                }
//                unordered_map<vector<int>,unsigned long>::const_iterator got = WholeMap.find (row);
//                
//                if ( got != WholeMap.end() ){
//                    LeafMap.erase(got);
//                    cout<<"a parent is found"<<endl;
//                }
            }
************************************/////////////////////////////////////////////
        }
        
        // Check if key  exists in the map
        if (LeafMap.count(keys[i]) > 0  && count == get_Dim()){// all direction has son, then delete
            //if count equal to or more than Dim, donot remove from leaf map
             //           cout<<"LeafMap delete key= "<<LeafMap[keys[i]]<<endl;
            LeafMap.erase(keys[i]);
        }
        else{
            //            cout<<"LeafMap delete error. key not exist. key= "<<LeafMap[keys[i]]<<endl;
        }
        

    }
}

void MeshMaker::refine1(vector<int> adds, vector<int> removes){
    
    //first, remove nodes from 2 maps
    //set value = -2
    //make sure not remove all nodes added last time
    //Otherwise, have to move a node from wholemap to leafmap.a little difficult
    vector<vector<int> > keys_remove = findks(removes);
    
    for(auto &key : keys_remove){
        // Check if key  exists in the map
        if (LeafMap.count(key) <= 0)
        {
            //            cout<<"leafmap  delete error. no key= "<<LeafMap[key]<<endl;
        }
        else{
            LeafMap[key] = -2;
            //            cout<<"leafmap  value changed.  key= -2"<<endl;
        }
    }
    
    for(auto &key : keys_remove){
        // Check if key  exists in the map
        if (WholeMap.count(key) <= 0)
        {
            //            cout<<"WholeMap delete error. no key= "<<WholeMap[key]<<endl;
        }
        else{
            WholeMap[key] = -2;
            //            cout<<"WholeMap value changed.  key= -2"<<endl;
        }
    }
    
    
    //second, add node to 2maps
    //remove a node from leaf node, because it is not a leaf anymore
    vector<vector<int> > keys = findks(adds);
    
    for(int i=0; i<keys.size(); ++i){
        
        //Level
        
        vector<vector<int>> newkey_L;
        vector<int> tmpL;
        for (int j=0; j<get_Dim(); ++j){
            tmpL.push_back(keys[i][j]);
        }
        for (int j=0; j<get_Dim(); ++j){
            newkey_L.push_back(tmpL);
        }
        
        for (int j=0; j<get_Dim(); ++j){
            newkey_L[j][j] = newkey_L[j][j] + 1;
        }
        
        //Cell
        
        vector<vector<int>> newkey_C;
        vector<int> tmpC;
        for (int j=get_Dim(); j<2*get_Dim(); ++j){
            tmpC.push_back(keys[i][j]);
        }
        for (int j=0; j<2*get_Dim(); ++j){
            newkey_C.push_back(tmpC);
        }
        
        for (int j=0; j<get_Dim(); ++j){
            newkey_C[2*j][j] = 2*newkey_C[2*j][j];
            if(keys[i][j] == 0){
                newkey_C[2*j+1][j] = 2*newkey_C[2*j+1][j];
            }
            else{
                newkey_C[2*j+1][j] = 2*newkey_C[2*j+1][j] + 1;
            }
        }
        
        //Level+Cell
        vector<vector<int>> newkey_LC;
        for(int j=0; j<newkey_L.size(); ++j){
            vector<int> row_LC;
            row_LC.reserve(newkey_L[0].size() + newkey_C[0].size());
            row_LC.insert( row_LC.end(), newkey_L[j].begin(), newkey_L[j].end() );
            row_LC.insert( row_LC.end(), newkey_C[2*j].begin(), newkey_C[2*j].end() );
            newkey_LC.push_back(row_LC);
            
            vector<int> row_LC2;
            row_LC2.reserve(newkey_L[0].size() + newkey_C[0].size());
            row_LC2.insert( row_LC2.end(), newkey_L[j].begin(), newkey_L[j].end() );
            row_LC2.insert( row_LC2.end(), newkey_C[2*j+1].begin(), newkey_C[2*j+1].end() );
            newkey_LC.push_back(row_LC2);
            
        }
        
        // Level + Cell + Degree
        vector<vector<int> > newkey_LCD;
        newkey_LCD = Compute_LevelCellDegree(newkey_LC);
        
        
        for(auto &key : newkey_LCD){
            // Check if key  exists in the map
            if (WholeMap.count(key) <= 0)
            {
                WholeMap[key] = WholeMap.size();
                LeafMap[key] = WholeMap[key];
                //                cout<<"LeafMap  new key= "<<WholeMap[key]<<endl;
                //                cout<<"WholeMap new key= "<<LeafMap[key]<<endl;
                
            }
            else{
                //                cout<<"WholeMap add error. key already exists.  key= "<<WholeMap[key]<<endl;
                //                cout<<"LeafMap  add error. key already exists.  key= "<<LeafMap[key]<<endl;
                
            }
            //search parents of the new node in leaf map
            //move them to the big map
            
            /********************************************************
             delete elements that are not leaf any more from leaf map
             vector<vector<int>> parents = findParents(key);
             for(auto &row : parents){
             
             
             
             if (LeafMap.count(row) > 0)
             {
             LeafMap.erase(row);
             //            cout<<"leafmap  delete error. no key= "<<LeafMap[key]<<endl;
             }
             //                unordered_map<vector<int>,unsigned long>::const_iterator got = WholeMap.find (row);
             //
             //                if ( got != WholeMap.end() ){
             //                    LeafMap.erase(got);
             //                    cout<<"a parent is found"<<endl;
             //                }
             }
             ************************************/////////////////////////////////////////////
        }
        
        // Check if key  exists in the map
        if (LeafMap.count(keys[i]) > 0){
            //            cout<<"LeafMap delete key= "<<LeafMap[keys[i]]<<endl;
            LeafMap.erase(keys[i]);
        }
        else{
            //            cout<<"LeafMap delete error. key not exist. key= "<<LeafMap[keys[i]]<<endl;
        }
        
        
    }
}

//compute level
//input:Dim,Np
//combinations of k numbers sums less than n
vector<vector<int> > MeshMaker::Compute_Level(int k,//Dim
                                              int n,//Np
                                              int isLeaf) {//if it it a leaf element
    vector<vector<int> > result;
    vector<int> sol;
    if (isLeaf == 0) {
        CombinationSumN0(k, n, result, sol);//sum == n
    }
    else if (isLeaf == 1){
        CombinationSumN1(k, n, result, sol);//sum <= n
    }
    else
        cout<<"error"<<endl;
    return result;
}

//compute cell from level
vector<vector<int> > MeshMaker::Compute_LevelCell(vector<vector<int> > keys_L){
    vector<vector<int> > keys_LC;
    
    for (const auto &key_L :keys_L) {
        //for every level, for every number, generate possible values
        //from 0 to pow(2,num-1)
        vector<vector<int> > key_LC;
        vector<int> com;
        vector<vector<int> > coms;
        for(const auto &num : key_L){//for every number
            //int tmp = ((num==0) ? 0 : pow(2,num-1)-1);
            com.resize(pow(2,max(0,num-1)));
            iota (begin(com), end(com), 0);
            coms.push_back(com);//add to coms
        }
        
        
        vector<vector<int> > key_C = combinationCell(coms);//cells for one level
        
        for(const auto &row : key_C){
            //one level and its cells
            vector<int> row_LC;
            row_LC.reserve(row.size() + key_L.size());
            row_LC.insert( row_LC.end(), key_L.begin(),key_L.end() );
            row_LC.insert( row_LC.end(), row.begin(), row.end());
            key_LC.push_back(row_LC);
        }
        
        //merge level and cells
        keys_LC.insert(keys_LC.end(), key_LC.begin(), key_LC.end());
    }
    
    return keys_LC;
}

//compute degrees from level&cell
vector<vector<int> > MeshMaker::Compute_LevelCellDegree(vector<vector<int> > keys_LC){
    
    vector<vector<int> > keys_LCD;
    
    for (const auto &key_LC :keys_LC) {//for every key_Level
        vector<vector<int> > key_LCD;
        vector<int> com;
        vector<vector<int> > coms;
        for (int i=0; i<Dim; ++i){
            com.resize(this->D + 1);
            iota (begin(com), end(com), 0);
            coms.push_back(com);
        }
        vector<vector<int> > key_D_full = combinationCell(coms);
        
        //delete sum k more than Dm
        vector<vector<int> > key_D;
        for(const auto & row : key_D_full){
            int sum_of_elems = 0;
            for (auto& n : row){
                sum_of_elems += n;
            }
            if(sum_of_elems <= this->Dm){
                key_D.push_back(row);
            }
        }
        
        for(const auto &row : key_D){
            vector<int> row_LCD;
            row_LCD.reserve(row.size() + key_LC.size());
            row_LCD.insert( row_LCD.end(), key_LC.begin(),key_LC.end() );
            row_LCD.insert( row_LCD.end(), row.begin(), row.end());
            key_LCD.push_back(row_LCD);
        }
        
        keys_LCD.insert(keys_LCD.end(), key_LCD.begin(), key_LCD.end());
    }
    return keys_LCD;
    
}

//k numbers sum to s     (s<=n)
//k:from 0 to n
void MeshMaker::CombinationSumN0( int k,
                                 int n,
                                 vector<vector<int> > &res,
                                 vector<int> &sol) {
    if( k<0 || n<0 ){
        return;
    }
    else if(k==0 && n>=0){
        res.push_back(sol);
        return;
    }
    else{
        for(int i=0; i<=n; ++i){
            sol.push_back(i);
            CombinationSumN0(k-1, n-i, res, sol);
            sol.pop_back();
        }
    }
}

//k numbers sum to n
//k:from 1 to n
void MeshMaker::CombinationSumN1( int k,
                                 int n,
                                 vector<vector<int> > &result,
                                 vector<int> &sol) {
    if( k<0 || n<0 ){
        return;
    }
    else if(k==0 && n==0){
        result.push_back(sol);
        return;
    }
    else{
        for(int i=0; i<=n; ++i){//if i from 1, there is no level 0 in  all leaf map level
            sol.push_back(i);
            CombinationSumN1(k-1, n-i, result, sol);
            sol.pop_back();
        }
    }
}



//for one level, compute combinations of cells
vector<vector<int> > MeshMaker::combinationCell(const vector<vector<int> > keys_Level) {
    vector<vector<int> > keys_Cell;
    vector<int>  sol;
    cart_product(keys_Cell, sol, keys_Level.begin(), keys_Level.end());
    return keys_Cell;
}

//Descartes product of multi vectors
void MeshMaker::cart_product(
                             vector<vector<int> >& rvvi,  // final result
                             vector<int>&  rvi,   // current result
                             vector<vector<int> >::const_iterator me, // current input
                             vector<vector<int> >::const_iterator end) // final input
{
    if(me == end) {
        // terminal condition of the recursion. We no longer have
        // any input vectors to manipulate. Add the current result (rvi)
        // to the total set of results (rvvvi).
        rvvi.push_back(rvi);
        return;
    }
    
    // need an easy name for my vector-of-ints
    const vector<int>& mevi = *me;
    for(vector<int>::const_iterator it = mevi.begin();
        it != mevi.end();
        it++) {
        // final rvi will look like "a, b, c, ME, d, e, f"
        // At the moment, rvi already has "a, b, c"
        rvi.push_back(*it);  // add ME
        cart_product(rvvi, rvi, me+1, end); //add "d, e, f"
        rvi.pop_back(); // clean ME off for next round
    }
}

void MeshMaker::PrintKeys(){
    cout<<"Dim= "<<Dim<<", "<<endl;
    cout<<"Np = "<<Np<<", "<<endl;
    cout<<"D  = "<<D<<", "<<endl;
    cout<<"Dm = "<<Dm<<", "<<endl;
    
    cout<<"keys"<<endl;
    for ( const auto &row : keys_LCD )
    {
        for ( const auto &s : row )
            cout << s << ' ';
        cout << endl;
    }
}

void MeshMaker::PrintWholeMap(){
    
    //sort by key
    vector<pair<vector<int>, int> > elems(WholeMap.begin(), WholeMap.end());
    sort(elems.begin(), elems.end());
    
    ofstream out;
    out.open("big map.dat");
    
    
//    out<<"-------------------------------"<<endl;
    
//    out<<"------Whole Map-------"<<endl;
    for(const auto &row:elems){
        vector<int> keys = row.first;
        
        for(auto &key : keys)
        {
            out << " " << key;
        }
        
        out << "  "<< row.second << endl;
    }
//    out<<"-----WholeMap: size:   "<<WholeMap.size()<<endl;
//    out<<"-------------------------------"<<endl;
    out<<endl<<endl;
    out.close();
}
void MeshMaker::PrintLeafMap(){
    
    //sort by key
    vector<pair<vector<int>, int> > elems(LeafMap.begin(), LeafMap.end());
    sort(elems.begin(), elems.end());
    
    //    //sort by value
    //    vector<pair<vector<int>, int>> elems(map.begin(), map.end());
    //    sort(elems.begin(), elems.end(), comp);
    //
    ofstream out;
    out.open("small map.dat");
    
    
    //    out<<"-------------------------------"<<endl;
    
    //    out<<"------Whole Map-------"<<endl;
    for(const auto &row:elems){
        vector<int> keys = row.first;
        
        for(auto &key : keys)
        {
            out << " " << key;
        }
        
        out << "  "<< row.second << endl;
    }
    //    out<<"-----WholeMap: size:   "<<WholeMap.size()<<endl;
    //    out<<"-------------------------------"<<endl;
    out<<endl<<endl;
    out.close();
}

//find keys in leafmap
vector<vector<int> > MeshMaker:: findks(vector<int> values){
    //find elements which need to do refinement
    //the elements are in leaf map now
    
    vector<vector<int> > keys;
    
    for (int i=0; i<values.size(); ++i) {
        unordered_map<vector<int>, unsigned long>::const_iterator it;
        for (it = LeafMap.begin(); it != LeafMap.end(); ++it) {
            if(it->second == values[i]){
                keys.push_back(it->first);
            }
        }
    }
    return keys;
}
//find key.only support find a key once
vector<int>  MeshMaker:: findk(int value){
    vector<int>  key;
    unordered_map<vector<int>, unsigned long>::const_iterator it;
    for (it = LeafMap.begin(); it != LeafMap.end(); ++it) {
        if (it->second == value) {
            key = it->first;
        }
        else{
//            cout<<"LeafMap value not found"<<value<<endl;
        }
    }
    return key;
}
int MeshMaker:: get_bigmapsize(){
    return WholeMap.size();
}

vector<vector<int>> MeshMaker:: findParents(vector<int> key){
    
    //Level
    
    vector<vector<int>> newkey_L;
    vector<int> tmpL;
    for (int j=0; j<get_Dim(); ++j){
        tmpL.push_back(key[j]);
    }
    for (int j=0; j<get_Dim(); ++j){
        newkey_L.push_back(tmpL);
    }
    
    for (int j=0; j<get_Dim(); ++j){
        if(newkey_L[j][j] != 0){
            newkey_L[j][j] = newkey_L[j][j] - 1;
        }
    }
    
    //Cell
    
    vector<vector<int>> newkey_C;
    vector<int> tmpC;
    for (int j=get_Dim(); j<2*get_Dim(); ++j){
        tmpC.push_back(key[j]);
    }
    for (int j=0; j<get_Dim(); ++j){
        newkey_C.push_back(tmpC);
    }
    
    for (int j=0; j<get_Dim(); ++j){
        newkey_C[j][j] = int(newkey_C[j][j] / 2);
    }

    
    //Level+Cell
    vector<vector<int>> newkey_LC;
    for(int j=0; j<newkey_L.size(); ++j){
        vector<int> row_LC;
        row_LC.reserve(newkey_L[0].size() + newkey_C[0].size());
        row_LC.insert( row_LC.end(), newkey_L[j].begin(), newkey_L[j].end() );
        row_LC.insert( row_LC.end(), newkey_C[j].begin(), newkey_C[j].end() );
        newkey_LC.push_back(row_LC);
        
    }
    
    // Level + Cell + Degree
    vector<vector<int> > newkey_LCD;
    newkey_LCD = Compute_LevelCellDegree(newkey_LC);
    
      return newkey_LCD;
}
