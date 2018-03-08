#include "mex.h"
#include "class_handle.hpp"
#include "MeshMaker.hpp"
using namespace mexutils;
// A set of utility functions are provided in class_handle.hpp
// in the namespace mexutils. These can be to convert between
// some matlab and std data types, and ease various tasks for
// the mex interface


// wrapper class to convert matlab arguments to appropriate arguments
// for wrapped class
class MeshMaker_wrapper
{
public:
    void findvs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        std::vector<int> nallowed;
        std::vector<mwSize> index;
        index.push_back(0);
        index.push_back(0);
        
        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        mxNumericArrayWrapper keys = mxnthargmatrix (nrhs, prhs, 1, 2);
        
        mwSize nrows = keys.getRows ();
        mwSize ncols = keys.getColumns ();
        
        vector<vector<int >> vkeys;
        for(mwSize i=0; i<nrows; ++i){
            vector<int> vkey;
            index[0] = i;
            for(mwSize j = 0; j<ncols; ++j){
                index[1] = j;
                vkey.push_back(keys.getDoubleValue(index));
            }
            vkeys.push_back(vkey);
        }
        vector<int> vvalues = mesh.findvs(vkeys);
   
        mxSetLHS (vvalues, 1, nlhs, plhs);
        
    }
   
    void findk(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {

        std::vector<int> nallowed;
        
        // one arg, the node number
        nallowed.push_back (1);
        
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        int value = int (mxnthargscalar (nrhs, prhs, 1, 2));
        
        vector<int> key = mesh.findk(value);
        
        if(key.size() == 0){
            mxSetLHS (-1, 1, nlhs, plhs);
            return;
        }
        else{
            mxSetLHS (key, 1, nlhs, plhs);
        }
        

/*find many keys from many values
 //////////////////////
        std::vector<int> nallowed;
        std::vector<mwSize> index;
        index.push_back(0);
        index.push_back(0);
        
        nallowed.push_back (1);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        mxNumericArrayWrapper values = mxnthargmatrix (nrhs, prhs, 1, 2);
        
        mwSize nrows = values.getRows ();
        mwSize ncols = values.getColumns ();
        
        vector<int > vvalues;
        
        index[0] = 0;
        for (mwSize j = 0; j<ncols; ++j){
            index[1] = j;
            vvalues.push_back(values.getDoubleValue(index));
        }
        
        vector<vector<int>> vkeys = mesh.findk(vvalues);
        
        nrows = vkeys.size();
        ncols = vkeys[0].size();
        
        plhs[0] = mxCreateDoubleMatrix(nrows,ncols,mxREAL);
        
        mxNumericArrayWrapper out = mxNumericArrayWrapper (plhs[0]);
        
        index[0] = 0;
        index[1] = 0;
        for(mwSize i= 0; i<nrows; ++i){
            index[0] = i;
            for (mwSize j = 0; j<ncols; ++j) {
                index[1] = j;
                out.setDoubleValue(index, vkeys[i][j]);
            }
        }
////////////////////////////////////////////////////*/
    }
    
    void refine(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
      
        std::vector<int> nallowed;
        std::vector<mwSize> index;
        index.push_back(0);
        index.push_back(0);
        mwSize nrows;
        mwSize ncols;
        
        nallowed.push_back (4);
        int nargin = mxnarginchk (nrhs, nallowed, 2);
        
        
        mxNumericArrayWrapper adds = mxnthargmatrix (nrhs, prhs, 1, 2);

        nrows = adds.getRows ();
        ncols = adds.getColumns ();
        
        vector<int > vvalues_add;
        
        index[0] = 0;
        for (mwSize j = 0; j<ncols; ++j){
            index[1] = j;
            vvalues_add.push_back(adds.getDoubleValue(index));
        }
        
        mxNumericArrayWrapper dims = mxnthargmatrix (nrhs, prhs, 2, 2);
        nrows = dims.getRows ();
        ncols = dims.getColumns ();
        
        vector<vector<int >> vvdims;
        
        index[0] = 0;
        index[1] = 0;
        for (mwSize i = 0; i<nrows; ++i){
            index[0] = i;
            vector<int> tmp;
            tmp.clear();
            for (mwSize j=0; j<ncols; ++j){
                index[1] = j;
                tmp.push_back(dims.getDoubleValue(index));
            }
            vvdims.push_back(tmp);
        }
  //      for(int i=0;i<vvdims.size(); ++i){
  //          for(int j=0; j<ncols; ++j){
  //              cout<<vvdims[i][j]<<" ";
  //          }
  //          cout<<endl;
  //      }
        
        
        mxNumericArrayWrapper removes = mxnthargmatrix (nrhs, prhs, 3, 2);
        nrows = removes.getRows ();
        ncols = removes.getColumns ();
        
        vector<int > vvalues_rem;
        
        index[0] = 0;
        for (mwSize j = 0; j<ncols; ++j){
            index[1] = j;
            vvalues_rem.push_back(removes.getDoubleValue(index));
        }
        
        int Nmax = int (mxnthargscalar (nrhs, prhs, 4, 2));
        
        mesh.refine(vvalues_add, vvdims, vvalues_rem, Nmax);
    }
    void init(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        int Dim,Np,D,Dm;
        Dim = mxGetScalar(prhs[2]);
        Np = mxGetScalar(prhs[3]);
        D = mxGetScalar(prhs[4]);
        Dm = mxGetScalar(prhs[5]);
        mesh.init(Dim,Np,D,Dm);
    }
    void print(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        mesh.PrintWholeMap();
        mesh.PrintLeafMap();
    }
    void allv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        vector<int> vvalues;
        for(const auto & row:mesh.get_WholeMap()){
            vvalues.push_back(row.second);
        }
        mxSetLHS (vvalues, 1, nlhs, plhs);
    }
    void leafv(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        vector<int> vvalues;
        for(const auto & row:mesh.get_LeafMap()){
            vvalues.push_back(row.second);
        }
        mxSetLHS (vvalues, 1, nlhs, plhs);
    }
    void size(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
    {
        mxSetLHS (mesh.get_bigmapsize(), 1, nlhs, plhs);
    }
private:
    
    // instance of the wrapped c++ class
    MeshMaker mesh;
    
};

// mexfunction defintion, all the interction with the class is done through
// this function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // use the macros provided by class_handle.hpp to create the interface
    // these macros assume the wrapper class can be constructed without arguments.
    // Note that the instance of the wrapper class is declared with 'new' and
    // created on the heap
    BEGIN_MEX_CLASS_WRAPPER(MeshMaker_wrapper)
    REGISTER_CLASS_METHOD(MeshMaker_wrapper,findvs)
    REGISTER_CLASS_METHOD(MeshMaker_wrapper,findk)
    REGISTER_CLASS_METHOD(MeshMaker_wrapper,refine)
    REGISTER_CLASS_METHOD(MeshMaker_wrapper,init)
    REGISTER_CLASS_METHOD(MeshMaker_wrapper,print)
    REGISTER_CLASS_METHOD(MeshMaker_wrapper,leafv)
    REGISTER_CLASS_METHOD(MeshMaker_wrapper,allv)
    REGISTER_CLASS_METHOD(MeshMaker_wrapper,size)

    END_MEX_CLASS_WRAPPER(MeshMaker_wrapper)
}



