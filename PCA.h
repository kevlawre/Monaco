#ifndef PCA_H
#define PCA_H
#include <vector>

using namespace std;

struct eigen_unit{
    vector<double> eigenvector; 
    vector<double> data; 
    double eigenvalue;
    eigen_unit(vector<double> ev, double eval){ 
        eigenvector = ev; 
        eigenvalue = eval; 
    }
    struct by_eigenvalue { 
        bool operator()(eigen_unit const &a, eigen_unit const &b){ 
            return a.eigenvalue < b.eigenvalue;
        }
    };
};

struct component_data {
    vector<eigen_unit> comp; //components
    vector<vector<double> > original_data;   //original data 
};

struct pivot{
    int I; //Pivot column
    int J; //Pivot row
};

struct jacobi_transform{
    vector<vector<double> > orig_data; // Saved copy of original data
    vector<vector<double> > eigenvectors; //OUTPUT Eigenvectors
    vector<vector<double> > covariance_matrix; //Working Covariance Matrix
    vector<double> eigenvalues; //OUTPUT Eigenvalues
    double thresh; //Working Threshold
    int it_count; //Working Iteration Count
    int it_max; //INPUT Const M
    vector<double> n_ev; //Transfer vector (donor)
    vector<double> o_ev; //Intermediary vector (acceptor)
    pivot p; //Working pivot
    
};

bool check_convergence(jacobi_transform &jf);
void jacobi_rotation(jacobi_transform &jf);
void eigenvector_init(jacobi_transform &jf);
jacobi_transform eigenvalue_calc(vector<vector<double> > c, int MAX_ITERATIONS);


#endif // PCA
