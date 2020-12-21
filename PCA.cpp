#include <vector> 
#include <functional> 
#include <algorithm> 
#include <math.h> 
#include <fstream>
#include <sstream> 
#include <string> 
#include <numeric>
#include <iostream>
#include <iomanip>
#include "PCA.h"
#include "../lib/error_logs.h"
#include "Monaco.h"
#include "../lib/statistics.hpp"

using namespace std; 

bool check_convergence(jacobi_transform &jf){ 
    double order = jf.covariance_matrix.size();  //Order of covariance matrix
    double thold = MODULI(jf.covariance_matrix);  //Modulus of covariance matrix
    thold = 0.2 * pow(thold/order, 2);  //Rutishauser's threshold formula
    jf.thresh = thold;  //Set threshold. 
    return thold == 0.0 ? false : true; //Return boolean. 
}

void init_identity(jacobi_transform &jf){
    int order = jf.covariance_matrix.size(); //Order of covariance matrix
    vector<vector<double> > seigen = IDENTITY_MATRIX(order); //Create identity matrix
    jf.eigenvectors = seigen; //Set identity matrix
    return; 
}

void load_components(jacobi_transform jf, vector<eigen_unit> &e_units){
    int order = jf.eigenvectors.size(); //Order of eigenvector matrix
    for(unsigned int i = 0; i < order; i++){
        eigen_unit ev(jf.eigenvectors[i], jf.eigenvalues[i]); //Zip Eigenvalue and Eigenvector
        e_units.push_back(ev); //accumulate zipped unit in reference vector
    } 
    return; 
}

bool pivot_check(jacobi_transform &jf, double &s){
    s = 100.0 * fabs(jf.covariance_matrix[jf.p.I][jf.p.J]);//epsilon value
    double epi = s + fabs(jf.eigenvalues[jf.p.I]); //i + epsilon
    double epj = s + fabs(jf.eigenvalues[jf.p.J]);//j + epsilon
    
    if(jf.it_count > 3 && epi == fabs(jf.eigenvalues[jf.p.I]) && epj == fabs(jf.eigenvalues[jf.p.J])){ 
        jf.covariance_matrix[jf.p.I][jf.p.J] = 0.0; //Annihilate tiny off-diagonal elements
        return false; 
    } 
    return true; 
}

void jacobi_rotation(jacobi_transform &jf){
    int order = jf.covariance_matrix.size(); 
    
    for(int i = 0; i < order; i++){
        for(int j = (i+1); j < order; j++){
            double epsilon = 0; //Declare an epsilon value
            jf.p.I = i; // set column of current sweep
            jf.p.J = j;  // set row of current sweep
            
            if (pivot_check(jf, epsilon) == true){
                 
                double delta = jf.eigenvalues[j] - jf.eigenvalues[i]; //difference of diagonal elements.
                double theta = 0.5 * (delta) / jf.covariance_matrix[i][j]; //angle of rotation 
                
                double tangent = ((fabs(delta) + epsilon) == fabs(delta)) ?
                    jf.covariance_matrix[i][j]/delta : theta < 0 ?
                    -1.0/(fabs(theta) + sqrt(1 + pow(theta, 2))) :
                    1.0/(fabs(theta) + sqrt(1 + pow(theta, 2))); //Tangent of theta
                
                double cosine = 1.0 / sqrt ( 1.0 + pow(tangent, 2)); //cosine of theta
                double sine = tangent * cosine; // sine of theta
                
                double tau = sine / ( 1.0 + cosine ); //tan of (theta/2)
                delta = tangent * jf.covariance_matrix[i][j]; //update delta
                
                //accumulate changes to diagonal elements
                jf.n_ev[i] = (double) jf.n_ev[i] - delta; 
                jf.n_ev[j] = (double) jf.n_ev[j] + delta; 
                jf.eigenvalues[i] = (double) jf.eigenvalues[i] - delta;
                jf.eigenvalues[j] = (double) jf.eigenvalues[j] + delta; 
                
                jf.covariance_matrix[i][j] = 0.0; //anihilate pivot from covariance matrix. 
                
                //perform the rotation
                for (int k = 0; k < i; k++){
                    double g = jf.covariance_matrix[k][i];
                    double h = jf.covariance_matrix[k][j];
                    jf.covariance_matrix[k][i] = g - sine * (h + g * tau);
                    jf.covariance_matrix[k][j] = h + sine * (g - h * tau);
                }
                for (int k = i + 1; k < j; k++){
                    double g = jf.covariance_matrix[i][k];
                    double h = jf.covariance_matrix[k][j];
                    jf.covariance_matrix[i][k] = g - sine * (h + g * tau);
                    jf.covariance_matrix[k][j] = h + sine * (g - h * tau);
                }
                for (int k = j + 1; k < order; k++){
                    double g = jf.covariance_matrix[i][k];
                    double h = jf.covariance_matrix[j][k];
                    jf.covariance_matrix[i][k] = g - sine * (h + g * tau);
                    jf.covariance_matrix[j][k] = h + sine * (g - h * tau);
                }
                //accumulate information in eigenvectors matrix. 
                for (int k = 0; k < order; k++){
                    double g = jf.eigenvectors[k][i];
                    double h = jf.eigenvectors[k][j];
                    jf.eigenvectors[k][i] = g + sine * (h - g * tau);
                    jf.eigenvectors[k][j] = h - sine * (g + h * tau);
                }
            }  
        }  
    }
    for (int i = 0; i < order; i++ ){
        jf.o_ev[i] = jf.o_ev[i] + jf.n_ev[i]; //aggregate the corrections
        jf.eigenvalues[i] = jf.o_ev[i]; //transfer the corrections to eigenvalues
        jf.n_ev[i] = 0; //reset the vector. 
    }
}



jacobi_transform eigenvalue_calc(vector<vector<double> > c, int MAX_ITERATIONS){
    int ITERATION_COUNT = 0;//Counter for iterations 
    int order = c[0].size(); //Number of dimensions in the dataset
    vector<double> nev(order, 0.0); // initialize an empty vector to hold corrections
    
    jacobi_transform jf; //create a jacobi_transformation structure
    jf.it_max = MAX_ITERATIONS; //Set max value for iterations. 
    jf.orig_data = c; //save a copy of original data
    jf.covariance_matrix = COVARIANCE_MATRIX(c); //Calculate covariance matrix. 
    init_identity(jf);//initialize identity matrix. 
    jf.eigenvalues = DIAGONAL(jf.covariance_matrix); //set initial eigenvalues
    jf.o_ev = jf.eigenvalues; //keep a copy of initial eigenvalues
    jf.n_ev = nev; //set the correction vector. 
    
    do{
        ITERATION_COUNT++;//Count the iterations
        jf.it_count = ITERATION_COUNT; //set the iteration count
        if(check_convergence(jf) == false){break;} //Check for convergence
        if(ITERATION_COUNT > 3){jf.thresh = 0;} //Check threshold
        jacobi_rotation(jf); //perform the rotation
        
    }while(ITERATION_COUNT < MAX_ITERATIONS); 
    
    return jf; 
}

component_data PCA(vector<vector<double> > Features){
    int max_it = Features.size()*Features[0].size(); 
    jacobi_transform jf = eigenvalue_calc(Features, max_it); 
    int order = jf.covariance_matrix[0].size();
    vector<eigen_unit> components; 
    load_components(jf, components); 
    EigenData_Dump(jf, 1);
    cout<<"MONACO RESULTS. . ."<<endl;
    component_data cd; 
    cd.comp = components; 
    cd.original_data = jf.orig_data;
    Monaco(cd); 
    return cd; 
}


