#include <vector> 
#include <functional> 
#include <algorithm> 
#include <math.h> 
#include <fstream>
#include <sstream> 
#include <string> 
#include <time.h>
#include <numeric>
#include <iostream>
#include <limits> //DBL_MAX; 
#include <iomanip>
#include "PCA.h"
#include "../lib/error_logs.h"
#include "Monaco.h"
#include "../lib/statistics.hpp"

using namespace std; 
const int MAX_SIM_COUNT = 999; 

void Kaiser_Guttman(component_data cd, vector<int> &c){
    sort(cd.comp.begin(), cd.comp.end(), eigen_unit::by_eigenvalue()); //Sort by Eigenvalue
    int csize = cd.comp.size(); //Order of component matrix. 
    
    for(int i = 0; i < csize; i++){
        if(cd.comp[i].eigenvalue>=1.0){ //"If greater than most"
             c[i] = 1; //Reflect in referenced vector. 
        }
    }
    return;   
}

void Anderson_Braak(component_data cd, vector<int> &c){
    int order = cd.comp.size(); //Order of component matrix. 
    int SIM_COUNT = 0; //Simulation counter
    vector<int> randoms(order, 1);//Initial vector to hold randoms 
    //Initialize with 1's so to account for observed Eigenvalue
    int Y = cd.original_data.size(); //Original data size. 
    cd.original_data = TRANSPOSE(cd.original_data); //Transpose original data. 

    while(SIM_COUNT < MAX_SIM_COUNT){ //MAX_SIM_COUNT == 999
        
        for(int i = 0; i < order; i++){
            //srand(time(0)); in Monaco Function. 
            random_shuffle(cd.original_data[i].begin(), cd.original_data[i].end()); 
        }
        vector<vector<double> > sim = TRANSPOSE(cd.original_data);
        //Transpose the original data back. 
        jacobi_transform jt = eigenvalue_calc(sim, Y*order);  //Calculate eigenvalues. 
        
        for(int i = 0; i < order; i++){
            
            if(jt.eigenvalues[i] > cd.comp[i].eigenvalue){ //if calculated > observed
                randoms[i] = randoms[i]+1; //Iterate the count of randoms
            }
        }
        SIM_COUNT++; //Iterate the simulation count. 
    }
    for(int i = 0; i < order; i++){
        if(randoms[i]<=50){ //Estimated P-Value limit. 
            c[i] = 1; //Reflect in referenced vector. 
        }
    }
    return; 
}

void Monte_Carlo(component_data cd, vector<int> &c, int flag){
    sort(cd.comp.begin(), cd.comp.end(), eigen_unit::by_eigenvalue()); //Sort eigenvalues
    int order = cd.comp.size(); //Order of component matrix.
    int Y = cd.original_data.size(); //Original data size. 
    int SIM_COUNT = 0; //Simulation counter
    vector<vector<double> > simulations; //Vector to hold results from simulations
    
    while(SIM_COUNT < MAX_SIM_COUNT){ //MAX_SIM_COUNT == 999
        vector<vector<double> > sim; // 2D Vector to hold simulation 
        
        for(int i = 0; i < Y; i++){
            vector<double> simrow; //Row of sim 
            
            for(int j = 0; j < order; j++){
                double r = rand(); //srand(time(null)) in Monaco Function
                simrow.push_back(r); //Aggregate random values. 
            }
            sim.push_back(simrow); //Aggregate simulation rows. 
        }
        jacobi_transform jt = eigenvalue_calc(sim ,order*Y); //Calculate eigenvalues for simulation
        sort(jt.eigenvalues.begin(), jt.eigenvalues.end()); //Sort simulated eigenvalues. 
        simulations.push_back(jt.eigenvalues); //Accumulate data. 
        SIM_COUNT++; //Iterate SIM_COUNT
    }
    simulations = TRANSPOSE(simulations); //Transpose simulations. 
    
    for(int i = 0; i < order; i++){
        double av = MEAN(simulations[i]); //Mean of simulations. 
        double sd = STANDARD_DEVIATION(simulations[i]); //Standard deviations of simulations
        
        if(flag == 0 && cd.comp[i].eigenvalue > (av+1.96*sd/sqrt(Y))){ //Calculate confidence interval. 
            c[i] = 1; //Reflect in referenced vector. 
        }
        else if(flag == 1 && cd.comp[i].eigenvalue > av){ //For Average Monte-Carlo
            c[i] = 1;  //Reflect in referenced vector. 
        }
    }
    return; 
}


void Minimum_Average_Partial(component_data cd, vector<int> &c){ 
    double MAP = numeric_limits<double>::max(); //Maximum numerical limit of double 
    int MAP_index = 0; //MAP index tracker. 
    cd.original_data = TRANSPOSE(cd.original_data); //Transpose original data. 
    int order = cd.original_data.size(); //Order of component matrix. 
    
    for(int j = 0; j < order; j++){
            cd.comp[j].data = cd.original_data[j]; //Set data properties. 
    }
    sort(cd.comp.begin(), cd.comp.end(), eigen_unit::by_eigenvalue()); //Sort by eigenvalue
    vector<vector<double> > od; //Original "sorted" data vector.  
    for(int j = 0; j < order; j++){
        od.push_back(cd.comp[j].data); //Fill in with original data in sorted form. 
    }
    od = TRANSPOSE(od);  //Transpose sorted original data.
    vector<vector<double> > cm = COVARIANCE_MATRIX(od); //Calculate covariance matrix.
    INVERSE(cm); //Invert covariance matrix. 

    for(int i = 0; i < order-2; i++){
        vector<vector<double> > crmn; //Partial covariance matrix. 
        int diff = order - (i+1); //Number of partialed out components. 
        for(int j = 0; j < diff; j++){
            vector<double> crmnr; //Partial Correlation Row. 
            for(int k = 0; k < diff; k++){
                double rjk = cm[j][k]; //Inverted position j,k
                double rjj =  cm[j][j]; //Inverted position j, j
                double rkk = cm[k][k]; //Inverted position k, ,
                double r = -1*rjk/sqrt(rjj*rkk); //Partial Correlation Coefficient
                crmnr.push_back(r); //Aggregate Partial Correlations
            }
            crmn.push_back(crmnr); //Aggregate Partial Correlation Rows. 
        }
        double porder = crmn.size(); //Order of Partial Correlation Matrix
        vector<double> ap; //Average partial 
        double counter = 0.0; //Counter. 
        for(int j = 0; j < porder; j++){
            for(int k = 0; k < porder; k++){
                if(j == k){continue;} //Skip diagonal values. 
                ap.push_back(pow(crmn[j][k], 2)); //Aggregate off-diagonal values. 
                counter++; //Iterate Counter. 
            } 
        }
        double dap = accumulate(ap.begin(), ap.end(), 0.0)/counter; //Calculate Average. 
        if(dap < MAP){ //Check if minimum. 
            MAP = dap; //Set minimum value. 
            MAP_index = i; //Set minimum index. 
        }
    }
    for(int i = 0; i < MAP_index; i++){
        c[i] = 1; //Reflect in reference vector. 
    }
    return; 
}

void Monaco(component_data cd){
    srand(time(0)); //Seed Random with time. 
    int order = cd.comp.size(); 
    data_fingerprint dt(order);  
    Kaiser_Guttman(cd, dt.kg); 
    cout<<"..............................................................................."<<endl; 
    cout<<"KAISER-GUTTMAN (Industry standard): "<<accumulate(dt.kg.begin(), dt.kg.end(), 0.0)<<endl; 
    Anderson_Braak(cd, dt.ran); 
    cout<<"ANDERSON-BRAAK PERMUTATION TEST (P-Score via shuffled simulations): "<<accumulate(dt.ran.begin(), dt.ran.end(), 0.0)<<endl;
    Minimum_Average_Partial(cd, dt.map);
    cout<<"VELICER'S MINIMUM AVERAGE PARTIAL (Optimized for minimal residual correlation): "<<accumulate(dt.map.begin(), dt.map.end(), 0.0)<<endl;
    Monte_Carlo(cd, dt.mcm, 0); 
    cout<<"HORN'S PARALLEL ANALYSIS (Confidence intervals via Monte-Carlo Simulations): "<<accumulate(dt.mcm.begin(), dt.mcm.end(), 0.0)<<endl;
    Monte_Carlo(cd, dt.mca, 1); 
    cout<<"AVERAGE MONTE-CARLO (Average comparision via Monte-Carlo Simulations): "<<accumulate(dt.mca.begin(), dt.mca.end(), 0.0)<<endl; 
    
}










