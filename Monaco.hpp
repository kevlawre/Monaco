#ifndef MONACO_H
#define MONACO_H
#include "PCA.hpp"
#include <vector> 

using namespace std;

struct data_fingerprint{
    vector<int> kg; 
    vector<int> ran; 
    vector<int> mcm; 
    vector<int> mca; 
    vector<int> map; 
    
    data_fingerprint(int order){
        vector<int> z(order, 0); 
        kg = z; 
        ran = z; 
        mcm = z; 
        mca = z; 
        map = z; 
    }
};

void Monaco(component_data cd);
#endif // MONACO
