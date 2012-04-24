/* runTPM.cpp: 

    This file contains a main function to run TPMs. 

*/ 

// Imports: 
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib> 
#include <vector>
#include "helper.h"
using namespace std; 

int main(int argc, char* argv[]){
    /* main for runTPM.cpp */

    int seed=0; 
    int K_=1, N_=10, L_=5; 
    int **x, **w, *h, tau; 
    int updateRule_=1; 

    int nSample_=10000, nReset_=100000; 

    for (int i=0; i<argc; i++){
        if (string(argv[i])=="-uR"){
            updateRule_ = atoi(argv[i+1]); 
        }
        if (string(argv[i])=="-s"){
            seed = atoi(argv[i+1]); 
        }
    }
    //     if (string(argv[i])=="-sC"){
    //         stopCount = atoi(argv[i+1]); 
    //     }
    //     if (string(argv[i])=="-mC"){
    //         maxCount = atoi(argv[i+1]); 
    //     }
    //     if (string(argv[i])=="-sC_"){
    //         stepCount = atoi(argv[i+1]); 
    //     }
    //     if (string(argv[i])=="-tC"){
    //         topCount = atoi(argv[i+1]); 
    //     }
    //     if (string(argv[i])=="-sC__"){
    //         saveCount = atoi(argv[i+1]); 
    //     }
    // }


    // Defining a TPM: 
    // TPM::PTPM M1, M2; 
    srand(seed); 
    TPM::pTPM M1(nSample_, nReset_, K_, N_, L_, updateRule_); 
    TPM::TPM M2(K_, N_, L_, updateRule_); 

    x = new int*[K_]; 
    float ***pW_; 
    pW_ = new float**[K_]; 
    for (int i=0; i<K_; i++){
        x[i] = new int[N_];
        pW_[i]  = new float*[N_]; 
        for (int j=0; j<N_; j++){
            pW_[i][j] = new float[2*L_+1]; 
            x[i][j] = 1; 
            for (int k=0; k<2*L_+1; k++){
                pW_[i][j][k] = 1.; 
            }
        }
    }

    for (int i=0; i<400; i++){
        x = h::generateX(K_, N_); 
        M1.setX(x); 
        M2.setX(x); 
        M2.compute(); 

        M1.setTargetTau(M2.getTau()); 
        M1.MC_updatePW(); 
    }

    w = M2.getW(); 
    M1.computeAvW(); 
    float **avW; 
    avW = M1.getAvW(); 

    for (int i=0; i<M1.getK(); i++){
        for (int j=0; j<M1.getN(); j++){
            for (int k=0; k<2*M1.getL()+1; k++){
                cout << M1.getPW()[i][j][k] << ":"; 
            }
            cout << "  "; 
        }
        cout << endl; 
    }

    int **mPW; 
    mPW = M1.getMostProbW(); 
    for (int i=0; i<K_;i++){
        for (int j=0; j<N_; j++){
            cout << w[i][j] << ":" << avW[i][j] << ":" << mPW[i][j] << "  "; 
        }
        cout << endl; 
    }



    return 0; 

}