//
// Created by Nurlan Akashayev on 11/29/15.
//

// Imports:
#include <iostream>

#include <vector>
#include "helper.h"
#include "geneticAttackInfrastructure.h"

using namespace std;

int main(int argc, char* argv[]) {
    /* main for runGeneticAttack.cpp */

    int seed=0;
    int K_=4, N_=100, L_=50;
    int **x, **w, *h, tau;
    int updateRule_=1;
    int TPMCount_ = 100;
    int M_ = 100;

    int nSample_=10000, nReset_=100000;

    for (int i=0; i<argc; i++){
        if (string(argv[i])=="-uR"){
            updateRule_ = atoi(argv[i+1]);
        }
        if (string(argv[i])=="-s"){
            seed = atoi(argv[i+1]);
        }
    }
    //	 if (string(argv[i])=="-sC"){
    //		 stopCount = atoi(argv[i+1]);
    //	 }
    //	 if (string(argv[i])=="-mC"){
    //		 maxCount = atoi(argv[i+1]);
    //	 }
    //	 if (string(argv[i])=="-sC_"){
    //		 stepCount = atoi(argv[i+1]);
    //	 }
    //	 if (string(argv[i])=="-tC"){
    //		 topCount = atoi(argv[i+1]);
    //	 }
    //	 if (string(argv[i])=="-sC__"){
    //		 saveCount = atoi(argv[i+1]);
    //	 }
    // }


    // Defining a TPM:
    // TPM::P_TPM ME, MA;
    srand(seed);
    TPM::TPM M1(K_, N_, L_, updateRule_), M2(K_, N_, L_, updateRule_);
    TPM2::GA GA(TPMCount_, K_, N_, L_, updateRule_, M_);

    /* syncTPM function:

			This function implements the synchronization of two TPMs.

		Inputs:
			>> M1, M2: TPMs to be synchronized.

		*/
//
//    // Declare/Init variables:
//    int K=M1.getK(), N=M1.getN(), L=M1.getL();
//    if ((K!=M2.getK()) or (N!=M2.getN()) or (L!=M2.getL())){
//        cout << "Different set of parameters in synchronizing machines. Not possible to synchronize: " << endl;
//        cout << "\t" << K << "  " << M2.getK() << endl;
//        cout << "\t" << N << "  " << M2.getN() << endl;
//        cout << "\t" << L << "  " << M2.getL() << endl;
//    }

    float rho=0;
    float ga_max_rho = 0;
    int **generatedX, **w1, **w2;

    int t=0;
    while (rho<1){
        t++;

        generatedX = h::generateX(K_, N_, generatedX);
        M1.setX(generatedX);
        M2.setX(generatedX);
        GA.setX(generatedX);

        if (M1.getTau()==M2.getTau()) {
            GA.adjust(M1.getTau());
            M1.updateW();
            M2.updateW();
        }

        w1 = M1.getW();
        w2 = M2.getW();
        rho = h::computeOverlapFull(K_, N_, L_, w1, w2);
        cout << t << "  " << rho << endl;

        ga_max_rho = h::computeOverlapFull(K_, N_, L_, w1, GA.getTPMs());
        cout << t << "  " << ga_max_rho << "<-- GA"<< endl;
    }

    return 0;

}
