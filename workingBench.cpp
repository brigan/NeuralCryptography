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
	int K_=2, N_=10, L_=50; 
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
	// TPM::PTPM ME, MA; 
	srand(seed); 
	TPM::pTPM ME(nSample_, nReset_, K_, N_, L_, updateRule_); 
	TPM::TPM MA(K_, N_, L_, updateRule_), MB(K_, N_, L_, updateRule_); 

	int **w1, **w2, **w3; 
	float ***pW1; 
	float rho12, rho12_, rho23; 

	// Full dynamic attack: 
	int t=0; 
	rho23 = 0;  
	rho12_ = 0; 
	while ((rho23<1 and rho23>-1) and (rho12_<1 and rho12_>-1)){
		t++; 
		x = h::generateX(K_, N_, x); 
		MA.setX(x); // setX() includes compute();
		MB.setX(x); 
		ME.setX(x); 
		ME.setTargetTau(MA.getTau()); 
		ME.MC_updatePW(); 
		if (MA.getTau()==MB.getTau()){
			ME.learn_updatePW(); 
			MA.updateW(); 
			MB.updateW(); 
		}

		pW1 = ME.getPW(); 
		w1 = ME.getMostProbW(); 
		w2 = MA.getW(); 
		w3 = MB.getW(); 
		rho23 = h::computeOverlapFull(K_, N_, L_, w2, w3); 
		rho12 = h::computeOverlapFull(K_, N_, L_, w2, pW1); 
		rho12_ = h::computeOverlapFull(K_, N_, L_, w2, w1); 
		cout << t << "  " << rho23 << "  " << rho12 << "  " << rho12_ << endl; 
		if (t==10){
			ME.preferMPW(-1,-1);
		}
		// if (t==50){
		// 	ME.preferMPW(0,8);
		// 	pW1 = ME.getPW(); 
		// 	for (int k=0; k<(2*ME.getL()+1); k++){
		// 		cout << k-ME.getL() << "  " << pW1[0][8][k] << endl; 
		// 	}
		// 	return 0; 
		// }
	}

	return 0; 

}