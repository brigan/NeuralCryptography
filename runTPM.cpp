/*

	Copyright (C) 2012, 2013, 2015 Luis F Seoane. 

		Contact: luis.seoane@upf.edu, brigan@gmail.com


	This file is part of NeuralCryptography.

    NeuralCryptography is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    NeuralCryptography is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with NeuralCryptography.  If not, see <http://www.gnu.org/licenses/>.

*/ 



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
	int K_=2, N_=100, L_=50; 
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
	// TPM::P_TPM ME, MA; 
	srand(seed); 
	TPM::P_TPM ME(nSample_, nReset_, K_, N_, L_, updateRule_); 
	TPM::TPM MA(K_, N_, L_, updateRule_), MB(K_, N_, L_, updateRule_); 

	int **w1, **w2, **w3; 
	float ***pW1; 
	float rho12, rho12_, rho23; 


	// // Just synchronize: 
	// h::synchTPM(MA, MB); 

	// // Statich attach: 
	// for (int it=0; it<400; it++){
	//	 x = h::generateX(K_, N_, x); 
	//	 MA.setX(x); 
	//	 ME.setX(x); 
	//	 ME.setTargetTau(MA.getTau()); 
	//	 ME.MC_updatePW(); 

	//	 pW1 = ME.getPW(); 
	//	 w2 = MA.getW(); 
	//	 rho12 = h::computeOverlapFull(K_, N_, L_, w2, pW1); 
	//	 cout << it << "  " << rho12 << endl; 
	// }

	// Full dynamic attack: 
	int t=0; 
	rho23 = 0;  
	rho12_ = 0; 
	while ((rho23<1 and rho23>-1) and (rho12_<1 and rho12_>-1)){
		t++; 
		x = h::generateX(K_, N_, x); 
		MA.setX(x); 
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
		if (t==100){
			pW1 = ME.getPW(); 
			for (int k=0; k<(2*ME.getL()+1); k++){
				cout << k-ME.getL() << "  " << pW1[0][0][k] << endl; 
			}
			return 0; 
		}
	}

	return 0; 

}