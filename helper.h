/* helper.h: 

    In this file it are implemented functions to run TPMs. 

*/ 

#define _helper_h_

// Imports: 
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib> 
#include <vector>
#include "TPM.h"
using namespace std; 


namespace h{

	int **copyW(int K_, int N_, int **&w_){
		/* copyW function: 

			This function copies a weights. 

		Inputs: 
			>> K_: number of hidden units. 
			>> N_: number of inputs per hidden unit. 
			>> w_: weights to be copied. 

		Returns:
			<< newW: copied weights. 

		*/ 

		// Declare/Init variables:
		int **newW; 
		newW = new int*[K_]; 
		for (int i=0; i<K_; i++){
			newW[i] = new int[N_]; 
			for (int j=0; j<N_; j++){
				newW[i][j] = w_[i][j]; 
			}
		}

		return newW; 
	}

	float ***copyPW(int K_, int N_, int L_, float ***&pW_){
		/* copyPW function: 

			This function copies a pWeights. 

		Inputs: 
			>> K_: number of hidden units. 
			>> N_: number of inputs per hidden unit. 
			>> L_: depth of the weights. 
			>> pW_: pWeights to be copied. 

		Returns:
			<< newPW: copied weights. 

		*/ 

		// Declare/Init variables:
		float ***newW; 
		newW = new float**[K_]; 
		for (int i=0; i<K_; i++){
			newW[i] = new float*[N_]; 
			for (int j=0; j<N_; j++){
				newW[i][j] = new float[2*L_+1]; 
				for (int k=0; k<2*L_+1; k++){
					newW[i][j][k] = pW_[i][j][k]; 
				}
			}
		}

		return newW; 
	}

	int **generateX(int K, int N, int **&x){
		/* generateX function: 

			This function generates a new random input for a TPM. 

		Inputs: 
			>> K: number of hidden units. 
			>> N: number of inputs to each hidden unit. 

		Returns: 
			<< x: random array of arrays enclosing an input to the TPM. 

		*/

		// Decalre/Init variables: 
		x = new int*[K]; 
		for (int i=0; i<K; i++){
			x[i] = new int[N]; 
			for (int j=0; j<N; j++){
				if (rand()%2){
					x[i][j] = 1; 
				}
				else{
					x[i][j] = -1; 
				}
			}
		}

		return x; 
	}

	float computeOverlapSingle(int N, int L, int *&w1, int *&w2){
		/* computeOverlapSingle function: 

			This function computes the overlap of two weight vectors. 

		Inputs: 
			>> N: number of input at hidden unit. 
			>> L: depth of the weights. 
			>> w1, w2: vectors from which the overlap will be calculated. 

		Returns: 
			<< rho: overlap of the vectors. 

		*/ 

		// Decalre/Init variables: 
		float rho, Q1=0, Q2=0, R; 

		for (int i=0; i<N; i++){
			Q1 += (w1[i]*w1[i])/N; 
			Q2 += (w2[i]*w2[i])/N; 
			R += (w1[i]*w2[i])/N; 
		}
		rho = R/sqrt(Q1*Q2); 

		return rho; 
	}

	float **computeOverlapMatrix(int N, int L, int *&w1, int *&w2){
		/* computeOverlapMatrix function: 

			This function computes the overlap matrix for the two vectors provided. 

		Inputs: 
			>> N: number of inputs to the hidden units. 
			>> L: depth of the weights. 
			>> w1, w2: weights which overlap is sought. 

		Returns: 
			<< f: overlap matrix of the vectors. 

		*/ 

		// Declare/Init variables: 
		float **f; 
		f = new float*[2*L+1];
		for (int i=0; i<2*L+1; i++){
			f[i] = new float[2*L+1]; 
			for (int j=0; j<2*L+1; j++){
				f[i][j] = 0; 
			}
		}

		for (int i=0; i<N; i++){
			f[w1[i]+L][w2[i]+L] += 1./N; 
		}

		return f; 
	}

	float **computeOverlapMatrix(int N, int L, int *&w1, float **&pW2){
		/* computeOverlapMatrix function: 

			This function overloads the previous one and the following. It is a mixed solution
			between the two: it computes an overlap Matrix with a weight and a pWeight.

		Inputs: 
			>> N: number of inputs per hidden unit. 
			>> L: depth of the weights. 
			>> w1: weight. 
			>> pW2: pWeight. 

		*/ 

		// Declare/Init variables: 
		float **f; 
		f = new float*[2*L+1]; 
		for (int i=0; i<2*L+1; i++){
			f[i] = new float[2*L+1]; 
			for (int j=0; j<2*L+1; j++){
				f[i][j]=0; 
			}
		}

		for (int i=0; i<N; i++){
			for (int j2=0; j2<2*L+1; j2++){
				f[w1[i]+L][j2] += pW2[i][j2]; 
			}
		}

		return f; 
	}

	float **computeOverlapMatrix(int N, int L, float **&pW1, float **&pW2){
		/* computeOverlapMatrix function: 

			This function overloads the previous one and its aim is the same as before but now
			pWeights are provided.

		Inputs: 
			>> N: number of inputs per hidden unit. 
			>> L: depth of the weights. 
			>> pW1, pW2: pWeights. 

		*/ 

		// Declare/Init variables: 
		float **f; 
		f = new float*[2*L+1]; 
		for (int i=0; i<2*L+1; i++){
			f[i] = new float[2*L+1]; 
			for (int j=0; j<2*L+1; j++){
				f[i][j]=0; 
			}
		}

		for (int i=0; i<N; i++){
			for (int j1=0; j1<2*L+1; j1++){
				for (int j2=0; j2<2*L+1; j2++){
					f[j1][j2] += pW1[i][j1]*pW2[i][j2]; 
				}
			}
		}

		return f; 
	}

	float computeOverlapFromMatrix(int N, int L, float **&f){
		/* computeOverlapFromMatrix function: 

			This function computes the overlap of two vectors from their overlap matrix. 

		Inputs: 
			>> N: numberof inputs per hidden unit. 
			>> L: depth of the weights. 
			>> f: overlap matrix. 

		*/ 

		// Declare/Init variables: 
		float rho, R=0, Q1=0, Q2=0; 

		for (int i=0; i<2*L+1; i++){
			for (int j=0; j<2*L+1; j++){
				Q1 += (i-L)*(i-L)*f[i][j]; 
				Q2 += (j-L)*(j-L)*f[i][j]; 
				R += (i-L)*(j-L)*f[i][j]; 
			}
		}
		rho = R/sqrt(Q1*Q2); 

		return rho; 
	}

	float computeOverlapFull(int K, int N, int L, int **&w1, int **&w2){
		/* computeOverlap function: 

			This function computes the overlap between the weights of two TPM. 

		Inputs: 
			>> K: number of hidden units. 
			>> N: number of inputs at hidden unit. 
			>> L: depth of the weights. 
			>> w1, w2: the weights of which overlap must be calculated. 

		Returns: 
			<< rho: overlap between the weights. 

		*/ 

		//Decalre/Init variables: 
		float rho=0; 

		// Loop to compute rho: 
		for (int i=0; i<K; i++){
			float **f = computeOverlapMatrix(N, L, w1[i], w2[i]); 
			rho += computeOverlapFromMatrix(N, L, f)/K; 
		}

		return rho; 
	}

	float computeOverlapFull(int K, int N, int L, int **&w1, float ***&pW2){
		/* computeOverlap function: 

			This function is a mixed solution between the previous and the following ones. 

		Inputs: 
			>> K: number of hidden units. 
			>> N: number of inputs at hidden unit. 
			>> L: depth of the weights. 
			>> w1: weight. 
			>> pW2: pWeight. 

		Returns: 
			<< rho: overlap between the weights. 

		*/ 

		//Decalre/Init variables: 
		float rho=0; 

		// Loop to compute rho: 
		for (int i=0; i<K; i++){
			float **f = computeOverlapMatrix(N, L, w1[i], pW2[i]); 
			rho += computeOverlapFromMatrix(N, L, f)/K; 
		}

		return rho; 
	}

	float computeOverlapFull(int K, int N, int L, float ***&pW1, float ***&pW2){
		/* computeOverlap function: 

			This function overwrites the previous one and computes the overlap from pWeights. 

		Inputs: 
			>> K: number of hidden units. 
			>> N: number of inputs at hidden unit. 
			>> L: depth of the weights. 
			>> pW1, pW2: the pWeights of which overlap must be calculated. 

		Returns: 
			<< rho: overlap between the weights. 

		*/ 

		//Decalre/Init variables: 
		float rho=0; 

		// Loop to compute rho: 
		for (int i=0; i<K; i++){
			float **f = computeOverlapMatrix(N, L, pW1[i], pW2[i]); 
			rho += computeOverlapFromMatrix(N, L, f)/K; 
		}

		return rho; 
	}

	void synchTPM(TPM::TPM M1, TPM::TPM M2){
		/* syncTPM function: 

			This function implements the synchronization of two TPMs. 

		Inputs: 
			>> M1, M2: TPMs to be synchronized. 

		*/ 

		// Declare/Init variables: 
		int K=M1.getK(), N=M1.getN(), L=M1.getL(); 
		if ((K!=M2.getK()) or (N!=M2.getN()) or (L!=M2.getL())){
			cout << "Different set of parameters in synchronizing machines. Not possible to synchronize: " << endl; 
			cout << "\t" << K << "  " << M2.getK() << endl; 
			cout << "\t" << N << "  " << M2.getN() << endl; 
			cout << "\t" << L << "  " << M2.getL() << endl; 
		}

		float rho=0; 
		int **x, **w1, **w2; 

		int t=0; 
		while (rho<1){
			t++; 
			
			x = generateX(K, N, x); 
			M1.setX(x); 
			M2.setX(x); 
			if (M1.getTau()==M2.getTau()){
				M1.updateW(); 
				M2.updateW(); 
			}

			w1 = M1.getW(); 
			w2 = M2.getW(); 
			rho = computeOverlapFull(K, N, L, w1, w2); 
			cout << t << "  " << rho << endl; 
		}

		return; 
	}

	float ***extendW(int K_, int N_, int L_, int **w){
		/* extendW function: 

			This function extends a weights instance into a pW one. 

		Inputs: 
			>> w: weights. 

		Returns: 
			<< pW: pWeights build from w. 

		*/ 

		// Declare/Init variables: 
		float ***pW; 
		pW = new float**[K_]; 
		for (int i=0; i<K_; i++){
			pW[i] = new float*[N_]; 
			for (int j=0; j<N_; j++){
				pW[i][j] = new float[2*L_+1]; 
				for (int k=0; k<2*L_+1; k++){
					pW[i][j][k] = 0.; 
				}
				pW[i][j][w[i][j]+L_] = 1.; 
			}
		}

		return pW; 
	}
}


