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


/* TPM.h: 

	In this file it is implemented the class of TPM objects. TPM stands for Tree Parity Machines,
	which is currently (2015) the most widely used means for Neural Cryptography. As of today, it
	has not been reported any successful attack on TPM-based cryptography.

	In this file it is also implemented a class P_TPM whose name stands for probabilistic Tree
	Parity Machine. This class aims at an attack on Tree Parity Machines based on a probabilistic
	inference of its weights. This line of attack is inspired by [1], where a probabilistic attack
	is implemented that breaks the Permutation Parity Machine (PPM)-based cryptography. The class
	P_TPM allows us to implement similar attacks on TPMs. Some such attacks have been tried and the
	results are promising, but they indicate that more research is necessary to finish breaking the
	TPM-based cryptography. 

	If you plan to do research on breaking the TPM or any newer machines, please contact me, I might
	be interested! Please, contact me and cite this software if you use it to implement your
	research on Neural Cryptography. (Thanks beforehand!)

	Code with a python implementation of TPMs and PPMs is also available upon request, but it has
	not been liberated because it is old and would require a review. Also, python is notably slower
	than c++... Notably, among the software available upon request is that used to implement the
	attacks in [1].

    [1] Seoane LF and Ruttor A, Successful attack on permutation-parity-machine-based neural
cryptography. Physical Review E, 85(2), 025101 (2012).

*/

#ifndef _TPM_h_
#define _TPM_h_


// Imports: 
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib> 
#include <vector>
using namespace std; 


namespace TPM{


	///////////////////////////////
	// 
	//  TPM Class: 

	class TPM {
			int K, N, L; 
			int **x, **w, *h, tau; 
			int updateRule; 
		public: 
			// Init functions:
			TPM(int, int, int, int); 
			void initW(); 

			// Get functions: 
			int getK(), getN(), getL(); 
			int **getX(), **getW(), *getH(), getTau(); 
			int getUpdateRule(); 

			// Set functions: 
			void setK(int), setN(int), setL(int); 
			void setX(int **), setW(int **); 
			void setUpdateRule(int); 

			// Compute functions: 
			void computeH(), computeTau(), compute(); 
			void updateW();
			void setPrecomputedResults(int **, int *, int);

	}; 


		
}
#endif