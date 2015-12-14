//
// Created by Nurlan Akashayev on 11/30/15.
//
#ifndef NEURALCRYPTOGRAPHY_GENETICATTACKINFRASTRUCTURE_H
#define NEURALCRYPTOGRAPHY_GENETICATTACKINFRASTRUCTURE_H

#include <list>
#include "TPM.h"

using namespace std;

namespace TPM2 {
    class GA {
        int TPMCount; // Total number of TPMs TODO: rename to initialTPMCount
        int M; // Limit which determines GA's behavior
        int K, N, L;
        int updateRule;
        int **x;
        //TPM *tpms;
        list<TPM::TPM> tpms;

    public:

        GA(int, int, int, int, int, int);

        // Get functions:
        int getTPMCount();
        list<TPM::TPM> *getTPMs();

        int getM();
        int getK(), getN(), getL();
        int getUpdateRule();
        int **getX();

        // Compute functions:
        void compute();
        void updateW();

        void adjust(int);

        // Set functions:
        void setX(int **);
    private:
        // Init functions:
        void initGA();

        void addTPMs(int);
        void deleteUnsuccessfulTPMs(int);

        // Set functions:
        void setTPMCount(int);
        void setM(int);
        void setK(int), setN(int), setL(int);
        void setUpdateRule(int);
        void setXForTPMs(int **);
    };
}
#endif