//
// Created by Nurlan Akashayev on 11/30/15.
//

#include <list>
#include "TPM.h"

using namespace std;

#ifndef NEURALCRYPTOGRAPHY_GENETICATTACKINFRASTRUCTURE_H
#define NEURALCRYPTOGRAPHY_GENETICATTACKINFRASTRUCTURE_H

#endif //NEURALCRYPTOGRAPHY_GENETICATTACKINFRASTRUCTURE_H

namespace TPM {
    class GA {
        int TPMCount; // Total number of TPMs
        int M; // Limit which determines GA's behavior
        int K, N, L;
        int updateRule;
        int **x;
        //TPM *tpms;
        list<TPM> tpms;

    public:

        GA(int, int, int, int, int, int);

        // Get functions:
        int getTPMCount();
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

    GA::GA(int TPMCount_, int K_, int N_, int L_, int updateRule_, int M_) {
        setTPMCount(TPMCount_);
        setK(K_);
        setN(N_);
        setL(L_);
        setUpdateRule(updateRule_);
        setM(M_);
        initGA();
    }

    void GA::initGA() {
        for (int i = 0; i < getTPMCount(); i++) {
            tpms.push_back(TPM::TPM(getK(), getL(), getN(), getUpdateRule()));
        }
    }

    int GA::getTPMCount() {
        return TPMCount;
    }

    int GA::getK() {
        return K;
    }

    int GA::getN() {
        return N;
    }

    int GA::getL() {
        return L;
    }

    int GA::getUpdateRule() {
        return updateRule;
    }


    int **GA::getX() {
        return x;
    }

    int GA::getM() {
        return M;
    }

    void GA::compute() {
        for (list<TPM>::iterator it = tpms.begin(); it != tpms.end(); it++) {
            it->compute();
        }
    }

    void GA::updateW() {

    }

    void GA::setX(int **x_) {
        x = x_;
    }

    void GA::setTPMCount(int TPMCount_) {
        TPMCount = TPMCount_;
    }

    void GA::setK(int K_) {
        K = K_;
    }

    void GA::setN(int N_) {
        N = N_;
    }

    void GA::setL(int L_) {
        L = L_;
    }

    void GA::setUpdateRule(int updateRule_) {
        updateRule = updateRule_;
    }

    void GA::setM(int M_) {
        M = M_;
    }

    void GA::setXForTPMs(int **x_) {
        for (list<TPM>::iterator it = tpms.begin(); it != tpms.end(); it++) {
            it->setX(x_);
        }
    }

    void GA::adjust(int tau) {
        if (getTPMCount() < getM()) {
            // for each TPM produce 4 TPMs with tau equal to input parameter
            addTPMs(tau);
            // it makes no sense to compute due to determined state
            // then we prepare X
            setXForTPMs(getX());
            // and just update weight of TPMs
            updateW();
        } else {
            setXForTPMs(getX());
            // compute
            compute();
            // cut off TPMs with tau not equal to input parameters.
            deleteUnsuccessfulTPMs(tau);
            // then update weights of remaining TPMs
            updateW();
        }
    }

    void GA::addTPMs(int tau) {
        list<TPM> newTpms;

        list<TPM>::iterator i = tpms.begin();
        while (i != tpms.end()) {

        }

        tpms.merge(newTpms);
    }

    void GA::deleteUnsuccessfulTPMs(int tau) {
        list<TPM>::iterator i = tpms.begin();
        while (i != tpms.end()) {
            if (i->getTau() != tau) {
                i = tpms.erase(i++);
            } else {
                i++;
            }
        }
    }

}