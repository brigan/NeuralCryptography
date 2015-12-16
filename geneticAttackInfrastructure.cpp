#include "geneticAttackInfrastructure.h"
#include <iterator>

namespace TPM2 {
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
            tpms.push_back(TPM::TPM(getK(), getN(), getL(), getUpdateRule()));
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
        for (list<TPM::TPM>::iterator it = tpms.begin(); it != tpms.end(); it++) {
            it->compute();
        }
    }

    void GA::updateW() {
        for (list<TPM::TPM>::iterator it = tpms.begin(); it != tpms.end(); it++) {
            it->updateW();
        }
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
        for (list<TPM::TPM>::iterator it = tpms.begin(); it != tpms.end(); it++) {
            it->setX(x_);
        }
    }

    void GA::adjust(int tau) {
        if (tpms.size() < getM()) {
            // for each TPM produce 4 TPMs with tau equal to input parameter
            addTPMs(tau);
            // update weight of TPMs
            updateW();
        } else {
            setXForTPMs(getX());
            // compute
            //compute();
            // cut off TPMs with tau not equal to input parameters.
            deleteUnsuccessfulTPMs(tau);
            // then update weights of remaining TPMs
            updateW();
        }
    }

    void GA::addTPMs(int tau) {
        list<TPM::TPM> newTpms;

        int *pIn = new int[2] {1, -1};
        int pInSize = 2;
        int pOutSize = K;

        list<TPM::TPM>::iterator i = tpms.begin();
        while (i != tpms.end()) {

            int count = 0;

            int* pOut = new int[pOutSize + 1]; // строка из K символов плюс 1 символ для терминального 0
            pOut[pOutSize] = 0;                  // помещаем 0 в конец строки
            int *stack = new int[(pOutSize-1) * 2],  // стек псевдорекурсии, глубина рекурсии K - 1
                    *pTop = stack,            // вершина стека
                    k = 0,                    // переменные цикла
                    n = 0,
                    j = 0;
            for (;;)                      // цикл псевдорекурсии
            {
                while(n < pInSize)
                {
                    pOut[k] = pIn[n++];
                    if (k == (pOutSize-1)) {
                        // have to find only 4 variants
                        if (count > 3) {
                            pTop = stack;
                            break;
                        }

                        int computedTau = 1;
                        for (int z = 0; z < K; z++) {
                            computedTau *= pOut[z];
                        }

                        if (computedTau == tau) {
                            // get new TPM with consistent tau and outputs of perceptrons
                            TPM::TPM newTPM = TPM::TPM(i->getK(), getN(), getL(), getUpdateRule());
                            newTPM.setX(getX());

                            // override results
                            // TODO: move to helper.h
                            int **newW = new int *[K];
                            for (int z = 0; z < K; z++) {
                                newW[z] = new int [N];
                                for (int jz = 0; jz < N; jz++) {
                                    newW[z][jz] = i->getW()[z][jz];
                                }
                            }

                            // TODO: move to helper.h
                            int* newH = new int[pOutSize + 1];
                            for (int z = 0; z < pOutSize + 1; z++) {
                                newH[z] = pOut[z];
                            }
                            newTPM.setPrecomputedResults(newW, newH, tau);
                            newTpms.push_back(newTPM);
                            count++;
                        }
                    } else
                    {
                        if (n < pInSize)
                        {
                            *pTop++ = k;          // сохраняем k и n в стеке
                            *pTop++ = n;
                        }
                        k++;                    // псевдорекурсивный вызов
                        n = 0;
                    }
                }
                if (pTop == stack)          // стек пуст, конец цикла
                    break;

                n = *(--pTop);              // выталкиваем k и n из стека
                k = *(--pTop);
            }
            delete[] pOut;
            delete[] stack;

            i = tpms.erase(i++);
        }

        i = newTpms.begin();
        while (i != newTpms.end()) {
            tpms.push_back(*i);
            newTpms.erase(i++);
        }
    }

    void GA::deleteUnsuccessfulTPMs(int tau) {
        list<TPM::TPM>::iterator i = tpms.begin();
        while (i != tpms.end()) {
            if (i->getTau() != tau) {
                i = tpms.erase(i++);
            } else {
                i++;
            }
        }
    }

    list<TPM::TPM> *GA::getTPMs() {
        return &tpms;
    }
}