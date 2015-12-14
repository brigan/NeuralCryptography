//
// Created by Nurlan Akashayev on 12/14/15.
//

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <vector>
#include "TPM.h"
using namespace std;


namespace TPM {

///////////////////////////////
//
//  TPM Functions:

///////////////////////////////
//
//  Init Functions

    TPM::TPM(int K_, int N_, int L_, int updateRule_) {
        /* init function for class TPM:

            This function initializes objects of class TPM. Weights are initialized randomly.

        Inputs:
            >> N_: number of inputs per hidden unit.
            >> K_: number of hidden units.
            >> L_: depth of the weights.
            >> updateRule: rule to update the weights:
                > 1: Hebbian learning.
                > 2: Anti-Hebbian learning.
                > 3: Random walk.

        */

        setK(K_);
        setN(N_);
        setL(L_);
        initW();
        setUpdateRule(updateRule_);


        return;
    }

    void TPM::initW() {
        /* initW function:

            This function initializes the weights of the TPM.

        */

        for (int i = 0; i < K; i++) {
            for (int j = 0; j < N; j++) {
                w[i][j] = rand() % (2 * L + 1) - L;
            }
        }
    }



///////////////////////////////
//
//  Get Functions:

    int TPM::getK() {
        /* getK function:

            This function returns the number of hidden units of the TPM.

        Returns:
            << K: number of hidden units.

        */

        return K;
    }

    int TPM::getN() {
        /* getN function:

            This function returns the number of inputs per hidden unit.

        Returns:
            << N: inputs per hidden unit.

        */

        return N;
    }

    int TPM::getL() {
        /* getL function:

            This function returns the depth of the weights.

        Returns:
            << L: depth of the weights.

        */

        return L;
    }

    int **TPM::getX() {
        /* getX function:

            This function returns the current input of the TPM.

        Returns:
            << x: current input.

        */

        return x;
    }

    int **TPM::getW() {
        /* getW function:

            This function returns the current weights.

        Returns:
            << w: weights of the TPM.

        */

        return w;
    }

    int *TPM::getH() {
        /* getK function:

            This function returns the current states of the hidden units.

        Returns:
            << h: array with the states of the hidden units.

        */

        return h;
    }

    int TPM::getTau() {
        /* getTau function:

            This function returns the output of the TPM.

        Returns:
            << tau: output of the TPM.

        */

        return tau;
    }

    int TPM::getUpdateRule() {
        /* getUpdateRule function:

            This function returns the update rule of the TPM as an integer.

        Returns:
            << updateRule: rule to update the weights encoded in an integer:
                > 1: Hebbian learning.
                > 2: Anti-Hebbian learning.
                > 3: Random walk update.

        */

        return updateRule;
    }


///////////////////////////////
//
//  Set Functions:

    void TPM::setK(int K_) {
        /* setK function:

            This function sets up the number of hidden unit. It also allocates memory for x, w and h.

        Inputs:
            >> K_: number of hidden unit.

        */

        K = K_;
        x = new int *[K];
        w = new int *[K];
        h = new int[K];

        return;
    }

    void TPM::setN(int N_) {
        /* setN function:

            This function sets up the number of inputs per hidden unit. It also allocates memory
        for x and w.

        Inputs:
            >> N_: number of input bits per hidden unit.

        */

        N = N_;
        for (int i = 0; i < K; i++) {
            x[i] = new int[N];
            w[i] = new int[N];
        }

        return;
    }

    void TPM::setL(int L_) {
        /* setL function:

            This function sets up the depth of the weights.

        Inputs:
            >> L_: depth of the weights.

        */

        L = L_;

        return;
    }

    void TPM::setX(int **x_) {
        /* setN function:

            This function sets up the input of the TPM.

        Inputs:
            >> x_: input of the TPM.

        */

        x = x_;
        compute();

        return;
    }

    void TPM::setW(int **w_) {
        /* setW function:

            This function sets up the weights of the TPM.

        Inputs:
            >> w_: weights of the TPM.

        */

        w = w_;

        return;
    }

    void TPM::setUpdateRule(int updateRule_) {
        /* setUpdateRule function:

            This function sets the rule used by the TPM to implement the update of the weights.

        */

        updateRule = updateRule_;

        return;
    }


///////////////////////////////
//
//  Compute Functions:

    void TPM::computeH() {
        /* computeH function:

            This function computes the activation levels of the hidden units with the current
            weights and inputs. The activation levels will be integers between -N*L and N*L,
            although what is needed for the TPM to compute are just the signs of these
            activation functions. It might become handy in the future to have this info and from
            there it is straight forward to get a sign function.

        */

        for (int i = 0; i < K; i++) {
            h[i] = 0;
            for (int j = 0; j < N; j++) {
                h[i] += x[i][j] * w[i][j];
            }
            if (h[i] == 0) h[i] = -1;    // Avoids having tau=0!!
        }

        return;
    }

    void TPM::computeTau() {
        /* computeTau function:

            This function computes the output of the TPM.

        */

        int hProc = 1;
        for (int i = 0; i < K; i++) {
            hProc *= h[i];
        }
        if (hProc <= 0) {
            tau = -1;
        }
        else {
            tau = 1;
        }

        return;
    }

    void TPM::compute() {
        /* compute function:

            This function performs a whole TPM computation by calling computeH() and
            computeTau().

        */

        computeH();
        computeTau();

        return;
    }

    void TPM::updateW() {
        /* updateW function:

            This function updates the weights of the TPM. The rule to update the weights is
            stored as an integer in updateRule:
                >> updateRule=1: Hebbian learning.
                >> updateRule=2: Hebbian learning.
                >> updateRule=3: Hebbian learning.

        */

        // updateFactor implements the different learning rules:
        int updateFactor;
        if (updateRule == 1) {
            updateFactor = tau;
        }
        else if (updateRule == 2) {
            updateFactor = -tau;
        }
        else if (updateRule == 3) {
            updateFactor = 1;
        }

        // Apply learning:
        for (int i = 0; i < K; i++) {
            if (tau * h[i] > 0) {
                for (int j = 0; j < N; j++) {
                    w[i][j] += updateFactor * x[i][j];
                    if (w[i][j] > L) w[i][j] = L;
                    if (w[i][j] < -L) w[i][j] = -L;
                }
            }
        }


        return;
    }

    void TPM::setPrecomputedResults(int **w_, int *h_, int tau_) {
        w = w_;
        h = h_;
        tau = tau_;
    }



///////////////////////////////
//
//  P_TPM Class:

    class P_TPM : public TPM {
        int nSample, nReset, targetTau;
        float ***pW, *pH, pTau;
        float **avW, **sig2W;
    public:
        // Init functions:
        P_TPM(int, int, int, int, int, int);

        void initPW(), initAsymPW(), normPW(int, int);

        // Get functions:
        float ***getPW(), *getPH(), getPTau();

        float **getAvW(), **getSig2W();

        int getNSample(), getNReset(), getTargetTau();

        int **getMostProbW();

        // Set functions:
        void setNSample(int), setNReset(int), setTargetTau(int);

        void setPW(float ***&);

        // Compute functions:
        void computeAvW(), computeSig2W();

        void computePH(), computePTau(), pCompute(bool fComputePTau);

        void drawW(int **&), MC_updatePW(), learn_updatePW();

        void resetPW(), preferMPW(int, int);

    };


///////////////////////////////
//
//  P_TPM Functions:

///////////////////////////////
//
//  Init Functions:

    P_TPM::P_TPM(int nSample_, int nReset_, int K_, int N_, int L_, int updateRule_) : TPM(K_, N_, L_, updateRule_) {
        /* Constructor function for class P_TPM:

            This constructor should implement a few things more (allocating and init pW) apart
            from calling super's constructor.

        Inputs:
            >> nSample_: number of MC samples.
            >> nReset_: number of failed samples after which weights of a unit are reset.
            >> K_: number of hidden units.
            >> N_: number of inputs per hidden unit.
            >> L_: depth of the weights.
            >> updateRule_: learning rule.

        */
        TPM(K_, N_, L_, updateRule_);
        setNSample(nSample_);
        setNReset(nReset_);
        // initPW();
        initAsymPW();
        pH = new float[getK()];

        return;
    }

    void P_TPM::initPW() {
        /* initPW function:

            This function initializes the probabilistic weights of the machine.
        */
        pW = new float **[getK()];
        avW = new float *[getK()];
        sig2W = new float *[getK()];
        for (int i = 0; i < getK(); i++) {
            pW[i] = new float *[getN()];
            avW[i] = new float[getN()];
            sig2W[i] = new float[getN()];
            for (int j = 0; j < getN(); j++) {
                pW[i][j] = new float[2 * getL() + 1];
                for (int k = 0; k < 2 * getL() + 1; k++) {
                    pW[i][j][k] = 1. / (2 * getL() + 1);
                }
            }
        }


        return;
    }

    void P_TPM::initAsymPW() {
        /* initAsymmPW function:

            This function initializes the pWeights such that an asymmetry is present since the
        very beginning. I hope this helps accelerate the convergence.

        */

        pW = new float **[getK()];
        avW = new float *[getK()];
        sig2W = new float *[getK()];
        for (int i = 0; i < getK(); i++) {
            pW[i] = new float *[getN()];
            avW[i] = new float[getN()];
            sig2W[i] = new float[getN()];
            for (int j = 0; j < getN(); j++) {
                pW[i][j] = new float[2 * getL() + 1];
                float sum = 0;
                for (int k = 0; k < 2 * getL() + 1; k++) {
                    pW[i][j][k] = float(rand());
                }
            }
        }
        normPW(-1, -1);

        return;
    }

    void P_TPM::normPW(int i = -1, int j = -1) {
        /* normPW function:

            This function normalizes the probability distribution associated with pW j of unit
            i. If negative indexes i,j are provided, the normalization is carried out for all
            the probability distributions.

        Inputs:
            >> i=-1,j=1: indexes of the pW to normalize. If they are negative, all the pW are
            normalized.

        */

        float sum = 0.;
        if (j < 0 and i < 0) {
            for (int i = 0; i < getK(); i++) {
                for (int j = 0; j < getN(); j++) {
                    sum = 0.;
                    for (int k = 0; k < 2 * getL() + 1; k++) {
                        sum += pW[i][j][k];
                    }
                    for (int k = 0; k < 2 * getL() + 1; k++) {
                        pW[i][j][k] /= sum;
                    }
                }
            }
        }
        else {
            sum = 0.;
            for (int k = 0; k < 2 * getL() + 1; k++) {
                sum += pW[i][j][k];
            }
            for (int k = 0; k < 2 * getL() + 1; k++) {
                pW[i][j][k] /= sum;
            }
        }

        return;
    }


///////////////////////////////
//
//  Get Functions:

    float ***P_TPM::getPW() {
        /* getPW function:

            This function returns the probability that the P_TPM has got each of the possible weights.

        Returns:
            >> pW: probability that the machine has got each of the possible weights.

        */

        return pW;
    }

    float *P_TPM::getPH() {
        /* getPH function:

            This function returns the probability that each of the hidden units is activated (+).
        In this case it is just stored the probability that the sign of the hidden unit is + or -
        as it would take a lot to compute the probability of the hidden state taking each of the
        (2*L+1)*N possible values.

        Returns:
            >> pH: probability that each of the hidden units is activated (+).

        */

        return pH;
    }

    float P_TPM::getPTau() {
        /* getPTau function:

            This function returns the probability that the whole P_TPM is activated.

        Returns:
            >> pTau: probability that the whole P_TPM is activated.

        */

        return pTau;
    }

    float **P_TPM::getAvW() {
        /* getAvW function:

            This function returns the average weights given the current pW.

        Returns:
            << avW: average weights.

        */

        return avW;
    }

    float **P_TPM::getSig2W() {
        /* getSig2W function:

            This function returns the squared average of deviations from the mean of the weights.

        Returns:
            << sig2W: squared average of deviations from the mean of the weights.

        */

        return sig2W;
    }

    int P_TPM::getNSample() {
        /* getNSample function:

            This function returns the number of MC samples.

        Returns:
            << nSample: number of MC samples.

        */

        return nSample;
    }

    int P_TPM::getNReset() {
        /* getNReset function:

            This function returns the number of failed samples after which weights are reset.

        Returns:
            << nReset: number of MC samples.

        */

        return nReset;
    }

    int P_TPM::getTargetTau() {
        /* getTargetTau function:

            This function returns the target output.

        Returns:
            << targetTau: current target tau.

        */

        return targetTau;
    }

    int **P_TPM::getMostProbW() {
        /* getMostProbW function:

            This function returns the most probable weight.

        Returns:
            >> mPW: most probable weight, meaning that with a larger probability of happening.

        */

        // Declare/Init variables:
        int **mPW;
        mPW = new int *[getK()];

        for (int i = 0; i < getK(); i++) {
            mPW[i] = new int[getN()];
            for (int j = 0; j < getN(); j++) {
                float maxP = -1;
                for (int k = 0; k < 2 * getL() + 1; k++) {
                    if (pW[i][j][k] > maxP) {
                        mPW[i][j] = k - getL();
                        maxP = pW[i][j][k];
                    }
                }
            }
        }

        return mPW;
    }


///////////////////////////////
//
//  Set Functions:

    void P_TPM::setNSample(int nSample_) {
        /* setNSample function:

            This function sets up the number of MC samples.

        */

        nSample = nSample_;

        return;
    }

    void P_TPM::setNReset(int nReset_) {
        /* setNReset function:

            This function sets the number of failed samples after which a pW instance is
            reset.

        */

        nReset = nReset_;

        return;
    }

    void P_TPM::setTargetTau(int targetTau_) {
        /* setTargetTau function:

            This function sets a target Tau.

        Inputs:
            >> targetTau_: desired target Tau.

        */

        targetTau = targetTau_;

        return;
    }

    void P_TPM::setPW(float ***&pW_) {
        /* setPW function:

            This function sets a new probability distribution of the weights.

        Inputs:
            >> pW_: new probability distribution of the weights.

        */

        for (int i = 0; i < getK(); i++) {
            for (int j = 0; j < getN(); j++) {
                for (int k = 0; k < 2 * getL() + 1; k++) {
                    pW[i][j][k] = pW_[i][j][k];
                }
            }
        }

        delete[] pW_;

        return;
    }


///////////////////////////////
//
//  Compute Functions:

    void P_TPM::computeAvW() {
        /* computeAvW function:

            This function computes the average weights.

        */

        for (int i = 0; i < getK(); i++) {
            for (int j = 0; j < getN(); j++) {
                avW[i][j] = 0;
                for (int k = 0; k < 2 * getL() + 2; k++) {
                    avW[i][j] += (k - getL()) * pW[i][j][k];
                }
            }
        }

        return;
    }

    void P_TPM::computeSig2W() {
        /* computesig2W function:

            This function computes the squared average deviation of the weights.

        */

        for (int i = 0; i < getK(); i++) {
            for (int j = 0; j < getN(); j++) {
                sig2W[i][j] = 0.;
                for (int k = 0; k < 2 * getL() + 1; k++) {
                    sig2W[i][j] += pW[i][j][k] * float(pow(k - float(getL()) - avW[i][j], 2));
                }
            }
        }

        return;
    }

    void P_TPM::computePH() {
        /* computePH function:

            This function computes pH given the current probability distribution of the weights.
            To implement this computation, the distribution of the weights is approximated by a
            Gaussian function and the sign of the average is inverted whenever the associated
            input is negative.

        */

        float avH, sig2H;
        for (int i = 0; i < getK(); i++) {

            // Compute averages and deviations:
            avH = 0;
            sig2H = 0;
            for (int j = 0; j < getN(); j++) {
                avH += getX()[i][j] * avW[i][j];
                sig2H += sig2W[i][j];
            }

            // Map into a N(0,1) distribution problem:
            float z = -avH / sqrt(sig2H);
            pH[i] = (1 / 2 + erfc(z / sqrt(2)) / 2);
        }

        return;
    }

    void P_TPM::computePTau() {
        /* computePTau function:

            This function computes the probability that tau=1 with the current settings. Because
        this computation can become more difficult for large K, it is only explicitly implemented
        for K=1, K=2, and K=3. For larger K it could be approximated as a product of Gaussian
        vairables, but this might not be neccessary as TPM with K>3 do not gain so much capabilities
        for cryptography.

        */

        if (getK() == 1) {
            pTau = pH[0];
        }
        else if (getK() == 2) {
            pTau = (pH[0] * pH[1]) + ((1 - pH[0]) * (1 - pH[1]));
        }
        else if (getK() == 3) {
            pTau = (pH[0] * pH[1] * pH[2]) + ((1 - pH[0]) * (1 - pH[1]) * pH[2]) +
                   ((1 - pH[0]) * pH[1] * (1 - pH[2])) + (pH[0] * (1 - pH[1]) * (1 - pH[2]));
        }


        return;
    }

    void P_TPM::pCompute(bool fComputePTau = false) {
        /* pCompute function:

            This function computes the probabilistic variables.

        */

        computeAvW();
        computeSig2W();
        computePH();
        if (fComputePTau) computePTau();

        return;
    }

    void P_TPM::drawW(int **&rW) {
        /* drawW function:

            This function draws random weights from the current probability distribution.

        Inputs:
            >> rW: array where the random weights will be stored. This works as a place to
            store the output.

        */

        // Declare/Init variables:
        rW = new int *[getK()];
        float rnd;

        for (int i = 0; i < getK(); i++) {
            rW[i] = new int[getN()];
            for (int j = 0; j < getN(); j++) {
                rnd = double(rand()) / RAND_MAX;
                rW[i][j] = -getL();
                for (int k = 0; k < 2 * getL() + 1; k++) {
                    if (rnd < pW[i][j][k]) break;
                    rW[i][j]++;
                    rnd -= pW[i][j][k];
                }
            }
        }

        return;
    }


    void P_TPM::MC_updatePW() {
        /* MC_updateW function:

            This function implements a Monte-Carlo sampling of the weights space according to
            their current probability distributions.

        */

        float ***averaged;                    // In this variable it is stored the new pW!
        averaged = new float **[getK()];
        for (int i = 0; i < getK(); i++) {
            averaged[i] = new float *[getN()];
            for (int j = 0; j < getN(); j++) {
                averaged[i][j] = new float[2 * getL() + 1];
                for (int k = 0; k < 2 * getL() + 1; k++) {
                    averaged[i][j][k] = 0.;
                }
            }
        }

        int sampleCount = 0, failedCount = 0;
        int **candidateW;
        while (sampleCount < nSample) {
            drawW(candidateW);
            setW(candidateW);
            compute();
            if (getTau() == targetTau) {
                sampleCount++;
                failedCount = 0;
                for (int i = 0; i < getK(); i++) {
                    for (int j = 0; j < getN(); j++) {
                        averaged[i][j][getW()[i][j] + getL()] += 1. / nSample;
                    }
                }
            }
            else {
                failedCount++;
                if (failedCount == nReset) {
                    resetPW();
                    failedCount = 0;
                }
            }
            for (int i = 0; i < getK(); i++) delete[] candidateW[i];
            // delete[] candidateW;
        }

        setPW(averaged);
        setW(getMostProbW());

        return;
    }

    void P_TPM::learn_updatePW() {
        /* learn_updatePW function:

            This function updates the pWeights of the P_TPM. The rule to update the weights is stored
        as an integer in updateRule:
                >> updateRule=1: Hebbian learning.
                >> updateRule=2: Hebbian learning.
                >> updateRule=3: Hebbian learning.

            The update for pW is tricky:

                >> Update is dictated by the TPM since we wish to update whenever the synchronizing
        machines do so! It could be implemented an update rule proportional to p(tau=targetTau), but
        in this case we would be acting as if the synchronizing TPMs were taking into account what
        the P_TPM says, which is not the case.
                >> The product: updateFactor*x[i][j] determines the direction of the update. The
        pWeights are shifted accordingly.
                >> pH[i] is taken into account indeed. Unlike with p(tau=targetTau), in which both
        TPM update notwithstanding what P_TPM says; the weights of a unit which doesn't match the
        output is not updated, even for the TPMs. We can not guess which units are updated in each.
        The best we can do is approximate this by the probability that each unit is active. This
        defines the flux transferred.

        */

        // pCompute: This step is necessary to calculate pH.
        int updateFactor;
        pCompute();

        // Factor factor implements the different learning rules:
        if (getUpdateRule() == 1) {
            updateFactor = targetTau;
        }
        else if (getUpdateRule() == 2) {
            updateFactor = -targetTau;
        }
        else if (getUpdateRule() == 3) {
            updateFactor = 1;
        }

        // Apply learning:
        for (int i = 0; i < getK(); i++) {
            float flux = pH[i];
            if (targetTau < 0) flux = 1 - pH[i];
            for (int j = 0; j < getN(); j++) {
                if (updateFactor * getX()[i][j] > 0) {
                    float temp = 0;
                    for (int k = 0; k < 2 * getL(); k++) {
                        float temp_ = flux * pW[i][j][k];
                        pW[i][j][k] += temp - flux * pW[i][j][k];
                        temp = temp_;
                    }
                    pW[i][j][2 * getL()] += temp;
                }
                else {
                    float temp = 0;
                    for (int k = 2 * getL(); k > 0; k--) {
                        float temp_ = flux * pW[i][j][k];
                        pW[i][j][k] += temp - flux * pW[i][j][k];
                        temp = temp_;
                    }
                    pW[i][j][0] += temp;
                }
            }
        }

        return;
    }

    void P_TPM::resetPW() {
        /* resetPW function:

            This function resets the weights corresponding to the input with lower sig2W. This
            means that it resets those pW which have already converged more. The idea behind
            this is that if a reset is needed it is because a wrong collapse exists.

            ACHTUNG!! Preliminary simulations show that this function destroys too many
            information. It must be rewritten!!

        */

        computeAvW();
        computeSig2W();

        int ind1 = 0, ind2 = 0;
        float minSig2 = sig2W[0][0];
        for (int i = 0; i < getK(); i++) {
            for (int j = 0; j < getN(); j++) {
                // ACHTUNG!! if (sig2W[i][j]>minSig2){
                if (sig2W[i][j] < minSig2) {
                    minSig2 = sig2W[i][j];
                    ind1 = i, ind2 = j;
                }
            }
        }

        // for (int k=0; k<2*getL()+1; k++){
        //	 pW[ind1][ind2][k] = 1./(2*getL()+1);
        // }
        float sum = 0;
        for (int k = 0; k < 2 * getL() + 1; k++) {
            pW[ind1][ind2][k] = float(rand());
        }
        normPW(ind1, ind2);

        return;
    }

    void P_TPM::preferMPW(int i = -1, int j = -1) {
        /* preferMPW function:

            This function reshapes the probability distribution to approach a Gaussian bell
            centered in the most probable weight. If i and j are lower than 0 then all pW are so
            reshaped. Else, the transformation affects only (i,j).

        Inputs:
            >> i=-1,j=-1: indexes indicating which units and weights are affected by the
            reshaping.

        */

        float z1, z2;
        computeAvW();
        computeSig2W();
        int **mPW;
        mPW = getMostProbW();

        if (i < 0 and j < 0) {
            for (int i = 0; i < getK(); i++) {
                for (int j = 0; j < getN(); j++) {
                    z1 = (-100 * getL() - avW[i][j]) / sqrt(sig2W[i][j]);
                    z2 = (-getL() + 0.5 - avW[i][j]) / sqrt(sig2W[i][j]);
                    pW[i][j][0] = erfc(z1 / sqrt(2)) / 2 - erfc(z2 / sqrt(2)) / 2;
                    for (int k = 1; k < 2 * getL(); k++) {
                        z1 = (k - getL() - 0.5 - avW[i][j]) / sqrt(sig2W[i][j]);
                        z2 = (k - getL() + 0.5 - avW[i][j]) / sqrt(sig2W[i][j]);
                        pW[i][j][k] = erfc(z1 / sqrt(2)) / 2 - erfc(z2 / sqrt(2)) / 2;
                    }
                    z1 = (getL() - 0.5 - avW[i][j]) / sqrt(sig2W[i][j]);
                    pW[i][j][2 * getL()] = 1. / 2 + erfc(z1 / sqrt(2)) / 2;
                }
            }
            normPW(-1, -1);
        }
        else {
            z1 = (-getL() - 0.5 - mPW[i][j]) / sqrt(sig2W[i][j]);
            z2 = (-getL() + 0.5 - mPW[i][j]) / sqrt(sig2W[i][j]);
            pW[i][j][0] = erfc(z1 / sqrt(2)) / 2 - erfc(z2 / sqrt(2)) / 2;
            for (int k = 1; k < 2 * getL(); k++) {
                z1 = (k - getL() - 0.5 - mPW[i][j]) / sqrt(sig2W[i][j]);
                z2 = (k - getL() + 0.5 - mPW[i][j]) / sqrt(sig2W[i][j]);
                pW[i][j][k] = erfc(z1 / sqrt(2)) / 2 - erfc(z2 / sqrt(2)) / 2;
            }
            z1 = (getL() - 0.5 - mPW[i][j]) / sqrt(sig2W[i][j]);
            z2 = (getL() + 0.5 - mPW[i][j]) / sqrt(sig2W[i][j]);
            pW[i][j][2 * getL()] = erfc(z1 / sqrt(2)) / 2 - erfc(z2 / sqrt(2)) / 2;
            normPW(i, j);
        }

        return;
    }
}