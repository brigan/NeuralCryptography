//
// Created by Nurlan Akashayev on 12/12/15.
//

#include <list>
#include "gtest/gtest.h"
#include "../../TPM.h"



// pIn - входной массив
// N - размер входного массива
// K - количество элементов в размещении

void PermutationWithRepetition(const int* pIn, int N, int K)
{
    int* pOut = new int[K + 1]; // строка из K символов плюс 1 символ для терминального 0
    pOut[K] = 0;                  // помещаем 0 в конец строки
    K--;
    int *stack = new int[K * 2],  // стек псевдорекурсии, глубина рекурсии K - 1
            *pTop = stack,            // вершина стека
            k = 0,                    // переменные цикла
            n = 0,
            j = 0;
    for (;;)                      // цикл псевдорекурсии
    {
        while(n < N)
        {
            pOut[k] = pIn[n++];
            if (k == K) {
                printf("%02d. ", ++j);
                for (int i = 0; i < k+1; i++) {
                    printf("%2d ", pOut[i]);
                }
                printf("\n");
            } else
            {
                if (n < N)
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
}

TEST(genetic_attack_check, tpm_deletion) {

    std::list<TPM::TPM> tpms;

    for (int i = 0; i < 5; i++) {
        TPM::TPM tpm = TPM::TPM(i, i, i, i);
        tpms.push_back(tpm);
    }

    std::list<TPM::TPM>::iterator i = tpms.begin();

    while (i != tpms.end()) {
        if (i->getK() == 3) {
            i = tpms.erase(i);
        } else {
            i++;
        }
    }

    EXPECT_EQ(4, tpms.size());
}

TEST(genetic_attack_check, permutation_with_repetition) {
    int charSet[2] = {-1, 1};
    int* ptr = charSet;

    PermutationWithRepetition(ptr, 2, 3);

    EXPECT_EQ(1, 1);
}