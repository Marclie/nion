/*
Copyright 2022 Marcus Dante Liebenthal

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
#if TEST_LARGE_DEGREE  == ON
#ifndef NION_HARDERTEST_HPP
#define NION_HARDERTEST_HPP






#include <fstream>
#include <chrono>
#include <random>
#include "Sedenion.h"
#include "nion.hpp"
#include <iostream>

template<typename T>
void printSpeedupSedenion(const T& niontime, const T& othertime) {
    // ANSI escape codes for colors
    const std::string red = "\033[0;31m";
    const std::string green = "\033[0;32m";
    const std::string reset = "\033[0m";

    T speedup = static_cast<T>(niontime) / static_cast<T>(othertime);
    if (speedup < 1)
        std::cout << "nion is " << green << 1 / speedup << " times FASTER" << reset << std::endl;
    else
        std::cout << "nion is " << red << speedup << " times SLOWER" << reset << std::endl;
}

template<typename T>
T getMAEsedenion(nion<T> nion, Sedenion<T> compare){
    int degree = nion.degree;
    T mae = 0;
    for (int i = 0; i < degree; i++) {
        mae += std::pow(nion[i] - compare[i], 2);
    }

    return std::sqrt(mae);
}

template <typename T>
void SedenionComparison(int trials){

    std::cout << "\n\n#### Comparing Sedenion with nion: " << trials << " trials ####\n" << std::endl;
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());

    // timers for norms
    long double normNionTimer = 0;
    long double normSedenionTimer = 0;
    T MRE_Norm = 0;
    T MAE_Norm = 0;
    T MAX_Norm = -1;
    T maxNormNion;
    T maxNormSedenion;
    Sedenion<T> max1value;
    Sedenion<T> max2value;

    // timers for addition
    long double addNionTimer = 0;
    long double addSedenionTimer = 0;
    T MRE_Add = 0;
    T MAE_Add= 0;
    T MAX_Add = -1;
    nion<T> maxAddNion;
    Sedenion<T> maxAddSedenion;
    Sedenion<T> max1valueAdd;
    Sedenion<T> max2valueAdd;

    // timers for conjugate
    long double conjNionTimer = 0;
    long double conjSedenionTimer = 0;
    T MRE_Conj = 0;
    T MAE_Conj= 0;
    T MAX_Conj = -1;
    nion<T> maxConjNion;
    Sedenion<T> maxConjSedenion;
    Sedenion<T> max1valueConj;
    Sedenion<T> max2valueConj;

    // timers for multiplication
    long double mulNionTimer = 0;
    long double mulSedenionTimer = 0;
    T MRE_Mul = 0;
    T MAE_Mul= 0;
    T MAX_Mul = -1;
    nion<T> maxMulNion;
    Sedenion<T> maxMulSedenion;
    Sedenion<T> max1valueMul;
    Sedenion<T> max2valueMul;

    // timers for division
    long double divNionTimer = 0;
    long double divSedenionTimer = 0;
    T MRE_Div = 0;
    T MAE_Div= 0;
    T MAX_Div = -1;
    nion<T> maxDivNion;
    Sedenion<T> maxDivSedenion;
    Sedenion<T> max1valueDiv;
    Sedenion<T> max2valueDiv;

    // timers for power
    long double powNionTimer = 0;
    long double powSedenionTimer = 0;
    T MRE_Pow = 0;
    T MAE_Pow= 0;
    T MAX_Pow = -1;
    nion<T> maxPowNion;
    Sedenion<T> maxPowSedenion;
    Sedenion<T> max1valuePow;
    Sedenion<T> max2valuePow;

    // timers for exponential
    long double expNionTimer = 0;
    long double expSedenionTimer = 0;
    T MRE_Exp = 0;
    T MAE_Exp= 0;
    T MAX_Exp = -1;
    nion<T> maxExpNion;
    Sedenion<T> maxExpSedenion;
    Sedenion<T> max1valueExp;
    Sedenion<T> max2valueExp;

    // timers for sine
    long double sinNionTimer = 0;
    long double sinSedenionTimer = 0;
    T MRE_Sin = 0;
    T MAE_Sin= 0;
    T MAX_Sin = -1;
    nion<T> maxSinNion;
    Sedenion<T> maxSinSedenion;
    Sedenion<T> max1valueSin;
    Sedenion<T> max2valueSin;

    // timers for cosine
    long double cosNionTimer = 0;
    long double cosSedenionTimer = 0;
    T MRE_Cos = 0;
    T MAE_Cos= 0;
    T MAX_Cos = -1;
    nion<T> maxCosNion;
    Sedenion<T> maxCosSedenion;
    Sedenion<T> max1valueCos;
    Sedenion<T> max2valueCos;

    // timers for tangent
    long double tanNionTimer = 0;
    long double tanSedenionTimer = 0;
    T MRE_Tan = 0;
    T MAE_Tan= 0;
    T MAX_Tan = -1;
    nion<T> maxTanNion;
    Sedenion<T> maxTanSedenion;
    Sedenion<T> max1valueTan;
    Sedenion<T> max2valueTan;

    // timers for sinh
    long double sinhNionTimer = 0;
    long double sinhSedenionTimer = 0;
    T MRE_Sinh = 0;
    T MAE_Sinh= 0;
    T MAX_Sinh = -1;
    nion<T> maxSinhNion;
    Sedenion<T> maxSinhSedenion;
    Sedenion<T> max1valueSinh;
    Sedenion<T> max2valueSinh;

    // timers for cosh
    long double coshNionTimer = 0;
    long double coshSedenionTimer = 0;
    T MRE_Cosh = 0;
    T MAE_Cosh= 0;
    T MAX_Cosh = -1;
    nion<T> maxCoshNion;
    Sedenion<T> maxCoshSedenion;
    Sedenion<T> max1valueCosh;
    Sedenion<T> max2valueCosh;

    // timers for tanh
    long double tanhNionTimer = 0;
    long double tanhSedenionTimer = 0;
    T MRE_Tanh = 0;
    T MAE_Tanh= 0;
    T MAX_Tanh = -1;
    nion<T> maxTanhNion;
    Sedenion<T> maxTanhSedenion;
    Sedenion<T> max1valueTanh;
    Sedenion<T> max2valueTanh;

    auto startNion = std::chrono::high_resolution_clock::now();
    auto endNion = std::chrono::high_resolution_clock::now();
    auto startSedenion = std::chrono::high_resolution_clock::now();
    auto endSedenion = std::chrono::high_resolution_clock::now();

    nion < T > nionResult;
    Sedenion<T> sedenionResult;
    T diff;

    T *vals1 = new T[16];
    T *vals2 = new T[16];
    std::uniform_real_distribution<T> distribution(-5, 5);
    for (int i = 0; i < trials; ++i) {

        //generate random quaternion numbers
        for (int j = 0; j < 16; ++j) {
            vals1[j] = distribution(generator);
            vals2[j] = distribution(generator);
        }

        nion<T> nion1(vals1, 16);
        nion<T> nion2(vals2, 16);
        Sedenion<T> sedenion1(vals1[0], vals1[1], vals1[2], vals1[3], vals1[4], vals1[5], vals1[6], vals1[7],
                              vals1[8], vals1[9], vals1[10], vals1[11], vals1[12], vals1[13], vals1[14], vals1[15]);
        Sedenion<T> sedenion2(vals2[0], vals2[1], vals2[2], vals2[3], vals2[4], vals2[5], vals2[6], vals2[7],
                              vals2[8], vals2[9], vals2[10], vals2[11], vals2[12], vals2[13], vals2[14], vals2[15]);

        // norm
        {
            startNion = std::chrono::high_resolution_clock::now();
            T nionNorm = nion1.abs();
            endNion = std::chrono::high_resolution_clock::now();
            normNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            T sedenionNorm = norm(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            normSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = fabs(nionNorm - sedenionNorm);
            if (diff > MAX_Norm) {
                MAX_Norm = diff;
                maxNormNion = nionNorm;
                maxNormSedenion = sedenionNorm;
                max1value = sedenion1;
            }
            MAE_Norm += diff;
            MRE_Norm += diff / sedenionNorm;
        }

        // conjugate
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1.conj();
            endNion = std::chrono::high_resolution_clock::now();
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = conj(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            conjSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Conj) {
                MAX_Conj = diff;
                maxConjNion = nionResult;
                maxConjSedenion = sedenionResult;
                max1valueConj = sedenion1;
            }
            MAE_Conj += diff;
            MRE_Conj += diff / norm(sedenionResult);
        }

        // add
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 + nion2;
            endNion = std::chrono::high_resolution_clock::now();
            addNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = sedenion1 + sedenion2;
            endSedenion = std::chrono::high_resolution_clock::now();
            addSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Add) {
                MAX_Add = diff;
                maxAddNion = nionResult;
                maxAddSedenion = sedenionResult;
                max1valueAdd = sedenion1;
                max2valueAdd = sedenion2;
            }
            MAE_Add += diff;
            MRE_Add += diff / norm(sedenionResult);
        }

        // multiply
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 * nion2;
            endNion = std::chrono::high_resolution_clock::now();
            mulNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = sedenion1 * sedenion2;
            endSedenion = std::chrono::high_resolution_clock::now();
            mulSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Mul) {
                MAX_Mul = diff;
                maxMulNion = nionResult;
                maxMulSedenion = sedenionResult;
                max1valueMul = sedenion1;
                max2valueMul = sedenion2;
            }
            MAE_Mul += diff;
            MRE_Mul += diff / norm(sedenionResult);
        }

        // divide
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 / nion2;
            endNion = std::chrono::high_resolution_clock::now();
            divNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = sedenion1 / sedenion2;
            endSedenion = std::chrono::high_resolution_clock::now();
            divSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Div) {
                MAX_Div = diff;
                maxDivNion = nionResult;
                maxDivSedenion = sedenionResult;
                max1valueDiv = sedenion1;
                max2valueDiv = sedenion2;
            }
            MAE_Div += diff;
            MRE_Div += diff / norm(sedenionResult);
        }

        // exponential
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = exp(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            expNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = exp(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            expSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Exp) {
                MAX_Exp = diff;
                maxExpNion = nionResult;
                maxExpSedenion = sedenionResult;
                max1valueExp = sedenion1;
            }
            MAE_Exp += diff;
            MRE_Exp += diff / norm(sedenionResult);
        }

        // pow
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = pow(nion1, 3);
            endNion = std::chrono::high_resolution_clock::now();
            powNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = pow(sedenion1, static_cast<T>(3));
            endSedenion = std::chrono::high_resolution_clock::now();
            powSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Pow) {
                MAX_Pow = diff;
                maxPowNion = nionResult;
                maxPowSedenion = sedenionResult;
                max1valuePow = sedenion1;
                max2valuePow = sedenion2;
            }
            MAE_Pow += diff;
            MRE_Pow += diff / norm(sedenionResult);
        }

        // sine
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sin(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = sin(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            sinSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Sin) {
                MAX_Sin = diff;
                maxSinNion = nionResult;
                maxSinSedenion = sedenionResult;
                max1valueSin = sedenion1;
            }
            MAE_Sin += diff;
            MRE_Sin += diff / norm(sedenionResult);
        }

        // cosine
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cos(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            cosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = cos(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            cosSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Cos) {
                MAX_Cos = diff;
                maxCosNion = nionResult;
                maxCosSedenion = sedenionResult;
                max1valueCos = sedenion1;
            }
            MAE_Cos += diff;
            MRE_Cos += diff / norm(sedenionResult);
        }

        // tangent
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tan(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = tan(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            tanSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Tan) {
                MAX_Tan = diff;
                maxTanNion = nionResult;
                maxTanSedenion = sedenionResult;
                max1valueTan = sedenion1;
            }
            MAE_Tan += diff;
            MRE_Tan += diff / norm(sedenionResult);
        }

        // hyperbolic sine
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sinh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = sinh(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            sinhSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Sinh) {
                MAX_Sinh = diff;
                maxSinhNion = nionResult;
                maxSinhSedenion = sedenionResult;
                max1valueSinh = sedenion1;
            }
            MAE_Sinh += diff;
            MRE_Sinh += diff / norm(sedenionResult);
        }

        // hyperbolic cosine
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cosh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            coshNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = cosh(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            coshSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Cosh) {
                MAX_Cosh = diff;
                maxCoshNion = nionResult;
                maxCoshSedenion = sedenionResult;
                max1valueCosh = sedenion1;
            }
            MAE_Cosh += diff;
            MRE_Cosh += diff / norm(sedenionResult);
        }

        // hyperbolic tangent
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tanh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenionResult = tanh(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            tanhSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            diff = getMAEsedenion(nionResult, sedenionResult);
            if (diff > MAX_Tanh) {
                MAX_Tanh = diff;
                maxTanhNion = nionResult;
                maxTanhSedenion = sedenionResult;
                max1valueTanh = sedenion1;
            }
            MAE_Tanh += diff;
            MRE_Tanh += diff / norm(sedenionResult);
        }

    }

    delete[] vals1;
    delete[] vals2;

/*** Norm ***/

    std::cout << "\n\n ---> Norm <---" << std::endl;
    std::cout << "Average norm time for nion: " << normNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average norm time for sedenion: " << normSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(normNionTimer, normSedenionTimer);
    std::cout << "Average norm error for nion: " << MAE_Norm / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Norm / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum norm error for nion: " << MAX_Norm << std::endl;
    std::cout << "nion: " << maxNormNion << "\nsedenion: " << maxNormSedenion << std::endl;
    std::cout << "input: " << max1value << std::endl;

/*** Conjugation ***/

    std::cout << "\n\n ---> Conjugation <---" << std::endl;
    std::cout << "Average conjugation time for nion: " << conjNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average conjugation time for sedenion: " << conjSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(conjNionTimer, conjSedenionTimer);
    std::cout << "Average conjugation error for nion: " << MAE_Conj / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Conj / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum conjugation error for nion: " << MAX_Conj << std::endl;
    std::cout << "nion: " << maxConjNion << "\nsedenion: " << maxConjSedenion << std::endl;
    std::cout << "input: " << max1valueConj << std::endl;

/*** Addition ***/

    std::cout << "\n\n ---> Addition <---" << std::endl;
    std::cout << "Average addition time for nion: " << addNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average addition time for sedenion: " << addSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(addNionTimer, addSedenionTimer);
    std::cout << "Average addition error for nion: " << MAE_Add / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Add / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum addition error for nion: " << MAX_Add << std::endl;
    std::cout << "nion: " << maxAddNion << "\nsedenion: " << maxAddSedenion << std::endl;
    std::cout << "input1: " << max1valueAdd << std::endl;
    std::cout << "input2: " << max2valueAdd << std::endl;

/*** Multiplication ***/

    std::cout << "\n\n ---> Multiplication <---" << std::endl;
    std::cout << "Average multiplication time for nion: " << mulNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average multiplication time for sedenion: " << mulSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(mulNionTimer, mulSedenionTimer);
    std::cout << "Average multiplication error for nion: " << MAE_Mul / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Mul / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum multiplication error for nion: " << MAX_Mul << std::endl;
    std::cout << "nion: " << maxMulNion << "\nsedenion: " << maxMulSedenion << std::endl;
    std::cout << "input1: " << max1valueMul << std::endl;
    std::cout << "input2: " << max2valueMul << std::endl;

/*** Division ***/

    std::cout << "\n\n ---> Division <---" << std::endl;
    std::cout << "Average division time for nion: " << divNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average division time for sedenion: " << divSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(divNionTimer, divSedenionTimer);
    std::cout << "Average division error for nion: " << MAE_Div / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Div / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum division error for nion: " << MAX_Div << std::endl;
    std::cout << "nion: " << maxDivNion << "\nsedenion: " << maxDivSedenion << std::endl;
    std::cout << "input1: " << max1valueDiv << std::endl;
    std::cout << "input2: " << max2valueDiv << std::endl;

/*** Exponential ***/

    std::cout << "\n\n ---> Exponential <---" << std::endl;
    std::cout << "Average exponential time for nion: " << expNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average exponential time for sedenion: " << expSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(expNionTimer, expSedenionTimer);
    std::cout << "Average exponential error for nion: " << MAE_Exp / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Exp / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum exponential error for nion: " << MAX_Exp << std::endl;
    std::cout << "nion: " << maxExpNion << "\nsedenion: " << maxExpSedenion << std::endl;
    std::cout << "input: " << max1valueExp << std::endl;

/*** Power ***/

    std::cout << "\n\n ---> Power <---" << std::endl;
    std::cout << "Average power time for nion: " << powNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average power time for sedenion: " << powSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(powNionTimer, powSedenionTimer);
    std::cout << "Average power error for nion: " << MAE_Pow / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Pow / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum power error for nion: " << MAX_Pow << std::endl;
    std::cout << "nion: " << maxPowNion << "\nsedenion: " << maxPowSedenion << std::endl;
    std::cout << "input: " << max1valuePow << std::endl;

/*** Sine ***/

    std::cout << "\n\n ---> Sine <---" << std::endl;
    std::cout << "Average sine time for nion: " << sinNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average sine time for sedenion: " << sinSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(sinNionTimer, sinSedenionTimer);
    std::cout << "Average sine error for nion: " << MAE_Sin / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Sin / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum sine error for nion: " << MAX_Sin << std::endl;
    std::cout << "nion: " << maxSinNion << "\nsedenion: " << maxSinSedenion << std::endl;
    std::cout << "input: " << max1valueSin << std::endl;

/*** Cosine ***/

    std::cout << "\n\n ---> Cosine <---" << std::endl;
    std::cout << "Average cosine time for nion: " << cosNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average cosine time for sedenion: " << cosSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(cosNionTimer, cosSedenionTimer);
    std::cout << "Average cosine error for nion: " << MAE_Cos / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Cos / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum cosine error for nion: " << MAX_Cos << std::endl;
    std::cout << "nion: " << maxCosNion << "\nsedenion: " << maxCosSedenion << std::endl;
    std::cout << "input: " << max1valueCos << std::endl;

/*** Tangent ***/

    std::cout << "\n\n ---> Tangent <---" << std::endl;
    std::cout << "Average tangent time for nion: " << tanNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average tangent time for sedenion: " << tanSedenionTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupSedenion(tanNionTimer, tanSedenionTimer);
    std::cout << "Average tangent error for nion: " << MAE_Tan / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Tan / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum tangent error for nion: " << MAX_Tan << std::endl;
    std::cout << "nion: " << maxTanNion << "\nsedenion: " << maxTanSedenion << std::endl;
    std::cout << "input: " << max1valueTan << std::endl;

/*** Hyperbolic Sine ***/

    std::cout << "\n\n ---> Hyperbolic Sine <---" << std::endl;
    std::cout << "Average hyperbolic sine time for nion: " << sinhNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average hyperbolic sine time for sedenion: " << sinhSedenionTimer / static_cast<T>(trials) << " ns"
              << std::endl;
    printSpeedupSedenion(sinhNionTimer, sinhSedenionTimer);
    std::cout << "Average hyperbolic sine error for nion: " << MAE_Sinh / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Sinh / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum hyperbolic sine error for nion: " << MAX_Sinh << std::endl;
    std::cout << "nion: " << maxSinhNion << "\nsedenion: " << maxSinhSedenion << std::endl;
    std::cout << "input: " << max1valueSinh << std::endl;

/*** Hyperbolic Cosine ***/

    std::cout << "\n\n ---> Hyperbolic Cosine <---" << std::endl;
    std::cout << "Average hyperbolic cosine time for nion: " << coshNionTimer / static_cast<T>(trials) << " ns"
              << std::endl;
    std::cout << "Average hyperbolic cosine time for sedenion: " << coshSedenionTimer / static_cast<T>(trials) << " ns"
              << std::endl;
    printSpeedupSedenion(coshNionTimer, coshSedenionTimer);
    std::cout << "Average hyperbolic cosine error for nion: " << MAE_Cosh / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Cosh / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum hyperbolic cosine error for nion: " << MAX_Cosh << std::endl;
    std::cout << "nion: " << maxCoshNion << "\nsedenion: " << maxCoshSedenion << std::endl;
    std::cout << "input: " << max1valueCosh << std::endl;

/*** Hyperbolic Tangent ***/

    std::cout << "\n\n ---> Hyperbolic Tangent <---" << std::endl;
    std::cout << "Average hyperbolic tangent time for nion: " << tanhNionTimer / static_cast<T>(trials) << " ns"
              << std::endl;
    std::cout << "Average hyperbolic tangent time for sedenion: " << tanhSedenionTimer / static_cast<T>(trials) << " ns"
              << std::endl;
    printSpeedupSedenion(tanhNionTimer, tanhSedenionTimer);
    std::cout << "Average hyperbolic tangent error for nion: " << MAE_Tanh / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and sedenion: " << MRE_Tanh / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum hyperbolic tangent error for nion: " << MAX_Tanh << std::endl;
    std::cout << "nion: " << maxTanhNion << "\nsedenion: " << maxTanhSedenion << std::endl;
    std::cout << "input: " << max1valueTanh << std::endl;

}


#endif //NION_HARDERTEST_HPP
#endif