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

#ifdef TEST_LARGE_DEGREE
#ifndef NION_TRIGINTADUONIONTEST_HPP
#define NION_TRIGINTADUONIONTEST_HPP

#include <fstream>
#include <chrono>
#include <random>
#include "Trigintaduonion.h"
#include "../nion.hpp"
#include <iostream>


template<typename T>
void printSpeedupTrigintaduonion(const T& niontime, const T& othertime) {
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
T getMAEtrigintaduonion(nion<T,32> nion, Trigintaduonion<T> compare){
    int degree = nion.size();
    T mae = 0;
    for (int i = 0; i < degree; i++) {
        mae += std::pow(nion[i] - compare[i], 2);
    }

    return std::sqrt(mae);
}

template<typename T>
void TrigintaduonionComparison(int trials) {
    std::cout << "\n\n#### Comparing Trigintaduonion with nion: " << trials << " trials ####\n" << std::endl;
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());

    // timers for norms
    T normNionTimer = 0;
    T normTrigintaduonionTimer = 0;
    T MAE_Norm = 0;
    T MRE_Norm = 0;
    T MAX_Norm = -1;
    T maxNormNion = 0;
    T maxNormTrigintaduonion = 0;
    Trigintaduonion<T> max1valueNorm;


    // timers for addition
    T addNionTimer = 0;
    T addTrigintaduonionTimer = 0;
    T MAE_Add = 0;
    T MRE_Add = 0;
    T MAX_Add = -1;
    nion<T,32> maxAddNion;
    Trigintaduonion<T> maxAddTrigintaduonion;
    Trigintaduonion<T> max1valueAdd;
    Trigintaduonion<T> max2valueAdd;

    // timers for conjugate
    T conjNionTimer = 0;
    T conjTrigintaduonionTimer = 0;
    T MAE_Conj = 0;
    T MRE_Conj = 0;
    T MAX_Conj = -1;
    nion<T,32> maxConjNion;
    Trigintaduonion<T> maxConjTrigintaduonion;
    Trigintaduonion<T> max1valueConj;
    Trigintaduonion<T> max2valueConj;

    // timers for multiplication
    T mulNionTimer = 0;
    T mulTrigintaduonionTimer = 0;
    T MAE_Mul = 0;
    T MRE_Mul = 0;
    T MAX_Mul = -1;
    nion<T,32> maxMulNion;
    Trigintaduonion<T> maxMulTrigintaduonion;
    Trigintaduonion<T> max1valueMul;
    Trigintaduonion<T> max2valueMul;

    // timers for division
    T divNionTimer = 0;
    T divTrigintaduonionTimer = 0;
    T MAE_Div = 0;
    T MRE_Div = 0;
    T MAX_Div = -1;
    nion<T,32> maxDivNion;
    Trigintaduonion<T> maxDivTrigintaduonion;
    Trigintaduonion<T> max1valueDiv;
    Trigintaduonion<T> max2valueDiv;

    // timers for power
    T powNionTimer = 0;
    T powTrigintaduonionTimer = 0;
    T MAE_Pow = 0;
    T MRE_Pow = 0;
    T MAX_Pow = -1;
    nion<T,32> maxPowNion;
    Trigintaduonion<T> maxPowTrigintaduonion;
    Trigintaduonion<T> max1valuePow;
    Trigintaduonion<T> max2valuePow;

    // timers for square
    T SqNionTimer = 0;
    T SqTrigintaduonionTimer = 0;
    T MAE_Sq = 0;
    T MRE_Sq = 0;
    T MAX_Sq = -1;
    nion<T,32> maxSqNion;
    Trigintaduonion<T> maxSqTrigintaduonion;
    Trigintaduonion<T> max1valueSq;
    Trigintaduonion<T> max2valueSq;

    // timers for square root
    T SqrtNionTimer = 0;
    T SqrtTrigintaduonionTimer = 0;
    T MAE_Sqrt = 0;
    T MRE_Sqrt = 0;
    T MAX_Sqrt = -1;
    nion<T,32> maxSqrtNion;
    Trigintaduonion<T> maxSqrtTrigintaduonion;
    Trigintaduonion<T> max1valueSqrt;
    Trigintaduonion<T> max2valueSqrt;

    // timers for exponential
    T expNionTimer = 0;
    T expTrigintaduonionTimer = 0;
    T MAE_Exp = 0;
    T MRE_Exp = 0;
    T MAX_Exp = -1;
    nion<T,32> maxExpNion;
    Trigintaduonion<T> maxExpTrigintaduonion;
    Trigintaduonion<T> max1valueExp;
    Trigintaduonion<T> max2valueExp;

    // timers for logarithm
    T logNionTimer = 0;
    T logTrigintaduonionTimer = 0;
    T MAE_Log = 0;
    T MRE_Log = 0;
    T MAX_Log = -1;
    nion<T,32> maxLogNion;
    Trigintaduonion<T> maxLogTrigintaduonion;
    Trigintaduonion<T> max1valueLog;
    Trigintaduonion<T> max2valueLog;

    // timers for sin
    T sinNionTimer = 0;
    T sinTrigintaduonionTimer = 0;
    T MAE_Sin = 0;
    T MRE_Sin = 0;
    T MAX_Sin = -1;
    nion<T,32> maxSinNion;
    Trigintaduonion<T> maxSinTrigintaduonion;
    Trigintaduonion<T> max1valueSin;
    Trigintaduonion<T> max2valueSin;

    // timers for asin
    T asinNionTimer = 0;
    T asinTrigintaduonionTimer = 0;
    T MAE_Asin = 0;
    T MRE_Asin = 0;
    T MAX_Asin = -1;
    nion<T,32> maxAsinNion;
    Trigintaduonion<T> maxAsinTrigintaduonion;
    Trigintaduonion<T> max1valueAsin;
    Trigintaduonion<T> max2valueAsin;

    // timers for cos
    T cosNionTimer = 0;
    T cosTrigintaduonionTimer = 0;
    T MAE_Cos = 0;
    T MRE_Cos = 0;
    T MAX_Cos = -1;
    nion<T,32> maxCosNion;
    Trigintaduonion<T> maxCosTrigintaduonion;
    Trigintaduonion<T> max1valueCos;
    Trigintaduonion<T> max2valueCos;

    // timers for acos
    T acosNionTimer = 0;
    T acosTrigintaduonionTimer = 0;
    T MAE_Acos = 0;
    T MRE_Acos = 0;
    T MAX_Acos = -1;
    nion<T,32> maxAcosNion;
    Trigintaduonion<T> maxAcosTrigintaduonion;
    Trigintaduonion<T> max1valueAcos;
    Trigintaduonion<T> max2valueAcos;

    // timers for tan
    T tanNionTimer = 0;
    T tanTrigintaduonionTimer = 0;
    T MAE_Tan = 0;
    T MRE_Tan = 0;
    T MAX_Tan = -1;
    nion<T,32> maxTanNion;
    Trigintaduonion<T> maxTanTrigintaduonion;
    Trigintaduonion<T> max1valueTan;
    Trigintaduonion<T> max2valueTan;

    // timers for atan
    T atanNionTimer = 0;
    T atanTrigintaduonionTimer = 0;
    T MAE_Atan = 0;
    T MRE_Atan = 0;
    T MAX_Atan = -1;
    nion<T,32> maxAtanNion;
    Trigintaduonion<T> maxAtanTrigintaduonion;
    Trigintaduonion<T> max1valueAtan;
    Trigintaduonion<T> max2valueAtan;

    // timers for sinh
    T sinhNionTimer = 0;
    T sinhTrigintaduonionTimer = 0;
    T MAE_Sinh = 0;
    T MRE_Sinh = 0;
    T MAX_Sinh = -1;
    nion<T,32> maxSinhNion;
    Trigintaduonion<T> maxSinhTrigintaduonion;
    Trigintaduonion<T> max1valueSinh;
    Trigintaduonion<T> max2valueSinh;

    // timers for asinh
    T asinhNionTimer = 0;
    T asinhTrigintaduonionTimer = 0;
    T MAE_Asinh = 0;
    T MRE_Asinh = 0;
    T MAX_Asinh = -1;
    nion<T,32> maxAsinhNion;
    Trigintaduonion<T> maxAsinhTrigintaduonion;
    Trigintaduonion<T> max1valueAsinh;
    Trigintaduonion<T> max2valueAsinh;

    // timers for cosh
    T coshNionTimer = 0;
    T coshTrigintaduonionTimer = 0;
    T MAE_Cosh = 0;
    T MRE_Cosh = 0;
    T MAX_Cosh = -1;
    nion<T,32> maxCoshNion;
    Trigintaduonion<T> maxCoshTrigintaduonion;
    Trigintaduonion<T> max1valueCosh;
    Trigintaduonion<T> max2valueCosh;

    // timers for acosh
    T acoshNionTimer = 0;
    T acoshTrigintaduonionTimer = 0;
    T MAE_Acosh = 0;
    T MRE_Acosh = 0;
    T MAX_Acosh = -1;
    nion<T,32> maxAcoshNion;
    Trigintaduonion<T> maxAcoshTrigintaduonion;
    Trigintaduonion<T> max1valueAcosh;
    Trigintaduonion<T> max2valueAcosh;

    // timers for tanh
    T tanhNionTimer = 0;
    T tanhTrigintaduonionTimer = 0;
    T MAE_Tanh = 0;
    T MRE_Tanh = 0;
    T MAX_Tanh = -1;
    nion<T,32> maxTanhNion;
    Trigintaduonion<T> maxTanhTrigintaduonion;
    Trigintaduonion<T> max1valueTanh;
    Trigintaduonion<T> max2valueTanh;

    // timers for atanh
    T atanhNionTimer = 0;
    T atanhTrigintaduonionTimer = 0;
    T MAE_Atanh = 0;
    T MRE_Atanh = 0;
    T MAX_Atanh = -1;
    nion<T,32> maxAtanhNion;
    Trigintaduonion<T> maxAtanhTrigintaduonion;
    Trigintaduonion<T> max1valueAtanh;
    Trigintaduonion<T> max2valueAtanh;


    auto startNion = std::chrono::high_resolution_clock::now();
    auto endNion = std::chrono::high_resolution_clock::now();
    auto startTrigintaduonion = std::chrono::high_resolution_clock::now();
    auto endTrigintaduonion = std::chrono::high_resolution_clock::now();

    nion<T,32> nionResult;
    Trigintaduonion<T> trigintaduonion_result;
    T diff;

    std::uniform_real_distribution<T> distribution(-5.0, 5.0);
    T* val1 = new T[32];
    T* val2 = new T[32];
    for (int i = 0; i < trials; ++i) {
        // get random Sedenion number
        for (int j = 0; j < 32; ++j) {
            val1[j] = distribution(generator);
            val2[j] = distribution(generator);
        }
        nion<T,32> nion1(val1, 32);
        nion<T,32> nion2(val2, 32);

        Trigintaduonion<T> trigintaduonion1(val1[0], val1[1], val1[2], val1[3], val1[4], val1[5], val1[6], val1[7], val1[8], val1[9], val1[10], val1[11], val1[12], val1[13], val1[14], val1[15], val1[16], val1[17], val1[18], val1[19], val1[20], val1[21], val1[22], val1[23], val1[24], val1[25], val1[26], val1[27], val1[28], val1[29], val1[30], val1[31]);
        Trigintaduonion<T> trigintaduonion2(val2[0], val2[1], val2[2], val2[3], val2[4], val2[5], val2[6], val2[7], val2[8], val2[9], val2[10], val2[11], val2[12], val2[13], val2[14], val2[15], val2[16], val2[17], val2[18], val2[19], val2[20], val2[21], val2[22], val2[23], val2[24], val2[25], val2[26], val2[27], val2[28], val2[29], val2[30], val2[31]);


        // norm
        {
            // evaluate nion norm, and time
            startNion = std::chrono::high_resolution_clock::now();
            T nionNorm = nion1.abs();
            endNion = std::chrono::high_resolution_clock::now();
            normNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion norm, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            T stdNorm = norm(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            normTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ normdifference
            diff = std::abs(nionNorm - stdNorm);
            MRE_Norm +=  diff / norm(trigintaduonion_result);
            diff /= stdNorm;
            MAE_Norm += diff;

            // get max difference between nion and Trigintaduonion norms
            if (diff > MAX_Norm) {
                MAX_Norm = diff;
                maxNormNion = nionNorm;
                maxNormTrigintaduonion = stdNorm;
                max1valueNorm = trigintaduonion1;
            }
        }

        /// addition
        {
            // evaluate nion addition, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 + nion2;
            endNion = std::chrono::high_resolution_clock::now();
            addNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion addition, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = trigintaduonion1 + trigintaduonion2;
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            addTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ additiondifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Add +=  diff / norm(trigintaduonion_result);
            
            MAE_Add += diff;

            // get max difference between nion and Trigintaduonion addition
            if (diff > MAX_Add) {
                MAX_Add = diff;
                maxAddNion = nionResult;
                maxAddTrigintaduonion = trigintaduonion_result;
                max1valueAdd = trigintaduonion1;
                max2valueAdd = trigintaduonion2;
            }
        }

        /// conjugate
        {
            // evaluate nion conjugate, and time
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();
            nionResult = nion2.conj();
            endNion = std::chrono::high_resolution_clock::now();
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = conj(trigintaduonion2);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            conjTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ conjugatedifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Conj +=  diff / norm(trigintaduonion_result);
            
            MAE_Conj += diff;

            // get max difference between nion and Trigintaduonion conjugate
            if (diff > MAX_Conj) {
                MAX_Conj = diff;
                maxConjNion = nionResult;
                maxConjTrigintaduonion = trigintaduonion_result;
                max1valueConj = trigintaduonion1;
                max2valueConj = trigintaduonion2;
            }

        }

        /// multiplication
        {
            // evaluate nion multiplication, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 * nion2;
            endNion = std::chrono::high_resolution_clock::now();
            mulNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion multiplication, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = trigintaduonion1 * trigintaduonion2;
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            mulTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ multiplicationdifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Mul +=  diff / norm(trigintaduonion_result);
            
            MAE_Mul += diff;

            // get max difference between nion and Trigintaduonion multiplication
            if (diff > MAX_Mul) {
                MAX_Mul = diff;
                maxMulNion = nionResult;
                maxMulTrigintaduonion = trigintaduonion_result;
                max1valueMul = trigintaduonion1;
                max2valueMul = trigintaduonion2;
            }

        }

        /// division
        {
            // evaluate nion division, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 / nion2;
            endNion = std::chrono::high_resolution_clock::now();
            divNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion division, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = trigintaduonion1 / trigintaduonion2;
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            divTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ divisiondifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Div +=  diff / norm(trigintaduonion_result);
            
            MAE_Div += diff;

            // get max difference between nion and Trigintaduonion division
            if (diff > MAX_Div) {
                MAX_Div = diff;
                maxDivNion = nionResult;
                maxDivTrigintaduonion = trigintaduonion_result;
                max1valueDiv = trigintaduonion1;
                max2valueDiv = trigintaduonion2;
            }
        }

        /// power
        {
            // evaluate nion power, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = pow(nion1, static_cast<T>(5));
            endNion = std::chrono::high_resolution_clock::now();
            powNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion power, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = pow(trigintaduonion1, static_cast<T>(5));
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            powTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ powerdifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Pow +=  diff / norm(trigintaduonion_result);
            
            MAE_Pow += diff;

            // get max difference between nion and Trigintaduonion power
            if (diff > MAX_Pow) {
                MAX_Pow = diff;
                maxPowNion = nionResult;
                maxPowTrigintaduonion = trigintaduonion_result;
                max1valuePow = trigintaduonion1;
                max2valuePow = trigintaduonion2;
            }
        }

        /// square
        {
            // evaluate nion square, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sqr(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            SqNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion square, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = sqr(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            SqTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ squaredifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Sq +=  diff / norm(trigintaduonion_result);
            
            MAE_Sq += diff;

            // get max difference between nion and Trigintaduonion square
            if (diff > MAX_Sq) {
                MAX_Sq = diff;
                maxSqNion = nionResult;
                maxSqTrigintaduonion = trigintaduonion_result;
                max1valueSq = trigintaduonion1;
                max2valueSq = trigintaduonion2;
            }
        }

        /// square root
        {
            // evaluate nion square, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = pow(nion1, 0.5);
            endNion = std::chrono::high_resolution_clock::now();
            SqrtNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion square, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = pow(trigintaduonion1, 0.5l);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            SqrtTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ squaredifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Sqrt +=  diff / norm(trigintaduonion_result);
            
            MAE_Sqrt += diff;

            // get max difference between nion and Trigintaduonion square
            if (diff > MAX_Sqrt) {
                MAX_Sqrt = diff;
                maxSqrtNion = nionResult;
                maxSqrtTrigintaduonion = trigintaduonion_result;
                max1valueSqrt = trigintaduonion1;
                max2valueSqrt = trigintaduonion2;
            }
        }

        /// exponential
        {
            // evaluate nion exponential, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = exp(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            expNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion exponential, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = exp(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            expTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ exponentialdifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Exp +=  diff / norm(trigintaduonion_result);
            
            MAE_Exp += diff;

            // get max difference between nion and Trigintaduonion exponential
            if (diff > MAX_Exp) {
                MAX_Exp = diff;
                maxExpNion = nionResult;
                maxExpTrigintaduonion = trigintaduonion_result;
                max1valueExp = trigintaduonion1;
                max2valueExp = trigintaduonion2;
            }
        }

        /// logarithm
        {
            // evaluate nion logarithm, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = log(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            logNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion logarithm, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = log(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            logTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ logarithmdifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Log +=  diff / norm(trigintaduonion_result);
            
            MAE_Log += diff;

            // get max difference between nion and Trigintaduonion logarithm
            if (diff > MAX_Log) {
                MAX_Log = diff;
                maxLogNion = nionResult;
                maxLogTrigintaduonion = trigintaduonion_result;
                max1valueLog = trigintaduonion1;
                max2valueLog = trigintaduonion2;
            }
        }

        /// sine
        {
            // evaluate nion sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sin(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion sine, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = sin(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            sinTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ sinedifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Sin +=  diff / norm(trigintaduonion_result);
            
            MAE_Sin += diff;

            // get max difference between nion and Trigintaduonion sine
            if (diff > MAX_Sin) {
                MAX_Sin = diff;
                maxSinNion = nionResult;
                maxSinTrigintaduonion = trigintaduonion_result;
                max1valueSin = trigintaduonion1;
                max2valueSin = trigintaduonion2;
            }
        }

        /// asine
        {
            // evaluate nion asine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = asin(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            asinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion asine, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = asin(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            asinTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ asinedifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Asin +=  diff / norm(trigintaduonion_result);
            
            MAE_Asin += diff;

            // get max difference between nion and Trigintaduonion asine
            if (diff > MAX_Asin) {
                MAX_Asin = diff;
                maxAsinNion = nionResult;
                maxAsinTrigintaduonion = trigintaduonion_result;
                max1valueAsin = trigintaduonion1;
                max2valueAsin = trigintaduonion2;
            }
        }

        /// cosine
        {
            // evaluate nion cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cos(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            cosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion cosine, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = cos(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            cosTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ cosinedifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Cos +=  diff / norm(trigintaduonion_result);
            
            MAE_Cos += diff;

            // get max difference between nion and Trigintaduonion cosine
            if (diff > MAX_Cos) {
                MAX_Cos = diff;
                maxCosNion = nionResult;
                maxCosTrigintaduonion = trigintaduonion_result;
                max1valueCos = trigintaduonion1;
                max2valueCos = trigintaduonion2;
            }
        }

        /// acosine
        {
            // evaluate nion acosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = acos(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            acosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion acosine, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = acos(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            acosTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ acosinedifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Acos +=  diff / norm(trigintaduonion_result);
            
            MAE_Acos += diff;

            // get max difference between nion and Trigintaduonion acosine
            if (diff > MAX_Acos) {
                MAX_Acos = diff;
                maxAcosNion = nionResult;
                maxAcosTrigintaduonion = trigintaduonion_result;
                max1valueAcos = trigintaduonion1;
                max2valueAcos = trigintaduonion2;
            }
        }

        /// tangent
        {
            // evaluate nion tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tan(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion tangent, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = tan(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            tanTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ tangentdifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Tan +=  diff / norm(trigintaduonion_result);
            
            MAE_Tan += diff;

            // get max difference between nion and Trigintaduonion tangent
            if (diff > MAX_Tan) {
                MAX_Tan = diff;
                maxTanNion = nionResult;
                maxTanTrigintaduonion = trigintaduonion_result;
                max1valueTan = trigintaduonion1;
                max2valueTan = trigintaduonion2;
            }
        }

        /// atan
        {
            // evaluate nion atan, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = atan(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            atanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion atan, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = atan(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            atanTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ atandifference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Atan +=  diff / norm(trigintaduonion_result);
            
            MAE_Atan += diff;

            // get max difference between nion and Trigintaduonion atan
            if (diff > MAX_Atan) {
                MAX_Atan = diff;
                maxAtanNion = nionResult;
                maxAtanTrigintaduonion = trigintaduonion_result;
                max1valueAtan = trigintaduonion1;
                max2valueAtan = trigintaduonion2;
            }
        }

        /// hyperbolic sine
        {
            // evaluate nion hyperbolic sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sinh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion hyperbolic sine, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = sinh(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            sinhTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ hyperbolicsine difference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Sinh +=  diff / norm(trigintaduonion_result);
            
            MAE_Sinh += diff;

            // get max difference between nion and Trigintaduonion hyperbolic sine
            if (diff > MAX_Sinh) {
                MAX_Sinh = diff;
                maxSinhNion = nionResult;
                maxSinhTrigintaduonion = trigintaduonion_result;
                max1valueSinh = trigintaduonion1;
                max2valueSinh = trigintaduonion2;
            }
        }

        /// hyperbolic cosine
        {
            // evaluate nion hyperbolic cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cosh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            coshNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion hyperbolic cosine, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = cosh(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            coshTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ hyperboliccosine difference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Cosh +=  diff / norm(trigintaduonion_result);
            
            MAE_Cosh += diff;

            // get max difference between nion and Trigintaduonion hyperbolic cosine
            if (diff > MAX_Cosh) {
                MAX_Cosh = diff;
                maxCoshNion = nionResult;
                maxCoshTrigintaduonion = trigintaduonion_result;
                max1valueCosh = trigintaduonion1;
                max2valueCosh = trigintaduonion2;
            }
        }

        /// hyperbolic tangent
        {
            // evaluate nion hyperbolic tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tanh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion hyperbolic tangent, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = tanh(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            tanhTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ hyperbolictangent difference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Tanh +=  diff / norm(trigintaduonion_result);
            
            MAE_Tanh += diff;

            // get max difference between nion and Trigintaduonion hyperbolic tangent
            if (diff > MAX_Tanh) {
                MAX_Tanh = diff;
                maxTanhNion = nionResult;
                maxTanhTrigintaduonion = trigintaduonion_result;
                max1valueTanh = trigintaduonion1;
                max2valueTanh = trigintaduonion2;
            }
        }

        /// hyperbolic arc sine
        {
            // evaluate nion hyperbolic arc sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = asinh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            asinhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion hyperbolic arc sine, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = asinh(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            asinhTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ hyperbolicarc sine difference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Asinh +=  diff / norm(trigintaduonion_result);
            
            MAE_Asinh += diff;

            // get max difference between nion and Trigintaduonion hyperbolic arc sine
            if (diff > MAX_Asinh) {
                MAX_Asinh = diff;
                maxAsinhNion = nionResult;
                maxAsinhTrigintaduonion = trigintaduonion_result;
                max1valueAsinh = trigintaduonion1;
                max2valueAsinh = trigintaduonion2;
            }
        }

        /// hyperbolic arc cosine
        {
            // evaluate nion hyperbolic arc cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = acosh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            acoshNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion hyperbolic arc cosine, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = acosh(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            acoshTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ hyperbolicarc cosine difference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Acosh +=  diff / norm(trigintaduonion_result);
            
            MAE_Acosh += diff;

            // get max difference between nion and Trigintaduonion hyperbolic arc cosine
            if (diff > MAX_Acosh) {
                MAX_Acosh = diff;
                maxAcoshNion = nionResult;
                maxAcoshTrigintaduonion = trigintaduonion_result;
                max1valueAcosh = trigintaduonion1;
                max2valueAcosh = trigintaduonion2;
            }
        }

        /// hyperbolic arc tangent
        {
            // evaluate nion hyperbolic arc tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = atanh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            atanhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Trigintaduonion hyperbolic arc tangent, and time
            startTrigintaduonion = std::chrono::high_resolution_clock::now();
            trigintaduonion_result = atanh(trigintaduonion1);
            endTrigintaduonion = std::chrono::high_resolution_clock::now();
            atanhTrigintaduonionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endTrigintaduonion - startTrigintaduonion).count();

            // get difference between nion and Trigintaduonion norms. Add to MAE_ hyperbolicarc tangent difference
            diff = getMAEtrigintaduonion<T>(nionResult, trigintaduonion_result);
            MRE_Atanh +=  diff / norm(trigintaduonion_result);
            
            MAE_Atanh += diff;

            // get max difference between nion and Trigintaduonion hyperbolic arc tangent
            if (diff > MAX_Atanh) {
                MAX_Atanh = diff;
                maxAtanhNion = nionResult;
                maxAtanhTrigintaduonion = trigintaduonion_result;
                max1valueAtanh = trigintaduonion1;
                max2valueAtanh = trigintaduonion2;
            }
        }
    }
    delete[] val1;
    delete[] val2;
    
    T trialfp = static_cast<T>(trials);

    /*** norm ***/
    std::cout << "----> Norm <---- " << std::endl;
    std::cout << "Average norm time for nion: " << normNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average norm time for Trigintaduonion: " << normTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(normNionTimer, normTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Norm / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Norm / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Norm << std::endl;
    std::cout << "nion: " << maxNormNion << "\nTrigintaduonion: " << maxNormTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueNorm << std::endl;


    /*** addition ***/
    std::cout << "\n\n---> Addition <--- " << std::endl;
    std::cout << "Average addition time for nion: " << addNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average addition time for Trigintaduonion: " << addTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(addNionTimer, addTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Add / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Add / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Add << std::endl;
    std::cout << "nion: " << maxAddNion << "\nTrigintaduonion: " << maxAddTrigintaduonion << std::endl;
    std::cout << "input1: " << max1valueAdd << "\ninput2: " << max2valueAdd << std::endl;


    /*** conjugate ***/
    std::cout << "\n\n---> Conjugate <--- " << std::endl;
    std::cout << "Average conjugate time for nion: " << conjNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average conjugate time for Trigintaduonion: " << conjTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(conjNionTimer, conjTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Conj / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Conj / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Conj << std::endl;
    std::cout << "nion: " << maxConjNion << "\nTrigintaduonion: " << maxConjTrigintaduonion << std::endl;
    std::cout << "input: " << max2valueConj << std::endl;


    /*** multiplication ***/
    std::cout << "\n\n---> Multiplication <--- " << std::endl;
    std::cout << "Average multiplication time for nion: " << mulNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average multiplication time for Trigintaduonion: " << mulTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(mulNionTimer, mulTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Mul / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Mul / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Mul << std::endl;
    std::cout << "nion: " << maxMulNion << "\nTrigintaduonion: " << maxMulTrigintaduonion << std::endl;
    std::cout << "input1: " << max1valueMul << "\ninput2: " << max2valueMul << std::endl;


    /*** division ***/
    std::cout << "\n\n---> Division <--- " << std::endl;
    std::cout << "Average division time for nion: " << divNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average division time for Trigintaduonion: " << divTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(divNionTimer, divTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Div / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Div / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Div << std::endl;
    std::cout << "nion: " << maxDivNion << "\nTrigintaduonion: " << maxDivTrigintaduonion << std::endl;
    std::cout << "input1: " << max1valueDiv << "\ninput2: " << max2valueDiv << std::endl;


    /*** power ***/
    std::cout << "\n\n---> Power <--- " << std::endl;
    std::cout << "Average power time for nion: " << powNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average power time for Trigintaduonion: " << powTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(powNionTimer, powTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Pow / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Pow / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Pow << std::endl;
    std::cout << "nion: " << maxPowNion << "\nTrigintaduonion: " << maxPowTrigintaduonion << std::endl;
    std::cout << "input1: " << max1valuePow << "\ninput2: " << max2valuePow << std::endl;

    /*** square ***/
    std::cout << "\n\n---> Square <--- " << std::endl;
    std::cout << "Average square time for nion: " << SqNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average square time for Trigintaduonion: " << SqTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(SqNionTimer, SqTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Sq / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Sq / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Sq << std::endl;
    std::cout << "nion: " << maxSqNion << "\nTrigintaduonion: " << maxSqTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueSq << std::endl;

    /*** square Root ***/
    std::cout << "\n\n---> Square root <--- " << std::endl;
    std::cout << "Average square time for nion: " << SqrtNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average square time for Trigintaduonion: " << SqrtTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(SqrtNionTimer, SqrtTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Sqrt / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Sqrt / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Sqrt << std::endl;
    std::cout << "nion: " << maxSqrtNion << "\nTrigintaduonion: " << maxSqrtTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueSqrt << std::endl;


    /*** exponential ***/
    std::cout << "\n\n---> Exponential <--- " << std::endl;
    std::cout << "Average exponential time for nion: " << expNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average exponential time for Trigintaduonion: " << expTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(expNionTimer, expTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Exp / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Exp / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Exp << std::endl;
    std::cout << "nion: " << maxExpNion << "\nTrigintaduonion: " << maxExpTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueExp << std::endl;


    /*** logarithm ***/
    std::cout << "\n\n---> Logarithm <--- " << std::endl;
    std::cout << "Average logarithm time for nion: " << logNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average logarithm time for Trigintaduonion: " << logTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(logNionTimer, logTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Log / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Log / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Log << std::endl;
    std::cout << "nion: " << maxLogNion << "\nTrigintaduonion: " << maxLogTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueLog << std::endl;


    /*** sine ***/
    std::cout << "\n\n---> Sine <--- " << std::endl;
    std::cout << "Average sin time for nion: " << sinNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average sin time for Trigintaduonion: " << sinTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(sinNionTimer, sinTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Sin / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Sin / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Sin << std::endl;
    std::cout << "nion: " << maxSinNion << "\nTrigintaduonion: " << maxSinTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueSin << std::endl;


    /***  inverse sine ***/
    std::cout << "\n\n---> Inverse Sine <--- " << std::endl;
    std::cout << "Average asin time for nion: " << asinNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average asin time for Trigintaduonion: " << asinTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(asinNionTimer, asinTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Asin / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Asin / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Asin << std::endl;
    std::cout << "nion: " << maxAsinNion << "\nTrigintaduonion: " << maxAsinTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueAsin << std::endl;


    /*** cosine ***/
    std::cout << "\n\n---> Cosine <--- " << std::endl;
    std::cout << "Average cos time for nion: " << cosNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average cos time for Trigintaduonion: " << cosTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(cosNionTimer, cosTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Cos / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Cos / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Cos << std::endl;
    std::cout << "nion: " << maxCosNion << "\nTrigintaduonion: " << maxCosTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueCos << std::endl;


    /***  inverse cosine ***/
    std::cout << "\n\n---> Inverse Cosine <--- " << std::endl;
    std::cout << "Average acos time for nion: " << acosNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average acos time for Trigintaduonion: " << acosTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(acosNionTimer, acosTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Acos / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Acos / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Acos << std::endl;
    std::cout << "nion: " << maxAcosNion << "\nTrigintaduonion: " << maxAcosTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueAcos << std::endl;


    /*** tangent ***/
    std::cout << "\n\n---> Tangent <--- " << std::endl;
    std::cout << "Average tan time for nion: " << tanNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average tan time for Trigintaduonion: " << tanTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(tanNionTimer, tanTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Tan / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Tan / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Tan << std::endl;
    std::cout << "nion: " << maxTanNion << "\nTrigintaduonion: " << maxTanTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueTan << std::endl;


    /***  inverse tangent ***/
    std::cout << "\n\n---> Inverse Tangent <--- " << std::endl;
    std::cout << "Average atan time for nion: " << atanNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average atan time for Trigintaduonion: " << atanTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(atanNionTimer, atanTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Atan / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Atan / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Atan << std::endl;
    std::cout << "nion: " << maxAtanNion << "\nTrigintaduonion: " << maxAtanTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueAtan << std::endl;


    /*** hyperbolic sine ***/
    std::cout << "\n\n---> Hyperbolic Sine <--- " << std::endl;
    std::cout << "Average sinh time for nion: " << sinhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average sinh time for Trigintaduonion: " << sinhTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(sinhNionTimer, sinhTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Sinh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Sinh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Sinh << std::endl;
    std::cout << "nion: " << maxSinhNion << "\nTrigintaduonion: " << maxSinhTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueSinh << std::endl;


    /*** inverse hyperbolic sine ***/
    std::cout << "\n\n---> Inverse Hyperbolic Sine <--- " << std::endl;
    std::cout << "Average asinh time for nion: " << asinhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average asinh time for Trigintaduonion: " << asinhTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(asinhNionTimer, asinhTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Asinh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Asinh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Asinh << std::endl;
    std::cout << "nion: " << maxAsinhNion << "\nTrigintaduonion: " << maxAsinhTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueAsinh << std::endl;


    /*** hyperbolic cosine ***/
    std::cout << "\n\n---> Hyperbolic Cosine <--- " << std::endl;
    std::cout << "Average cosh time for nion: " << coshNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average cosh time for Trigintaduonion: " << coshTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(coshNionTimer, coshTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Cosh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Cosh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Cosh << std::endl;
    std::cout << "nion: " << maxCoshNion << "\nTrigintaduonion: " << maxCoshTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueCosh << std::endl;


    /*** inverse hyperbolic cosine ***/
    std::cout << "\n\n---> Inverse Hyperbolic Cosine <--- " << std::endl;
    std::cout << "Average acosh time for nion: " << acoshNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average acosh time for Trigintaduonion: " << acoshTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(acoshNionTimer, acoshTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Acosh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Acosh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Acosh << std::endl;
    std::cout << "nion: " << maxAcoshNion << "\nTrigintaduonion: " << maxAcoshTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueAcosh << std::endl;


    /*** hyperbolic tangent ***/
    std::cout << "\n\n---> Hyperbolic Tangent <--- " << std::endl;
    std::cout << "Average tanh time for nion: " << tanhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average tanh time for Trigintaduonion: " << tanhTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(tanhNionTimer, tanhTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Tanh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Tanh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Tanh << std::endl;
    std::cout << "nion: " << maxTanhNion << "\nTrigintaduonion: " << maxTanhTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueTanh << std::endl;


    /*** inverse hyperbolic tangent ***/
    std::cout << "\n\n---> Inverse Hyperbolic Tangent <--- " << std::endl;
    std::cout << "Average atanh time for nion: " << atanhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average atanh time for Trigintaduonion: " << atanhTrigintaduonionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(atanhNionTimer, atanhTrigintaduonionTimer);
    std::cout << "Average difference between nion and Trigintaduonion: " << MRE_Atanh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Trigintaduonion: " << MAE_Atanh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Trigintaduonion: " << MAX_Atanh << std::endl;
    std::cout << "nion: " << maxAtanhNion << "\nTrigintaduonion: " << maxAtanhTrigintaduonion << std::endl;
    std::cout << "input: " << max1valueAtanh << std::endl;
}

#endif //NION_TRIGINTADUONIONTEST_HPP
#endif
