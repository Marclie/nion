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
#include "../nion.hpp"
#include <iostream>


using Nion::nion;

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
template<typename T>
void SedenionComparison(int trials) {
    std::cout << "\n\n#### Comparing Sedenion with nion: " << trials << " trials ####\n" << std::endl;
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());

    // timers for norms
    T normNionTimer = 0;
    T normSedenionTimer = 0;
    T MAE_Norm = 0;
    T MRE_Norm = 0;
    T MAX_Norm = -1;
    T maxNormNion = 0;
    T maxNormSedenion = 0;
    Sedenion<T> max1valueNorm;


    // timers for addition
    T addNionTimer = 0;
    T addSedenionTimer = 0;
    T MAE_Add = 0;
    T MRE_Add = 0;
    T MAX_Add = -1;
    nion<T> maxAddNion;
    Sedenion<T> maxAddSedenion;
    Sedenion<T> max1valueAdd;
    Sedenion<T> max2valueAdd;

    // timers for conjugate
    T conjNionTimer = 0;
    T conjSedenionTimer = 0;
    T MAE_Conj = 0;
    T MRE_Conj = 0;
    T MAX_Conj = -1;
    nion<T> maxConjNion;
    Sedenion<T> maxConjSedenion;
    Sedenion<T> max1valueConj;
    Sedenion<T> max2valueConj;

    // timers for multiplication
    T mulNionTimer = 0;
    T mulSedenionTimer = 0;
    T MAE_Mul = 0;
    T MRE_Mul = 0;
    T MAX_Mul = -1;
    nion<T> maxMulNion;
    Sedenion<T> maxMulSedenion;
    Sedenion<T> max1valueMul;
    Sedenion<T> max2valueMul;

    // timers for division
    T divNionTimer = 0;
    T divSedenionTimer = 0;
    T MAE_Div = 0;
    T MRE_Div = 0;
    T MAX_Div = -1;
    nion<T> maxDivNion;
    Sedenion<T> maxDivSedenion;
    Sedenion<T> max1valueDiv;
    Sedenion<T> max2valueDiv;

    // timers for power
    T powNionTimer = 0;
    T powSedenionTimer = 0;
    T MAE_Pow = 0;
    T MRE_Pow = 0;
    T MAX_Pow = -1;
    nion<T> maxPowNion;
    Sedenion<T> maxPowSedenion;
    Sedenion<T> max1valuePow;
    Sedenion<T> max2valuePow;

    // timers for exponential
    T expNionTimer = 0;
    T expSedenionTimer = 0;
    T MAE_Exp = 0;
    T MRE_Exp = 0;
    T MAX_Exp = -1;
    nion<T> maxExpNion;
    Sedenion<T> maxExpSedenion;
    Sedenion<T> max1valueExp;
    Sedenion<T> max2valueExp;

    // timers for logarithm
    T logNionTimer = 0;
    T logSedenionTimer = 0;
    T MAE_Log = 0;
    T MRE_Log = 0;
    T MAX_Log = -1;
    nion<T> maxLogNion;
    Sedenion<T> maxLogSedenion;
    Sedenion<T> max1valueLog;
    Sedenion<T> max2valueLog;

    // timers for sin
    T sinNionTimer = 0;
    T sinSedenionTimer = 0;
    T MAE_Sin = 0;
    T MRE_Sin = 0;
    T MAX_Sin = -1;
    nion<T> maxSinNion;
    Sedenion<T> maxSinSedenion;
    Sedenion<T> max1valueSin;
    Sedenion<T> max2valueSin;

    // timers for asin
    T asinNionTimer = 0;
    T asinSedenionTimer = 0;
    T MAE_Asin = 0;
    T MRE_Asin = 0;
    T MAX_Asin = -1;
    nion<T> maxAsinNion;
    Sedenion<T> maxAsinSedenion;
    Sedenion<T> max1valueAsin;
    Sedenion<T> max2valueAsin;

    // timers for cos
    T cosNionTimer = 0;
    T cosSedenionTimer = 0;
    T MAE_Cos = 0;
    T MRE_Cos = 0;
    T MAX_Cos = -1;
    nion<T> maxCosNion;
    Sedenion<T> maxCosSedenion;
    Sedenion<T> max1valueCos;
    Sedenion<T> max2valueCos;

    // timers for acos
    T acosNionTimer = 0;
    T acosSedenionTimer = 0;
    T MAE_Acos = 0;
    T MRE_Acos = 0;
    T MAX_Acos = -1;
    nion<T> maxAcosNion;
    Sedenion<T> maxAcosSedenion;
    Sedenion<T> max1valueAcos;
    Sedenion<T> max2valueAcos;

    // timers for tan
    T tanNionTimer = 0;
    T tanSedenionTimer = 0;
    T MAE_Tan = 0;
    T MRE_Tan = 0;
    T MAX_Tan = -1;
    nion<T> maxTanNion;
    Sedenion<T> maxTanSedenion;
    Sedenion<T> max1valueTan;
    Sedenion<T> max2valueTan;

    // timers for atan
    T atanNionTimer = 0;
    T atanSedenionTimer = 0;
    T MAE_Atan = 0;
    T MRE_Atan = 0;
    T MAX_Atan = -1;
    nion<T> maxAtanNion;
    Sedenion<T> maxAtanSedenion;
    Sedenion<T> max1valueAtan;
    Sedenion<T> max2valueAtan;

    // timers for sinh
    T sinhNionTimer = 0;
    T sinhSedenionTimer = 0;
    T MAE_Sinh = 0;
    T MRE_Sinh = 0;
    T MAX_Sinh = -1;
    nion<T> maxSinhNion;
    Sedenion<T> maxSinhSedenion;
    Sedenion<T> max1valueSinh;
    Sedenion<T> max2valueSinh;

    // timers for asinh
    T asinhNionTimer = 0;
    T asinhSedenionTimer = 0;
    T MAE_Asinh = 0;
    T MRE_Asinh = 0;
    T MAX_Asinh = -1;
    nion<T> maxAsinhNion;
    Sedenion<T> maxAsinhSedenion;
    Sedenion<T> max1valueAsinh;
    Sedenion<T> max2valueAsinh;

    // timers for cosh
    T coshNionTimer = 0;
    T coshSedenionTimer = 0;
    T MAE_Cosh = 0;
    T MRE_Cosh = 0;
    T MAX_Cosh = -1;
    nion<T> maxCoshNion;
    Sedenion<T> maxCoshSedenion;
    Sedenion<T> max1valueCosh;
    Sedenion<T> max2valueCosh;

    // timers for acosh
    T acoshNionTimer = 0;
    T acoshSedenionTimer = 0;
    T MAE_Acosh = 0;
    T MRE_Acosh = 0;
    T MAX_Acosh = -1;
    nion<T> maxAcoshNion;
    Sedenion<T> maxAcoshSedenion;
    Sedenion<T> max1valueAcosh;
    Sedenion<T> max2valueAcosh;

    // timers for tanh
    T tanhNionTimer = 0;
    T tanhSedenionTimer = 0;
    T MAE_Tanh = 0;
    T MRE_Tanh = 0;
    T MAX_Tanh = -1;
    nion<T> maxTanhNion;
    Sedenion<T> maxTanhSedenion;
    Sedenion<T> max1valueTanh;
    Sedenion<T> max2valueTanh;

    // timers for atanh
    T atanhNionTimer = 0;
    T atanhSedenionTimer = 0;
    T MAE_Atanh = 0;
    T MRE_Atanh = 0;
    T MAX_Atanh = -1;
    nion<T> maxAtanhNion;
    Sedenion<T> maxAtanhSedenion;
    Sedenion<T> max1valueAtanh;
    Sedenion<T> max2valueAtanh;


    auto startNion = std::chrono::high_resolution_clock::now();
    auto endNion = std::chrono::high_resolution_clock::now();
    auto startSedenion = std::chrono::high_resolution_clock::now();
    auto endSedenion = std::chrono::high_resolution_clock::now();

    nion<T> nionResult;
    Sedenion<T> sedenion_result;
    T diff;

    std::uniform_real_distribution<T> distribution(-10.0, 10.0);
    T* val1 = new T[16];
    T* val2 = new T[16];
    for (int i = 0; i < trials; ++i) {
        // get random Sedenion number
        for (int j = 0; j < 16; ++j) {
            val1[j] = distribution(generator);
            val2[j] = distribution(generator);
        }
        nion<T> nion1(val1, 16);
        nion<T> nion2(val2, 16);

        Sedenion<T> sedenion1(val1[0], val1[1], val1[2], val1[3], val1[4], val1[5], val1[6], val1[7], val1[8], val1[9], val1[10], val1[11], val1[12], val1[13], val1[14], val1[15]);
        Sedenion<T> sedenion2(val2[0], val2[1], val2[2], val2[3], val2[4], val2[5], val2[6], val2[7], val2[8], val2[9], val2[10], val2[11], val2[12], val2[13], val2[14], val2[15]);


        // norm
        {
            // evaluate nion norm, and time
            startNion = std::chrono::high_resolution_clock::now();
            T nionNorm = nion1.abs();
            endNion = std::chrono::high_resolution_clock::now();
            normNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion norm, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            T stdNorm = norm(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            normSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ normdifference
            diff = std::abs(nionNorm - stdNorm);
            MRE_Norm += diff;
            diff /= stdNorm;
            MAE_Norm += diff;

            // get max difference between nion and Sedenion norms
            if (diff > MAX_Norm) {
                MAX_Norm = diff;
                maxNormNion = nionNorm;
                maxNormSedenion = stdNorm;
                max1valueNorm = sedenion1;
            }
        }

        /// addition
        {
            // evaluate nion addition, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 + nion2;
            endNion = std::chrono::high_resolution_clock::now();
            addNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion addition, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = sedenion1 + sedenion2;
            endSedenion = std::chrono::high_resolution_clock::now();
            addSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ additiondifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Add += diff;
            diff /= norm(sedenion_result);
            MAE_Add += diff;

            // get max difference between nion and Sedenion addition
            if (diff > MAX_Add) {
                MAX_Add = diff;
                maxAddNion = nionResult;
                maxAddSedenion = sedenion_result;
                max1valueAdd = sedenion1;
                max2valueAdd = sedenion2;
            }
        }

        /// conjugate
        {
            // evaluate nion conjugate, and time
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();
            nionResult = nion2.conj();
            endNion = std::chrono::high_resolution_clock::now();
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = conj(sedenion2);
            endSedenion = std::chrono::high_resolution_clock::now();
            conjSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ conjugatedifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Conj += diff;
            diff /= norm(sedenion_result);
            MAE_Conj += diff;

            // get max difference between nion and Sedenion conjugate
            if (diff > MAX_Conj) {
                MAX_Conj = diff;
                maxConjNion = nionResult;
                maxConjSedenion = sedenion_result;
                max1valueConj = sedenion1;
                max2valueConj = sedenion2;
            }

        }

        /// multiplication
        {
            // evaluate nion multiplication, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 * nion2;
            endNion = std::chrono::high_resolution_clock::now();
            mulNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion multiplication, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = sedenion1 * sedenion2;
            endSedenion = std::chrono::high_resolution_clock::now();
            mulSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ multiplicationdifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Mul += diff;
            diff /= norm(sedenion_result);
            MAE_Mul += diff;

            // get max difference between nion and Sedenion multiplication
            if (diff > MAX_Mul) {
                MAX_Mul = diff;
                maxMulNion = nionResult;
                maxMulSedenion = sedenion_result;
                max1valueMul = sedenion1;
                max2valueMul = sedenion2;
            }

        }

        /// division
        {
            // evaluate nion division, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 / nion2;
            endNion = std::chrono::high_resolution_clock::now();
            divNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion division, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = sedenion1 / sedenion2;
            endSedenion = std::chrono::high_resolution_clock::now();
            divSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ divisiondifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Div += diff;
            diff /= norm(sedenion_result);
            MAE_Div += diff;

            // get max difference between nion and Sedenion division
            if (diff > MAX_Div) {
                MAX_Div = diff;
                maxDivNion = nionResult;
                maxDivSedenion = sedenion_result;
                max1valueDiv = sedenion1;
                max2valueDiv = sedenion2;
            }
        }

        /// power
        {
            // evaluate nion power, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = pow(nion1, static_cast<T>(5));
            endNion = std::chrono::high_resolution_clock::now();
            powNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion power, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = pow(sedenion1, static_cast<T>(5));
            endSedenion = std::chrono::high_resolution_clock::now();
            powSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ powerdifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Pow += diff;
            diff /= norm(sedenion_result);
            MAE_Pow += diff;

            // get max difference between nion and Sedenion power
            if (diff > MAX_Pow) {
                MAX_Pow = diff;
                maxPowNion = nionResult;
                maxPowSedenion = sedenion_result;
                max1valuePow = sedenion1;
                max2valuePow = sedenion2;
            }
        }

        /// exponential
        {
            // evaluate nion exponential, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = exp(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            expNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion exponential, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = exp(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            expSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ exponentialdifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Exp += diff;
            diff /= norm(sedenion_result);
            MAE_Exp += diff;

            // get max difference between nion and Sedenion exponential
            if (diff > MAX_Exp) {
                MAX_Exp = diff;
                maxExpNion = nionResult;
                maxExpSedenion = sedenion_result;
                max1valueExp = sedenion1;
                max2valueExp = sedenion2;
            }
        }

        /// logarithm
        {
            // evaluate nion logarithm, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = log(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            logNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion logarithm, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = log(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            logSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ logarithmdifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Log += diff;
            diff /= norm(sedenion_result);
            MAE_Log += diff;

            // get max difference between nion and Sedenion logarithm
            if (diff > MAX_Log) {
                MAX_Log = diff;
                maxLogNion = nionResult;
                maxLogSedenion = sedenion_result;
                max1valueLog = sedenion1;
                max2valueLog = sedenion2;
            }
        }

        /// sine
        {
            // evaluate nion sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sin(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion sine, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = sin(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            sinSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ sinedifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Sin += diff;
            diff /= norm(sedenion_result);
            MAE_Sin += diff;

            // get max difference between nion and Sedenion sine
            if (diff > MAX_Sin) {
                MAX_Sin = diff;
                maxSinNion = nionResult;
                maxSinSedenion = sedenion_result;
                max1valueSin = sedenion1;
                max2valueSin = sedenion2;
            }
        }

        /// asine
        {
            // evaluate nion asine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = asin(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            asinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion asine, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = asin(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            asinSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ asinedifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Asin += diff;
            diff /= norm(sedenion_result);
            MAE_Asin += diff;

            // get max difference between nion and Sedenion asine
            if (diff > MAX_Asin) {
                MAX_Asin = diff;
                maxAsinNion = nionResult;
                maxAsinSedenion = sedenion_result;
                max1valueAsin = sedenion1;
                max2valueAsin = sedenion2;
            }
        }

        /// cosine
        {
            // evaluate nion cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cos(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            cosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion cosine, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = cos(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            cosSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ cosinedifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Cos += diff;
            diff /= norm(sedenion_result);
            MAE_Cos += diff;

            // get max difference between nion and Sedenion cosine
            if (diff > MAX_Cos) {
                MAX_Cos = diff;
                maxCosNion = nionResult;
                maxCosSedenion = sedenion_result;
                max1valueCos = sedenion1;
                max2valueCos = sedenion2;
            }
        }

        /// acosine
        {
            // evaluate nion acosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = acos(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            acosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion acosine, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = acos(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            acosSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ acosinedifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Acos += diff;
            diff /= norm(sedenion_result);
            MAE_Acos += diff;

            // get max difference between nion and Sedenion acosine
            if (diff > MAX_Acos) {
                MAX_Acos = diff;
                maxAcosNion = nionResult;
                maxAcosSedenion = sedenion_result;
                max1valueAcos = sedenion1;
                max2valueAcos = sedenion2;
            }
        }

        /// tangent
        {
            // evaluate nion tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tan(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion tangent, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = tan(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            tanSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ tangentdifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Tan += diff;
            diff /= norm(sedenion_result);
            MAE_Tan += diff;

            // get max difference between nion and Sedenion tangent
            if (diff > MAX_Tan) {
                MAX_Tan = diff;
                maxTanNion = nionResult;
                maxTanSedenion = sedenion_result;
                max1valueTan = sedenion1;
                max2valueTan = sedenion2;
            }
        }

        /// atan
        {
            // evaluate nion atan, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = atan(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            atanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion atan, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = atan(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            atanSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ atandifference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Atan += diff;
            diff /= norm(sedenion_result);
            MAE_Atan += diff;

            // get max difference between nion and Sedenion atan
            if (diff > MAX_Atan) {
                MAX_Atan = diff;
                maxAtanNion = nionResult;
                maxAtanSedenion = sedenion_result;
                max1valueAtan = sedenion1;
                max2valueAtan = sedenion2;
            }
        }

        /// hyperbolic sine
        {
            // evaluate nion hyperbolic sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sinh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion hyperbolic sine, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = sinh(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            sinhSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ hyperbolicsine difference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Sinh += diff;
            diff /= norm(sedenion_result);
            MAE_Sinh += diff;

            // get max difference between nion and Sedenion hyperbolic sine
            if (diff > MAX_Sinh) {
                MAX_Sinh = diff;
                maxSinhNion = nionResult;
                maxSinhSedenion = sedenion_result;
                max1valueSinh = sedenion1;
                max2valueSinh = sedenion2;
            }
        }

        /// hyperbolic cosine
        {
            // evaluate nion hyperbolic cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cosh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            coshNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion hyperbolic cosine, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = cosh(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            coshSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ hyperboliccosine difference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Cosh += diff;
            diff /= norm(sedenion_result);
            MAE_Cosh += diff;

            // get max difference between nion and Sedenion hyperbolic cosine
            if (diff > MAX_Cosh) {
                MAX_Cosh = diff;
                maxCoshNion = nionResult;
                maxCoshSedenion = sedenion_result;
                max1valueCosh = sedenion1;
                max2valueCosh = sedenion2;
            }
        }

        /// hyperbolic tangent
        {
            // evaluate nion hyperbolic tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tanh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion hyperbolic tangent, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = tanh(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            tanhSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ hyperbolictangent difference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Tanh += diff;
            diff /= norm(sedenion_result);
            MAE_Tanh += diff;

            // get max difference between nion and Sedenion hyperbolic tangent
            if (diff > MAX_Tanh) {
                MAX_Tanh = diff;
                maxTanhNion = nionResult;
                maxTanhSedenion = sedenion_result;
                max1valueTanh = sedenion1;
                max2valueTanh = sedenion2;
            }
        }

        /// hyperbolic arc sine
        {
            // evaluate nion hyperbolic arc sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = asinh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            asinhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion hyperbolic arc sine, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = asinh(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            asinhSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ hyperbolicarc sine difference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Asinh += diff;
            diff /= norm(sedenion_result);
            MAE_Asinh += diff;

            // get max difference between nion and Sedenion hyperbolic arc sine
            if (diff > MAX_Asinh) {
                MAX_Asinh = diff;
                maxAsinhNion = nionResult;
                maxAsinhSedenion = sedenion_result;
                max1valueAsinh = sedenion1;
                max2valueAsinh = sedenion2;
            }
        }

        /// hyperbolic arc cosine
        {
            // evaluate nion hyperbolic arc cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = acosh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            acoshNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion hyperbolic arc cosine, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = acosh(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            acoshSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ hyperbolicarc cosine difference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Acosh += diff;
            diff /= norm(sedenion_result);
            MAE_Acosh += diff;

            // get max difference between nion and Sedenion hyperbolic arc cosine
            if (diff > MAX_Acosh) {
                MAX_Acosh = diff;
                maxAcoshNion = nionResult;
                maxAcoshSedenion = sedenion_result;
                max1valueAcosh = sedenion1;
                max2valueAcosh = sedenion2;
            }
        }

        /// hyperbolic arc tangent
        {
            // evaluate nion hyperbolic arc tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = atanh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            atanhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate Sedenion hyperbolic arc tangent, and time
            startSedenion = std::chrono::high_resolution_clock::now();
            sedenion_result = atanh(sedenion1);
            endSedenion = std::chrono::high_resolution_clock::now();
            atanhSedenionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endSedenion - startSedenion).count();

            // get difference between nion and Sedenion norms. Add to MAE_ hyperbolicarc tangent difference
            diff = getMAEsedenion<T>(nionResult, sedenion_result);
            MRE_Atanh += diff;
            diff /= norm(sedenion_result);
            MAE_Atanh += diff;

            // get max difference between nion and Sedenion hyperbolic arc tangent
            if (diff > MAX_Atanh) {
                MAX_Atanh = diff;
                maxAtanhNion = nionResult;
                maxAtanhSedenion = sedenion_result;
                max1valueAtanh = sedenion1;
                max2valueAtanh = sedenion2;
            }
        }
    }
    delete[] val1;
    delete[] val2;

    T trialfp = static_cast<T>(trials);

    /*** norm ***/
    std::cout << "----> Norm <---- " << std::endl;
    std::cout << "Average norm time for nion: " << normNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average norm time for Sedenion: " << normSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(normNionTimer, normSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Norm / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Norm / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Norm << std::endl;
    std::cout << "nion: " << maxNormNion << "\nSedenion: " << maxNormSedenion << std::endl;
    std::cout << "input: " << max1valueNorm << std::endl;


    /*** addition ***/
    std::cout << "\n\n---> Addition <--- " << std::endl;
    std::cout << "Average addition time for nion: " << addNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average addition time for Sedenion: " << addSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(addNionTimer, addSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Add / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Add / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Add << std::endl;
    std::cout << "nion: " << maxAddNion << "\nSedenion: " << maxAddSedenion << std::endl;
    std::cout << "input1: " << max1valueAdd << "\ninput2: " << max2valueAdd << std::endl;


    /*** conjugate ***/
    std::cout << "\n\n---> Conjugate <--- " << std::endl;
    std::cout << "Average conjugate time for nion: " << conjNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average conjugate time for Sedenion: " << conjSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(conjNionTimer, conjSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Conj / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Conj / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Conj << std::endl;
    std::cout << "nion: " << maxConjNion << "\nSedenion: " << maxConjSedenion << std::endl;
    std::cout << "input: " << max2valueConj << std::endl;


    /*** multiplication ***/
    std::cout << "\n\n---> Multiplication <--- " << std::endl;
    std::cout << "Average multiplication time for nion: " << mulNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average multiplication time for Sedenion: " << mulSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(mulNionTimer, mulSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Mul / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Mul / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Mul << std::endl;
    std::cout << "nion: " << maxMulNion << "\nSedenion: " << maxMulSedenion << std::endl;
    std::cout << "input1: " << max1valueMul << "\ninput2: " << max2valueMul << std::endl;


    /*** division ***/
    std::cout << "\n\n---> Division <--- " << std::endl;
    std::cout << "Average division time for nion: " << divNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average division time for Sedenion: " << divSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(divNionTimer, divSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Div / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Div / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Div << std::endl;
    std::cout << "nion: " << maxDivNion << "\nSedenion: " << maxDivSedenion << std::endl;
    std::cout << "input1: " << max1valueDiv << "\ninput2: " << max2valueDiv << std::endl;


    /*** power ***/
    std::cout << "\n\n---> Power <--- " << std::endl;
    std::cout << "Average power time for nion: " << powNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average power time for Sedenion: " << powSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(powNionTimer, powSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Pow / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Pow / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Pow << std::endl;
    std::cout << "nion: " << maxPowNion << "\nSedenion: " << maxPowSedenion << std::endl;
    std::cout << "input1: " << max1valuePow << "\ninput2: " << max2valuePow << std::endl;


    /*** exponential ***/
    std::cout << "\n\n---> Exponential <--- " << std::endl;
    std::cout << "Average exponential time for nion: " << expNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average exponential time for Sedenion: " << expSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(expNionTimer, expSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Exp / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Exp / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Exp << std::endl;
    std::cout << "nion: " << maxExpNion << "\nSedenion: " << maxExpSedenion << std::endl;
    std::cout << "input: " << max1valueExp << std::endl;


    /*** logarithm ***/
    std::cout << "\n\n---> Logarithm <--- " << std::endl;
    std::cout << "Average logarithm time for nion: " << logNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average logarithm time for Sedenion: " << logSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(logNionTimer, logSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Log / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Log / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Log << std::endl;
    std::cout << "nion: " << maxLogNion << "\nSedenion: " << maxLogSedenion << std::endl;
    std::cout << "input: " << max1valueLog << std::endl;


    /*** sine ***/
    std::cout << "\n\n---> Sine <--- " << std::endl;
    std::cout << "Average sin time for nion: " << sinNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average sin time for Sedenion: " << sinSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(sinNionTimer, sinSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Sin / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Sin / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Sin << std::endl;
    std::cout << "nion: " << maxSinNion << "\nSedenion: " << maxSinSedenion << std::endl;
    std::cout << "input: " << max1valueSin << std::endl;


    /***  inverse sine ***/
    std::cout << "\n\n---> Inverse Sine <--- " << std::endl;
    std::cout << "Average asin time for nion: " << asinNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average asin time for Sedenion: " << asinSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(asinNionTimer, asinSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Asin / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Asin / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Asin << std::endl;
    std::cout << "nion: " << maxAsinNion << "\nSedenion: " << maxAsinSedenion << std::endl;
    std::cout << "input: " << max1valueAsin << std::endl;


    /*** cosine ***/
    std::cout << "\n\n---> Cosine <--- " << std::endl;
    std::cout << "Average cos time for nion: " << cosNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average cos time for Sedenion: " << cosSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(cosNionTimer, cosSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Cos / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Cos / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Cos << std::endl;
    std::cout << "nion: " << maxCosNion << "\nSedenion: " << maxCosSedenion << std::endl;
    std::cout << "input: " << max1valueCos << std::endl;


    /***  inverse cosine ***/
    std::cout << "\n\n---> Inverse Cosine <--- " << std::endl;
    std::cout << "Average acos time for nion: " << acosNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average acos time for Sedenion: " << acosSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(acosNionTimer, acosSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Acos / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Acos / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Acos << std::endl;
    std::cout << "nion: " << maxAcosNion << "\nSedenion: " << maxAcosSedenion << std::endl;
    std::cout << "input: " << max1valueAcos << std::endl;


    /*** tangent ***/
    std::cout << "\n\n---> Tangent <--- " << std::endl;
    std::cout << "Average tan time for nion: " << tanNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average tan time for Sedenion: " << tanSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(tanNionTimer, tanSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Tan / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Tan / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Tan << std::endl;
    std::cout << "nion: " << maxTanNion << "\nSedenion: " << maxTanSedenion << std::endl;
    std::cout << "input: " << max1valueTan << std::endl;


    /***  inverse tangent ***/
    std::cout << "\n\n---> Inverse Tangent <--- " << std::endl;
    std::cout << "Average atan time for nion: " << atanNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average atan time for Sedenion: " << atanSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(atanNionTimer, atanSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Atan / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Atan / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Atan << std::endl;
    std::cout << "nion: " << maxAtanNion << "\nSedenion: " << maxAtanSedenion << std::endl;
    std::cout << "input: " << max1valueAtan << std::endl;


    /*** hyperbolic sine ***/
    std::cout << "\n\n---> Hyperbolic Sine <--- " << std::endl;
    std::cout << "Average sinh time for nion: " << sinhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average sinh time for Sedenion: " << sinhSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(sinhNionTimer, sinhSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Sinh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Sinh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Sinh << std::endl;
    std::cout << "nion: " << maxSinhNion << "\nSedenion: " << maxSinhSedenion << std::endl;
    std::cout << "input: " << max1valueSinh << std::endl;


    /*** inverse hyperbolic sine ***/
    std::cout << "\n\n---> Inverse Hyperbolic Sine <--- " << std::endl;
    std::cout << "Average asinh time for nion: " << asinhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average asinh time for Sedenion: " << asinhSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(asinhNionTimer, asinhSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Asinh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Asinh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Asinh << std::endl;
    std::cout << "nion: " << maxAsinhNion << "\nSedenion: " << maxAsinhSedenion << std::endl;
    std::cout << "input: " << max1valueAsinh << std::endl;


    /*** hyperbolic cosine ***/
    std::cout << "\n\n---> Hyperbolic Cosine <--- " << std::endl;
    std::cout << "Average cosh time for nion: " << coshNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average cosh time for Sedenion: " << coshSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(coshNionTimer, coshSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Cosh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Cosh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Cosh << std::endl;
    std::cout << "nion: " << maxCoshNion << "\nSedenion: " << maxCoshSedenion << std::endl;
    std::cout << "input: " << max1valueCosh << std::endl;


    /*** inverse hyperbolic cosine ***/
    std::cout << "\n\n---> Inverse Hyperbolic Cosine <--- " << std::endl;
    std::cout << "Average acosh time for nion: " << acoshNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average acosh time for Sedenion: " << acoshSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(acoshNionTimer, acoshSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Acosh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Acosh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Acosh << std::endl;
    std::cout << "nion: " << maxAcoshNion << "\nSedenion: " << maxAcoshSedenion << std::endl;
    std::cout << "input: " << max1valueAcosh << std::endl;


    /*** hyperbolic tangent ***/
    std::cout << "\n\n---> Hyperbolic Tangent <--- " << std::endl;
    std::cout << "Average tanh time for nion: " << tanhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average tanh time for Sedenion: " << tanhSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(tanhNionTimer, tanhSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Tanh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Tanh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Tanh << std::endl;
    std::cout << "nion: " << maxTanhNion << "\nSedenion: " << maxTanhSedenion << std::endl;
    std::cout << "input: " << max1valueTanh << std::endl;


    /*** inverse hyperbolic tangent ***/
    std::cout << "\n\n---> Inverse Hyperbolic Tangent <--- " << std::endl;
    std::cout << "Average atanh time for nion: " << atanhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average atanh time for Sedenion: " << atanhSedenionTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(atanhNionTimer, atanhSedenionTimer);
    std::cout << "Average difference between nion and Sedenion: " << MRE_Atanh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and Sedenion: " << MAE_Atanh / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and Sedenion: " << MAX_Atanh << std::endl;
    std::cout << "nion: " << maxAtanhNion << "\nSedenion: " << maxAtanhSedenion << std::endl;
    std::cout << "input: " << max1valueAtanh << std::endl;
}

#endif //NION_HARDERTEST_HPP
#endif