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

#ifndef NION_BOOSTQUATERNIONTEST_HPP
#define NION_BOOSTQUATERNIONTEST_HPP

#include <fstream>
#include <boost/math/quaternion.hpp>
#include <chrono>
#include <random>
#include "nion.hpp"
#include <iostream>

using Nion::nion;

template<typename T>
void printSpeedupQuaternion(const T& niontime, const T& othertime) {
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
T getMAEquaternion(nion<T> nion, boost::math::quaternion<T> compare){
    int degree = nion.degree;
    T mae = 0;
    mae += std::pow(nion[0] - compare.R_component_1(), 2);
    mae += std::pow(nion[1] - compare.R_component_2(), 2);
    mae += std::pow(nion[2] - compare.R_component_3(), 2);
    mae += std::pow(nion[3] - compare.R_component_4(), 2);
    return std::sqrt(mae);
}

template<typename T>
void boostQuaternionComparison(int trials) {
    std::cout << "\n\n#### Comparing boost::math::quaternion with nion: " << trials << " trials ####\n" << std::endl;
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());

    // timers for norms
    T normNionTimer = 0;
    T normBoostTimer = 0;
    T MRE_Norm = 0;
    T MAE_Norm = 0;
    T MAX_Norm = -1;
    T maxNormNion;
    T maxNormBoost;
    boost::math::quaternion<T> max1value;
    boost::math::quaternion<T> max2value;

    // timers for addition
    T addNionTimer = 0;
    T addBoostTimer = 0;
    T MRE_Add = 0;
    T MAE_Add= 0;
    T MAX_Add = -1;
    nion<T> maxAddNion;
    boost::math::quaternion<T> maxAddBoost;
    boost::math::quaternion<T> max1valueAdd;
    boost::math::quaternion<T> max2valueAdd;

    // timers for conjugate
    T conjNionTimer = 0;
    T conjBoostTimer = 0;
    T MRE_Conj = 0;
    T MAE_Conj= 0;
    T MAX_Conj = -1;
    nion<T> maxConjNion;
    boost::math::quaternion<T> maxConjBoost;
    boost::math::quaternion<T> max1valueConj;
    boost::math::quaternion<T> max2valueConj;

    // timers for multiplication
    T mulNionTimer = 0;
    T mulBoostTimer = 0;
    T MRE_Mul = 0;
    T MAE_Mul= 0;
    T MAX_Mul = -1;
    nion<T> maxMulNion;
    boost::math::quaternion<T> maxMulBoost;
    boost::math::quaternion<T> max1valueMul;
    boost::math::quaternion<T> max2valueMul;

    // timers for division
    T divNionTimer = 0;
    T divBoostTimer = 0;
    T MRE_Div = 0;
    T MAE_Div= 0;
    T MAX_Div = -1;
    nion<T> maxDivNion;
    boost::math::quaternion<T> maxDivBoost;
    boost::math::quaternion<T> max1valueDiv;
    boost::math::quaternion<T> max2valueDiv;

    // timers for power
    T powNionTimer = 0;
    T powBoostTimer = 0;
    T MRE_Pow = 0;
    T MAE_Pow= 0;
    T MAX_Pow = -1;
    nion<T> maxPowNion;
    boost::math::quaternion<T> maxPowBoost;
    boost::math::quaternion<T> max1valuePow;
    boost::math::quaternion<T> max2valuePow;

    // timers for exponential
    T expNionTimer = 0;
    T expBoostTimer = 0;
    T MRE_Exp = 0;
    T MAE_Exp= 0;
    T MAX_Exp = -1;
    nion<T> maxExpNion;
    boost::math::quaternion<T> maxExpBoost;
    boost::math::quaternion<T> max1valueExp;
    boost::math::quaternion<T> max2valueExp;

    // timers for sine
    T sinNionTimer = 0;
    T sinBoostTimer = 0;
    T MRE_Sin = 0;
    T MAE_Sin= 0;
    T MAX_Sin = -1;
    nion<T> maxSinNion;
    boost::math::quaternion<T> maxSinBoost;
    boost::math::quaternion<T> max1valueSin;
    boost::math::quaternion<T> max2valueSin;

    // timers for cosine
    T cosNionTimer = 0;
    T cosBoostTimer = 0;
    T MRE_Cos = 0;
    T MAE_Cos= 0;
    T MAX_Cos = -1;
    nion<T> maxCosNion;
    boost::math::quaternion<T> maxCosBoost;
    boost::math::quaternion<T> max1valueCos;
    boost::math::quaternion<T> max2valueCos;

    // timers for tangent
    T tanNionTimer = 0;
    T tanBoostTimer = 0;
    T MRE_Tan = 0;
    T MAE_Tan= 0;
    T MAX_Tan = -1;
    nion<T> maxTanNion;
    boost::math::quaternion<T> maxTanBoost;
    boost::math::quaternion<T> max1valueTan;
    boost::math::quaternion<T> max2valueTan;

    // timers for sinh
    T sinhNionTimer = 0;
    T sinhBoostTimer = 0;
    T MRE_Sinh = 0;
    T MAE_Sinh= 0;
    T MAX_Sinh = -1;
    nion<T> maxSinhNion;
    boost::math::quaternion<T> maxSinhBoost;
    boost::math::quaternion<T> max1valueSinh;
    boost::math::quaternion<T> max2valueSinh;

    // timers for cosh
    T coshNionTimer = 0;
    T coshBoostTimer = 0;
    T MRE_Cosh = 0;
    T MAE_Cosh= 0;
    T MAX_Cosh = -1;
    nion<T> maxCoshNion;
    boost::math::quaternion<T> maxCoshBoost;
    boost::math::quaternion<T> max1valueCosh;
    boost::math::quaternion<T> max2valueCosh;

    // timers for tanh
    T tanhNionTimer = 0;
    T tanhBoostTimer = 0;
    T MRE_Tanh = 0;
    T MAE_Tanh= 0;
    T MAX_Tanh = -1;
    nion<T> maxTanhNion;
    boost::math::quaternion<T> maxTanhBoost;
    boost::math::quaternion<T> max1valueTanh;
    boost::math::quaternion<T> max2valueTanh;

    auto startNion = std::chrono::high_resolution_clock::now();
    auto endNion = std::chrono::high_resolution_clock::now();
    auto startBoost = std::chrono::high_resolution_clock::now();
    auto endBoost = std::chrono::high_resolution_clock::now();

    nion < T > nionResult;
    boost::math::quaternion<T> boostResult;
    T diff;

    T *vals1 = new T[4];
    T *vals2 = new T[4];
    T nionNorm;
    T boostNorm;
    std::uniform_real_distribution<T> distribution(-10, 10);
    for (int i = 0; i < trials; ++i) {

        //generate random quaternion numbers
        for (int j = 0; j < 4; ++j) {

            vals1[j] = distribution(generator);
            vals2[j] = distribution(generator);
        }

        nion < T > nion1(vals1, 4);
        nion < T > nion2(vals2, 4);
        boost::math::quaternion<T> boost1(vals1[0], vals1[1], vals1[2], vals1[3]);
        boost::math::quaternion<T> boost2(vals2[0], vals2[1], vals2[2], vals2[3]);

        // norm
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = abs(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            normNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::norm(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            normBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            nionNorm = abs(nionResult);
            boostNorm = boost::math::norm(boostResult);

            diff = fabs(nionNorm - boostNorm);
            MRE_Norm += diff;
            MAE_Norm += diff / boostNorm;

            if (diff > MAX_Norm) {
                MAX_Norm = diff;
                maxNormNion = nionNorm;
                maxNormBoost = boostNorm;
                max1value = boost1;
            }

        }

        // addition
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 + nion2;
            endNion = std::chrono::high_resolution_clock::now();
            addNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost1 + boost2;
            endBoost = std::chrono::high_resolution_clock::now();
            addBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Add += diff;
            MAE_Add+= diff / boostNorm;

            if (diff > MAX_Add) {
                MAX_Add = diff;
                maxAddNion = nionResult;
                maxAddBoost = boostResult;
                max1valueAdd = boost1;
                max2valueAdd = boost2;
            }
        }

        // conjugate
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = conj(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::conj(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            conjBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Conj += diff;
            MAE_Conj+= diff / boostNorm;

            if (diff > MAX_Conj) {
                MAX_Conj = diff;
               maxConjNion = nionResult;
               maxConjBoost = boostResult;
                max1valueConj = boost1;
            }
        }

        // multiplication
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 * nion2;
            endNion = std::chrono::high_resolution_clock::now();
            mulNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost1 * boost2;
            endBoost = std::chrono::high_resolution_clock::now();
            mulBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Mul += diff;
            MAE_Mul+= diff / boostNorm;

            if (diff > MAX_Mul) {
                MAX_Mul = diff;
               maxMulNion = nionResult;
               maxMulBoost = boostResult;
                max1valueMul = boost1;
                max2valueMul = boost2;
            }
        }

        // division
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 / nion2;
            endNion = std::chrono::high_resolution_clock::now();
            divNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost1 / boost2;
            endBoost = std::chrono::high_resolution_clock::now();
            divBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Div += diff;
            MAE_Div+= diff / boostNorm;

            if (diff > MAX_Div) {
                MAX_Div = diff;
               maxDivNion = nionResult;
               maxDivBoost = boostResult;
                max1valueDiv = boost1;
                max2valueDiv = boost2;
            }
        }

        // exp
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = exp(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            expNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::exp(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            expBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Exp += diff;
            MAE_Exp+= diff / boostNorm;

            if (diff > MAX_Exp) {
                MAX_Exp = diff;
               maxExpNion = nionResult;
               maxExpBoost = boostResult;
                max1valueExp = boost1;
            }
        }

        // pow
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = pow(nion1, 3);
            endNion = std::chrono::high_resolution_clock::now();
            powNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::pow(boost1, 3);
            endBoost = std::chrono::high_resolution_clock::now();
            powBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Pow += diff;
            MAE_Pow+= diff / boostNorm;

            if (diff > MAX_Pow) {
                MAX_Pow = diff;
               maxPowNion = nionResult;
               maxPowBoost = boostResult;
                max1valuePow = boost1;
            }
        }

        // sin
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sin(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::sin(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            sinBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Sin += diff;
            MAE_Sin+= diff / boostNorm;

            if (diff > MAX_Sin) {
                MAX_Sin = diff;
               maxSinNion = nionResult;
               maxSinBoost = boostResult;
                max1valueSin = boost1;
            }
        }

        // cos
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cos(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            cosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::cos(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            cosBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Cos += diff;
            MAE_Cos+= diff / boostNorm;

            if (diff > MAX_Cos) {
                MAX_Cos = diff;
               maxCosNion = nionResult;
               maxCosBoost = boostResult;
                max1valueCos = boost1;
            }
        }

        // tan
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tan(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::tan(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            tanBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Tan += diff;
            MAE_Tan+= diff / boostNorm;

            if (diff > MAX_Tan) {
                MAX_Tan = diff;
               maxTanNion = nionResult;
               maxTanBoost = boostResult;
                max1valueTan = boost1;
            }
        }

        // sinh
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sinh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::sinh(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            sinhBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Sinh += diff;
            MAE_Sinh+= diff / boostNorm;

            if (diff > MAX_Sinh) {
                MAX_Sinh = diff;
               maxSinhNion = nionResult;
               maxSinhBoost = boostResult;
                max1valueSinh = boost1;
            }
        }

        // cosh
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cosh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            coshNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::cosh(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            coshBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Cosh += diff;
            MAE_Cosh+= diff / boostNorm;

            if (diff > MAX_Cosh) {
                MAX_Cosh = diff;
               maxCoshNion = nionResult;
               maxCoshBoost = boostResult;
                max1valueCosh = boost1;
            }
        }

        // tanh
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tanh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::tanh(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            tanhBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEquaternion<T>(nionResult, boostResult);
            MRE_Tanh += diff;
            MAE_Tanh+= diff / boostNorm;

            if (diff > MAX_Tanh) {
                MAX_Tanh = diff;
               maxTanhNion = nionResult;
               maxTanhBoost = boostResult;
                max1valueTanh = boost1;
            }
        }
    }

    delete[] vals1;
    delete[] vals2;


    /*** Norm ***/
    std::cout << "\n\n ---> Norm <---" << std::endl;
    std::cout << "Average norm time for nion: " << normNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average norm time for quaternion: " << normBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(normNionTimer, normBoostTimer);
    std::cout << "Average norm error for nion: " << MRE_Norm / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Norm / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum norm error between nion and quaternion: " << MAX_Norm << std::endl;
    std::cout << "nion: " << maxNormNion << "\nquaternion: " << maxNormBoost << std::endl;
    std::cout << "input: " << max1value << std::endl;


    /*** Conjugation ***/
    std::cout << "\n\n ---> Conjugation <---" << std::endl;
    std::cout << "Average conjugate time for nion: " << conjNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average conjugate time for quaternion: " << conjBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(conjNionTimer, conjBoostTimer);
    std::cout << "Average conjugate error for nion: " << MRE_Conj / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Conj/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum conjugate error between nion and quaternion: " << MAX_Conj << std::endl;
    std::cout << "nion: " << maxConjNion << "\nquaternion: " << maxConjBoost << std::endl;
    std::cout << "input: " << max1valueConj << std::endl;


    /*** Addition ***/
    std::cout << "\n\n ---> Addition <---" << std::endl;
    std::cout << "Average addition time for nion: " << addNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average addition time for quaternion: " << addBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(addNionTimer, addBoostTimer);
    std::cout << "Average addition error for nion: " << MRE_Add / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Add/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum addition error between nion and quaternion: " << MAX_Add << std::endl;
    std::cout << "nion: " << maxAddNion << "\nquaternion: " << maxAddBoost << std::endl;
    std::cout << "input1: " << max1valueAdd << std::endl;
    std::cout << "input2: " << max2valueAdd << std::endl;


    /*** Multiplication ***/
    std::cout << "\n\n ---> Multiplication <---" << std::endl;
    std::cout << "Average multiplication time for nion: " << mulNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average multiplication time for quaternion: " << mulBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(mulNionTimer, mulBoostTimer);
    std::cout << "Average multiplication error for nion: " << MRE_Mul / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Mul/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum multiplication error between nion and quaternion: " << MAX_Mul << std::endl;
    std::cout << "nion: " << maxMulNion << "\nquaternion: " << maxMulBoost << std::endl;
    std::cout << "input1: " << max1valueMul << std::endl;
    std::cout << "input2: " << max2valueMul << std::endl;


    /*** Division ***/
    std::cout << "\n\n ---> Division <---" << std::endl;
    std::cout << "Average division time for nion: " << divNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average division time for quaternion: " << divBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(divNionTimer, divBoostTimer);
    std::cout << "Average division error for nion: " << MRE_Div / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Div/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum division error between nion and quaternion: " << MAX_Div << std::endl;
    std::cout << "nion: " << maxDivNion << "\nquaternion: " << maxDivBoost << std::endl;
    std::cout << "input1: " << max1valueDiv << std::endl;
    std::cout << "input2: " << max2valueDiv << std::endl;


    /*** Exponential ***/
    std::cout << "\n\n ---> Exponential <---" << std::endl;
    std::cout << "Average exp time for nion: " << expNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average exp time for quaternion: " << expBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(expNionTimer, expBoostTimer);
    std::cout << "Average exp error for nion: " << MRE_Exp / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Exp/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum exp error between nion and quaternion: " << MAX_Exp << std::endl;
    std::cout << "nion: " << maxExpNion << "\nquaternion: " << maxExpBoost << std::endl;
    std::cout << "input: " << max1valueExp << std::endl;


    /*** Power ***/
    std::cout << "\n\n ---> Power <---" << std::endl;
    std::cout << "Average pow time for nion: " << powNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average pow time for quaternion: " << powBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(powNionTimer, powBoostTimer);
    std::cout << "Average pow error for nion: " << MRE_Pow / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Pow/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum pow error between nion and quaternion: " << MAX_Pow << std::endl;
    std::cout << "nion: " << maxPowNion << "\nquaternion: " << maxPowBoost << std::endl;
    std::cout << "input: " << max1valuePow << std::endl;


    /*** Sine ***/
    std::cout << "\n\n ---> Sine <---" << std::endl;
    std::cout << "Average sin time for nion: " << sinNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average sin time for quaternion: " << sinBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(sinNionTimer, sinBoostTimer);
    std::cout << "Average sin error for nion: " << MRE_Sin / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Sin/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum sin error between nion and quaternion: " << MAX_Sin << std::endl;
    std::cout << "nion: " << maxSinNion << "\nquaternion: " << maxSinBoost << std::endl;
    std::cout << "input: " << max1valueSin << std::endl;


    /*** Cosine ***/
    std::cout << "\n\n ---> Cosine <---" << std::endl;
    std::cout << "Average cos time for nion: " << cosNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average cos time for quaternion: " << cosBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(cosNionTimer, cosBoostTimer);
    std::cout << "Average cos error for nion: " << MRE_Cos / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Cos/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum cos error between nion and quaternion: " << MAX_Cos << std::endl;
    std::cout << "nion: " << maxCosNion << "\nquaternion: " << maxCosBoost << std::endl;
    std::cout << "input: " << max1valueCos << std::endl;


    /*** Tangent ***/
    std::cout << "\n\n ---> Tangent <---" << std::endl;
    std::cout << "Average tan time for nion: " << tanNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average tan time for quaternion: " << tanBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(tanNionTimer, tanBoostTimer);
    std::cout << "Average tan error for nion: " << MRE_Tan / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Tan/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum tan error between nion and quaternion: " << MAX_Tan << std::endl;
    std::cout << "nion: " << maxTanNion << "\nquaternion: " << maxTanBoost << std::endl;
    std::cout << "input: " << max1valueTan << std::endl;


    /*** Hyperbolic Sine ***/
    std::cout << "\n\n ---> Hyperbolic Sine <---" << std::endl;
    std::cout << "Average sinh time for nion: " << sinhNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average sinh time for quaternion: " << sinhBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(sinhNionTimer, sinhBoostTimer);
    std::cout << "Average sinh error for nion: " << MRE_Sinh / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Sinh/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum sinh error between nion and quaternion: " << MAX_Sinh << std::endl;
    std::cout << "nion: " << maxSinhNion << "\nquaternion: " << maxSinhBoost << std::endl;
    std::cout << "input: " << max1valueSinh << std::endl;


    /*** Hyperbolic Cosine ***/
    std::cout << "\n\n ---> Hyperbolic Cosine <---" << std::endl;
    std::cout << "Average cosh time for nion: " << coshNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average cosh time for quaternion: " << coshBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(coshNionTimer, coshBoostTimer);
    std::cout << "Average cosh error for nion: " << MRE_Cosh / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Cosh/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum cosh error between nion and quaternion: " << MAX_Cosh << std::endl;
    std::cout << "nion: " << maxCoshNion << "\nquaternion: " << maxCoshBoost << std::endl;
    std::cout << "input: " << max1valueCosh << std::endl;


    /*** Hyperbolic Tangent ***/
    std::cout << "\n\n ---> Hyperbolic Tangent <---" << std::endl;
    std::cout << "Average tanh time for nion: " << tanhNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average tanh time for quaternion: " << tanhBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupQuaternion(tanhNionTimer, tanhBoostTimer);
    std::cout << "Average tanh error for nion: " << MRE_Tanh / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and quaternion: " << MAE_Tanh/ static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum tanh error between nion and quaternion: " << MAX_Tanh << std::endl;
    std::cout << "nion: " << maxTanhNion << "\nquaternion: " << maxTanhBoost << std::endl;
    std::cout << "input: " << max1valueTanh << std::endl;
}

#endif //NION_BOOSTQUATERNIONTEST_HPP
