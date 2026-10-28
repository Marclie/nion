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

#ifndef NION_BOOSTOCTONIONTEST_HPP
#define NION_BOOSTOCTONIONTEST_HPP

#include <fstream>
#include <chrono>
#include <random>
#include <boost/math/octonion.hpp>
#include "../nion.hpp"
#include <iostream>


template<typename T>
void printSpeedupOctonion(const T& niontime, const T& othertime) {
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
T getMAEoctonion(nion<T,8> nion, boost::math::octonion<T> compare){
    int degree = nion.size();
    T mae = 0;
    mae += std::pow(nion[0] - compare.R_component_1(), 2);
    mae += std::pow(nion[1] - compare.R_component_2(), 2);
    mae += std::pow(nion[2] - compare.R_component_3(), 2);
    mae += std::pow(nion[3] - compare.R_component_4(), 2);
    mae += std::pow(nion[4] - compare.R_component_5(), 2);
    mae += std::pow(nion[5] - compare.R_component_6(), 2);
    mae += std::pow(nion[6] - compare.R_component_7(), 2);
    mae += std::pow(nion[7] - compare.R_component_8(), 2);
    return std::sqrt(mae);
}

template <typename T>
void boostOctonionComparison(int trials){

    std::cout << "\n\n#### Comparing boost::math::octonion with nion: " << trials << " trials ####\n" << std::endl;
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());

    // timers for norms
    T normNionTimer = 0;
    T normBoostTimer = 0;
    T MRE_Norm = 0;
    T MAE_Norm = 0;
    T MAX_Norm = -1;
    T maxNormNion;
    T maxNormBoost;
    boost::math::octonion<T> max1value;
    boost::math::octonion<T> max2value;

    // timers for addition
    T addNionTimer = 0;
    T addBoostTimer = 0;
    T MRE_Add = 0;
    T MAE_Add= 0;
    T MAX_Add = -1;
    nion<T,8> maxAddNion;
    boost::math::octonion<T> maxAddBoost;
    boost::math::octonion<T> max1valueAdd;
    boost::math::octonion<T> max2valueAdd;

    // timers for conjugate
    T conjNionTimer = 0;
    T conjBoostTimer = 0;
    T MRE_Conj = 0;
    T MAE_Conj= 0;
    T MAX_Conj = -1;
    nion<T,8> maxConjNion;
    boost::math::octonion<T> maxConjBoost;
    boost::math::octonion<T> max1valueConj;
    boost::math::octonion<T> max2valueConj;

    // timers for multiplication
    T mulNionTimer = 0;
    T mulBoostTimer = 0;
    T MRE_Mul = 0;
    T MAE_Mul= 0;
    T MAX_Mul = -1;
    nion<T,8> maxMulNion;
    boost::math::octonion<T> maxMulBoost;
    boost::math::octonion<T> max1valueMul;
    boost::math::octonion<T> max2valueMul;

    // timers for division
    T divNionTimer = 0;
    T divBoostTimer = 0;
    T MRE_Div = 0;
    T MAE_Div= 0;
    T MAX_Div = -1;
    nion<T,8> maxDivNion;
    boost::math::octonion<T> maxDivBoost;
    boost::math::octonion<T> max1valueDiv;
    boost::math::octonion<T> max2valueDiv;

    // timers for power
    T powNionTimer = 0;
    T powBoostTimer = 0;
    T MRE_Pow = 0;
    T MAE_Pow= 0;
    T MAX_Pow = -1;
    nion<T,8> maxPowNion;
    boost::math::octonion<T> maxPowBoost;
    boost::math::octonion<T> max1valuePow;
    boost::math::octonion<T> max2valuePow;

    // timers for exponential
    T expNionTimer = 0;
    T expBoostTimer = 0;
    T MRE_Exp = 0;
    T MAE_Exp= 0;
    T MAX_Exp = -1;
    nion<T,8> maxExpNion;
    boost::math::octonion<T> maxExpBoost;
    boost::math::octonion<T> max1valueExp;
    boost::math::octonion<T> max2valueExp;

    // timers for sine
    T sinNionTimer = 0;
    T sinBoostTimer = 0;
    T MRE_Sin = 0;
    T MAE_Sin= 0;
    T MAX_Sin = -1;
    nion<T,8> maxSinNion;
    boost::math::octonion<T> maxSinBoost;
    boost::math::octonion<T> max1valueSin;
    boost::math::octonion<T> max2valueSin;

    // timers for cosine
    T cosNionTimer = 0;
    T cosBoostTimer = 0;
    T MRE_Cos = 0;
    T MAE_Cos= 0;
    T MAX_Cos = -1;
    nion<T,8> maxCosNion;
    boost::math::octonion<T> maxCosBoost;
    boost::math::octonion<T> max1valueCos;
    boost::math::octonion<T> max2valueCos;

    // timers for tangent
    T tanNionTimer = 0;
    T tanBoostTimer = 0;
    T MRE_Tan = 0;
    T MAE_Tan= 0;
    T MAX_Tan = -1;
    nion<T,8> maxTanNion;
    boost::math::octonion<T> maxTanBoost;
    boost::math::octonion<T> max1valueTan;
    boost::math::octonion<T> max2valueTan;

    // timers for sinh
    T sinhNionTimer = 0;
    T sinhBoostTimer = 0;
    T MRE_Sinh = 0;
    T MAE_Sinh= 0;
    T MAX_Sinh = -1;
    nion<T,8> maxSinhNion;
    boost::math::octonion<T> maxSinhBoost;
    boost::math::octonion<T> max1valueSinh;
    boost::math::octonion<T> max2valueSinh;

    // timers for cosh
    T coshNionTimer = 0;
    T coshBoostTimer = 0;
    T MRE_Cosh = 0;
    T MAE_Cosh= 0;
    T MAX_Cosh = -1;
    nion<T,8> maxCoshNion;
    boost::math::octonion<T> maxCoshBoost;
    boost::math::octonion<T> max1valueCosh;
    boost::math::octonion<T> max2valueCosh;

    // timers for tanh
    T tanhNionTimer = 0;
    T tanhBoostTimer = 0;
    T MRE_Tanh = 0;
    T MAE_Tanh= 0;
    T MAX_Tanh = -1;
    nion<T,8> maxTanhNion;
    boost::math::octonion<T> maxTanhBoost;
    boost::math::octonion<T> max1valueTanh;
    boost::math::octonion<T> max2valueTanh;

    auto startNion = std::chrono::high_resolution_clock::now();
    auto endNion = std::chrono::high_resolution_clock::now();
    auto startBoost = std::chrono::high_resolution_clock::now();
    auto endBoost = std::chrono::high_resolution_clock::now();

    nion <T,8> nionResult;
    boost::math::octonion<T> boostResult;
    T diff;

    T *vals1 = new T[8];
    T *vals2 = new T[8];
    std::uniform_real_distribution<T> distribution(-5, 5);
    for (int i = 0; i < trials; ++i) {

        //generate random quaternion numbers
        for (int j = 0; j < 8; ++j) {
            vals1[j] = distribution(generator);
            vals2[j] = distribution(generator);
        }

        nion<T,8> nion1(vals1, 8);
        nion<T,8> nion2(vals2, 8);
        boost::math::octonion<T> boost1(vals1[0], vals1[1], vals1[2], vals1[3], vals1[4], vals1[5], vals1[6], vals1[7]);
        boost::math::octonion<T> boost2(vals2[0], vals2[1], vals2[2], vals2[3], vals2[4], vals2[5], vals2[6], vals2[7]);

        // norm
        {
            startNion = std::chrono::high_resolution_clock::now();
            T nionNorm = nion1.abs();
            endNion = std::chrono::high_resolution_clock::now();
            normNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            T boostNorm = boost::math::norm(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            normBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = fabs(nionNorm - boostNorm);
            if (diff > MAX_Norm) {
                MAX_Norm = diff;
                maxNormNion = nionNorm;
                maxNormBoost = boostNorm;
                max1value = boost1;
            }
            MAE_Norm += diff;
            MRE_Norm += diff / boostNorm;
        }

        // conjugate
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1.conj();
            endNion = std::chrono::high_resolution_clock::now();
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::conj(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            conjBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Conj) {
                MAX_Conj = diff;
                maxConjNion = nionResult;
                maxConjBoost = boostResult;
                max1valueConj = boost1;
            }
            MAE_Conj += diff;
            MRE_Conj += diff / norm(boostResult);
        }

        // add
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 + nion2;
            endNion = std::chrono::high_resolution_clock::now();
            addNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost1 + boost2;
            endBoost = std::chrono::high_resolution_clock::now();
            addBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Add) {
                MAX_Add = diff;
                maxAddNion = nionResult;
                maxAddBoost = boostResult;
                max1valueAdd = boost1;
                max2valueAdd = boost2;
            }
            MAE_Add += diff;
            MRE_Add += diff / norm(boostResult);
        }

        // multiply
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 * nion2;
            endNion = std::chrono::high_resolution_clock::now();
            mulNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost1 * boost2;
            endBoost = std::chrono::high_resolution_clock::now();
            mulBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Mul) {
                MAX_Mul = diff;
                maxMulNion = nionResult;
                maxMulBoost = boostResult;
                max1valueMul = boost1;
                max2valueMul = boost2;
            }
            MAE_Mul += diff;
            MRE_Mul += diff / norm(boostResult);
        }

        // divide
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion1 / nion2;
            endNion = std::chrono::high_resolution_clock::now();
            divNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost1 / boost2;
            endBoost = std::chrono::high_resolution_clock::now();
            divBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Div) {
                MAX_Div = diff;
                maxDivNion = nionResult;
                maxDivBoost = boostResult;
                max1valueDiv = boost1;
                max2valueDiv = boost2;
            }
            MAE_Div += diff;
            MRE_Div += diff / norm(boostResult);
        }

        // exponential
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = exp(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            expNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::exp(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            expBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Exp) {
                MAX_Exp = diff;
                maxExpNion = nionResult;
                maxExpBoost = boostResult;
                max1valueExp = boost1;
            }
            MAE_Exp += diff;
            MRE_Exp += diff / norm(boostResult);
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

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Pow) {
                MAX_Pow = diff;
                maxPowNion = nionResult;
                maxPowBoost = boostResult;
                max1valuePow = boost1;
                max2valuePow = boost2;
            }
            MAE_Pow += diff;
            MRE_Pow += diff / norm(boostResult);
        }

        // sine
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sin(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::sin(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            sinBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Sin) {
                MAX_Sin = diff;
                maxSinNion = nionResult;
                maxSinBoost = boostResult;
                max1valueSin = boost1;
            }
            MAE_Sin += diff;
            MRE_Sin += diff / norm(boostResult);
        }

        // cosine
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cos(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            cosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::cos(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            cosBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Cos) {
                MAX_Cos = diff;
                maxCosNion = nionResult;
                maxCosBoost = boostResult;
                max1valueCos = boost1;
            }
            MAE_Cos += diff;
            MRE_Cos += diff / norm(boostResult);
        }

        // tangent
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tan(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::tan(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            tanBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Tan) {
                MAX_Tan = diff;
                maxTanNion = nionResult;
                maxTanBoost = boostResult;
                max1valueTan = boost1;
            }
            MAE_Tan += diff;
            MRE_Tan += diff / norm(boostResult);
        }

        // hyperbolic sine
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sinh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            sinhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::sinh(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            sinhBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Sinh) {
                MAX_Sinh = diff;
                maxSinhNion = nionResult;
                maxSinhBoost = boostResult;
                max1valueSinh = boost1;
            }
            MAE_Sinh += diff;
            MRE_Sinh += diff / norm(boostResult);
        }

        // hyperbolic cosine
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cosh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            coshNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::cosh(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            coshBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Cosh) {
                MAX_Cosh = diff;
                maxCoshNion = nionResult;
                maxCoshBoost = boostResult;
                max1valueCosh = boost1;
            }
            MAE_Cosh += diff;
            MRE_Cosh += diff / norm(boostResult);
        }

        // hyperbolic tangent
        {
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tanh(nion1);
            endNion = std::chrono::high_resolution_clock::now();
            tanhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startBoost = std::chrono::high_resolution_clock::now();
            boostResult = boost::math::tanh(boost1);
            endBoost = std::chrono::high_resolution_clock::now();
            tanhBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

            diff = getMAEoctonion(nionResult, boostResult);
            if (diff > MAX_Tanh) {
                MAX_Tanh = diff;
                maxTanhNion = nionResult;
                maxTanhBoost = boostResult;
                max1valueTanh = boost1;
            }
            MAE_Tanh += diff;
            MRE_Tanh += diff / norm(boostResult);
        }

    }

    delete[] vals1;
    delete[] vals2;

    /*** Norm ***/
    std::cout << "\n\n ---> Norm <---" << std::endl;
    std::cout << "Average norm time for nion: " << normNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average norm time for octonion: " << normBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(normNionTimer, normBoostTimer);
    std::cout << "Average norm error for nion: " << MAE_Norm / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Norm / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum norm error for nion: " << MAX_Norm << std::endl;
    std::cout << "nion: " << maxNormNion << "\noctonion: " << maxNormBoost << std::endl;
    std::cout << "input: " << max1value << std::endl;

    /*** Conjugation ***/
    std::cout << "\n\n ---> Conjugation <---" << std::endl;
    std::cout << "Average conjugation time for nion: " << conjNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average conjugation time for octonion: " << conjBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(conjNionTimer, conjBoostTimer);
    std::cout << "Average conjugation error for nion: " << MAE_Conj / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Conj / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum conjugation error for nion: " << MAX_Conj << std::endl;
    std::cout << "nion: " << maxConjNion << "\noctonion: " << maxConjBoost << std::endl;
    std::cout << "input: " << max1valueConj << std::endl;

    /*** Addition ***/
    std::cout << "\n\n ---> Addition <---" << std::endl;
    std::cout << "Average addition time for nion: " << addNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average addition time for octonion: " << addBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(addNionTimer, addBoostTimer);
    std::cout << "Average addition error for nion: " << MAE_Add / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Add / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum addition error for nion: " << MAX_Add << std::endl;
    std::cout << "nion: " << maxAddNion << "\noctonion: " << maxAddBoost << std::endl;
    std::cout << "input1: " << max1valueAdd << std::endl;
    std::cout << "input2: " << max2valueAdd << std::endl;

    /*** Multiplication ***/
    std::cout << "\n\n ---> Multiplication <---" << std::endl;
    std::cout << "Average multiplication time for nion: " << mulNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average multiplication time for octonion: " << mulBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(mulNionTimer, mulBoostTimer);
    std::cout << "Average multiplication error for nion: " << MAE_Mul / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Mul / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum multiplication error for nion: " << MAX_Mul << std::endl;
    std::cout << "nion: " << maxMulNion << "\noctonion: " << maxMulBoost << std::endl;
    std::cout << "input1: " << max1valueMul << std::endl;
    std::cout << "input2: " << max2valueMul << std::endl;

    /*** Division ***/
    std::cout << "\n\n ---> Division <---" << std::endl;
    std::cout << "Average division time for nion: " << divNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average division time for octonion: " << divBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(divNionTimer, divBoostTimer);
    std::cout << "Average division error for nion: " << MAE_Div / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Div / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum division error for nion: " << MAX_Div << std::endl;
    std::cout << "nion: " << maxDivNion << "\noctonion: " << maxDivBoost << std::endl;
    std::cout << "input1: " << max1valueDiv << std::endl;
    std::cout << "input2: " << max2valueDiv << std::endl;

    /*** Exponential ***/
    std::cout << "\n\n ---> Exponential <---" << std::endl;
    std::cout << "Average exponential time for nion: " << expNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average exponential time for octonion: " << expBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(expNionTimer, expBoostTimer);
    std::cout << "Average exponential error for nion: " << MAE_Exp / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Exp / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum exponential error for nion: " << MAX_Exp << std::endl;
    std::cout << "nion: " << maxExpNion << "\noctonion: " << maxExpBoost << std::endl;
    std::cout << "input: " << max1valueExp << std::endl;

    /*** Power ***/
    std::cout << "\n\n ---> Power <---" << std::endl;
    std::cout << "Average power time for nion: " << powNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average power time for octonion: " << powBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(powNionTimer, powBoostTimer);
    std::cout << "Average power error for nion: " << MAE_Pow / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Pow / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum power error for nion: " << MAX_Pow << std::endl;
    std::cout << "nion: " << maxPowNion << "\noctonion: " << maxPowBoost << std::endl;
    std::cout << "input: " << max1valuePow << std::endl;

    /*** Sine ***/
    std::cout << "\n\n ---> Sine <---" << std::endl;
    std::cout << "Average sine time for nion: " << sinNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average sine time for octonion: " << sinBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(sinNionTimer, sinBoostTimer);
    std::cout << "Average sine error for nion: " << MAE_Sin / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Sin / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum sine error for nion: " << MAX_Sin << std::endl;
    std::cout << "nion: " << maxSinNion << "\noctonion: " << maxSinBoost << std::endl;
    std::cout << "input: " << max1valueSin << std::endl;

    /*** Cosine ***/
    std::cout << "\n\n ---> Cosine <---" << std::endl;
    std::cout << "Average cosine time for nion: " << cosNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average cosine time for octonion: " << cosBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(cosNionTimer, cosBoostTimer);
    std::cout << "Average cosine error for nion: " << MAE_Cos / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Cos / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum cosine error for nion: " << MAX_Cos << std::endl;
    std::cout << "nion: " << maxCosNion << "\noctonion: " << maxCosBoost << std::endl;
    std::cout << "input: " << max1valueCos << std::endl;

    /*** Tangent ***/
    std::cout << "\n\n ---> Tangent <---" << std::endl;
    std::cout << "Average tangent time for nion: " << tanNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average tangent time for octonion: " << tanBoostTimer / static_cast<T>(trials) << " ns" << std::endl;
    printSpeedupOctonion(tanNionTimer, tanBoostTimer);
    std::cout << "Average tangent error for nion: " << MAE_Tan / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Tan / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum tangent error for nion: " << MAX_Tan << std::endl;
    std::cout << "nion: " << maxTanNion << "\noctonion: " << maxTanBoost << std::endl;
    std::cout << "input: " << max1valueTan << std::endl;

    /*** Hyperbolic Sine ***/
    std::cout << "\n\n ---> Hyperbolic Sine <---" << std::endl;
    std::cout << "Average hyperbolic sine time for nion: " << sinhNionTimer / static_cast<T>(trials) << " ns" << std::endl;
    std::cout << "Average hyperbolic sine time for octonion: " << sinhBoostTimer / static_cast<T>(trials) << " ns"
              << std::endl;
    printSpeedupOctonion(sinhNionTimer, sinhBoostTimer);
    std::cout << "Average hyperbolic sine error for nion: " << MAE_Sinh / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Sinh / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum hyperbolic sine error for nion: " << MAX_Sinh << std::endl;
    std::cout << "nion: " << maxSinhNion << "\noctonion: " << maxSinhBoost << std::endl;
    std::cout << "input: " << max1valueSinh << std::endl;

    /*** Hyperbolic Cosine ***/
    std::cout << "\n\n ---> Hyperbolic Cosine <---" << std::endl;
    std::cout << "Average hyperbolic cosine time for nion: " << coshNionTimer / static_cast<T>(trials) << " ns"
              << std::endl;
    std::cout << "Average hyperbolic cosine time for octonion: " << coshBoostTimer / static_cast<T>(trials) << " ns"
                << std::endl;
    printSpeedupOctonion(coshNionTimer, coshBoostTimer);
    std::cout << "Average hyperbolic cosine error for nion: " << MAE_Cosh / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Cosh / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum hyperbolic cosine error for nion: " << MAX_Cosh << std::endl;
    std::cout << "nion: " << maxCoshNion << "\noctonion: " << maxCoshBoost << std::endl;
    std::cout << "input: " << max1valueCosh << std::endl;

    /*** Hyperbolic Tangent ***/
    std::cout << "\n\n ---> Hyperbolic Tangent <---" << std::endl;
    std::cout << "Average hyperbolic tangent time for nion: " << tanhNionTimer / static_cast<T>(trials) << " ns"
              << std::endl;
    std::cout << "Average hyperbolic tangent time for octonion: " << tanhBoostTimer / static_cast<T>(trials) << " ns"
              << std::endl;
    printSpeedupOctonion(tanhNionTimer, tanhBoostTimer);
    std::cout << "Average hyperbolic tangent error for nion: " << MAE_Tanh / static_cast<T>(trials) << std::endl;
    std::cout << "Average relative difference between nion and octonion: " << MRE_Tanh / static_cast<T>(trials)
              << std::endl;
    std::cout << "\nMaximum hyperbolic tangent error for nion: " << MAX_Tanh << std::endl;
    std::cout << "nion: " << maxTanhNion << "\noctonion: " << maxTanhBoost << std::endl;
    std::cout << "input: " << max1valueTanh << std::endl;



    std::cout << "\n\n###### Comparing sine series expansion for octonions to function: 100 terms ######" << std::endl;
    nion<T,8> n2({1, 2, 3, 4, 5, 6, 7, 8});
    boost::math::octonion<T> b2(1, 2, 3, 4, 5, 6, 7, 8);

    std::cout << "input: " << b2 << std::endl;

    nion<T,8> sinNion({0, 0, 0, 0, 0, 0, 0, 0});
    boost::math::octonion<T> sinBoost(0, 0, 0, 0, 0, 0, 0, 0);
    std::cout << "sin = sum_{n=0}^{\\infty} (-1)^(2n+1) * (x^(2n+1))/(1+2n)!\n\n" << std::endl;

    for (int i = 0; i < 100; ++i) {
        nion<T,8> nionTerm = pow(n2, 2 * i + 1);
        boost::math::octonion<T> boostTerm = pow(b2, 2 * i + 1);
        for (int j = 0; j < (1 + 2 * i); ++j) {
            nionTerm /= (T)((1 + 2 * i)-j);
            boostTerm /= (T)((1 + 2 * i)-j);
        }
        if (i % 2 == 0) {
            sinNion += nionTerm;
            sinBoost += boostTerm;
        } else {
            sinNion -= nionTerm;
            sinBoost -= boostTerm;
        }
    }

    std::cout << "sin(nion series):\t" << sinNion << std::endl;
    std::cout << "sin(boost series):\t" << sinBoost << std::endl;
    std::cout << "\nsin(nion function):\t" << sin(n2) << std::endl;
    std::cout << "sin(boost function):\t" << sin(b2) << std::endl;


    std::cout << "\n\n###### Comparing hyperbolic sine series expansion for octonions to function: 100 terms ######" << std::endl;
    nion<T,8> n({1, 2, 3, 4, 5, 6, 7, 8});
    boost::math::octonion<T> b(1, 2, 3, 4, 5, 6, 7, 8);

    std::cout << "input: " << b << std::endl;

    nion<T,8> sinhNion({0, 0, 0, 0, 0, 0, 0, 0});
    boost::math::octonion<T> sinhBoost(0, 0, 0, 0, 0, 0, 0, 0);
    std::cout << "sinh = sum_{n=0}^{\\infty} (x^(2n+1))/(1+2n)!\n\n" << std::endl;
    for (int i = 0; i < 100; ++i) {
        nion<T,8> nionTerm = pow(n, 2 * i + 1);
        boost::math::octonion<T> boostTerm = pow(b, 2 * i + 1);
        for (int j = 0; j < (1 + 2 * i); ++j) {
            nionTerm /= (T)((1 + 2 * i)-j);
            boostTerm /= (T)((1 + 2 * i)-j);
        }
        sinhNion += nionTerm;
        sinhBoost += boostTerm;
    }

    std::cout << "sinh(nion series):\t" << sinhNion << std::endl;
    std::cout << "sinh(boost series):\t" << sinhBoost << std::endl;
    std::cout << "\nsinh(nion function):\t" << sinh(n) << std::endl;
    std::cout << "sinh(boost function):\t" << sinh(b) << std::endl;

}

#endif //NION_BOOSTOCTONIONTEST_HPP
