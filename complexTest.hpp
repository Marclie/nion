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

#ifndef NION_COMPLEXTEST_HPP
#define NION_COMPLEXTEST_HPP

#include <fstream>
#include <chrono>
#include <random>
#include <complex>
#include "nion.hpp"
#include <iostream>

template<typename T>
void printSpeedupComplex(const T& niontime, const T& othertime) {
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
T getMAEcomplex(nion<T> nion, std::complex<T> compare){
    int degree = nion.degree;
    T mae = 0;
    mae += std::pow(nion[0] - compare.real(), 2);
    mae += std::pow(nion[1] - compare.imag(), 2);
    return std::sqrt(mae);
}

template<typename T>
void stdComplexComparison(int trials) {
    std::cout << "\n\n#### Comparing std::complex with nion: " << trials << " trials ####\n" << std::endl;
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());

    // timers for norms
    long double normNionTimer = 0;
    long double normStdTimer = 0;
    T MAE_Norm = 0;
    T MRE_Norm = 0;
    T MAX_Norm = -1;
    T maxNormNion = 0;
    T maxNormStd = 0;
    std::complex<T> max1valueNorm;


    // timers for addition
    long double addNionTimer = 0;
    long double addStdTimer = 0;
    T MAE_Add= 0;
    T MRE_Add = 0;
    T MAX_Add = -1;
    nion<T> maxAddNion;
    std::complex<T> maxAddStd;
    std::complex<T> max1valueAdd;
    std::complex<T> max2valueAdd;

    // timers for conjugate
    long double conjNionTimer = 0;
    long double conjStdTimer = 0;
    T MAE_Conj= 0;
    T MRE_Conj = 0;
    T MAX_Conj = -1;
    nion<T> maxConjNion;
    std::complex<T> maxConjStd;
    std::complex<T> max1valueConj;
    std::complex<T> max2valueConj;

    // timers for multiplication
    long double mulNionTimer = 0;
    long double mulStdTimer = 0;
    T MAE_Mul= 0;
    T MRE_Mul = 0;
    T MAX_Mul = -1;
    nion<T> maxMulNion;
    std::complex<T> maxMulStd;
    std::complex<T> max1valueMul;
    std::complex<T> max2valueMul;

    // timers for division
    long double divNionTimer = 0;
    long double divStdTimer = 0;
    T MAE_Div= 0;
    T MRE_Div = 0;
    T MAX_Div = -1;
    nion<T> maxDivNion;
    std::complex<T> maxDivStd;
    std::complex<T> max1valueDiv;
    std::complex<T> max2valueDiv;

    // timers for power
    long double powNionTimer = 0;
    long double powStdTimer = 0;
    T MAE_Pow= 0;
    T MRE_Pow = 0;
    T MAX_Pow = -1;
    nion<T> maxPowNion;
    std::complex<T> maxPowStd;
    std::complex<T> max1valuePow;
    std::complex<T> max2valuePow;

    // timers for exponential
    long double expNionTimer = 0;
    long double expStdTimer = 0;
    T MAE_Exp= 0;
    T MRE_Exp = 0;
    T MAX_Exp = -1;
    nion<T> maxExpNion;
    std::complex<T> maxExpStd;
    std::complex<T> max1valueExp;
    std::complex<T> max2valueExp;

    // timers for logarithm
    long double logNionTimer = 0;
    long double logStdTimer = 0;
    T MAE_Log= 0;
    T MRE_Log = 0;
    T MAX_Log = -1;
    nion<T> maxLogNion;
    std::complex<T> maxLogStd;
    std::complex<T> max1valueLog;
    std::complex<T> max2valueLog;

    // timers for sin
    long double sinNionTimer = 0;
    long double sinStdTimer = 0;
    T MAE_Sin= 0;
    T MRE_Sin = 0;
    T MAX_Sin = -1;
    nion<T> maxSinNion;
    std::complex<T> maxSinStd;
    std::complex<T> max1valueSin;
    std::complex<T> max2valueSin;

    // timers for asin
    long double asinNionTimer = 0;
    long double asinStdTimer = 0;
    T MAE_Asin= 0;
    T MRE_Asin = 0;
    T MAX_Asin = -1;
    nion<T> maxAsinNion;
    std::complex<T> maxAsinStd;
    std::complex<T> max1valueAsin;
    std::complex<T> max2valueAsin;

    // timers for cos
    long double cosNionTimer = 0;
    long double cosStdTimer = 0;
    T MAE_Cos= 0;
    T MRE_Cos = 0;
    T MAX_Cos = -1;
    nion<T> maxCosNion;
    std::complex<T> maxCosStd;
    std::complex<T> max1valueCos;
    std::complex<T> max2valueCos;

    // timers for acos
    long double acosNionTimer = 0;
    long double acosStdTimer = 0;
    T MAE_Acos= 0;
    T MRE_Acos = 0;
    T MAX_Acos = -1;
    nion<T> maxAcosNion;
    std::complex<T> maxAcosStd;
    std::complex<T> max1valueAcos;
    std::complex<T> max2valueAcos;

    // timers for tan
    long double tanNionTimer = 0;
    long double tanStdTimer = 0;
    T MAE_Tan= 0;
    T MRE_Tan = 0;
    T MAX_Tan = -1;
    nion<T> maxTanNion;
    std::complex<T> maxTanStd;
    std::complex<T> max1valueTan;
    std::complex<T> max2valueTan;

    // timers for atan
    long double atanNionTimer = 0;
    long double atanStdTimer = 0;
    T MAE_Atan= 0;
    T MRE_Atan = 0;
    T MAX_Atan = -1;
    nion<T> maxAtanNion;
    std::complex<T> maxAtanStd;
    std::complex<T> max1valueAtan;
    std::complex<T> max2valueAtan;

    // timers for sinh
    long double sinhNionTimer = 0;
    long double sinhStdTimer = 0;
    T MAE_Sinh= 0;
    T MRE_Sinh = 0;
    T MAX_Sinh = -1;
    nion<T> maxSinhNion;
    std::complex<T> maxSinhStd;
    std::complex<T> max1valueSinh;
    std::complex<T> max2valueSinh;

    // timers for asinh
    long double asinhNionTimer = 0;
    long double asinhStdTimer = 0;
    T MAE_Asinh= 0;
    T MRE_Asinh = 0;
    T MAX_Asinh = -1;
    nion<T> maxAsinhNion;
    std::complex<T> maxAsinhStd;
    std::complex<T> max1valueAsinh;
    std::complex<T> max2valueAsinh;

    // timers for cosh
    long double coshNionTimer = 0;
    long double coshStdTimer = 0;
    T MAE_Cosh= 0;
    T MRE_Cosh = 0;
    T MAX_Cosh = -1;
    nion<T> maxCoshNion;
    std::complex<T> maxCoshStd;
    std::complex<T> max1valueCosh;
    std::complex<T> max2valueCosh;

    // timers for acosh
    long double acoshNionTimer = 0;
    long double acoshStdTimer = 0;
    T MAE_Acosh= 0;
    T MRE_Acosh = 0;
    T MAX_Acosh = -1;
    nion<T> maxAcoshNion;
    std::complex<T> maxAcoshStd;
    std::complex<T> max1valueAcosh;
    std::complex<T> max2valueAcosh;

    // timers for tanh
    long double tanhNionTimer = 0;
    long double tanhStdTimer = 0;
    T MAE_Tanh= 0;
    T MRE_Tanh = 0;
    T MAX_Tanh = -1;
    nion<T> maxTanhNion;
    std::complex<T> maxTanhStd;
    std::complex<T> max1valueTanh;
    std::complex<T> max2valueTanh;

    // timers for atanh
    long double atanhNionTimer = 0;
    long double atanhStdTimer = 0;
    T MAE_Atanh= 0;
    T MRE_Atanh = 0;
    T MAX_Atanh = -1;
    nion<T> maxAtanhNion;
    std::complex<T> maxAtanhStd;
    std::complex<T> max1valueAtanh;
    std::complex<T> max2valueAtanh;


    auto startNion = std::chrono::high_resolution_clock::now();
    auto endNion = std::chrono::high_resolution_clock::now();
    auto startStd = std::chrono::high_resolution_clock::now();
    auto endStd = std::chrono::high_resolution_clock::now();

    nion < T > nionResult;
    std::complex<T> stdResult;
    T diff;

    std::uniform_real_distribution<T> distribution(-10.0, 10.0);
    for (int i = 0; i < trials; ++i) {
        // get random std::complex number
        std::complex<T> complex1(distribution(generator), distribution(generator));
        std::complex<T> complex2(distribution(generator), distribution(generator));

        // assign std::complex to nion
        nion < T > nion_complex1(complex1);
        nion < T > nion_complex2(complex2);


        // norm
        {
            // evaluate nion norm, and time
            startNion = std::chrono::high_resolution_clock::now();
            T nionNorm = nion_complex1.abs();
            endNion = std::chrono::high_resolution_clock::now();
            normNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex norm, and time
            startStd = std::chrono::high_resolution_clock::now();
            T stdNorm = norm(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            normStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ normdifference
            diff = std::abs(nionNorm - stdNorm);
            MRE_Norm += diff;
            diff /= stdNorm;
            MAE_Norm += diff;

            // get max difference between nion and std::complex norms
            if (diff > MAX_Norm) {
                MAX_Norm = diff;
                maxNormNion = nionNorm;
                maxNormStd = stdNorm;
                max1valueNorm = complex1;
            }
        }

        /// addition
        {
            // evaluate nion addition, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion_complex1 + nion_complex2;
            endNion = std::chrono::high_resolution_clock::now();
            addNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex addition, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = complex1 + complex2;
            endStd = std::chrono::high_resolution_clock::now();
            addStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ additiondifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Add += diff;
            diff /= norm(stdResult);
            MAE_Add+= diff;

            // get max difference between nion and std::complex addition
            if (diff > MAX_Add) {
                MAX_Add = diff;
                maxAddNion = nionResult;
                maxAddStd = stdResult;
                max1valueAdd = complex1;
                max2valueAdd = complex2;
            }
        }

        /// conjugate
        {
            // evaluate nion conjugate, and time
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();
            nionResult = nion_complex2.conj();
            endNion = std::chrono::high_resolution_clock::now();
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            startStd = std::chrono::high_resolution_clock::now();
            stdResult = std::conj(complex2);
            endStd = std::chrono::high_resolution_clock::now();
            conjStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ conjugatedifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Conj += diff;
            diff /= norm(stdResult);
            MAE_Conj+= diff;

            // get max difference between nion and std::complex conjugate
            if (diff > MAX_Conj) {
                MAX_Conj = diff;
                maxConjNion = nionResult;
                maxConjStd = stdResult;
                max1valueConj = complex1;
                max2valueConj = complex2;
            }

        }

        /// multiplication
        {
            // evaluate nion multiplication, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion_complex1 * nion_complex2;
            endNion = std::chrono::high_resolution_clock::now();
            mulNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex multiplication, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = complex1 * complex2;
            endStd = std::chrono::high_resolution_clock::now();
            mulStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ multiplicationdifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Mul += diff;
            diff /= norm(stdResult);
            MAE_Mul+= diff;

            // get max difference between nion and std::complex multiplication
            if (diff > MAX_Mul) {
                MAX_Mul = diff;
                maxMulNion = nionResult;
                maxMulStd = stdResult;
                max1valueMul = complex1;
                max2valueMul = complex2;
            }

        }

        /// division
        {
            // evaluate nion division, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion_complex1 / nion_complex2;
            endNion = std::chrono::high_resolution_clock::now();
            divNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex division, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = complex1 / complex2;
            endStd = std::chrono::high_resolution_clock::now();
            divStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ divisiondifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Div += diff;
            diff /= norm(stdResult);
            MAE_Div+= diff;

            // get max difference between nion and std::complex division
            if (diff > MAX_Div) {
                MAX_Div = diff;
                maxDivNion = nionResult;
                maxDivStd = stdResult;
                max1valueDiv = complex1;
                max2valueDiv = complex2;
            }
        }

        /// power
        {
            // evaluate nion power, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = pow(nion_complex1, 5);
            endNion = std::chrono::high_resolution_clock::now();
            powNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex power, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = std::pow(complex1, 5);
            endStd = std::chrono::high_resolution_clock::now();
            powStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ powerdifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Pow += diff;
            diff /= norm(stdResult);
            MAE_Pow+= diff;

            // get max difference between nion and std::complex power
            if (diff > MAX_Pow) {
                MAX_Pow = diff;
                maxPowNion = nionResult;
                maxPowStd = stdResult;
                max1valuePow = complex1;
                max2valuePow = complex2;
            }
        }

        /// exponential
        {
            // evaluate nion exponential, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = exp(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            expNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex exponential, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = std::exp(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            expStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ exponentialdifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Exp += diff;
            diff /= norm(stdResult);
            MAE_Exp+= diff;

            // get max difference between nion and std::complex exponential
            if (diff > MAX_Exp) {
                MAX_Exp = diff;
                maxExpNion = nionResult;
                maxExpStd = stdResult;
                max1valueExp = complex1;
                max2valueExp = complex2;
            }
        }

        /// logarithm
        {
            // evaluate nion logarithm, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = log(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            logNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex logarithm, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = std::log(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            logStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ logarithmdifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Log += diff;
            diff /= norm(stdResult);
            MAE_Log+= diff;

            // get max difference between nion and std::complex logarithm
            if (diff > MAX_Log) {
                MAX_Log = diff;
                maxLogNion = nionResult;
                maxLogStd = stdResult;
                max1valueLog = complex1;
                max2valueLog = complex2;
            }
        }

        /// sine
        {
            // evaluate nion sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sin(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            sinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex sine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = sin(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            sinStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ sinedifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Sin += diff;
            diff /= norm(stdResult);
            MAE_Sin+= diff;

            // get max difference between nion and std::complex sine
            if (diff > MAX_Sin) {
                MAX_Sin = diff;
                maxSinNion = nionResult;
                maxSinStd = stdResult;
                max1valueSin = complex1;
                max2valueSin = complex2;
            }
        }

        /// asine
        {
            // evaluate nion asine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = asin(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            asinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex asine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = asin(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            asinStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ asinedifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Asin += diff;
            diff /= norm(stdResult);
            MAE_Asin+= diff;

            // get max difference between nion and std::complex asine
            if (diff > MAX_Asin) {
                MAX_Asin = diff;
                maxAsinNion = nionResult;
                maxAsinStd = stdResult;
                max1valueAsin = complex1;
                max2valueAsin = complex2;
            }
        }

        /// cosine
        {
            // evaluate nion cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cos(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            cosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex cosine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = cos(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            cosStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ cosinedifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Cos += diff;
            diff /= norm(stdResult);
            MAE_Cos+= diff;

            // get max difference between nion and std::complex cosine
            if (diff > MAX_Cos) {
                MAX_Cos = diff;
                maxCosNion = nionResult;
                maxCosStd = stdResult;
                max1valueCos = complex1;
                max2valueCos = complex2;
            }
        }

        /// acosine
        {
            // evaluate nion acosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = acos(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            acosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex acosine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = acos(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            acosStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ acosinedifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Acos += diff;
            diff /= norm(stdResult);
            MAE_Acos+= diff;

            // get max difference between nion and std::complex acosine
            if (diff > MAX_Acos) {
                MAX_Acos = diff;
                maxAcosNion = nionResult;
                maxAcosStd = stdResult;
                max1valueAcos = complex1;
                max2valueAcos = complex2;
            }
        }

        /// tangent
        {
            // evaluate nion tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tan(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            tanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex tangent, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = tan(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            tanStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ tangentdifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Tan += diff;
            diff /= norm(stdResult);
            MAE_Tan+= diff;

            // get max difference between nion and std::complex tangent
            if (diff > MAX_Tan) {
                MAX_Tan = diff;
                maxTanNion = nionResult;
                maxTanStd = stdResult;
                max1valueTan = complex1;
                max2valueTan = complex2;
            }
        }

        /// atan
        {
            // evaluate nion atan, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = atan(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            atanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex atan, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = atan(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            atanStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ atandifference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Atan += diff;
            diff /= norm(stdResult);
            MAE_Atan+= diff;

            // get max difference between nion and std::complex atan
            if (diff > MAX_Atan) {
                MAX_Atan = diff;
                maxAtanNion = nionResult;
                maxAtanStd = stdResult;
                max1valueAtan = complex1;
                max2valueAtan = complex2;
            }
        }

        /// hyperbolic sine
        {
            // evaluate nion hyperbolic sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = sinh(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            sinhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex hyperbolic sine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = sinh(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            sinhStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ hyperbolicsine difference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Sinh += diff;
            diff /= norm(stdResult);
            MAE_Sinh+= diff;

            // get max difference between nion and std::complex hyperbolic sine
            if (diff > MAX_Sinh) {
                MAX_Sinh = diff;
                maxSinhNion = nionResult;
                maxSinhStd = stdResult;
                max1valueSinh = complex1;
                max2valueSinh = complex2;
            }
        }

        /// hyperbolic cosine
        {
            // evaluate nion hyperbolic cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = cosh(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            coshNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex hyperbolic cosine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = cosh(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            coshStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ hyperboliccosine difference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Cosh += diff;
            diff /= norm(stdResult);
            MAE_Cosh+= diff;

            // get max difference between nion and std::complex hyperbolic cosine
            if (diff > MAX_Cosh) {
                MAX_Cosh = diff;
                maxCoshNion = nionResult;
                maxCoshStd = stdResult;
                max1valueCosh = complex1;
                max2valueCosh = complex2;
            }
        }

        /// hyperbolic tangent
        {
            // evaluate nion hyperbolic tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = tanh(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            tanhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex hyperbolic tangent, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = tanh(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            tanhStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ hyperbolictangent difference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Tanh += diff;
            diff /= norm(stdResult);
            MAE_Tanh+= diff;

            // get max difference between nion and std::complex hyperbolic tangent
            if (diff > MAX_Tanh) {
                MAX_Tanh = diff;
                maxTanhNion = nionResult;
                maxTanhStd = stdResult;
                max1valueTanh = complex1;
                max2valueTanh = complex2;
            }
        }

        /// hyperbolic arc sine
        {
            // evaluate nion hyperbolic arc sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = asinh(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            asinhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex hyperbolic arc sine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = asinh(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            asinhStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ hyperbolicarc sine difference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Asinh += diff;
            diff /= norm(stdResult);
            MAE_Asinh+= diff;

            // get max difference between nion and std::complex hyperbolic arc sine
            if (diff > MAX_Asinh) {
                MAX_Asinh = diff;
                maxAsinhNion = nionResult;
                maxAsinhStd = stdResult;
                max1valueAsinh = complex1;
                max2valueAsinh = complex2;
            }
        }

        /// hyperbolic arc cosine
        {
            // evaluate nion hyperbolic arc cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = acosh(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            acoshNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex hyperbolic arc cosine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = acosh(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            acoshStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ hyperbolicarc cosine difference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Acosh += diff;
            diff /= norm(stdResult);
            MAE_Acosh+= diff;

            // get max difference between nion and std::complex hyperbolic arc cosine
            if (diff > MAX_Acosh) {
                MAX_Acosh = diff;
                maxAcoshNion = nionResult;
                maxAcoshStd = stdResult;
                max1valueAcosh = complex1;
                max2valueAcosh = complex2;
            }
        }

        /// hyperbolic arc tangent
        {
            // evaluate nion hyperbolic arc tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = atanh(nion_complex1);
            endNion = std::chrono::high_resolution_clock::now();
            atanhNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex hyperbolic arc tangent, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = atanh(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            atanhStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to MAE_ hyperbolicarc tangent difference
            diff = getMAEcomplex<long double>(nionResult, stdResult);
            MRE_Atanh += diff;
            diff /= norm(stdResult);
            MAE_Atanh+= diff;

            // get max difference between nion and std::complex hyperbolic arc tangent
            if (diff > MAX_Atanh) {
                MAX_Atanh = diff;
                maxAtanhNion = nionResult;
                maxAtanhStd = stdResult;
                max1valueAtanh = complex1;
                max2valueAtanh = complex2;
            }
        }
    }

    T trialfp = static_cast<T>(trials);

    /*** norm ***/
    std::cout << "----> Norm <---- " << std::endl;
    std::cout << "Average norm time for nion: " << normNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average norm time for std::complex: " << normStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(normNionTimer, normStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Norm / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Norm / trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Norm << std::endl;
    std::cout << "nion: " << maxNormNion << "\nstd::complex: " << maxNormStd << std::endl;
    std::cout << "input: " << max1valueNorm << std::endl;


    /*** addition ***/
    std::cout << "\n\n---> Addition <--- " << std::endl;
    std::cout << "Average addition time for nion: " << addNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average addition time for std::complex: " << addStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(addNionTimer, addStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Add / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Add/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Add << std::endl;
    std::cout << "nion: " << maxAddNion << "\nstd::complex: " << maxAddStd << std::endl;
    std::cout << "input1: " << max1valueAdd << "\ninput2: " << max2valueAdd << std::endl;


    /*** conjugate ***/
    std::cout << "\n\n---> Conjugate <--- " << std::endl;
    std::cout << "Average conjugate time for nion: " << conjNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average conjugate time for std::complex: " << conjStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(conjNionTimer, conjStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Conj / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Conj/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Conj << std::endl;
    std::cout << "nion: " << maxConjNion << "\nstd::complex: " << maxConjStd << std::endl;
    std::cout << "input: " << max2valueConj << std::endl;


    /*** multiplication ***/
    std::cout << "\n\n---> Multiplication <--- " << std::endl;
    std::cout << "Average multiplication time for nion: " << mulNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average multiplication time for std::complex: " << mulStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(mulNionTimer, mulStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Mul / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Mul/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Mul << std::endl;
    std::cout << "nion: " << maxMulNion << "\nstd::complex: " << maxMulStd << std::endl;
    std::cout << "input1: " << max1valueMul << "\ninput2: " << max2valueMul << std::endl;


    /*** division ***/
    std::cout << "\n\n---> Division <--- " << std::endl;
    std::cout << "Average division time for nion: " << divNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average division time for std::complex: " << divStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(divNionTimer, divStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Div / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Div/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Div << std::endl;
    std::cout << "nion: " << maxDivNion << "\nstd::complex: " << maxDivStd << std::endl;
    std::cout << "input1: " << max1valueDiv << "\ninput2: " << max2valueDiv << std::endl;


    /*** power ***/
    std::cout << "\n\n---> Power <--- " << std::endl;
    std::cout << "Average power time for nion: " << powNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average power time for std::complex: " << powStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(powNionTimer, powStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Pow / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Pow/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Pow << std::endl;
    std::cout << "nion: " << maxPowNion << "\nstd::complex: " << maxPowStd << std::endl;
    std::cout << "input1: " << max1valuePow << "\ninput2: " << max2valuePow << std::endl;


    /*** exponential ***/
    std::cout << "\n\n---> Exponential <--- " << std::endl;
    std::cout << "Average exponential time for nion: " << expNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average exponential time for std::complex: " << expStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(expNionTimer, expStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Exp / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Exp/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Exp << std::endl;
    std::cout << "nion: " << maxExpNion << "\nstd::complex: " << maxExpStd << std::endl;
    std::cout << "input: " << max1valueExp << std::endl;


    /*** logarithm ***/
    std::cout << "\n\n---> Logarithm <--- " << std::endl;
    std::cout << "Average logarithm time for nion: " << logNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average logarithm time for std::complex: " << logStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(logNionTimer, logStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Log / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Log/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Log << std::endl;
    std::cout << "nion: " << maxLogNion << "\nstd::complex: " << maxLogStd << std::endl;
    std::cout << "input: " << max1valueLog << std::endl;


    /*** sine ***/
    std::cout << "\n\n---> Sine <--- " << std::endl;
    std::cout << "Average sin time for nion: " << sinNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average sin time for std::complex: " << sinStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(sinNionTimer, sinStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Sin / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Sin/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Sin << std::endl;
    std::cout << "nion: " << maxSinNion << "\nstd::complex: " << maxSinStd << std::endl;
    std::cout << "input: " << max1valueSin << std::endl;


    /***  inverse sine ***/
    std::cout << "\n\n---> Inverse Sine <--- " << std::endl;
    std::cout << "Average asin time for nion: " << asinNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average asin time for std::complex: " << asinStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(asinNionTimer, asinStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Asin / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Asin/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Asin << std::endl;
    std::cout << "nion: " << maxAsinNion << "\nstd::complex: " << maxAsinStd << std::endl;
    std::cout << "input: " << max1valueAsin << std::endl;


    /*** cosine ***/
    std::cout << "\n\n---> Cosine <--- " << std::endl;
    std::cout << "Average cos time for nion: " << cosNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average cos time for std::complex: " << cosStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(cosNionTimer, cosStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Cos / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Cos/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Cos << std::endl;
    std::cout << "nion: " << maxCosNion << "\nstd::complex: " << maxCosStd << std::endl;
    std::cout << "input: " << max1valueCos << std::endl;


    /***  inverse cosine ***/
    std::cout << "\n\n---> Inverse Cosine <--- " << std::endl;
    std::cout << "Average acos time for nion: " << acosNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average acos time for std::complex: " << acosStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(acosNionTimer, acosStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Acos / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Acos/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Acos << std::endl;
    std::cout << "nion: " << maxAcosNion << "\nstd::complex: " << maxAcosStd << std::endl;
    std::cout << "input: " << max1valueAcos << std::endl;


    /*** tangent ***/
    std::cout << "\n\n---> Tangent <--- " << std::endl;
    std::cout << "Average tan time for nion: " << tanNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average tan time for std::complex: " << tanStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(tanNionTimer, tanStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Tan / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Tan/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Tan << std::endl;
    std::cout << "nion: " << maxTanNion << "\nstd::complex: " << maxTanStd << std::endl;
    std::cout << "input: " << max1valueTan << std::endl;


    /***  inverse tangent ***/
    std::cout << "\n\n---> Inverse Tangent <--- " << std::endl;
    std::cout << "Average atan time for nion: " << atanNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average atan time for std::complex: " << atanStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(atanNionTimer, atanStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Atan / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Atan/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Atan << std::endl;
    std::cout << "nion: " << maxAtanNion << "\nstd::complex: " << maxAtanStd << std::endl;
    std::cout << "input: " << max1valueAtan << std::endl;


    /*** hyperbolic sine ***/
    std::cout << "\n\n---> Hyperbolic Sine <--- " << std::endl;
    std::cout << "Average sinh time for nion: " << sinhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average sinh time for std::complex: " << sinhStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(sinhNionTimer, sinhStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Sinh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Sinh/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Sinh << std::endl;
    std::cout << "nion: " << maxSinhNion << "\nstd::complex: " << maxSinhStd << std::endl;
    std::cout << "input: " << max1valueSinh << std::endl;


    /*** inverse hyperbolic sine ***/
    std::cout << "\n\n---> Inverse Hyperbolic Sine <--- " << std::endl;
    std::cout << "Average asinh time for nion: " << asinhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average asinh time for std::complex: " << asinhStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(asinhNionTimer, asinhStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Asinh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Asinh/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Asinh << std::endl;
    std::cout << "nion: " << maxAsinhNion << "\nstd::complex: " << maxAsinhStd << std::endl;
    std::cout << "input: " << max1valueAsinh << std::endl;


    /*** hyperbolic cosine ***/
    std::cout << "\n\n---> Hyperbolic Cosine <--- " << std::endl;
    std::cout << "Average cosh time for nion: " << coshNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average cosh time for std::complex: " << coshStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(coshNionTimer, coshStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Cosh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Cosh/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Cosh << std::endl;
    std::cout << "nion: " << maxCoshNion << "\nstd::complex: " << maxCoshStd << std::endl;
    std::cout << "input: " << max1valueCosh << std::endl;


    /*** inverse hyperbolic cosine ***/
    std::cout << "\n\n---> Inverse Hyperbolic Cosine <--- " << std::endl;
    std::cout << "Average acosh time for nion: " << acoshNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average acosh time for std::complex: " << acoshStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(acoshNionTimer, acoshStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Acosh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Acosh/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Acosh << std::endl;
    std::cout << "nion: " << maxAcoshNion << "\nstd::complex: " << maxAcoshStd << std::endl;
    std::cout << "input: " << max1valueAcosh << std::endl;


    /*** hyperbolic tangent ***/
    std::cout << "\n\n---> Hyperbolic Tangent <--- " << std::endl;
    std::cout << "Average tanh time for nion: " << tanhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average tanh time for std::complex: " << tanhStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(tanhNionTimer, tanhStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Tanh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Tanh/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Tanh << std::endl;
    std::cout << "nion: " << maxTanhNion << "\nstd::complex: " << maxTanhStd << std::endl;
    std::cout << "input: " << max1valueTanh << std::endl;


    /*** inverse hyperbolic tangent ***/
    std::cout << "\n\n---> Inverse Hyperbolic Tangent <--- " << std::endl;
    std::cout << "Average atanh time for nion: " << atanhNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average atanh time for std::complex: " << atanhStdTimer / trialfp << " ns" << std::endl;
    printSpeedupComplex(atanhNionTimer, atanhStdTimer);
    std::cout << "Average difference between nion and std::complex: " << MRE_Atanh / trialfp << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << MAE_Atanh/ trialfp << std::endl;
    std::cout << "\nMaximum difference between nion and std::complex: " << MAX_Atanh << std::endl;
    std::cout << "nion: " << maxAtanhNion << "\nstd::complex: " << maxAtanhStd << std::endl;
    std::cout << "input: " << max1valueAtanh << std::endl;
}

#endif //NION_COMPLEXTEST_HPP
