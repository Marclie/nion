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

#include <iostream>
#include "nion.hpp"
#include <complex>
#include <random>
#include <chrono>
#include <boost/math/quaternion.hpp>

void stdQuaternionComparison(int trials);

void order2Test(){

    std::cout << "machine epsilon:\t" << nion<long double>::epsilon << std::endl;


    nion<long double> n1({42,69});
    std::cout  << "n1:\t"<< n1 <<  std::endl;

    nion<long double> n2({69,42});
    std::cout  << "n2:\t"<< n2 <<  std::endl;

    std::cout << "\nnegation" << std::endl;
    std::cout << "n1:\t" << -n1 << std::endl;
    std::cout << "n2:\t" << -n2 << std::endl;

    std::cout << "\nconjugates" << std::endl;
    std::cout  << "n1 conjugate:\t"<< n1.conj() <<  std::endl;
    std::cout  << "n2 conjugate:\t"<< n2.conj() <<  std::endl;

    std::cout << "\n inverse" << std::endl;
    std::cout  << "n1 inverse:\t"<< n1.inv() <<  std::endl; // go here: https://www.wolframalpha.com/input/?i=1%2F(42%2B69i)
    std::cout  << "n2 inverse:\t"<< n2.inv() <<  std::endl; // go here: https://www.wolframalpha.com/input/?i=1%2F(69%2B42i)

    std::cout << "\nnorms" << std::endl;
    std::cout  << "n1 norm:\t"<< n1.norm() <<  std::endl;
    std::cout  << "n2 norm:\t"<< n2.norm() <<  std::endl;

    std::cout << "\naddition" << std::endl;
    std::cout  << "n1 + n2:\t"<< n1 + n2 <<  std::endl;
    std::cout  << "n1 - n2:\t"<< n1 - n2 <<  std::endl;
    nion<long double> n3(2);
    n3 += n1;
    n3 -= n2;
    std::cout  << "n3 = 0 + n1 - n2:\t"<< n3 <<  std::endl;


    std::cout << "\nmultiplication" << std::endl;
    std::cout  << "n1 * n2:\t"<< (n1 * n2) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=N%5B%2842%2B69i%29*%2869%2B42i%29%5D
    std::cout  << "n1 / n1:\t"<< (n1 / n1) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%2F%2842%2B69i%29
    std::cout  << "n1 / n2:\t"<< (n1 / n2) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%2F%2869%2B42i%29
    nion<long double> n4(n1);
    n4 *= n2;
    std::cout  << "n4 = n1 * n2:\t"<< n4 <<  std::endl;
    n4 /= n2;
    std::cout  << "n4 = (n1 * n2) / n2:\t"<< n4 <<  std::endl;


    std::cout  << "dot product:\t"<< dot(n1, n2) <<  std::endl;


    std::cout << "\n#### ALGEBRA ####\n" << std::endl;


    std::cout << "\n\t power functions" << std::endl;
    std::cout << "pow(n1, 0):\t"<< pow(n1, 0) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E0
    std::cout << "pow(n1, 1):\t"<< pow(n1, 1) <<  std::endl;
    std::cout << "pow(n1, 2):\t"<< pow(n1, 2) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E2
    std::cout << "pow(n1, 3):\t"<< pow(n1, 3) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E3
    std::cout << "pow(n1, -1):\t"<< pow(n1, -1) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E-1
    std::cout << "sqrt(n1):\t"<< sqrt(n1) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=sqrt%2842%2B69i%29
    std::cout << "cbrt(n1):\t"<< cbrt(n1) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=cbrt%2842%2B69i%29

    std::cout << "\n\t test power functions" << std::endl;
    std::cout << "n1 * pow(n1, -1):\t"<< n1 * pow(n1, -1) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29*%2842%2B69i%29%5E-1
    std::cout << "pow(n1, 2) * pow(n1, -2):\t"<< pow(n1, 2) * pow(n1, -2) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E2*%2842%2B69i%29%5E-2
    std::cout << "pow(n1, 2) / pow(n1, 2):\t"<< pow(n1, 2) / pow(n1, 2) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E2/%2842%2B69i%29%5E2
    std::cout << "sqrt(pow(n1, 2)):\t"<< sqrt(pow(n1, 2)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=sqrt%28%2842%2B69i%29%5E2%29
    std::cout << "sqrt(n1) * sqrt(n1):\t"<< sqrt(n1) * sqrt(n1) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=sqrt%2842%2B69i%29*sqrt%2842%2B69i%29
    std::cout << "cbrt(n1) * cbrt(n1) * cbrt(n1):\t"<< cbrt(n1) * cbrt(n1) * cbrt(n1) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=cbrt%2842%2B69i%29*cbrt%2842%2B69i%29*cbrt%2842%2B69i%29
    std::cout << "pow(n1, n2):\t"<< pow(n1, n2) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E%2869%2B42i%29
    std::cout << "pow(n1, -n2):\t"<< pow(n1, -n2) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E-%2869%2B42i%29
    std::cout << "pow(n1, n2) * pow(n1, -n2):\t"<< pow(n1, n2) * pow(n1, -n2) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E%2869%2B42i%29*%2842%2B69i%29%5E-%2869%2B42i%29

    std::cout << "\n\t transcendental functions" << std::endl;
    nion<long double> expo = exp(n1); // go here: https://www.wolframalpha.com/input?i=exp%2842%2B69i%29
    nion<long double> ln = log(n1); // go here: https://www.wolframalpha.com/input?i=log%2842%2B69i%29
    std::cout  << "exp(n1):\t"<< expo <<  std::endl;
    std::cout  << "exp(std::complex)" << std::complex<long double>(42, 69) << ":\t"<< exp(std::complex<long double>(42, 69)) <<  std::endl;
    std::cout  << "log(n1):\t"<< ln <<  std::endl;

    std::cout << "\n\t test transcendental functions" << std::endl;
    std::cout  << "exp(log(n1)):\t"<< exp(ln) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=exp%28log%2842%2B69i%29%29
    std::cout  << "log(exp(n1)):\t"<< log(expo) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=log%28exp%2842%2B69i%29%29

    nion<long double>n5({1, -2});
    std::cout << "n5:\t"<< n5 <<  std::endl;
    std::cout << "\n\t trigonometric functions" << std::endl;
    std::cout << "sin(n5):\t"<< sin(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=sin%281-2i%29
    std::cout << "cos(n5):\t"<< cos(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=cos%281-2i%29
    std::cout << "tan(n5):\t"<< tan(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=tan%281-2i%29
    std::cout << "cot(n5):\t"<< cot(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=cot%281-2i%29
    std::cout << "sec(n5):\t"<< sec(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=sec%281-2i%29
    std::cout << "csc(n5):\t"<< csc(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=csc%281-2i%29

    std::cout << "\n\t inverse trigonometric functions" << std::endl;
    std::cout << "asin(n5):\t"<< asin(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=asin%281-2i%29
    std::cout << "acos(n5):\t"<< acos(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acos%281-2i%29
    std::cout << "atan(n5):\t"<< atan(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=atan%281-2i%29
    std::cout << "acot(n5):\t"<< acot(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acot%281-2i%29
    std::cout << "asec(n5):\t"<< asec(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=asec%281-2i%29
    std::cout << "acsc(n5):\t"<< acsc(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acsc%281-2i%29

    std::cout << "\n\t test inverse trigonometric functions" << std::endl;
    std::cout << "asin(sin(n5)):\t"<< asin(sin(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=asin%28sin%281-2i%29%29
    std::cout << "acos(cos(n5)):\t"<< acos(cos(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acos%28cos%281-2i%29%29
    std::cout << "atan(tan(n5)):\t"<< atan(tan(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=atan%28tan%281-2i%29%29
    std::cout << "acot(cot(n5)):\t"<< acot(cot(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acot%28cot%281-2i%29%29
    std::cout << "asec(sec(n5)):\t"<< asec(sec(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=asec%28sec%281-2i%29%29
    std::cout << "acsc(csc(n5)):\t"<< acsc(csc(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acsc%28csc%281-2i%29%29

    std::cout << "\n\t hyperbolic trigonometric functions" << std::endl;
    std::cout << "sinh(n1):\t"<< sinh(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=sinh%281-2i%29
    std::cout << "cosh(n1):\t"<< cosh(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=cosh%281-2i%29
    std::cout << "tanh(n1):\t"<< tanh(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=tanh%281-2i%29
    std::cout << "coth(n1):\t"<< coth(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=coth%281-2i%29
    std::cout << "sech(n1):\t"<< sech(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=sech%281-2i%29
    std::cout << "csch(n1):\t"<< csch(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=csch%281-2i%29

    std::cout << "\n\t inverse hyperbolic trigonometric functions" << std::endl;
    std::cout << "asinh(n5):\t"<< asinh(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=asinh%281-2i%29
    std::cout << "acosh(n5):\t"<< acosh(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acosh%281-2i%29
    std::cout << "atanh(n5):\t"<< atanh(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=atanh%281-2i%29
    std::cout << "acoth(n5):\t"<< acoth(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acoth%281-2i%29
    std::cout << "asech(n5):\t"<< asech(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=asech%281-2i%29
    std::cout << "acsch(n5):\t"<< acsch(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acsch%281-2i%29

    std::cout << "\n\t test inverse hyperbolic trigonometric functions" << std::endl;
    std::cout << "asinh(sinh(n5)):\t"<< asinh(sinh(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=asinh%28sinh%281-2i%29%29
    std::cout << "acosh(cosh(n5)):\t"<< acosh(cosh(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acosh%28cosh%281-2i%29%29
    std::cout << "atanh(tanh(n5)):\t"<< atanh(tanh(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=atanh%28tanh%281-2i%29%29
    std::cout << "acoth(coth(n5)):\t"<< acoth(coth(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acoth%28coth%281-2i%29%29
    std::cout << "asech(sech(n5)):\t"<< asech(sech(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=asech%28sech%281-2i%29%29
    std::cout << "acsch(csch(n5)):\t"<< acsch(csch(n5)) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=acsch%28csch%281-2i%29%29


    std::cout << "\n\t gamma function" << std::endl;
    std::cout << "gamma(n5):\t"<< gamma(n5) <<  std::endl; // go here: https://www.wolframalpha.com/input?i=gamma%281-2i%29
}

template <typename T>
void stdComplexComparison(int trials) {
    std::cout << "\nComparing std::complex with nion: " << trials << " trials\n" << std::endl;
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());

    // timers for norms
    long double normNionTimer = 0;
    long double normStdTimer = 0;
    T totalNorm = 0;

    // timers for addition
    long double addNionTimer = 0;
    long double addStdTimer = 0;
    T totalAdd = 0;

    // timers for conjugate
    long double conjNionTimer = 0;
    long double conjStdTimer = 0;
    T totalConj = 0;

    // timers for multiplication
    long double mulNionTimer = 0;
    long double mulStdTimer = 0;
    T totalMul = 0;

    // timers for division
    long double divNionTimer = 0;
    long double divStdTimer = 0;
    T totalDiv = 0;

    // timers for power
    long double powNionTimer = 0;
    long double powStdTimer = 0;
    T totalPow = 0;

    // timers for exponential
    long double expNionTimer = 0;
    long double expStdTimer = 0;
    T totalExp = 0;

    // timers for logarithm
    long double logNionTimer = 0;
    long double logStdTimer = 0;
    T totalLog = 0;

    // timers for atan(tan)
    long double tanNionTimer = 0;
    long double tanStdTimer = 0;
    T totalTan = 0;

    // timers for asin(sin)
    long double sinNionTimer = 0;
    long double sinStdTimer = 0;
    T totalSin = 0;

    // timers for acos(cos)
    long double cosNionTimer = 0;
    long double cosStdTimer = 0;
    T totalCos = 0;





    auto startNion = std::chrono::high_resolution_clock::now();
    auto endNion = std::chrono::high_resolution_clock::now();
    auto startStd = std::chrono::high_resolution_clock::now();
    auto endStd = std::chrono::high_resolution_clock::now();

    nion<T> nionResult;
    std::complex<T> stdResult;
    T diff;

    std::uniform_real_distribution<T> distribution(-100.0, 100.0);
    for (int i = 0; i < trials; ++i) {
        // get random std::complex number
        std::complex<T> complex1(distribution(generator), distribution(generator));
        std::complex<T> complex2(distribution(generator), distribution(generator));

        // assign std::complex to nion
        nion<T> nion_complex1(complex1);
        nion<T> nion_complex2(complex2);


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

            // get difference between nion and std::complex norms. Add to total norm difference
            diff = nionNorm - stdNorm;
            diff /= stdNorm;
            totalNorm += diff;
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

            // get difference between nion and std::complex norms. Add to total addition difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalAdd += diff;
        }

        /// conjugate
        {
            // evaluate nion conjugate, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = nion_complex1.conj();
            endNion = std::chrono::high_resolution_clock::now();

            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();
            nionResult = nion_complex2.conj();
            endNion = std::chrono::high_resolution_clock::now();
            conjNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex conjugate, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = std::conj(complex1);
            endStd = std::chrono::high_resolution_clock::now();
            conjStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            startStd = std::chrono::high_resolution_clock::now();
            stdResult = std::conj(complex2);
            endStd = std::chrono::high_resolution_clock::now();
            conjStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to total conjugate difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalConj += diff;

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

            // get difference between nion and std::complex norms. Add to total multiplication difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalMul += diff;
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

            // get difference between nion and std::complex norms. Add to total division difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalDiv += diff;
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

            // get difference between nion and std::complex norms. Add to total power difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalPow += diff;
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

            // get difference between nion and std::complex norms. Add to total exponential difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalExp += diff;
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

            // get difference between nion and std::complex norms. Add to total logarithm difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalLog += diff;
        }

        /// tangent
        {
            // evaluate nion tangent, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = atan(tan(nion_complex1));
            endNion = std::chrono::high_resolution_clock::now();
            tanNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex tangent, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = atan(std::tan(complex1));
            endStd = std::chrono::high_resolution_clock::now();
            tanStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to total tangent difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalTan += diff;
        }

        /// sine
        {
            // evaluate nion sine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = asin(sin(nion_complex1));
            endNion = std::chrono::high_resolution_clock::now();
            sinNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex sine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = asin(std::sin(complex1));
            endStd = std::chrono::high_resolution_clock::now();
            sinStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to total sine difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalSin += diff;
        }

        /// cosine
        {
            // evaluate nion cosine, and time
            startNion = std::chrono::high_resolution_clock::now();
            nionResult = acos(cos(nion_complex1));
            endNion = std::chrono::high_resolution_clock::now();
            cosNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

            // evaluate std::complex cosine, and time
            startStd = std::chrono::high_resolution_clock::now();
            stdResult = acos(std::cos(complex1));
            endStd = std::chrono::high_resolution_clock::now();
            cosStdTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endStd - startStd).count();

            // get difference between nion and std::complex norms. Add to total cosine difference
            diff = nionResult.abs() - norm(stdResult);
            diff /= norm(stdResult);
            totalCos += diff;
        }

    }

    T trialfp = static_cast<T>(trials);
    std::cout << "Average norm time for nion: " << addNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average norm time for std::complex: " << addStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalAdd / trialfp << std::endl;

    std::cout << "\nAverage addition time for nion: " << addNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average addition time for std::complex: " << addStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalAdd / trialfp << std::endl;

    std::cout << "\nAverage conjugate time for nion: " << conjNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average conjugate time for std::complex: " << conjStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalConj / trialfp << std::endl;

    std::cout << "\nAverage multiplication time for nion: " << mulNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average multiplication time for std::complex: " << mulStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalMul / trialfp << std::endl;

    std::cout << "\nAverage division time for nion: " << divNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average division time for std::complex: " << divStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalDiv / trialfp << std::endl;

    std::cout << "\nAverage power time for nion: " << powNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average power time for std::complex: " << powStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalPow / trialfp << std::endl;

    std::cout << "\nAverage exponential time for nion: " << expNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average exponential time for std::complex: " << expStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalExp / trialfp << std::endl;

    std::cout << "\nAverage logarithm time for nion: " << logNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average logarithm time for std::complex: " << logStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalLog / trialfp << std::endl;

    std::cout << "\nAverage atan(tan) time for nion: " << tanNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average atan(tan) time for std::complex: " << tanStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalTan / trialfp << std::endl;

    std::cout << "\nAverage asin(sin) time for nion: " << sinNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average asin(sin) time for std::complex: " << sinStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalSin / trialfp << std::endl;

    std::cout << "\nAverage acos(cos) time for nion: " << cosNionTimer / trialfp << " ns" << std::endl;
    std::cout << "Average acos(cos) time for std::complex: " << cosStdTimer / trialfp << " ns" << std::endl;
    std::cout << "Average relative difference between nion and std::complex: " << totalCos / trialfp << std::endl;

}
template <typename T>
void boostQuaternionComparison(int trials) {
    std::cout << "\n\n#### Comparing boost::math::quaternion with nion: " << trials << " trials ####\n" << std::endl;
    std::default_random_engine generator(std::chrono::steady_clock::now().time_since_epoch().count());

    // timers for norms
    long double normNionTimer = 0;
    long double normBoostTimer = 0;
    T totalNorm = 0;

    // timers for addition
    long double addNionTimer = 0;
    long double addBoostTimer = 0;
    T totalAdd = 0;

    // timers for conjugate
    long double conjNionTimer = 0;
    long double conjBoostTimer = 0;
    T totalConj = 0;

    // timers for multiplication
    long double mulNionTimer = 0;
    long double mulBoostTimer = 0;
    T totalMul = 0;

    // timers for division
    long double divNionTimer = 0;
    long double divBoostTimer = 0;
    T totalDiv = 0;

    // timers for power
    long double powNionTimer = 0;
    long double powBoostTimer = 0;
    T totalPow = 0;

    // timers for exponential
    long double expNionTimer = 0;
    long double expBoostTimer = 0;
    T totalExp = 0;

    // timers for logarithm
    long double logNionTimer = 0;
    long double logBoostTimer = 0;
    T totalLog = 0;

    // timers for atan(tan)
    long double tanNionTimer = 0;
    long double tanBoostTimer = 0;
    T totalTan = 0;

    // timers for asin(sin)
    long double sinNionTimer = 0;
    long double sinBoostTimer = 0;
    T totalSin = 0;

    // timers for acos(cos)
    long double cosNionTimer = 0;
    long double cosBoostTimer = 0;
    T totalCos = 0;

    auto startNion = std::chrono::high_resolution_clock::now();
    auto endNion = std::chrono::high_resolution_clock::now();
    auto startBoost = std::chrono::high_resolution_clock::now();
    auto endBoost = std::chrono::high_resolution_clock::now();

    nion<T> nionResult;
    boost::math::quaternion<T> boostResult;
    T diff;
    for (int i = 0; i < trials; ++i) {

        //generate random quaternion numbers
        T* vals1 = new T[4];
        T* vals2 = new T[4];
        for (int j = 0; j < 4; ++j) {
            std::uniform_real_distribution<T> distribution(-100, 100);
            vals1[j] = distribution(generator);
            vals2[j] = distribution(generator);
        }

        nion<T> nion1(vals1, 4);
        nion<T> nion2(vals2, 4);
        boost::math::quaternion<T> boost1(vals1[0], vals1[1], vals1[2], vals1[3]);
        boost::math::quaternion<T> boost2(vals2[0], vals2[1], vals2[2], vals2[3]);

        // norm
        startNion = std::chrono::high_resolution_clock::now();
        nionResult = nion1.abs();
        endNion = std::chrono::high_resolution_clock::now();
        normNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

        startBoost = std::chrono::high_resolution_clock::now();
        boostResult = boost::math::norm(boost1);
        endBoost = std::chrono::high_resolution_clock::now();
        normBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

        diff = std::abs(nionResult.norm() - boost::math::norm(boostResult));
        totalNorm += diff / nionResult.real();

        // addition
        startNion = std::chrono::high_resolution_clock::now();
        nionResult = nion1 + nion2;
        endNion = std::chrono::high_resolution_clock::now();
        addNionTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endNion - startNion).count();

        startBoost = std::chrono::high_resolution_clock::now();
        boostResult = boost1 + boost2;
        endBoost = std::chrono::high_resolution_clock::now();
        addBoostTimer += std::chrono::duration_cast<std::chrono::nanoseconds>(endBoost - startBoost).count();

        diff = std::abs(nionResult.abs() - boost::math::norm(boostResult));
        totalAdd += diff / nionResult.real();
    }

    std::cout << "Average norm time for nion: " << normNionTimer / trials << " ns" << std::endl;
    std::cout << "Average norm time for boost::math::quaternion: " << normBoostTimer / trials << " ns" << std::endl;
    std::cout << "Average relative difference between nion and boost::math::quaternion: " << totalNorm / trials << std::endl;

    std::cout << "\nAverage addition time for nion: " << addNionTimer / trials << " ns" << std::endl;
    std::cout << "Average addition time for boost::math::quaternion: " << addBoostTimer / trials << " ns" << std::endl;
    std::cout << "Average relative difference between nion and boost::math::quaternion: " << totalAdd / trials << std::endl;
}
int main() {

    std::cout.precision(17);

    int trials = 100000;
    std::cout << "nion complex number library" << std::endl;
    order2Test();
    stdComplexComparison<long double>(trials);
    boostQuaternionComparison<long double>(trials);
/*
    std::cout << "\n\t mixed order test\n\tThere may not be references for this behavior, but it works with my recursion" << std::endl;

    nion<long double>n6({1, 2});
    std::cout << "n6:\t"<< n6 <<  std::endl;

    nion<long double>n7({3, 4, 1, 2});
    std::cout << "n7:\t"<< n7 <<  std::endl;

    std::cout << "n6 + n7:\t"<< n6 + n7 <<  std::endl;
    std::cout << "n6 - n7:\t"<< n6 - n7 <<  std::endl;
    std::cout << "n6 * n7:\t"<< n6 * n7 <<  std::endl;
    std::cout << "n6 / n7:\t"<< n6 / n7 <<  std::endl;
    std::cout << "n6 * n7 / n6:\t" << n6 * n7 / n6 <<  std::endl;
    std::cout << "n7 * n6:\t"<< n7 * n6 <<  std::endl;
    std::cout << "n7 / n6:\t"<< n7 / n6 <<  std::endl;
    std::cout << "n7 * n6 / n7:\t" << n7 * n6 / n7 <<  std::endl;

    std::cout << "\n\t non 2 power order test\n\tThere may not be references for this behavior, but it works with my recursion" << std::endl;
    nion<long double> n8({1,2,3});
    std::cout << "n8:\t"<< n8 <<  std::endl;
    std::cout << "n8 + n8:\t"<< n8 + n8 <<  std::endl;
    std::cout << "n8 - n8:\t"<< n8 - n8 <<  std::endl;
    std::cout << "n8 * n8:\t"<< n8 * n8 <<  std::endl;
    std::cout << "n8 / n8:\t"<< n8 / n8 <<  std::endl;
    std::cout << "n8 * n8 / n8:\t" << n8 * n8 / n8 <<  std::endl;
    std::cout << "n8 / n8 * n8:\t" << n8 / n8 * n8 <<  std::endl;
    std::cout << "n8 + n7:\t"<< n8 + n7 <<  std::endl;
    std::cout << "n8 - n7:\t"<< n8 - n7 <<  std::endl;
    std::cout << "n8 * n7:\t"<< n8 * n7 <<  std::endl;
    std::cout << "n8 / n7:\t"<< n8 / n7 <<  std::endl;
    std::cout << "n8 * n7 / n8:\t" << n8 * n7 / n8 <<  std::endl;
    std::cout << "n7 * n8 / n7:\t" << n7 * n8 / n7 <<  std::endl;
*/

    return 0;
}


