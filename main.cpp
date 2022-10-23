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
#include "complexTest.hpp"
#include "quaternionTest.hpp"
#include "octonionTest.hpp"
#include <complex>


void order2Test() {

    std::cout << "\n####### Selected Examples #######" << std::endl;


    nion<long double> n1({42, 69});
    std::cout << "n1:\t" << n1 << std::endl;

    nion<long double> n2({69, 42});
    std::cout << "n2:\t" << n2 << std::endl;

    std::cout << "\nnegation" << std::endl;
    std::cout << "n1:\t" << -n1 << std::endl;
    std::cout << "n2:\t" << -n2 << std::endl;

    std::cout << "\nconjugates" << std::endl;
    std::cout << "n1 conjugate:\t" << n1.conj() << std::endl;
    std::cout << "n2 conjugate:\t" << n2.conj() << std::endl;

    std::cout << "\n inverse" << std::endl;
    std::cout << "n1 inverse:\t" << n1.inv()
              << std::endl; // go here: https://www.wolframalpha.com/input/?i=1%2F(42%2B69i)
    std::cout << "n2 inverse:\t" << n2.inv()
              << std::endl; // go here: https://www.wolframalpha.com/input/?i=1%2F(69%2B42i)

    std::cout << "\nnorms" << std::endl;
    std::cout << "n1 norm:\t" << n1.norm() << std::endl;
    std::cout << "n2 norm:\t" << n2.norm() << std::endl;

    std::cout << "\naddition" << std::endl;
    std::cout << "n1 + n2:\t" << n1 + n2 << std::endl;
    std::cout << "n1 - n2:\t" << n1 - n2 << std::endl;
    nion<long double> n3(2);
    n3 += n1;
    n3 -= n2;
    std::cout << "n3 = 0 + n1 - n2:\t" << n3 << std::endl;


    std::cout << "\nmultiplication" << std::endl;
    std::cout << "n1 * n2:\t" << (n1 * n2)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=N%5B%2842%2B69i%29*%2869%2B42i%29%5D
    std::cout << "n1 / n1:\t" << (n1 / n1)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%2F%2842%2B69i%29
    std::cout << "n1 / n2:\t" << (n1 / n2)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%2F%2869%2B42i%29
    nion<long double> n4(n1);
    n4 *= n2;
    std::cout << "n4 = n1 * n2:\t" << n4 << std::endl;
    n4 /= n2;
    std::cout << "n4 = (n1 * n2) / n2:\t" << n4 << std::endl;


    std::cout << "dot product:\t" << dot(n1, n2) << std::endl;


    std::cout << "\n#### ALGEBRA ####\n" << std::endl;


    std::cout << "\n\t power functions" << std::endl;
    std::cout << "pow(n1, 0):\t" << pow(n1, 0)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E0
    std::cout << "pow(n1, 1):\t" << pow(n1, 1) << std::endl;
    std::cout << "pow(n1, 2):\t" << pow(n1, 2)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E2
    std::cout << "pow(n1, 3):\t" << pow(n1, 3)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E3
    std::cout << "pow(n1, -1):\t" << pow(n1, -1)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E-1
    std::cout << "sqrt(n1):\t" << sqrt(n1)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=sqrt%2842%2B69i%29
    std::cout << "cbrt(n1):\t" << cbrt(n1)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=cbrt%2842%2B69i%29

    std::cout << "\n\t test power functions" << std::endl;
    std::cout << "n1 * pow(n1, -1):\t" << n1 * pow(n1, -1)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29*%2842%2B69i%29%5E-1
    std::cout << "pow(n1, 2) * pow(n1, -2):\t" << pow(n1, 2) * pow(n1, -2)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E2*%2842%2B69i%29%5E-2
    std::cout << "pow(n1, 2) / pow(n1, 2):\t" << pow(n1, 2) / pow(n1, 2)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E2/%2842%2B69i%29%5E2
    std::cout << "sqrt(pow(n1, 2)):\t" << sqrt(pow(n1, 2))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=sqrt%28%2842%2B69i%29%5E2%29
    std::cout << "sqrt(n1) * sqrt(n1):\t" << sqrt(n1) * sqrt(n1)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=sqrt%2842%2B69i%29*sqrt%2842%2B69i%29
    std::cout << "cbrt(n1) * cbrt(n1) * cbrt(n1):\t" << cbrt(n1) * cbrt(n1) * cbrt(n1)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=cbrt%2842%2B69i%29*cbrt%2842%2B69i%29*cbrt%2842%2B69i%29
    std::cout << "pow(n1, n2):\t" << pow(n1, n2)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E%2869%2B42i%29
    std::cout << "pow(n1, -n2):\t" << pow(n1, -n2)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E-%2869%2B42i%29
    std::cout << "pow(n1, n2) * pow(n1, -n2):\t" << pow(n1, n2) * pow(n1, -n2)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=%2842%2B69i%29%5E%2869%2B42i%29*%2842%2B69i%29%5E-%2869%2B42i%29

    std::cout << "\n\t transcendental functions" << std::endl;
    nion<long double> expo = exp(n1); // go here: https://www.wolframalpha.com/input?i=exp%2842%2B69i%29
    nion<long double> ln = log(n1); // go here: https://www.wolframalpha.com/input?i=log%2842%2B69i%29
    std::cout << "exp(n1):\t" << expo << std::endl;
    std::cout << "exp(std::complex)" << std::complex<long double>(42, 69) << ":\t"
              << exp(std::complex<long double>(42, 69)) << std::endl;
    std::cout << "log(n1):\t" << ln << std::endl;

    std::cout << "\n\t test transcendental functions" << std::endl;
    std::cout << "exp(log(n1)):\t" << exp(ln)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=exp%28log%2842%2B69i%29%29
    std::cout << "log(exp(n1)):\t" << log(expo)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=log%28exp%2842%2B69i%29%29

    nion<long double> n5({1, -2});
    std::cout << "n5:\t" << n5 << std::endl;
    std::cout << "\n\t trigonometric functions" << std::endl;
    std::cout << "sin(n5):\t" << sin(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=sin%281-2i%29
    std::cout << "cos(n5):\t" << cos(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=cos%281-2i%29
    std::cout << "tan(n5):\t" << tan(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=tan%281-2i%29
    std::cout << "cot(n5):\t" << cot(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=cot%281-2i%29
    std::cout << "sec(n5):\t" << sec(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=sec%281-2i%29
    std::cout << "csc(n5):\t" << csc(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=csc%281-2i%29

    std::cout << "\n\t inverse trigonometric functions" << std::endl;
    std::cout << "asin(n5):\t" << asin(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=asin%281-2i%29
    std::cout << "acos(n5):\t" << acos(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=acos%281-2i%29
    std::cout << "atan(n5):\t" << atan(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=atan%281-2i%29
    std::cout << "acot(n5):\t" << acot(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=acot%281-2i%29
    std::cout << "asec(n5):\t" << asec(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=asec%281-2i%29
    std::cout << "acsc(n5):\t" << acsc(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=acsc%281-2i%29

    std::cout << "\n\t test inverse trigonometric functions" << std::endl;
    std::cout << "asin(sin(n5)):\t" << asin(sin(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=asin%28sin%281-2i%29%29
    std::cout << "acos(cos(n5)):\t" << acos(cos(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=acos%28cos%281-2i%29%29
    std::cout << "atan(tan(n5)):\t" << atan(tan(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=atan%28tan%281-2i%29%29
    std::cout << "acot(cot(n5)):\t" << acot(cot(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=acot%28cot%281-2i%29%29
    std::cout << "asec(sec(n5)):\t" << asec(sec(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=asec%28sec%281-2i%29%29
    std::cout << "acsc(csc(n5)):\t" << acsc(csc(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=acsc%28csc%281-2i%29%29

    std::cout << "\n\t hyperbolic trigonometric functions" << std::endl;
    std::cout << "sinh(n1):\t" << sinh(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=sinh%281-2i%29
    std::cout << "cosh(n1):\t" << cosh(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=cosh%281-2i%29
    std::cout << "tanh(n1):\t" << tanh(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=tanh%281-2i%29
    std::cout << "coth(n1):\t" << coth(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=coth%281-2i%29
    std::cout << "sech(n1):\t" << sech(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=sech%281-2i%29
    std::cout << "csch(n1):\t" << csch(n5) << std::endl; // go here: https://www.wolframalpha.com/input?i=csch%281-2i%29

    std::cout << "\n\t inverse hyperbolic trigonometric functions" << std::endl;
    std::cout << "asinh(n5):\t" << asinh(n5)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=asinh%281-2i%29
    std::cout << "acosh(n5):\t" << acosh(n5)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=acosh%281-2i%29
    std::cout << "atanh(n5):\t" << atanh(n5)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=atanh%281-2i%29
    std::cout << "acoth(n5):\t" << acoth(n5)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=acoth%281-2i%29
    std::cout << "asech(n5):\t" << asech(n5)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=asech%281-2i%29
    std::cout << "acsch(n5):\t" << acsch(n5)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=acsch%281-2i%29

    std::cout << "\n\t test inverse hyperbolic trigonometric functions" << std::endl;
    std::cout << "asinh(sinh(n5)):\t" << asinh(sinh(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=asinh%28sinh%281-2i%29%29
    std::cout << "acosh(cosh(n5)):\t" << acosh(cosh(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=acosh%28cosh%281-2i%29%29
    std::cout << "atanh(tanh(n5)):\t" << atanh(tanh(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=atanh%28tanh%281-2i%29%29
    std::cout << "acoth(coth(n5)):\t" << acoth(coth(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=acoth%28coth%281-2i%29%29
    std::cout << "asech(sech(n5)):\t" << asech(sech(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=asech%28sech%281-2i%29%29
    std::cout << "acsch(csch(n5)):\t" << acsch(csch(n5))
              << std::endl; // go here: https://www.wolframalpha.com/input?i=acsch%28csch%281-2i%29%29


    std::cout << "\n\t gamma function" << std::endl;
    std::cout << "gamma(n5):\t" << gamma(n5)
              << std::endl; // go here: https://www.wolframalpha.com/input?i=gamma%281-2i%29
}

void mixedOrderTest() {
    std::cout << "\n####### Testing multiplications with non-powers of twos and different orders" << std::endl;

    nion<long double> n6({1, 2});
    std::cout << "n6:\t" << n6 << std::endl;

    nion<long double> n7({3, 4, 1, 2});
    std::cout << "n7:\t" << n7 << std::endl;

    std::cout << "n6 + n7:\t" << n6 + n7 << std::endl;
    std::cout << "n6 - n7:\t" << n6 - n7 << std::endl;
    std::cout << "n6 * n7:\t" << n6 * n7 << std::endl;
    std::cout << "n6 / n7:\t" << n6 / n7 << std::endl;
    std::cout << "n6 * n7 / n6:\t" << n6 * n7 / n6 << std::endl;
    std::cout << "n7 * n6:\t" << n7 * n6 << std::endl;
    std::cout << "n7 / n6:\t" << n7 / n6 << std::endl;
    std::cout << "n7 * n6 / n7:\t" << n7 * n6 / n7 << std::endl;

    std::cout << "\n\t non 2 power order test\n" << std::endl;
    nion<long double> n8({1, 2, 3});
    std::cout << "n8:\t" << n8 << std::endl;
    std::cout << "n8 + n8:\t" << n8 + n8 << std::endl;
    std::cout << "n8 - n8:\t" << n8 - n8 << std::endl;
    std::cout << "n8 * n8:\t" << n8 * n8 << std::endl;
    std::cout << "n8 / n8:\t" << n8 / n8 << std::endl;
    std::cout << "n8 * n8 / n8:\t" << n8 * n8 / n8 << std::endl;
    std::cout << "n8 / n8 * n8:\t" << n8 / n8 * n8 << std::endl;
    std::cout << "n8 + n7:\t" << n8 + n7 << std::endl;
    std::cout << "n8 - n7:\t" << n8 - n7 << std::endl;
    std::cout << "n8 * n7:\t" << n8 * n7 << std::endl;
    std::cout << "n8 / n7:\t" << n8 / n7 << std::endl;
    std::cout << "n8 * n7 / n8:\t" << n8 * n7 / n8 << std::endl;
    std::cout << "n7 * n8 / n7:\t" << n7 * n8 / n7 << std::endl;
}

void writeSinData(int N) {
    long double h = 1.0 / (N + 1);
    std::ofstream myfile("sin_data.txt");
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; ++j) {
            long double x = (10 * i) * h - 5;
            long double y = (4 * j) * h - 2;
            myfile << sin(nion < long double > ({ x, y }))[0] << " ";
        }
        myfile << std::endl;
    }
    myfile.close();
}

int main() {

    std::cout.precision(10);
//

    std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%% nion complex number library benchmarks %%%%%%%%%%%%%%%%%%%%%%%%%%\n\n" << std::endl;
    order2Test();
    mixedOrderTest();

    int trials = 100000000;
    stdComplexComparison<long double>(trials);
    boostQuaternionComparison<long double>(trials);
    boostOctonionComparison<long double>(trials);

    std::cout << "\n\n%%%%%%%%%%%%%%%%%%%%%%%%%% nion complex number library benchmarks %%%%%%%%%%%%%%%%%%%%%%%%%%\n" << std::endl;

    int N = 100;
//    writeSinData(100);

    return 0;
}

