Nion - A Library for Hyper-Complex Arithmetic
---------------------------------------------

<p align="center">
<img src="NION_LOGO.svg" height="540" alt="logo for nion library">
</p>

Nion is a library for hyper-complex arithmetic. It is written in C++ and
is designed to be portable and fast. It is a header-only library, so
there is no need to compile it. It is licensed under the Apache License.

Introduction
------------
The Nion is a generalization of complex numbers to higher degrees.
Nions form a [division algebra](https://en.wikipedia.org/wiki/Division_algebra) up to degree 3,
and more generally a [Subalgebra](https://en.wikipedia.org/wiki/Subalgebra) of the
[Hypercomplex numbers](https://en.wikipedia.org/wiki/Hypercomplex_number).

The Nion is defined as a [Cayley-Dickson construction].
A nion of degree **n** is defined as the pairing of two nions of degree **n-1**.
Each degree doubles the number of dimensions (i.e. **2<sup>n</sup>**) and removes symmetries from the algebra.

Properties
----------
- The [**Reals**](https://en.wikipedia.org/wiki/Real_number) are a Nion of degree **0**, which is a single real number.
- [**Complex**](https://en.wikipedia.org/wiki/Complex_number) numbers are a Nion of degree **1**, which is the pairing of two reals.
  - Complex numbers are no longer an [ordered field](https://en.wikipedia.org/wiki/Ordered_field).
- The [**Quaternions**](https://en.wikipedia.org/wiki/Quaternion) are a Nion of degree **2**, which is the pairing of two complex numbers.
  - Quaternions lose the [commutative property](https://en.wikipedia.org/wiki/Commutative_property) of multiplication.
- The [**Octonions**](https://en.wikipedia.org/wiki/Octonion) are a Nion of degree **3**, which is the pairing of two quaternions.
  - Octonions lose the [associative property](https://en.wikipedia.org/wiki/Associative_property) of multiplication.
- The [**Sedenions**](https://en.wikipedia.org/wiki/Sedenion) are a Nion of degree **4**, which is the pairing of two octonions.
- Sedenions (and higher order nions) lose the [alternative property](https://en.wikipedia.org/wiki/Alternative_algebra) of multiplication and can have [zero divisors](https://en.wikipedia.org/wiki/Zero_divisor).

Nions support the generalization of complex numbers to even higher arbitrary degrees.
While the higher order nions are not as well known as the lower order nions, they exhibit emergent properties that
have potential applications in physics, computer graphics, and other fields.

Nions have the following properties:

* They form a [Subalgebra](https://en.wikipedia.org/wiki/Subalgebra). That is, they have the following properties:
  * The nions form a [vector space](https://en.wikipedia.org/wiki/Vector_space).
  * The [identity element](https://en.wikipedia.org/wiki/Identity_element) belongs to the set of nions.
  * The product of two nions is a nion (i.e., the nions form a [ring](https://en.wikipedia.org/wiki/Ring_(mathematics))).
  * The [inverse element](https://en.wikipedia.org/wiki/Inverse_element) of a nion is a nion.
  * The product of a nion and its inverse is the identity.

Applications
------------
Nions can be used to [improve the accuracy and efficiency of machine learning](https://doi.org/10.1109%2FACCESS.2020.3014690),
and to [accomplish tasks that cannot be done with real or complex numbers](https://www.sciencedirect.com/science/article/abs/pii/0550321383905849?via%3Dihub).

License
-------
Nion is licensed under the Apache License. See the [LICENSE](./LICENSE) file for more information.



**Table of Contents**

- [Nion] - A Library for Hyper-Complex Arithmetic
- [Usage]
  - [Operations]
  - [Functions]
    - [Algebraic Functions]
    - [Trigonometric Functions]
- [Examples]
- [Contributing]
- [License]
- [How to Cite]
- [References]
- [More Information]


[Usage]: #usage
[Operations]: #operations
[Examples]: #examples
[Contributing]: #contributing
[License]: #license
[How to Cite]: #how-to-cite
[References]: #references
[More Information]: #more-information

## Usage

### Operations

* `nion = nion + nion`
* `nion = nion - nion`
* `nion = nion * nion`
* `nion = nion / nion`
* `nion = nion + scalar`
* `nion = nion - scalar`
* `nion = nion * scalar`
* `nion = nion / scalar`
* `nion++`
* `nion--`

### Functions
[Functions]: #functions

* `nion.conj`
* `nion.inv`
* `nion.abs`
* `nion.norm`
* `nion.real`
* `nion.imag`
* `nion.resize`
* `nion.is_real`

### Algebraic Functions
[Algebraic Functions]: #algebraic-functions
* `exp(nion)`
* `log(nion)`
* `pow(nion)`
* `sqr(nion)`
* `sqrt(nion)`
* `cbrt(nion)`


### Trigonometric / Hyperbolic Functions
[Trigonometric Functions]: #trigonometric-functions
* `sin(nion)`    `sinh(nion)`
* `cos(nion)`    `cosh(nion)`
* `tan(nion)`    `tanh(nion)`
* `cot(nion)`    `coth(nion)`
* `sec(nion)`    `sech(nion)`
* `csc(nion)`    `csch(nion)`

#### Inverse Trigonometric / Hyperbolic Functions
* `asin(nion)`   `asinh(nion)`
* `acos(nion)`   `acosh(nion)`
* `atan(nion)`   `atanh(nion)`
* `acot(nion)`   `acoth(nion)`
* `asec(nion)`   `asech(nion)`
* `acsc(nion)`   `acsch(nion)`


## Examples

### Operations
* Addition

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = {5, 6, 7, 8};
nion<float> c = a + b;
std::cout << "a + b = " << c << std::endl; // (6,8,10,12)
```

* Subtraction

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = {5, 6, 7, 8};
nion<float> c = a - b;
std::cout << "a - b = " << c << std::endl; // (-4,-4,-4,-4)
```

* Multiplication

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = {5, 6, 7, 8};
nion<float> c = a * b;
std::cout << "a * b = " << c << std::endl; // (-60,12,30,24)
```

* Division

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = {5, 6, 7, 8};
nion<float> c = a / b;
std::cout << "a / b = " << c << std::endl; // (0.402299,0.045977,0,0.091954)
```

* Increment

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = ++a;
std::cout << "++a = " << b << std::endl; // (2,2,3,4)
```

* Decrement

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = --a;
std::cout << "--a = " << b << std::endl; // (0,2,3,4)
```

### Nion Functions

* Conjugation

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = a.conj();
std::cout << "a.conj() = " << b << std::endl; // (1,-2,-3,-4)
```

* Inversion

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = a.inv();
std::cout << "a.inv() = " << b << std::endl; // (0.0333333,-0.0666667,-0.1,-0.133333)
```

* Absolute Value

```
nion<float> a = {1, 2, 3, 4};
float b = a.abs();
std::cout << "a.abs() = " << b << std::endl; // 30
```

* Norm

```
nion<float> a = {1, 2, 3, 4};
float b = a.norm();
std::cout << "a.norm() = " << b << std::endl; // 5.47723
```

* Real Component

```
nion<float> a = {1, 2, 3, 4};
float b = a.real();
std::cout << "a.real() = " << b << std::endl; // 1
```

* Imaginary Component

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = a.imag();
std::cout << "a.imag() = " << b << std::endl; // (0,2,3,4)
```

* Resize

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = a.resize(8);
std::cout << "a.resize(8) = " << b << std::endl; // (1,2,3,4,0,0,0,0)
```

* Is Real

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = {2, 0};
std::cout << "a.is_real() = " << a.is_real() << std::endl; // false
std::cout << "b.is_real() = " << b.is_real() << std::endl; // true
```

### Algebraic Functions
* Exponential

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = exp(a);
std::cout << "exp(a) = " << b << std::endl; // (1.69392,-0.78956,-1.18434,-1.57912)
```

* Logarithm

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = log(a);
std::cout << "log(a) = " << b << std::endl; // (1.7006,0.51519,0.772785,1.03038)
```

* Power

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = pow(a, -2);
std::cout << "pow(a, -2) = " << b << std::endl; // (-0.0311111,-0.00444444,-0.00666667,-0.00888889)
```

### Trigonometric Functions
* Sine

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = sin(a);
std::cout << "sin(a) = " << b << std::endl; // (91.7837,21.8865,32.8297,43.773)
```

* Cosine

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = cos(a);
std::cout << "cos(a) = " << b << std::endl; // (58.9336,-34.0862,-51.1293,-68.1724)
```

* Tangent

```
nion<float> a = {1, 2, 3, 4};
nion<float> b = tan(a);
std::cout << "tan(a) = " << b << std::endl; // (3.82662e-05,0.371397,0.557096,0.742794)
```


## Contributing

Please feel free to contribute.

## License

This project is licensed under the Apache License, Version 2.0 - see the [LICENSE.md] file for details

[LICENSE.md]: ./LICENSE.md

## How to cite


When using Nion for research projects, please cite:

```
@misc{NionAlgebra,
    author       = {Marcus Dante Liebenthal},
    title        = {{Nion}: A Library for Hyper-Complex Arithmetic},
    month        = {December},
    year         = {2022},
    url          = {https://github.com/Marclie/nion} 
}
```


## References

1. [Cayley-Dickson construction]
1. [Algebraic Properties of the Quaternions and Octonions]
1. [Division Algebra]
1. [Frobenius Algebra]
1. [Algebraic Group]

[Cayley-Dickson construction]: https://en.wikipedia.org/wiki/Cayley–Dickson_construction
[Algebraic Properties of the Quaternions and Octonions]: https://www.math.auckland.ac.nz/~observato/AlgPropQuat.pdf
[Division Algebra]: https://en.wikipedia.org/wiki/Division_algebra
[Frobenius Algebra]: https://en.wikipedia.org/wiki/Frobenius_algebra
[Algebraic Group]: https://en.wikipedia.org/wiki/Algebraic_group