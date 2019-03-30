# BernsteinFEM - Library Overview

## Summary

1. [Introduction](#1-Introduction)
2. [Folders](#2-Folders)
3. [Classes](#3-Classes)

## 1. Introduction
BernsteinFEM is a C++ library that has tools to make fast Finite Elements computations, such as computing Stiffness and Mass Matrices and Load Vectors, through the use of a new polynomial basis, the Bernstein-Bézier polynomials.

To learn how the Finite Element computations were formulated, read the other PDF docs in this folder (in portguese only).

## 2. Folders
You'll see ahead an overview of how the project's files are organized in each folder.
###### These folders are listed here in the order they are consumed by some main application
#### 2.1 Elements Folder
Contains definition of the computational (geometrically defined) elements that are implemented in this library.
##### Classes
```C++
enum class Element_t
```
Defines computational element types that can be used by the library. Defined in [Element_t.h](../Elements/Element_t.h).
```C++
template <Element_t El>
class Element
```
Defines elements of the types specified by the template parameter. Defined in [Elements.h](../Elements/Element.h)
```C++
template <Element_t El>
class Permutation
```
Defines objects that are used by the `Element` perform permutation on the vector's and matrix's position. Defined in [Permutations.h](../Elements/Permutations.h)

#### 2.2 Quadra Folder
Defines auxiliary functions to get the integration points and weights of the Gauss-Jacobi and Gauss-Legendre Quadratures.

#### 2.3 Moments Folder
Defines classes for objects that compute the Load Vector using the Bernstein Polynomials, also they are auxiliary or base classes for other classes that perform computation.
###### The Moments classes compute the load vector, and are called this way because this work was done based in the article written by M. Ainsworth, G. Andriamaro and O. Davydov, with the help and orientation of S. M. Gomes and P. R. Devloo.

##### Classes
```C++
template<Element_t El>
class BMoment
```
Base (abstract) class for all Moments classes. Holds most of the data and methods just to not repeat code a lot. Defined in [Moments.h](../Moments/Moments.h)

Every other class defined in this folder are in the same header file and implemented in its own `.cpp` file, they represent the computation of the Moments of the `Element_t` it references through inheritance of the `BMoment<El>` class.

#### 2.4 Mass Folder
Defines classes for objects that compute the Mass Matrix using the Bernstein Polynomial.
##### Classes
```C++
template<Element_t El>
class BMass
```
Base (abstract) class for all Mass classes. Holds most of the data and methods just to not repeat code through all of the classes that compute the Mass Matrix. Defined in [Mass.h](../Mass/MassM.h).

Just as in the Moments folder, every other class in this folder are in the same header file and implemented in its own `.cpp` file, they represent the computation of the Mass Matrix for each of its element types.

#### 2.4 Stiffness Folder
Defines classes for objects that compute the Stiffness Matrix using the Bernstein Polynomial.
##### Classes
```C++
template<Element_t El>
class BStiff
```
Base (abstract) class for all Stiffness classes. Holds most of the data and methods just to not repeat code through all of the classes that compute the Stiffness Matrix. Defined in [Stiffness.h](../Stiffness/StiffM.h).

Just as in the Moments folder, every other class in this folder are in the same header file and implemented in its own `.cpp` file, they represent the computation of the Stiffness Matrix for each of its element types.

#### 2.5 Derivatives Folder
Defines an auxiliary class for the computation of the computation of the Quadrilateral element Stiffness Matrix (and probably of the Cube element in the future too).

#### 2.6 Tests Folder
Contains `.cpp` files with each with a `main()` function to either test the funcionality of the library or to demonstrate its usage.

## 3. Classes

In this section we'll make an overview of what the classes in this library do.

The classes are divided into the folders listed above.

Before explaining what the classes do we shall se a special enumeration used a lot here
#### 3.0 Element_t enum
Defines which are the computational elements used in this library, each representing a geometrical element, they are as follow:
* LinearEl
* QuadrilateralEl
* TriangularEl
* CubeEl
* TetrahedronEl

#### 3.1 Element Class
The `Element` class is a _template_ class, taking as _template_ parameter a [`Element_t` enum](#3.0-Element_t-enum). It has its specialization for each of the enum parameters, each in its own `.cpp` file.


This class is responsible for defining geometrically its element, through the vertices coordinates, making space transformations to the _Master Element_, used for the computations (also obtaining the Jacobian), and positioning each of the element's _domain points_ in the matrices or vectors.

#### 3.2 Permutation Class
The `Permutation` class is another _template_ class, taking as _template_ parameter a [`Element_t` enum](#3.0-Element_t-enum). Also having its specialization for each of the enum parameters.

This class is responsible for preparing a permutation vector of indices for the positioning of each of the element's _domain points_ in the matrices or vectors. It is auxiliary in the `Element` class for this process.

#### 3.3 BMoment Classes
The `BMoment` class is another _template_ class, taking as _template_ parameter a [`Element_t` enum](#3.0-Element_t-enum). This time its specializiations are made through inheritance. Each of the derived classes has the name `BMoment` with a suffix indicating the kind of element it computes (these suffixes carry out to the other classes of computation):
* BMoment1D - LinearEl
* BMoment2DQuad - QuadrilateralEl
* BMoment2DTri - TriangularEl
* BMoment3DCube - CubeEl
* BMoment3DTetra - TetrahedronEl

This class is responsible for storing and computing the Bernstein polynomial _Moments_, which are equivalent to computing the Load Vector's values of an element. 

To obtain the Load Vector one should initialize an object of this class with its constructor, call the `loadFunction` method with the _load function_ and then get the Vector by calling the `computeMoments` method.

#### 3.4 BMass Classes
The `BMass` classes are just like the [`BMoment` classes](#3.3-BMoment-Classes), except that they are used to compute the Mass Matrix of an element. It holds the same suffixes for each type of element.
###### Remark: these classes inherit from the BMoment classes in a manner not usual in OOP.

#### 3.5 BStiff Classes
The `BStiff` classes are also just like the [`BMoment` classes](#3.3-BMoment-Classes), just now they are used to compute the Stiffness Matrix of an element, also holding the same naming suffixes.
###### Remark 1: these classes also inherit from the BMoment classes in a manner not usual in OOP.
###### Remark 2: the QuadrilateralEl (BStiff2DQuad) class, due to a difficulty on making the calculations, has been made in a auxiliary class at the Derivatives folder.

#### 3.6 BEval Classes
The `BEval` class is another _template_ class, taking as _template_ parameter a [`Element_t` enum](#3.0Element_t-enum). Its specializations are made through inheritance.

This class is responsible for evaluating the resulting polynomial, of the FE analysis, in the quadrature points.
The resulting polynomial is obtained by solving a _Linear System_, the resulting vector out of this _system_ is the input to these classes.
