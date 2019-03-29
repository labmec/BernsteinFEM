# BernsteinFEM - Library Overview

## 1. Introduction
BernsteinFEM is a C++ library that has tools to make fast Finite Elements computations, such as Stiffness and Mass Matrices and Load Vectors.

## 2. Organization
You'll see ahead an overview of the organization of the files and concepts here implemented
#### 2.1 Folders 
###### These folders are listed here in the order they are consumed by some main application
##### 2.1.1 Elements Folder
Contains definition of the computational (geometrically defined) elements that are implemented in this library.
##### Classes
```C++
enum class Element_t
```
Defines copmutational element types that can be used by the library. Defined in [Element_t.h](../Elements/Element_t.h).
```C++
template <Element_t El>
class Element
```
Defines elements of the types specified by the template parameter. Defined in [Elements.h](../Elements/Elements.h)
```C++
template <Element_t El>
class Permutation
```
Defines objects that are used by the `Element` perform permutation on the vector's and matrix's position. Defined in [Permutations.h](../Elements/Permutations.h)

##### 2.1.2 Quadra Folder
Defines auxiliary functions to get the integration points and weights of the Gauss-Jacobi and Gauss-Legendre Quadratures.

##### 2.1.3 Moments Folder
Defines classes for objects that compute the Load Vector using the Bernstein Polynomials, also they are auxiliary or base classes for other classes that perform computation.
###### The Moments classes compute the load vector, and are called this way because this work was done based in the article written by M. Ainsworth, G. Andriamaro and O. Davydov, with the help and orientation of S. M. Gomes and P. R. Devloo.

##### Classes
```C++
template<Element_t El>
class BMoment
```
Base (abstract) class for all Moments classes. Holds most of the data and methods just to not repeat code a lot. Defined in [Moments.h](../Moments/Moments.h)

Every other class defined in this folder are in the same header file and implemented in its own `.cpp` file, they represent the computation of the Moments of the `Element_t` it references through inheritance of the `BMoment<El>` class.

##### 2.1.4 Mass Folder
Defines classes for objects that compute the Mass Matrix using the Bernstein Polynomial.
##### Classes
```C++
template<Element_t El>
class BMass
```
Base (abstract) class for all Mass classes. Holds most of the data and methods just to not repeat code through all of the classes that compute the Mass Matrix. Defined in [Mass.h](../Mass/MassM.h).

Just as in the Moments folder, every other class in this folder are in the same header file and implemented in its own `.cpp` file, they represent the computation of the Mass Matrix for each of its element types.

##### 2.1.4 Stiffness Folder
Defines classes for objects that compute the Stiffness Matrix using the Bernstein Polynomial.
##### Classes
```C++
template<Element_t El>
class BStiff
```
Base (abstract) class for all Stiffness classes. Holds most of the data and methods just to not repeat code through all of the classes that compute the Stiffness Matrix. Defined in [Stiffness.h](../Stiffness/StiffM.h).

Just as in the Moments folder, every other class in this folder are in the same header file and implemented in its own `.cpp` file, they represent the computation of the Stiffness Matrix for each of its element types.

##### 2.1.5 Derivatives Folder
Defines an auxiliary class for the computation of the computation of the Quadrilateral element Stiffness Matrix (and probably of the Cube element in the future too).

##### 2.1.6 Tests Folder
Contains `.cpp` files with each with a `main()` function to either test the funcionality of the library or to demonstrate its usage
