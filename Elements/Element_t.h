/***************************************************************
 * Author: Lucas B. Andrade
 * 
 * This file contains the definition of the enumeration class
 * Element_t which indicates the kind of elements this library
 * can work with.
 * 
 * The types of elements are in respect to its geometrical
 * properties (e.g. LinearEl is the 1D line element)
 * 
 * This enumeration is used throughout the library as a template
 * parameter making specializations for generic classes
 ***************************************************************/

#pragma once

enum class Element_t
{
    LinearEl,
    QuadrilateralEl,
    TriangularEl,
    CubeEl,
    TetrahedronEl
    // Add whatever kind of element here, then specialize the template parameter in its own .cpp file
};
