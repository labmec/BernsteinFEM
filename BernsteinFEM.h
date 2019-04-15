/*************** BernstienFEM.h ***************
 * This header wraps all of the include files
 * needed to work with this library
 *********************************************/

#include <armadillo>

#include "Moments.h"
#include "MassM.h"
#include "StiffM.h"
#include "Evaluation.h"
#include "Derivatives.h"

namespace BernsteinFEM {
    // Typedefs to simplify using Load Vectors
    using BLoad1D = BMoment1D;
    using BLoad2DQuad = BMoment2DQuad;
    using BLoad2DTri = BMoment2DTri;
    using BLoad3DCube = BMoment3DCube;
    using BLoad3DTetra = BMoment3DTetra;

    namespace BernsteinFEM_1D
    {
    typedef BMoment1D BMoment;
    typedef BMass1D Mass;
    typedef BStiff1D Stiffness;
    } // namespace BernsteinFEM_1D

    namespace BernsteinFEM_2D_Quad
    {
    typedef BMoment2DQuad BMoment;
    typedef BMass2DQuad Mass;
    typedef BStiff2DQuad Stiffness;
    using namespace QuadD;
    } // namespace BernsteinFEM_2D_Quad

    namespace BernsteinFEM_2D_Tri
    {
    typedef BMoment2DTri BMoment;
    typedef BMass2DTri Mass;
    typedef BStiff2DTri Stiffness;
    } // namespace BernsteinFEM_2D_Tri
} // BernsteinFEM



