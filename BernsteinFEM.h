/*************** BernstienFEM.h ***************
 * This header file defines namespaces to be
 * able to work under the same identifiers
 * for classes work for different types
 * of elements
 * 
 * type the following before your code:
 * using namespace BernsteinFEM1D;
 *
 *********************************************/

#include <armadillo>

#include "Moments.h"
#include "MassM.h"
#include "StiffM.h"
#include "Elem.h"


namespace BernsteinFEM1D {
    typedef BMoment1D BMoment;
    typedef BMass1D Mass;
    typedef BStiff1D Stiffness;
    typedef BElement1D Element;
}

namespace BernsteinFEM2DQuad {
    typedef BMoment2DQuad BMoment;
    typedef BMass2DQuad Mass;
    typedef BStiff2DQuad Stiffness;
    typedef BElement2DQuad Element;
}

namespace BernsteinFEM2DTri {
    typedef BMoment2DTri BMoment;
    typedef BMass2DTri Mass;
    typedef BStiff2DTri Stiffness;
    typedef BElement2DTri Element;
}