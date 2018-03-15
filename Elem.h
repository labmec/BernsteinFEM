#ifndef ELEM_H
#define ELEM_H

#include "MassM.h"
#include "StiffM.h"

class BElement1D {
    int q;
    int n;
    double *BBVector;
    double *MassFval;
    double *StiffFval;
    BMass1D MassMat;
    BStiff1D StiffMat;
public:
};

class BElement2DTri {
    int q;
    int n;
    double *BBVector;
    double *MassFunctionValues;
    double *StiffFunctionValues;
    BMass2DTri MassMat;
    BStiff2DTri StiffMat;
};

class BElement2DQuad {
    int q;
    int n;
    double *BBVector;
    double *MassFunctionValues;
    double *StiffFunctionValues;
    BMass2DQuad MassMat;
    BStiff2DQuad StiffMat;
};

#endif