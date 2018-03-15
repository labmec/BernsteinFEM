#include "Moments.h"

#ifndef STIFFM_H
#define STIFFM_H

/*to implement yet
  About computing the Stiffness matrix
*/
class BStiff1D : public BMoment1D
{

};

class BStiff2DTri : public BMoment2DTri
{
    int qS;
    int nS;
public:
    // check about the 'q' value
    BStiff2DTri (int q, int n) : BMoment2DTri (q, 2*n - 2)
    {
        qS = q;
        nS = n;
    }

    BStiff2DTri (int q, int n, double T[][2]) : BMoment2DTri (q, 2*n - 2, T)
    {    
        qS = q;
        nS = n;
    }
};

class BStiff2DQuad : public BMoment2DQuad
{

};

class BStiff3D : public BMoment3D
{

};

#endif