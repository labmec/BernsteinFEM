#pragma once

/* Defines functions to get the quadrature points and weights */

#define MAX_QUADRA_ORDER 80

double jacobi_xi(int q, int i);
double jacobi_w(int q, int i);

double legendre_xi(int q, int i);
double legendre_w(int q, int i);

double jacobi2_xi(int q, int i);
double jacobi1_w(int q, int i);
