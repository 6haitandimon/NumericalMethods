#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include "gauss/gauss.h"

const int iter_max = 100;

typedef double(*pf)(double*, double*, double, double);

double f1(double *uk1, double *uk, double t, double Tau);
double f2(double *uk1, double *uk, double t, double Tau);


double Differential(pf f, double *uk1, double *uk, double t, double Tau, int n);
double* Neuton(double *yk_plus, double *yk, double tk, double Tau, int n);